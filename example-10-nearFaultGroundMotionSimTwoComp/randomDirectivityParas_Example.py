#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2026-07-06
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import shutil
import math
import pandas as pd
import time
import openseespy.opensees as ops
import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
import numpy as np
from scipy.stats import truncnorm, weibull_min
import plotly.graph_objects as go
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
########################################################################################################################
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
########################################################################################################################
############################----import auxiliary modules---#############################################################
from auxiliaryModules.mainMod import NearFaultGroundMotionSim_twoComp_specificPulseParas
########################################################################################################################
########################################################################################################################
###---constants

SIGMA_LN_ZTOR = 0.860
SIGMA_LN_LR = 0.401
SIGMA_LN_WR_SURFACE = 0.278
SIGMA_LN_WR_BURIED = 0.286
RHO_LR_WR_BURIED = 0.244

WEIBULL = {
    0: (0.626, 1.0 / 3.921),  # SS
    1: (0.692, 1.0 / 3.394),  # REV
    2: (0.692, 1.0 / 3.394),  # N
}

TABLE2_TARGET = {
    "LR": (7.9, 18.4),
    "WR": (5.2, 20.3),
    "ZTOR": (3.1, 9.3),
    "RRUP": (6.4, 10.8),
    "d": (1.7, 17.0),
    "f": (0.4, 16.8),
}
########################################################################################################################
########################################################################################################################
def wr_upper_bound(ZTOR, d_deg, WF, ZTOF=0.0):
    return max(WF - (ZTOR - ZTOF) / np.sin(np.radians(d_deg)), 0.0)
########################################################################################################################
########################################################################################################################
def ztor_upper_bound(Mw, F):
    ub = np.exp(mean_ln_ztor(Mw, F) + 1.75 * SIGMA_LN_ZTOR)
    return min(ub, 15.0) if F == 1 else min(ub, 10.0)
########################################################################################################################
########################################################################################################################
def mean_ln_ztor(Mw, F):
    if F == 1:
        return 6.362 - 0.789 * Mw  ###-REV
    else:
        return 9.361 - 1.400 * Mw ###-SS and N
########################################################################################################################
########################################################################################################################
def mean_ln_lr(Mw, F):
    return {0: -5.821,1: -6.106,2: -5.919}[F] + 1.402 * Mw
########################################################################################################################
########################################################################################################################
def mean_ln_wr(Mw, F, surface):
    if surface:
        return {0: 0.860,1: 1.591,2: 0.958}[F] + 0.263 * Mw ###-ZTOR=0
    else:
        return {0: -3.403,1: -3.211,2: -3.164}[F] + 0.907 * Mw ###-ZTOR>0
########################################################################################################################
########################################################################################################################
def prob_buried(Mw, F):
    """ZTOR>0"""
    if F == 1: ###-REV
        return 1.0 / (1.0 + np.exp(-17.220 + 2.293 * Mw))
    else: ###-SS and N
        return 1.0 / (1.0 + np.exp(-15.293 + 2.293 * Mw))
########################################################################################################################
########################################################################################################################
def simulate_rupture_geometry(Mw, F, LF, WF, ZTOF, d_deg, rng, nmax=10000):
    """"""
    p_buried = prob_buried(Mw, F)
    for _ in range(nmax):
        buried = rng.random() < p_buried
        if buried: ###-ZTOR>0
            mu = np.array([mean_ln_ztor(Mw, F),mean_ln_lr(Mw, F),mean_ln_wr(Mw, F, surface=False)])
            s_l = SIGMA_LN_LR
            s_w = SIGMA_LN_WR_BURIED
            cov = np.array([[SIGMA_LN_ZTOR ** 2, 0.0, 0.0],
                [0.0, s_l ** 2, RHO_LR_WR_BURIED * s_l * s_w],
                [0.0, RHO_LR_WR_BURIED * s_l * s_w, s_w ** 2]]) ###-convariance matrix
            lnZ, lnL, lnW = rng.multivariate_normal(mu, cov) ###-multi-Gaussian sampling
            ZTOR = float(np.exp(lnZ))
            LR = float(np.exp(lnL))
            WR = float(np.exp(lnW))
            if not (0.0 < ZTOR <= ztor_upper_bound(Mw, F)):
                continue
        else: ###-ZTOR==0
            ZTOR = 0.0
            LR = float(np.exp(rng.normal(mean_ln_lr(Mw, F), SIGMA_LN_LR)))
            WR = float(np.exp(rng.normal(mean_ln_wr(Mw, F, True), SIGMA_LN_WR_SURFACE)))
        if LR > LF:
            continue
        wr_max = wr_upper_bound(ZTOR, d_deg, WF, ZTOF)
        if not (0.0 < WR <= wr_max):
            continue
        return {"buried": buried, "ZTOR": ZTOR, "LR": LR, "WR": WR}
    raise RuntimeError("Rupture geometry sampling did not converge.")
########################################################################################################################
########################################################################################################################
def unit(v):
    """retur normalized vector"""
    n = np.linalg.norm(v)
    return np.zeros_like(v) if n <= 1e-12 else v / n
########################################################################################################################
########################################################################################################################
def fault_basis(top_pts, d_deg):
    """"""
    p0 = np.asarray(top_pts[0], dtype=float)
    p1 = np.asarray(top_pts[1], dtype=float)
    sv = p1 - p0 ###-strike vector
    sv[2] = 0.0
    sv = unit(sv) ###-normalized vector
    ph = unit(np.array([-sv[1], sv[0], 0.0], dtype=float)) ###-perpendicular to sv,horizontal plane
    dip = np.radians(d_deg) ###-dip angle radian
    dv = unit(ph * np.cos(dip) + np.array([0.0, 0.0, np.sin(dip)], dtype=float))
    return p0.copy(), sv, dv, ph
########################################################################################################################
########################################################################################################################
def simulate_rupture_plane(top_pts, LF, LR, WR, ZTOR, d_deg, ZTOF, rng):
    """"""
    origin, sv, dv, ph = fault_basis(top_pts, d_deg)
    dip_offset = (ZTOR - ZTOF) / np.sin(np.radians(d_deg)) ###-dip plane distance
    strike_offset = rng.uniform(0.0, max(LF - LR, 0.0)) ###-random ruputre start position
    top_left = origin + strike_offset * sv + dip_offset * dv ###-global cooridinate
    corners = np.array([top_left,top_left + LR * sv,top_left + LR * sv + WR * dv,top_left + WR * dv]) ###-(P0,P1,P2,P3)
    return {"corners": corners,"top_left": top_left,"sv": sv,"dv": dv,"ph": ph}
########################################################################################################################
########################################################################################################################
def sample_truncnorm(mu, sigma, lo, hi, rng):
    """"""
    if sigma <= 1e-12:
        return float(np.clip(mu, lo, hi))
    a = (lo - mu) / sigma ###-normalized a
    b = (hi - mu) / sigma
    return float(truncnorm.rvs(a, b, loc=mu, scale=sigma, random_state=rng))
########################################################################################################################
########################################################################################################################
def sample_truncweibull(shape, scale, lo, hi, rng, nmax=5000):
    """"""
    x = lo
    for _ in range(nmax):
        x = weibull_min.rvs(shape, scale=scale, random_state=rng)
        if lo <= x <= hi:
            return float(x)
    return float(np.clip(x, lo, hi)) ###-limited between low and high
########################################################################################################################
########################################################################################################################
def simulate_hypocenter(F, LR, WR, top_left, sv, dv, rng):
    """"""
    hs = sample_truncnorm(0.5 * LR, 0.23 * LR, 0.0, LR, rng) ###-strike direction
    shape, scale = WEIBULL[F]
    hd = WR * sample_truncweibull(shape, scale, 1e-6, 1.0 - 1e-6, rng)
    hyp = top_left + hs * sv + hd * dv
    return {"hs": hs, "hd": hd, "hyp": hyp}
########################################################################################################################
########################################################################################################################
def point_to_triangle_distance(p, a, b, c):
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = np.dot(ab, ap)
    d2 = np.dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0: ###
        return np.linalg.norm(ap) ###-neearst point is a
    bp = p - b
    d3 = np.dot(ab, bp)
    d4 = np.dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return np.linalg.norm(bp)
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        proj = a + (d1 / (d1 - d3)) * ab
        return np.linalg.norm(p - proj)
    cp = p - c
    d5 = np.dot(ab, cp)
    d6 = np.dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return np.linalg.norm(cp)
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        proj = a + (d2 / (d2 - d6)) * ac
        return np.linalg.norm(p - proj)
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        proj = b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b)
        return np.linalg.norm(p - proj)
    denom = 1.0 / (va + vb + vc)
    proj = a + ab * (vb * denom) + ac * (vc * denom)
    return np.linalg.norm(p - proj)
########################################################################################################################
########################################################################################################################
def calc_rrup(site, corners):
    """calculate Rrup"""
    a, b, c, d = corners
    return min(point_to_triangle_distance(site, a, b, c),point_to_triangle_distance(site, a, c, d))
########################################################################################################################
########################################################################################################################
def acute_angle_deg_strike(v1, v2, eps=1e-12):
    """
    返回两个三维向量水平投影后的锐角夹角，单位为度，范围 [0, 90]。
    """
    v1 = np.asarray(v1, dtype=float).copy()
    v2 = np.asarray(v2, dtype=float).copy()
    v1[2] = 0.0
    v2[2] = 0.0
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < eps or n2 < eps:
        return 0.0
    cosang = np.dot(v1, v2) / (n1 * n2)
    cosang = np.clip(cosang, -1.0, 1.0)
    # 取锐角
    ang = np.degrees(np.arccos(abs(cosang)))
    return float(min(ang, 90.0))
########################################################################################################################
########################################################################################################################
def closest_point_on_segment_strike(p, a, b, eps=1e-12):
    """"""
    p = np.asarray(p, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ab = b - a
    ab2 = np.dot(ab, ab)
    if ab2 < eps:
        return a.copy()
    t = np.dot(p - a, ab) / ab2
    t = np.clip(t, 0.0, 1.0)
    return a + t * ab
########################################################################################################################
########################################################################################################################
def calc_directivity_strike(site, hyp, corners, sv, eps=1e-12):
    """"""
    site = np.asarray(site, dtype=float)
    hyp = np.asarray(hyp, dtype=float)
    sv = np.asarray(sv, dtype=float)
    top_left = np.asarray(corners[0], dtype=float)
    top_right = np.asarray(corners[1], dtype=float)
    bottom_right = np.asarray(corners[2], dtype=float)
    bottom_left = np.asarray(corners[3], dtype=float)
    top_edge = top_right - top_left
    left_edge = bottom_left - top_left
    lr = np.linalg.norm(top_edge) ###
    if lr < eps:
        raise ValueError("Invalid fault geometry: top edge length is zero.")
    left_edge_len2 = np.dot(left_edge, left_edge)
    if left_edge_len2 < eps:
        raise ValueError("Invalid fault geometry: left edge length is zero.")
    alpha = np.dot(hyp - top_left, left_edge) / left_edge_len2
    alpha = np.clip(alpha, 0.0, 1.0)
    left_pt = top_left + alpha * left_edge
    right_pt = top_right + alpha * (bottom_right - top_right)
    left_g = left_pt.copy()
    right_g = right_pt.copy()
    epic = hyp.copy()
    site_g = site.copy()
    left_g[2] = 0.0
    right_g[2] = 0.0
    epic[2] = 0.0
    site_g[2] = 0.0
    # site 到投影线段的最近点
    q = closest_point_on_segment_strike(site_g, left_g, right_g, eps=eps)
    # s = 最近点到震中的距离
    s = float(np.linalg.norm(q - epic))
    # theta = site-epic 与投影线段方向的锐角
    seg_dir = right_g - left_g
    theta_deg = acute_angle_deg_strike(site_g - epic, seg_dir, eps=eps)
    return s, theta_deg
########################################################################################################################
########################################################################################################################
def acute_angle_deg_dip(v1, v2, eps=1e-12):
    """
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < eps or n2 < eps:
        return 0.
    cosang = np.dot(v1, v2) / (n1 * n2)
    cosang = np.clip(cosang, -1.0, 1.0)
    ang = np.degrees(np.arccos(cosang))  # 原始夹角，范围 [0, 180]
    if ang > 90.0:
        ang = 180.0 - ang
    return float(ang)
########################################################################################################################
########################################################################################################################
def closest_point_on_segment_dip(p, a, b, eps=1e-12):
    """"""
    p = np.asarray(p, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ab = b - a
    ab2 = np.dot(ab, ab)
    if ab2 < eps:
        return a.copy()
    t = np.dot(p - a, ab) / ab2
    t = np.clip(t, 0.0, 1.0)
    return a + t * ab
########################################################################################################################
########################################################################################################################
def project_to_x0_plane(point):
    """"""
    p = np.asarray(point, dtype=float).copy()
    p[0] = 0.0
    return p
########################################################################################################################
########################################################################################################################
def calc_directivity_dip(site, hyp, corners, eps=1e-12):
    """"""
    site = np.asarray(site, dtype=float)
    hyp = np.asarray(hyp, dtype=float)
    top_left = np.asarray(corners[0], dtype=float)
    bottom_left = np.asarray(corners[3], dtype=float)
    # 将断层侧边投影到 x=0 平面
    side_top = project_to_x0_plane(top_left)
    side_bottom = project_to_x0_plane(bottom_left)
    side_vec = side_bottom - side_top
    side_len2 = np.dot(side_vec, side_vec)
    if side_len2 < eps:
        raise ValueError("Invalid fault geometry: projected side edge length is zero.")
    # 将震源和场地投影到 x=0 平面
    hyp_p = project_to_x0_plane(hyp)
    site_p = project_to_x0_plane(site)
    # 震源投影点应在线段上；这里做一次投影夹取，避免数值误差
    t_hyp = np.dot(hyp_p - side_top, side_vec) / side_len2
    t_hyp = np.clip(t_hyp, 0.0, 1.0)
    hyp_on_seg = side_top + t_hyp * side_vec
    # site 投影点到投影线段的最近点
    q = closest_point_on_segment_dip(site_p, side_top, side_bottom, eps=eps)
    # d 为最近点与震源投影点的距离
    d = float(np.linalg.norm(q - hyp_on_seg))
    # pha 为 site 投影点与震源投影点连线和侧边投影线段的夹角，取锐角
    pha_deg = acute_angle_deg_dip(site_p - hyp_on_seg, side_vec, eps=eps)
    return d, pha_deg
########################################################################################################################
########################################################################################################################
def run_one_p1(Mw, F, top_pts, LF, WF, ZTOF, d_deg, site, rng):
    """"""
    ###-return {"buried": buried, "ZTOR": ZTOR, "LR": LR, "WR": WR}
    geom = simulate_rupture_geometry(Mw, F, LF, WF, ZTOF, d_deg, rng)
    plane = simulate_rupture_plane(top_pts=top_pts,LF=LF,LR=geom["LR"],WR=geom["WR"],ZTOR=geom["ZTOR"],d_deg=d_deg,
                                   ZTOF=ZTOF,rng=rng)
    hypo = simulate_hypocenter(F=F,LR=geom["LR"],WR=geom["WR"],top_left=plane["top_left"],sv=plane["sv"],dv=plane["dv"],
                               rng=rng)
    rrup = calc_rrup(site, plane["corners"])
    out = {"ZTOR": geom["ZTOR"],"LR": geom["LR"],"WR": geom["WR"],"RRUP": rrup,"corners": plane["corners"],
           "hyp": hypo["hyp"]}
    if F == 0:
        s_eff, u_deg = calc_directivity_strike(site, hypo["hyp"], plane["corners"], plane["sv"])
        out["sord"] = s_eff
        out["thetaorpha"] = u_deg
    else:
        d_eff, f_deg = calc_directivity_dip(site, hypo["hyp"], plane["corners"])
        out["sord"] = d_eff
        out["thetaorpha"] = f_deg
    return out
########################################################################################################################
########################################################################################################################
def run_p1_batch(Mw, F, top_pts, LF, WF, ZTOF, d_deg, site, n_sims=100, seed=42):
    """"""
    rng = np.random.default_rng(seed)
    return [run_one_p1(Mw, F, top_pts, LF, WF, ZTOF, d_deg, site, rng) for _ in range(n_sims)]
########################################################################################################################
########################################################################################################################
def summarize_results(results, F):
    """"""
    keys = ["LR", "WR", "ZTOR", "RRUP","sord","thetaorpha"]
    out = {}
    for k in keys:
        vals = np.array([r[k] for r in results], dtype=float)
        out[k] = {"min": float(vals.min()),"max": float(vals.max()),"mean": float(vals.mean()),"std": float(vals.std(ddof=0))}
    return out
########################################################################################################################
########################################################################################################################
def print_summary(summary, F):
    """"""
    units = {"LR": "km","WR": "km","ZTOR": "km","RRUP": "km","sord": "km","thetaorpha": "deg"}
    keys = ["LR", "WR", "ZTOR", "RRUP","sord","thetaorpha"]
    print("=" * 84)
    print(f"{'Parameter':<14}{'Min':>12}{'Max':>12}{'Mean':>12}{'Std':>12}")
    print("=" * 84)
    for k in keys:
        s = summary[k]
        print(f"{(k + ' (' + units[k] + ')'):<14}" f"{s['min']:>12.3f}{s['max']:>12.3f}{s['mean']:>12.3f}{s['std']:>12.3f}")
    print("=" * 84)
########################################################################################################################
########################################################################################################################
def fault_plane_polygon(top_pts, WF, d_deg):
    """"""
    p0 = np.asarray(top_pts[0], dtype=float)
    p1 = np.asarray(top_pts[1], dtype=float)
    _, _, dv, _ = fault_basis(top_pts, d_deg)
    p2 = p1 + WF * dv
    p3 = p0 + WF * dv
    return np.array([p0, p1, p2, p3])
########################################################################################################################
########################################################################################################################
def rectangle_edges_trace(corners, name=None, color="rgba(80,80,80,0.9)", width=4, dash=None, showlegend=False):
    """"""
    closed = np.vstack([corners, corners[0]])
    return go.Scatter3d(x=closed[:, 0],y=closed[:, 1],z=closed[:, 2],mode="lines",line=dict(color=color, width=width,
            dash=dash),name=name,showlegend=showlegend,hoverinfo="skip")
########################################################################################################################
########################################################################################################################
def rectangle_mesh_trace(corners, color="rgba(160,160,160,0.45)", name=None, showlegend=False):
    """"""
    x = corners[:, 0]
    y = corners[:, 1]
    z = corners[:, 2]
    i = [0, 0]
    j = [1, 2]
    k = [2, 3]
    return go.Mesh3d(x=x,y=y,z=z,i=i,j=j,k=k,color=color,opacity=0.45,flatshading=True,name=name,showlegend=showlegend,
                     hoverinfo="skip")
########################################################################################################################
########################################################################################################################
def plot_figure9_style_plotly(top_pts, WF, d_deg, site, results, n_show=12):
    """"""
    fig = go.Figure()
    fault_poly = fault_plane_polygon(top_pts, WF, d_deg)
    fig.add_trace(rectangle_edges_trace(fault_poly,name="Fault plane",color="black",width=6,dash="dash",showlegend=True))
    show_results = results[:min(n_show, len(results))]
    hyp_xyz = []
    for idx, r in enumerate(show_results):
        corners = r["corners"]
        hyp = r["hyp"]
        hyp_xyz.append(hyp)

        fig.add_trace(rectangle_mesh_trace(corners,color="rgba(150,150,150,0.55)",name="Rupture realization"
        if idx == 0 else None,showlegend=(idx == 0)))
        fig.add_trace(rectangle_edges_trace(corners,color="rgba(90,90,90,0.9)",width=3,showlegend=False))
    hyp_xyz = np.array(hyp_xyz)
    fig.add_trace(go.Scatter3d(x=hyp_xyz[:, 0],y=hyp_xyz[:, 1],z=hyp_xyz[:, 2],mode="markers",marker=dict(size=6,
        color="black", symbol="diamond"),name="Hypocenters"))
    fig.add_trace(go.Scatter3d(x=[site[0]],y=[site[1]],z=[site[2]],mode="markers",marker=dict(size=7, color="red",
        symbol="circle"),name="Site"))
    all_pts = [fault_poly, hyp_xyz, np.asarray([site])]
    for r in show_results:
        all_pts.append(r["corners"])
    all_pts = np.vstack(all_pts)
    xmin, ymin, zmin = all_pts.min(axis=0)
    xmax, ymax, zmax = all_pts.max(axis=0)
    span = max(xmax - xmin, ymax - ymin, zmax - zmin)
    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    zmid = 0.5 * (zmin + zmax)
    half = 0.55 * span if span > 0 else 1.0
    fig.update_layout(title="P1 example: selected rupture realizations and hypocenter locations",
        scene=dict(xaxis=dict(title="X along strike (km)", range=[xmid - half, xmid + half]),
            yaxis=dict(title="Y horizontal normal (km)", range=[ymid - half, ymid + half]),
            zaxis=dict(title="Depth z (km)",range=[zmid + half, zmid - half],),aspectmode="cube",
            camera=dict(eye=dict(x=1.6, y=-1.5, z=1.1))),legend=dict(x=0.01, y=0.99),margin=dict(l=0, r=0, t=50, b=0))
    return fig
########################################################################################################################
########################################################################################################################
import numpy as np
import matplotlib
matplotlib.use("TkAgg")   # 若本机无 Tk，可改为 "QtAgg"
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_fault_figure(result, F, top_pts, WF, d_deg, site, number,
                      figsize=(14, 6), savepath=None, show=True):
    """
    左图为三维示意图，右图为方向性投影图。

    约定：
    - 左图三维坐标显示为 (Y, X, Z)
    - 若 F == 0，右图为 Y-X 排布，且轴正向与左图三维显示一致
    - 若 F != 0，右图为 Y-Z 排布，且轴正向与左图三维显示一致
    """
    result = dict(result)
    corners = np.asarray(result["corners"], dtype=float)
    hyp = np.asarray(result["hyp"], dtype=float)
    site = np.asarray(site, dtype=float)
    top_pts = np.asarray(top_pts, dtype=float)

    ZTOR = float(result["ZTOR"])
    LR = float(result["LR"])
    WR = float(result["WR"])
    RRUP = float(result["RRUP"])
    sord = float(result["sord"])
    thetaorpha = float(result["thetaorpha"])

    d_rad = np.deg2rad(d_deg)

    def _symbol_label(ax, pt, text, color, xlim, ylim, fx=0.008, fy=0.008):
        xr = max(abs(xlim[1] - xlim[0]), 1.0)
        yr = max(abs(ylim[1] - ylim[0]), 1.0)
        ax.text(
            pt[0] + fx * xr,
            pt[1] + fy * yr,
            text,
            color=color,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.06", fc="white", ec="none", alpha=0.65),
            zorder=6,
        )

    def _add_info_box(ax, lines, loc=(0.98, 0.98)):
        text = "\n".join(lines)
        ax.text(
            loc[0],
            loc[1],
            text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
            family="monospace",
            linespacing=1.0,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.6", alpha=0.88),
            zorder=10,
        )

    # 1. 构造整个断层面的四个角点
    top_left = top_pts[0]
    top_right = top_pts[1]

    strike_vec = top_right - top_left
    fault_length = np.linalg.norm(strike_vec)
    if fault_length < 1e-12:
        raise ValueError("top_pts 的两个点不能重合。")

    strike_unit = strike_vec / fault_length

    dip_h = np.array([-strike_unit[1], strike_unit[0], 0.0], dtype=float)
    dip_h_norm = np.linalg.norm(dip_h)
    if dip_h_norm < 1e-12:
        raise ValueError("无法从 top_pts 计算走向方向。")
    dip_h = dip_h / dip_h_norm
    if dip_h[1] < 0:
        dip_h = -dip_h

    dip_unit = np.array([
        dip_h[0] * np.cos(d_rad),
        dip_h[1] * np.cos(d_rad),
        np.sin(d_rad)
    ], dtype=float)

    fault_bottom_left = top_left + WF * dip_unit
    fault_bottom_right = top_right + WF * dip_unit
    fault_corners = np.array([top_left, top_right, fault_bottom_right, fault_bottom_left])

    # 2. 统一几何范围
    all_pts_raw = np.vstack([fault_corners, corners, hyp[None, :], site[None, :]])
    xmin, xmax = all_pts_raw[:, 0].min(), all_pts_raw[:, 0].max()
    ymin, ymax = all_pts_raw[:, 1].min(), all_pts_raw[:, 1].max()
    zmin, zmax = all_pts_raw[:, 2].min(), all_pts_raw[:, 2].max()

    xrange_ = max(xmax - xmin, 1.0)
    yrange_ = max(ymax - ymin, 1.0)
    zrange_ = max(zmax - zmin, 1.0)

    xpad = 0.08 * xrange_
    ypad = 0.08 * yrange_
    zpad = 0.08 * zrange_

    xmin_p, xmax_p = xmin - xpad, xmax + xpad
    ymin_p, ymax_p = ymin - ypad, ymax + ypad
    zmin_p, zmax_p = zmin - zpad, zmax + zpad

    xmid = 0.5 * (xmin_p + xmax_p)
    ymid = 0.5 * (ymin_p + ymax_p)
    zmid = 0.5 * (zmin_p + zmax_p)
    radius = 0.5 * max(xmax_p - xmin_p, ymax_p - ymin_p, zmax_p - zmin_p)

    # 左图三维显示范围，显示顺序为 (Y, X, Z)
    x3d_lim = (ymid - radius, ymid + radius)
    y3d_lim = (xmid - radius, xmid + radius)
    z3d_lim = (zmid - radius, zmid + radius)

    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2)

    # 显式启用三维图鼠标交互：左键旋转，中键/右键依后端支持缩放
    ax1.mouse_init()

    # 3. 左图：三维图，交换 X/Y 轴显示为 (Y, X, Z)
    fault_corners_xy = fault_corners[:, [1, 0, 2]]
    corners_xy = corners[:, [1, 0, 2]]
    hyp_xy = hyp[[1, 0, 2]]
    site_xy = site[[1, 0, 2]]

    fault_poly = Poly3DCollection(
        [fault_corners_xy],
        facecolors=(0, 0, 0, 0.0),
        edgecolors="0.45",
        linewidths=1.5,
        linestyles="--"
    )
    ax1.add_collection3d(fault_poly)

    rup_poly = Poly3DCollection(
        [corners_xy],
        facecolors=(0.95, 0.76, 0.30, 0.45),
        edgecolors="#c17c00",
        linewidths=2.0
    )
    ax1.add_collection3d(rup_poly)

    fc = np.vstack([fault_corners_xy, fault_corners_xy[0]])
    rc = np.vstack([corners_xy, corners_xy[0]])
    ax1.plot(fc[:, 0], fc[:, 1], fc[:, 2], ls="--", lw=1.5, c="0.45")
    ax1.plot(rc[:, 0], rc[:, 1], rc[:, 2], ls="-", lw=2.0, c="#c17c00")

    # 震源与场地点严格按三维显示坐标 (Y, X, Z) 放置
    ax1.scatter(
        hyp_xy[0], hyp_xy[1], hyp_xy[2],
        c="red", marker="*", s=180, label="Hypocenter", zorder=5
    )
    ax1.scatter(
        site_xy[0], site_xy[1], site_xy[2],
        c="#1565c0", marker="o", s=60, label="Site", zorder=5
    )

    ax1.set_xlim(*x3d_lim)
    ax1.set_ylim(*y3d_lim)
    ax1.set_zlim(*z3d_lim)
    ax1.invert_zaxis()
    ax1.view_init(elev=24, azim=-56)
    ax1.set_box_aspect([1, 1, 1])

    ax1.set_xlabel("Y (km)")
    ax1.set_ylabel("X (km)")
    ax1.set_zlabel("Depth Z (km)")
    ax1.set_title(
        f"3D Fault Geometry\n"
        f"Fault: L = {fault_length:.2f} km, W = {WF:.2f} km; "
        f"Rupture: L = {LR:.2f} km, W = {WR:.2f} km\n"
        f"ZTOR = {ZTOR:.2f} km, RRUP = {RRUP:.2f} km\n"
        f"H (Y,X,Z)=({hyp_xy[0]:.2f},{hyp_xy[1]:.2f},{hyp_xy[2]:.2f}) km; "
        f"S (Y,X,Z)=({site_xy[0]:.2f},{site_xy[1]:.2f},{site_xy[2]:.2f}) km"
    )
    ax1.legend(loc="upper left")

    if F == 0:
        # 走滑断层：投影到 z=0 平面，按 Y-X 排布显示
        tl, tr, br, bl = corners

        su = tr - tl
        su = su / np.linalg.norm(su)

        du = bl - tl
        du = du / np.linalg.norm(du)

        rel = hyp - tl
        xd = np.dot(rel, du)

        p0 = tl + xd * du
        p1 = p0 + LR * su

        p0_xy = np.array([p0[0], p0[1]])
        p1_xy = np.array([p1[0], p1[1]])
        hyp_xy2 = np.array([hyp[0], hyp[1]])
        site_xy2 = np.array([site[0], site[1]])

        seg = p1_xy - p0_xy
        seg2 = np.dot(seg, seg)
        if seg2 < 1e-12:
            q_xy = p0_xy.copy()
        else:
            t = np.dot(site_xy2 - p0_xy, seg) / seg2
            t = np.clip(t, 0.0, 1.0)
            q_xy = p0_xy + t * seg

        # 右图显示为 (Y, X)
        p0p = np.array([p0_xy[1], p0_xy[0]])
        p1p = np.array([p1_xy[1], p1_xy[0]])
        hypp = np.array([hyp_xy2[1], hyp_xy2[0]])
        sitep = np.array([site_xy2[1], site_xy2[0]])
        q = np.array([q_xy[1], q_xy[0]])

        ax2.plot(
            [p0p[0], p1p[0]], [p0p[1], p1p[1]],
            c="black", lw=2.2,
            label="Projected rupture line: rupture trace in Y-X view"
        )
        ax2.scatter(
            hypp[0], hypp[1],
            c="black", marker="*", s=160,
            label="H: hypocenter projection",
            zorder=5
        )
        ax2.scatter(
            sitep[0], sitep[1],
            c="#1565c0", marker="o", s=55,
            label="S: site projection",
            zorder=5
        )
        ax2.scatter(
            q[0], q[1],
            c="#ef6c00", marker="D", s=45,
            label="Q: closest point on rupture trace",
            zorder=5
        )

        ax2.plot(
            [hypp[0], sitep[0]], [hypp[1], sitep[1]],
            c="#c62828", lw=1.8,
            label="H-S line: from hypocenter to site"
        )
        ax2.plot(
            [sitep[0], q[0]], [sitep[1], q[1]],
            c="#ef6c00", lw=1.2, ls="--",
            label="S-Q line: shortest link from site to rupture"
        )

        ax2.set_xlim(*x3d_lim)
        ax2.set_ylim(*y3d_lim)

        _symbol_label(ax2, hypp, "H", "black", x3d_lim, y3d_lim)
        _symbol_label(ax2, sitep, "S", "#1565c0", x3d_lim, y3d_lim)
        _symbol_label(ax2, q, "Q", "#ef6c00", x3d_lim, y3d_lim)

        info_lines = [
            f"H=({hypp[0]:.2f},{hypp[1]:.2f})",
            f"S=({sitep[0]:.2f},{sitep[1]:.2f})",
            f"Q=({q[0]:.2f},{q[1]:.2f})",
            f"sd={sord:.2f}",
            f"th={thetaorpha:.2f}deg",
        ]
        _add_info_box(ax2, info_lines, loc=(0.98, 0.98))

        ax2.set_xlabel("Y (km)")
        ax2.set_ylabel("X (km)")
        ax2.set_title("Plan View (Projection to z=0, Y-X)")
        fig.suptitle(f"Strike-Slip Fault_{number}", y=0.98, fontsize=14)

    else:
        # 倾滑断层：投影到 x=0 平面，用 (Y, Z) 图表示
        side_a = corners[0]
        side_b = corners[3]

        side_a_p = np.array([side_a[1], side_a[2]])
        side_b_p = np.array([side_b[1], side_b[2]])
        hyp_p = np.array([hyp[1], hyp[2]])
        site_p = np.array([site[1], site[2]])

        seg = side_b_p - side_a_p
        seg2 = np.dot(seg, seg)
        if seg2 < 1e-12:
            q = side_a_p.copy()
        else:
            t = np.dot(site_p - side_a_p, seg) / seg2
            t = np.clip(t, 0.0, 1.0)
            q = side_a_p + t * seg

        ax2.plot(
            [side_a_p[0], side_b_p[0]], [side_a_p[1], side_b_p[1]],
            c="black", lw=2.2,
            label="Projected rupture side: rupture edge in Y-Z view"
        )
        ax2.scatter(
            hyp_p[0], hyp_p[1],
            c="black", marker="*", s=160,
            label="H: hypocenter projection",
            zorder=5
        )
        ax2.scatter(
            site_p[0], site_p[1],
            c="#1565c0", marker="o", s=55,
            label="S: site projection",
            zorder=5
        )
        ax2.scatter(
            q[0], q[1],
            c="#ef6c00", marker="D", s=45,
            label="Q: closest point on rupture side",
            zorder=5
        )

        ax2.plot(
            [hyp_p[0], site_p[0]], [hyp_p[1], site_p[1]],
            c="#c62828", lw=1.8,
            label="H-S line: from hypocenter to site"
        )
        ax2.plot(
            [site_p[0], q[0]], [site_p[1], q[1]],
            c="#ef6c00", lw=1.2, ls="--",
            label="S-Q line: shortest link from site to rupture"
        )

        yz_xlim = x3d_lim
        yz_ylim = (z3d_lim[1], z3d_lim[0])

        ax2.set_xlim(*yz_xlim)
        ax2.set_ylim(*yz_ylim)

        _symbol_label(ax2, hyp_p, "H", "black", yz_xlim, yz_ylim)
        _symbol_label(ax2, site_p, "S", "#1565c0", yz_xlim, yz_ylim)
        _symbol_label(ax2, q, "Q", "#ef6c00", yz_xlim, yz_ylim)

        info_lines = [
            f"H=({hyp_p[0]:.2f},{hyp_p[1]:.2f})",
            f"S=({site_p[0]:.2f},{site_p[1]:.2f})",
            f"Q=({q[0]:.2f},{q[1]:.2f})",
            f"sd={sord:.2f}",
            f"th={thetaorpha:.2f}deg",
        ]
        _add_info_box(ax2, info_lines, loc=(0.98, 0.98))

        ax2.set_xlabel("Y (km)")
        ax2.set_ylabel("Depth Z (km)")
        ax2.set_title("Directivity View (Projection to x=0)")
        fig.suptitle(f"Dip-Slip Fault_{number}", y=0.98, fontsize=14)

    ax2.set_aspect("equal", adjustable="box")
    ax2.grid(True, ls=":", alpha=0.5)
    ax2.legend(loc="best", fontsize=8.5)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig, (ax1, ax2)
########################################################################################################################
if __name__ == "__main__":
    Mw = 5.5
    F = 1
    Vs30 = 390.0
    LF = 18.7
    WF = 20.4
    ZTOF = 3.1
    d_deg = 85.0

    top_pts = [
        [0.0, 0.0, ZTOF],
        [LF, 0.0, ZTOF]]

    site = np.array([9.0, -5, 0.0], dtype=float)
    #######################################
    # n_sims = 100
    # seed = None
    # results = run_p1_batch(Mw=Mw,F=F,top_pts=top_pts,LF=LF,WF=WF,ZTOF=ZTOF,d_deg=d_deg,site=site,n_sims=n_sims,seed=seed)
    # with open("simulationParas.pkl",'wb') as f:
    #     pickle.dump(results, f)
    ####################################################################################################################
    ####################################################################################################################
    # summary = summarize_results(results, F=F)
    # print("\nSimulation summary")
    # print_summary(summary, F=F)
    # ####################################################################################################################
    # print("\nFirst 10 realizations")
    # print("=" * 78)
    # print(f"{'#':<4}{'ZTOR':>10}{'LR':>10}{'WR':>10}{'RRUP':>10}{'sord':>10}{'thetaorpha':>10}")
    # print("=" * 78)
    # for i, r in enumerate(results[:10], start=1):
    #     print(
    #         f"{i:<4}{r['ZTOR']:>10.3f}{r['LR']:>10.3f}{r['WR']:>10.3f}"
    #         f"{r['RRUP']:>10.3f}{r['sord']:>10.3f}{r['thetaorpha']:>10.3f}")
    # print("=" * 78)
    # ##################################################################################################################
    ####################################################################################################################
    # with open("simulationParas.pkl", "rb") as f:
    #     results = pickle.load(f)
    # for i1 in range(100):
    #     simParas=results[i1]
    #     print(simParas)
    #     plot_fault_figure(result=simParas, F=F, top_pts=top_pts, WF=WF, d_deg=d_deg, site=site,number=i1+1)







