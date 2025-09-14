# hover_info_stable.py
# 保存并运行：python hover_info_stable.py
# 依赖：pyvista, pyvistaqt, PyQt5, vtk, numpy

import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter
import vtk
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel
from PyQt5.QtCore import Qt, QTimer

def create_demo_data():
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.5, 0.0],
        [6.5, -0.5, 0.2],
        [9.0, 1.2, 0.0],
        [12.0, -0.8, 0.0],
    ])
    segments = [(0, 1), (1, 2), (2, 3), (3, 4)]
    return nodes, segments

def build_polydata(nodes, segments):
    points_pd = pv.PolyData(nodes)
    lines_list = []
    for (i, j) in segments:
        lines_list.extend([2, i, j])
    lines = np.array(lines_list, dtype=np.int64)
    lines_pd = pv.PolyData(nodes, lines=lines)
    return points_pd, lines_pd

def _get_plotter_parent(plotter):
    parent = None
    if hasattr(plotter, "window"):
        try:
            w = plotter.window
            parent = w() if callable(w) else w
        except Exception:
            parent = None
    if parent is None and hasattr(plotter, "app"):
        try:
            parent = plotter.app.activeWindow()
        except Exception:
            parent = None
    return parent

def closest_distance_between_segments(p1, q1, p2, q2, eps=1e-12):
    u = q1 - p1
    v = q2 - p2
    w = p1 - p2
    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a * c - b * b
    if D < eps:
        s = 0.0
        if c > eps:
            t = e / c
            t = np.clip(t, 0.0, 1.0)
        else:
            t = 0.0
    else:
        s = (b * e - c * d) / D
        t = (a * e - b * d) / D
        s = np.clip(s, 0.0, 1.0)
        t = np.clip(t, 0.0, 1.0)
    pc = p1 + s * u
    qc = p2 + t * v
    dist = np.linalg.norm(pc - qc)
    return dist, s, t

class HoverDialog(QDialog):
    """简洁的提示框，重用并靠位置显示"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.Tool | Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        self.setAttribute(Qt.WA_ShowWithoutActivating)
        self.label = QLabel("")
        self.label.setWordWrap(True)
        layout = QVBoxLayout()
        layout.setContentsMargins(8, 6, 8, 6)
        layout.addWidget(self.label)
        self.setLayout(layout)
        self.setStyleSheet("""
            QDialog { background-color: rgba(255,255,240,230); border: 1px solid #999; }
            QLabel { color: #222; }
        """)
    def setText(self, text):
        self.label.setText(text)
        self.adjustSize()

def main():
    nodes, segments = create_demo_data()
    points_pd, lines_pd = build_polydata(nodes, segments)

    pv.set_plot_theme("document")
    plotter = BackgroundPlotter(window_size=(900, 600), title="稳定的悬停提示（延迟显示，防止闪烁）")

    plotter.add_mesh(points_pd, render_points_as_spheres=True, point_size=12, color="navy", name="nodes")
    plotter.add_mesh(lines_pd, color="black", line_width=4, name="segments")
    plotter.add_text("鼠标悬停在节点或线段上显示信息（延迟显示以防闪烁）",
                     position="upper_left", font_size=11, name="info")

    point_picker = vtk.vtkPointPicker()
    cell_picker = vtk.vtkCellPicker()
    try:
        cell_picker.SetTolerance(0.01)
    except Exception:
        pass

    bbox_diag = np.linalg.norm(nodes.max(axis=0) - nodes.min(axis=0))

    current_dialog = None
    current_target = None   # ('node', idx) or ('seg', idx) or None
    scheduled_target = None  # 同上
    last_mouse_pos = (0, 0)

    # 延迟参数（可调）
    hover_delay_ms = 400   # 鼠标停留超过此时间才显示提示
    hide_delay_ms = 200    # 离开后延迟关闭，避免短时抖动

    # 定时器句柄状态（用于避免重复启动）
    show_timer_active = {"val": False}
    hide_timer_active = {"val": False}

    def pick_segment_by_ray(x, y, renderer):
        try:
            renderer.SetDisplayPoint(x, y, 0.0)
            renderer.DisplayToWorld()
            wp_near = np.array(renderer.GetWorldPoint()[:3])
            renderer.SetDisplayPoint(x, y, 1.0)
            renderer.DisplayToWorld()
            wp_far = np.array(renderer.GetWorldPoint()[:3])
        except Exception:
            return -1

        ray_p0 = wp_near
        ray_p1 = wp_far

        best_idx = -1
        best_dist = float("inf")
        threshold = max(1e-6, bbox_diag * 0.03)

        for k, (i, j) in enumerate(segments):
            a = nodes[i]; b = nodes[j]
            dist, s, t = closest_distance_between_segments(a, b, ray_p0, ray_p1)
            if dist < best_dist:
                best_dist = dist
                best_idx = k

        if best_dist <= threshold:
            return best_idx
        return -1

    def schedule_show(target):
        """安排显示 target（tuple），如果已有计划且相同则不重复"""
        nonlocal scheduled_target
        if scheduled_target == target and show_timer_active["val"]:
            return
        scheduled_target = target
        show_timer_active["val"] = True
        # 单次延迟后调用 attempt_show_for_scheduled
        QTimer.singleShot(hover_delay_ms, attempt_show_for_scheduled)

    def cancel_scheduled_show():
        nonlocal scheduled_target
        scheduled_target = None
        show_timer_active["val"] = False

    def attempt_show_for_scheduled():
        """定时器回调：确认鼠标当前位置仍对应 scheduled_target 再显示"""
        nonlocal scheduled_target, current_dialog, current_target, show_timer_active
        show_timer_active["val"] = False
        if scheduled_target is None:
            return
        # 重新拾取当前位置，确保仍在目标上
        x, y = last_mouse_pos
        renderer = plotter.renderer
        picked_target = None
        try:
            if point_picker.Pick(x, y, 0.0, renderer):
                pid = point_picker.GetPointId()
                if pid is not None and pid >= 0:
                    picked_target = ("node", int(pid))
        except Exception:
            picked_target = None

        if picked_target is None:
            try:
                seg_idx = -1
                if cell_picker.Pick(x, y, 0.0, renderer):
                    cid = cell_picker.GetCellId()
                    if cid is not None and cid >= 0:
                        seg_idx = int(cid)
                if seg_idx < 0:
                    seg_idx = pick_segment_by_ray(x, y, renderer)
                if seg_idx is not None and seg_idx >= 0:
                    picked_target = ("seg", int(seg_idx))
            except Exception:
                picked_target = None

        # 只有当重新拾取的目标与计划目标一致时才真正显示
        if picked_target is not None and picked_target == scheduled_target:
            # 显示或更新 current_dialog
            if current_dialog is None:
                parent = _get_plotter_parent(plotter)
                current_dialog = HoverDialog(parent)
            ttype, tid = picked_target
            if ttype == "node":
                coord = nodes[tid]
                txt = f"节点索引: {tid}\n坐标: ({coord[0]:.6f}, {coord[1]:.6f}, {coord[2]:.6f})"
            else:
                i, j = segments[tid]
                p0 = nodes[i]; p1 = nodes[j]
                length = np.linalg.norm(p1 - p0)
                txt = f"线段索引: {tid}\n节点: {i} - {j}\n长度: {length:.6f}"
            # 设置文本并移动到鼠标附近
            try:
                current_dialog.setText(txt)
                # 将全局鼠标位置转为屏幕坐标（caller gives display coordinates relative to render window).
                # 使用 last_mouse_pos (renderer display coords) 与渲染器窗口的位置进行映射：
                # 简单方案：把 dialog 放在鼠标全局位置偏移处
                from PyQt5.QtGui import QCursor
                gp = QCursor.pos()
                current_dialog.move(gp.x() + 12, gp.y() + 18)
                current_dialog.show()
            except Exception:
                pass
            current_target = picked_target
            # 取消 hide 计划（如果有）
            hide_timer_active["val"] = False
            scheduled_target = None
        else:
            # 目标变化或鼠标已离开，撤销计划
            scheduled_target = None

    def schedule_hide():
        """安排隐藏当前对话框"""
        if hide_timer_active["val"]:
            return
        hide_timer_active["val"] = True
        QTimer.singleShot(hide_delay_ms, attempt_hide_if_still_needed)

    def cancel_hide():
        hide_timer_active["val"] = False

    def attempt_hide_if_still_needed():
        """定时器回调：如果当前鼠标下不再是 current_target，则隐藏"""
        nonlocal current_dialog, current_target, hide_timer_active
        hide_timer_active["val"] = False
        if current_target is None:
            return
        x, y = last_mouse_pos
        renderer = plotter.renderer
        picked_target = None
        try:
            if point_picker.Pick(x, y, 0.0, renderer):
                pid = point_picker.GetPointId()
                if pid is not None and pid >= 0:
                    picked_target = ("node", int(pid))
        except Exception:
            picked_target = None

        if picked_target is None:
            try:
                seg_idx = -1
                if cell_picker.Pick(x, y, 0.0, renderer):
                    cid = cell_picker.GetCellId()
                    if cid is not None and cid >= 0:
                        seg_idx = int(cid)
                if seg_idx < 0:
                    seg_idx = pick_segment_by_ray(x, y, renderer)
                if seg_idx is not None and seg_idx >= 0:
                    picked_target = ("seg", int(seg_idx))
            except Exception:
                picked_target = None

        if picked_target != current_target:
            try:
                if current_dialog is not None:
                    current_dialog.close()
                current_dialog = None
            except Exception:
                pass
            current_target = None

    def on_mouse_move(caller, event):
        nonlocal last_mouse_pos
        try:
            x, y = caller.GetEventPosition()
        except Exception:
            return
        last_mouse_pos = (x, y)

        renderer = plotter.renderer
        picked = None

        # 优先点拾取
        try:
            if point_picker.Pick(x, y, 0.0, renderer):
                pid = point_picker.GetPointId()
                if pid is not None and pid >= 0:
                    picked = ("node", int(pid))
        except Exception:
            picked = None

        # 线段拾取
        if picked is None:
            try:
                seg_idx = -1
                if cell_picker.Pick(x, y, 0.0, renderer):
                    cid = cell_picker.GetCellId()
                    if cid is not None and cid >= 0:
                        seg_idx = int(cid)
                if seg_idx < 0:
                    seg_idx = pick_segment_by_ray(x, y, renderer)
                if seg_idx is not None and seg_idx >= 0:
                    picked = ("seg", int(seg_idx))
            except Exception:
                picked = None

        # 决策逻辑：根据 picked、current_target、scheduled_target 调度显示或隐藏
        if picked is not None:
            # 如果当前已经显示且目标相同，更新位置和文本（平滑）
            if current_target == picked and current_dialog is not None:
                # 仅更新位置以跟随鼠标（不重绘文本）
                try:
                    from PyQt5.QtGui import QCursor
                    gp = QCursor.pos()
                    current_dialog.move(gp.x() + 12, gp.y() + 18)
                except Exception:
                    pass
                cancel_scheduled_show()
                cancel_hide()
                return
            # 如果计划显示的是同一目标，保持等待
            if scheduled_target == picked and show_timer_active["val"]:
                cancel_hide()
                return
            # 否则安排显示该目标
            cancel_hide()
            schedule_show(picked)
        else:
            # 未拾取到目标：取消任何计划显示；若当前对话框在显示，则安排隐藏
            cancel_scheduled_show()
            if current_dialog is not None:
                schedule_hide()

    # 注册鼠标移动回调
    plotter.iren.add_observer("MouseMoveEvent", on_mouse_move)

    # 进入 Qt 事件循环
    app = plotter.app
    if hasattr(app, "exec_"):
        app.exec_()
    else:
        app.exec()

if __name__ == "__main__":
    main()