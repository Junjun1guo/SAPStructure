########################################################################################################################
# -*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
#####Units: Length-mm, Force-N, mass-ton, Stress-Mpa, g=9810mm/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-09-06
########################################################################################################################
########################################################################################################################
# import necessary modules
import numpy as np
import math
import os
import matplotlib.pyplot as plt
########################################################################################################################
########################################################################################################################
def tz_sand(phiDegree,d,gammaEff,z,elelength):
    """
    --------------------------------------------------------------------------------------------------------------------
    Method to compute ultimate resistance, tult, and displacement at 50% mobilization of tult, z50, for use in t-z
    curves for cohesionless soil (sand)
    --------------------------------------------------------------------------------------------------------------------
    references
	    Mosher, R.L. (1984). “Load transfer criteria for numerical analysis of axial loaded piles in sand.” U.S. Army
	    Engineering and Waterways Experimental Station, Automatic Data Processing Center, Vicksburg, Miss.
	    Kulhawy, F.H. (1991). "Drilled shaft foundations." Foundation engineering handbook, 2nd Ed., Chap 14, H.-Y.
	    Fang ed., Van Nostrand Reinhold, New York
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        phiDeree(degree)-friction angle
	    d(m)-pile diameter
	    gammaEff(kN/m3)-the effective density
	    z(m)-the depth from the ground surface
	    elelength(m)-length of the computed pile element
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    # Compute tult based on tult = Ko*sigV*pi*d*tan(delta), where
    #   Ko    is coeff. of lateral earth pressure at rest,
    #         taken as Ko = 0.4
    #   delta is interface friction between soil and pile,
    #         taken as delta = 0.8*phi to be representative of a
    #         smooth precast concrete pile after Kulhawy (1991)
    delta = 0.8 * phiDegree * np.pi / float(180.0)
    ##vertical effective stress at current depth
    sigamaV = gammaEff * z
    # if z = 0 (ground surface) need to specify a small non-zero value of sigV
    if sigamaV == 0.0:
        sigamaV = 0.01
    tu = 0.4 * sigamaV * np.pi * d * np.tan(delta)
    # TzSimple1 material formulated with tult as force, not stress, multiply by tributary length of pile
    tult = tu * elelength
    # Mosher (1984) provides recommended initial tangents based on friction angle
    # values are in units of psf/in
    kf = [6000, 10000, 10000, 14000, 14000, 18000]
    fric = [28, 31, 32, 34, 35, 38]
    # determine kf for input value of phi, linear interpolation for intermediate values
    if phiDegree < fric[0]:
        k = kf[0]
    elif phiDegree > fric[5]:
        k = kf[5]
    else:
        for i in range(5):
            if (fric[i] <= phiDegree) and (phiDegree <= fric[i + 1]):
                k = (kf[i + 1] - kf[i]) / (fric[i + 1] - fric[i]) * (phiDegree - fric[i]) + kf[i]
    # need to convert kf to units of kN/m^3 1psf=0.04788 kn/m2, 1in=0.0254m
    kSIunits = k * 1.885
    # based on a t-z curve of the shape recommended by Mosher (1984), z50 = tult/kf
    # set z50 [expr $tult/$kSIunits] ;#wrong, it should be changed as follows
    z50 = tu / (np.pi * d * kSIunits)
    return tult, z50
########################################################################################################################
########################################################################################################################
def tz_clay(cu, d, gammaEff, z, elelength):
    """
    --------------------------------------------------------------------------------------------------------------------
    Method to compute ultimate resistance, tult, and displacement at 50% mobilization of tult, z50, for use in t-z
    curves for clay soil
    --------------------------------------------------------------------------------------------------------------------
    references
        S. Soltanieh, M.M. Memarpour?, F. Kilanehei.Performance assessment of bridge-soil-foundation system with irregular
        configuration considering ground motion directionality effect. Soil dynamics and Earthquake Engineering.118(2019).
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        cu(kPa)-undrained shear strength at depth Z
        d(m)-the pile diameter
        gammaEff(kN/m3)-submerged unit weight
        z(m)-the depth from the ground surface
        elelength(m)-length of the pile element
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    As = np.pi * d * elelength
    sigmaV = gammaEff * z
    alpha1 = 0.5 * ((cu / float(sigmaV)) ** (-0.5)) * (cu / float(sigmaV))
    alpha2 = 0.5 * ((cu / float(sigmaV)) ** (-0.25)) * (cu / float(sigmaV))
    if alpha1 <= 1:
        alpha = alpha1
    elif alpha2 > 1:
        alpha = alpha2
    tult = alpha * cu * As
    z50 = d * 0.0031
    return tult, z50
########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
    ####################################################################################################################
    pass
