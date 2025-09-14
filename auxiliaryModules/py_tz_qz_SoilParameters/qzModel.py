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
def qz_sand(phiDegree,d,gammaEff,z,G):
    """
    --------------------------------------------------------------------------------------------------------------------
    Method to compute ultimate resistance, tult, and displacement at 50% mobilization of tult, z50, for use in q-z
    curves for cohesionless soil (sand)
    --------------------------------------------------------------------------------------------------------------------
    references
        Meyerhof G.G. (1976). "Bearing capacity and settlement of pile foundations." J. Geotech. Eng. Div., ASCE, 102(3), 195-228.
        Vijayvergiya, V.N. (1977). “Load-movement characteristics of piles.” Proc., Ports 77 Conf., ASCE, New York.
        Kulhawy, F.H. ad Mayne, P.W. (1990). Manual on Estimating Soil Properties for Foundation Design. Electrical Power
        Research Institute. EPRI EL-6800, Project 1493-6 Final Report.
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        phiDeree(degree)-friction angle
	    d(m)-pile diameter
	    gammaEff(kN/m3)-the effective density
	    z(m)-the depth from the ground surface
	    G(kPa)-Shear modulus of soil at the pile tip
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    k0 = 1 - np.sin(phiDegree * np.pi / float(180.0))
    # ultimate tip pressure can be computed by qult = Nq*sigV after Meyerhof (1976)
    #  where Nq is a bearing capacity factor, phi is friction angle, and sigV is eff. overburden
    #  stress at the pile tip.
    phi = phiDegree * np.pi / float(180.0)
    sigmaV = gammaEff * z
    # rigidity index
    Ir = G / float(sigmaV * np.tan(phi))
    # bearing capacity factor
    Nq = (1 + 2 * k0) * (1.0 / float(3.0 - np.sin(phi))) * np.exp(np.pi / 2.0 - phi) \
         * ((np.tan(np.pi / 4.0 + phi / 2.0)) ** 2) * ((Ir) ** ((4 * np.sin(phi)) / float(3.0 * (1 + np.sin(phi)))))
    # tip resistance
    qu = Nq * sigmaV
    # QzSimple1 material formulated with qult as force, not stress, multiply by area of pile tip
    qult = qu * np.pi * d ** 2 / 4.0
    # the q-z curve of Vijayvergiya (1977) has the form, q(z) = qult*(z/zc)^(1/3)
    #  where zc is critical tip deflection given as ranging from 3-9% of the
    #  pile diameter at the tip.

    # assume zc is 5% of pile diameter
    # assume zc is 5% of pile diameter
    # based on Vijayvergiya (1977) curve, z50 = 0.125*zc
    zc = 0.05 * d
    # based on Vijayvergiya (1977) curve, z50 = 0.125*zc
    z50 = 0.125 * zc
    return qult, z50
########################################################################################################################
########################################################################################################################
def qz_clay(cu,d):
    """
    --------------------------------------------------------------------------------------------------------------------
    Method to compute ultimate resistance, tult, and displacement at 50% mobilization of tult, z50, for use in q-z
    curves for clay soil
    --------------------------------------------------------------------------------------------------------------------
    references
        S. Soltanieh, M.M. Memarpour?, F. Kilanehei.Performance assessment of bridge-soil-foundation system with irregular
        configuration considering ground motion directionality effect. Soil dynamics and Earthquake Engineering.118(2019).
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        cu(kPa)-the undrained shear strength
        d(m)-pile diameter
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    Ag = np.pi * d ** 2 / 4.0
    qult = 9 * cu * Ag
    z50 = 0.013 * d
    return qult, z50
########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
    ####################################################################################################################
    phiDegree = 30
    d = 2
    gammaEff = 19
    G = 1e5
    z = 50
    print(qz_sand(phiDegree, d, gammaEff, z, G))
