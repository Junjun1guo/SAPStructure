########################################################################################################################
# -*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
#####Units: Length-mm, Force-N, mass-ton, Stress-Mpa, g=9810mm/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-09-04
########################################################################################################################
###########################################################
#                                                         #
# Procedure to compute ultimate lateral resistance, p_u,  #
#  and displacement at 50% of lateral capacity, y50, for  #
#  p-y springs representing cohesionless soil.            #
#                                                         #
#   Created by:   Hyung-suk Shin                          #
#                 University of Washington                #
#   Modified by:  Chris McGann                            #
#                 Pedro Arduino                           #
#                 Peter Mackenzie-Helnwein                #
#                 University of Washington                #
#                 Xiaowei Wang,Tongji University          #
#                 Jingcheng Wang Fuzhou University        #
#                 Junjun Guo,Beijing Jiaotong University  #
#                                                         #
###########################################################
# references
#  American Petroleum Institute (API) (1987). Recommended Practice for Planning, Designing and
#   Constructing Fixed Offshore Platforms. API Recommended Practice 2A(RP-2A), Washington D.C,
#   17th edition.
#
# Brinch Hansen, J. (1961). ﷿0﷿3﷿﷿The ultimate resistance of rigid piles against transversal forces.﷿0﷿3﷿﷿
#  Bulletin No. 12, Geoteknisk Institute, Copenhagen, 59.
#
#  Boulanger, R. W., Kutter, B. L., Brandenberg, S. J., Singh, P., and Chang, D. (2003). Pile
#   Foundations in liquefied and laterally spreading ground during earthquakes: Centrifuge experiments
#   and analyses. Center for Geotechnical Modeling, University of California at Davis, Davis, CA.
#   Rep. UCD/CGM-03/01.
#
#  Reese, L.C. and Van Impe, W.F. (2001), Single Piles and Pile Groups Under Lateral Loading.
#    A.A. Balkema, Rotterdam, Netherlands.
########################################################################################################################
########################################################################################################################
# import necessary modules
import numpy as np
import math
import os
import matplotlib.pyplot as plt
########################################################################################################################
########################################################################################################################
def py_sand(pyDepth, gamma, phiDegree, b, pEleLength, puSwitch, ASwitch, kSwitch, gwtSwitch):
    """
    --------------------------------------------------------------------------------------------------------------------
    A method for calculating soil parameters of cohesiveless soil (sand)
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        pyDepth(m)-the depth of the soil (start from ground),deep soil corresponding to large pult
        gamma(kN/m3)-the effective density
        phiDegree(degrees)-friction angle
        b(m)-pile diameter
        pEleLength(m)-pile element size
        puSwitch(int)-select pult definition method for p-y curves (1 for API (default), 2 for Brinch Hansen)
        ASwitch(int)-values of pult depend of the loading coefficent, A, which can be defined in several ways.
                     (1-coefficient A for static loading obtained from Reese, 2011,2-coefficient A for cyclic loading
                     obtained from API, add by Lianxu 2019.04.05)
        kSwitch(int)-variation in coefficent of subgrade reaction with depth for p-y curves (1 for API linear variation
                     (default), 2 for modified API parabolic variation,3 for ATC method)
        gwtSwitch(int)-effect of ground water on subgrade reaction modulus for p-y curves (1 for above gwt, 2 for below gwt)
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    pi = 3.14159265358979
    phi = phiDegree * (pi / 180.0)
    zbRatio = pyDepth / b
    # obtain loading-type coefficient A (static loading )for given depth-to-diameter ratio zb
    #  ---> values are obtained from a figure [Reese,Single Piles and Pile Groups Under Lateral Loading, Figure3.24]
    #  and are therefore approximate
    if ASwitch == 1:
        zb = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0,
              2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4.0, 4.125,
              4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5.0]
        As = [2.846, 2.7105, 2.6242, 2.5257, 2.4271, 2.3409, 2.2546, 2.1437, 2.0575, 1.9589, 1.8973, 1.8111, 1.7372,
              1.6632, 1.5893, 1.5277, 1.4415, 1.3799, 1.3368, 1.269, 1.2074, 1.1581, 1.1211, 1.078, 1.0349, 1.0164,
              0.9979, 0.9733, 0.9610, 0.9487, 0.9363, 0.9117, 0.8994, 0.8994, 0.8871, 0.8871, 0.8809, 0.8809, 0.8809,
              0.8809, 0.8809]
        dataNum = 41
        ###---linear interpolation to define A for intermediate values of depth:diameter ratio
        for i in range(0, dataNum - 1):
            if (zb[i] <= zbRatio) and (zbRatio <= zb[i + 1]):
                A = (As[i + 1] - As[i]) / (zb[i + 1] - zb[i]) * (zbRatio - zb[i]) + As[i]
            elif zbRatio>=5:
                A = 0.88
    if ASwitch == 2:
        A = 0.9  # A=0.9 for cyclic loading [API]
    # -------------------------------------------------------------------------------------------------------------------------------
    #  2. Define ultimate lateral resistance, pult
    # -------------------------------------------------------------------------------------------------------------------------------
    # pult is defined per API recommendations (Reese and Van Impe, 2001 or API, 1987) for puSwitch = 1
    #  OR per the method of Brinch Hansen (1961) for puSwitch = 2
    # -------API recommended method------------------------------
    if puSwitch == 1:
        # define common terms
        alpha = 0.5 * phi
        beta = pi / 4.0 + phi / 2.0
        k0 = 0.4
        ka = (np.tan(pi / 4.0 - phi / 2.0)) ** 2
        # terms for Equation (3.44), Reese and Van Impe (2001)
        c1 = k0 * np.tan(phi) * np.sin(beta) / float(np.tan(beta - phi) * np.cos(alpha))
        c2 = np.tan(beta) / float(np.tan(beta - phi)) * np.tan(beta) * np.tan(alpha)
        c3 = k0 * np.tan(beta) * (np.tan(phi) * np.sin(beta) - np.tan(alpha))
        c4 = np.tan(beta) / float(np.tan(beta - phi)) - ka
        # terms for Equation (3.45), Reese and Van Impe (2001)
        c5 = ka * ((np.tan(beta)) ** 8 - 1.0)
        c6 = k0 * np.tan(phi) * (np.tan(beta)) ** 4
        ###---Equation (3.44), Reese and Van Impe (2001)
        pst = gamma * pyDepth * (pyDepth * (c1 + c2 + c3) + b * c4)
        ####---Equation (3.45), Reese and Van Impe (2001)
        psd = b * gamma * pyDepth * (c5 + c6)
        ###---pult is the lesser of pst and psd. At surface, an arbitrary value is defined
        if pst <= psd:
            if pyDepth == 0:
                pu = 0.01
            else:
                pu = A * pst
        else:
            pu = A * psd
        # PySimple1 material formulated with pult as a force, not force/length, multiply by trib. length
        pult = pu * pEleLength
    # ------- Brinch Hansen method ------------------------------
    if puSwitch == 2:
        # pressure at ground surface
        Kqo = (np.exp((0.5 * pi + phi) * np.tan(phi)) * np.cos(phi) * np.tan(pi / 4.0 + phi / 2.0) - np.exp(
            -(pi / 2.0 - phi) * np.tan(phi)) *
               np.cos(phi) * np.tan(pi / 4.0 - phi / 2.0))
        Kco = (1.0 / np.tan(phi)) * (
                    np.exp((pi / 2.0 + phi) * np.tan(phi)) * np.cos(phi) * np.tan(pi / 4.0 + phi / 2.0) - 1.0)
        # pressure at great depth
        dcinf = 1.58 + 4.09 * (np.tan(phi)) ** 4
        Nc = (1.0 / np.tan(phi)) * (np.exp(pi * np.tan(phi))) * ((np.tan(pi / 4.0 + phi / 2.0)) ** 2 - 1.0)
        Ko = 1.0 - np.sin(phi)
        Kcinf = Nc * dcinf
        Kqinf = Kcinf * Ko * np.tan(phi)
        ###---pressure at an arbitrary depth
        aq = (Kqo / (Kqinf - Kqo)) * (Ko * np.sin(phi) / np.sin(pi / 4.0 + phi / 2.0))
        KqD = (Kqo + Kqinf * aq * zbRatio) / (1.0 + aq * zbRatio)
        ###---ultimate lateral resistance
        if pyDepth == 0:
            pu = 0.01
        else:
            pu = gamma * pyDepth * KqD * b
        # PySimple1 material formulated with pult as a force, not force/length, multiply by trib. length
        pult = pu * pEleLength
    ##------------------------------------------------------------------------------------------------------------------
    #  3. Define displacement at 50% lateral capacity, y50
    # values of y50 depend of the coefficent of subgrade reaction, k, which can be defined in several ways.
    #  for gwtSwitch = 1, k reflects soil above the groundwater table
    #  for gwtSwitch = 2, k reflects soil below the groundwater table
    #  a linear variation of k with depth is defined for kSwitch = 1 after API (1987) (default)
    #  a parabolic variation of k with depth is defined for kSwitch = 2 after Boulanger et al. (2003)
    # ATC recommend method，kSwitch = 3  add by Lianxu 2019.04.05
    ##------------------------------------------------------------------------------------------------------------------
    # API (1987) recommended subgrade modulus for given friction angle, values obtained from figure (approximate)
    ph = [28.3, 29.5, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 50.0]
    # ---subgrade modulus above the water table
    if gwtSwitch == 1:
        # units of k are lb/in^3
        k = [10, 23, 45, 61, 80, 100, 120, 140, 160, 182, 215, 250, 275, 400]
    else:
        # subgrade modulus below the water table
        k = [10, 20, 33, 42, 50, 60, 70, 85, 95, 107, 122, 141, 155, 225]
    dataNum = 14
    # linear interpolation for values of phi not represented above
    for i in range(0, dataNum - 1):
        if (ph[i] <= phiDegree) and (phiDegree <= ph[i + 1]):
            khat = (k[i + 1] - k[i]) / (ph[i + 1] - ph[i]) * (phiDegree - ph[i]) + k[i]
    # change units from (lb/in^3) to (kN/m^3)
    k_SIunits = khat * 271.45
    # define parabolic distribution of k with depth if desired (i.e. lin_par switch == 2)
    sigV = pyDepth * gamma
    if sigV == 0:
        sigV = 0.01
    if kSwitch == 2:
        # Equation (5-16), Boulanger et al. (2003)
        cSigma = (50.0 / sigV) ** 0.5
        # Equation (5-15), Boulanger et al. (2003)
        k_SIunits = cSigma * k_SIunits
    if kSwitch == 3:
        # ATC-32(1996) recommended subgrade modulus for given friction angle, values obtained from figure (approximate)
        ph=[29,30,31,32,33,34,35,36,37,38,39,40,41]
        # subgrade modulus above the water table
        if gwtSwitch == 1:
            # units of k are kN/m^3
            k=[720,3200,4150,5250,6500,8000,9600,11400,13000,14600,16200,18000,20100]
        # subgrade modulus below the water table
        else:
            # units of k are kN/m^3
            k=[720,2100,2700,3300,4000,4650,5450,6300,7000,7800,8600,9500,10700]
        dataNum = 13
        for i in range(0, dataNum - 1):
            if (ph[i] <= phiDegree) and (phiDegree <= ph[i + 1]):
                k_SIunits=(k[i+1] - k[i]) / (ph[i+1] - ph[i]) * (phiDegree - ph[i]) + k[i]

    # define y50 based on pult and subgrade modulus k
    # based on API (1987) recommendations, p-y curves are described using tanh functions.
    #  tcl does not have the atanh function, so must define this specifically
    #  i.e.  atanh(x) = 1/2*ln((1+x)/(1-x)), |x| < 1
    # ------------------------------------------------------------------------------------------------------
    #        p=A*pu*tanh[k*H/(A*pu)*y]  API 2002 (Eq.6.8.7-1)
    #        y=0.5*A*pu/(k*H)*ln[(1+p/pu/A)/(1-p/pu/A)]
    #        NOTE: Ln In Opensees is expressed as Log, i.e. log(e)=log(2.71828)=1
    # ------------------------------------------------------------------------------------------------------
    # when half of full resistance has been mobilized, p(y50)/pult = 0.5
    x=0.5 /A
    atanh_value = 0.5 * np.log((1.0 + x) / (1.0 - x))
    # need to be careful at ground surface (don't want to divide by zero)
    if pyDepth == 0:
        pyDepth = 0.01
    # compute y50 (need to use pult in units of force/length, and also multiply the coeff. A)
    y50 = (A * pu) / (k_SIunits * pyDepth) * atanh_value
    return pult, y50
########################################################################################################################
########################################################################################################################
def py_softClay_freeWater(gamma,z,cu,d,eps50,J,pEleLength):
    """
    --------------------------------------------------------------------------------------------------------------------
    Calculating soil parameters of soft clay with free water
    --------------------------------------------------------------------------------------------------------------------
    Matlock H. Correlation for design of laterally loaded piles in soft clay. InOffshore technology conference
	1970 Apr 21. OnePetro.
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        gamma(kN/m3)-submerged unit weight
        z(m)-the depth from the ground surface
	    Cu(kN/m2)-the undrained shear strength at depth Z
	    d(m)-the pile diameter
	    eps50-as follows
	    J-is a dimensionless factor;Matlock found that the value of J is 0.5 for soft clay and 0.25 for medium clay,
	      which was determined experimentally. The value of J is often 0.5
	    pEleLength(m)-length of the pile element
    --------------------------------------------------------------------------------------------------------------------
    Representative values of Cu and ε50 for normally consolidated clays
    [Peck RB, Hanson WE, Thornburn TH. Foundation Engineering, 2nd Ed. 1974.]
    --------------------------------------------------------------------------------------------------------------------
    consistency of clay        undrained shear strength Cu(kN/m2)(kPa)        eps50
           soft                            <48                                0.020
          medium                          48-96                               0.010
          stiff                           96-192                              0.005
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    p1 = (3 + gamma * z / float(cu) + J * z / float(d)) * cu * d
    p2 = 9 * cu * d
    pult = min(p1, p2) * pEleLength

    y50 = 2.5 * eps50 * d
    return pult, y50
########################################################################################################################
########################################################################################################################
def py_stiffClay_freeWater(ca,gamma,z,d,cu,eps50,pEleLength):
    """
    --------------------------------------------------------------------------------------------------------------------
    Calculating soil parameters of stiff clay with free water
    --------------------------------------------------------------------------------------------------------------------
    Reese LC, Cox WR, Koop FD. Field testing and analysis of laterally loaded piles om stiff clay. InOffshore
	technology conference 1975 May 4. OnePetro.
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        ca(kpa)-The average shear strength should be computed from the shear strength of the soil to a depth of 5 pile
                diameters. It should be defined as half the total maximum principal stress difference in an
                unconsolidated undrained triaxial test.
        gamma(kN/m2)-submerged unit weight
        z(m)-the depth from the ground surface
        d(m)-the pile diameter
        cu(kPa)-the undrained shear strength at depth Z
        eps50-
        pEleLength(m)-length of the pile element

    --------------------------------------------------------------------------------------------------------------------
    Representative values of Cu and ε50 for normally consolidated clays
    [Peck RB, Hanson WE, Thornburn TH. Foundation Engineering, 2nd Ed. 1974.]
    --------------------------------------------------------------------------------------------------------------------
    Representative values of Ɛ50 for overconsolidated clays
    -------------------------------------------------
                 average undrained shear strength(kPa)                        eps50
                              50-100                                           0.007
                              100-200                                          0.005
                              300-400                                          0.004
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    pct = 2.0 * ca * d + gamma * d * z + 2.83 * ca * z
    pcd = 11.0 * cu * d
    pult = min(pct, pcd) * pEleLength
    y50 = eps50 * d
    return pult, y50

########################################################################################################################
########################################################################################################################
def py_stiffClay_withoutFreeWater(ca,gamma,z,d,cu,eps50,pEleLength,J=0.5):
    """
    --------------------------------------------------------------------------------------------------------------------
    Calculating soil parameters of stiff clay without free water
    --------------------------------------------------------------------------------------------------------------------
    Welch, R.C. & L.C. Reese. (1972). Laterally loaded behavior of drilled shafts. Research Report 35-65-89.
    Center for Highway Research. University of Texas, Austin.
    --------------------------------------------------------------------------------------------------------------------
    Inputs:
        ca(kpa)-The average shear strength should be computed from the shear strength of the soil to a depth of 5 pile
                diameters. It should be defined as half the total maximum principal stress difference in an
                unconsolidated undrained triaxial test.
        gamma(kN/m2)-submerged unit weight
        z(m)-the depth from the ground surface
        d(m)-the pile diameter
        cu(kPa)-the undrained shear strength at depth Z
        eps50-
        pEleLength(m)-length of the pile element
    --------------------------------------------------------------------------------------------------------------------
    Representative values of Cu and ε50 for normally consolidated clays
    [Peck RB, Hanson WE, Thornburn TH. Foundation Engineering, 2nd Ed. 1974.]
    --------------------------------------------------------------------------------------------------------------------
    Representative values of Ɛ50 for overconsolidated clays
    -------------------------------------------------
                 average undrained shear strength(kPa)                       eps50
                              50-100                                           0.007
                              100-200                                          0.005
                              300-400                                          0.004
    --------------------------------------------------------------------------------------------------------------------
    Outputs:
        pult(kN)
        y50(m)
    --------------------------------------------------------------------------------------------------------------------
    """
    p1 = (3 + gamma * z / float(ca) + J * z / float(d)) * ca * d
    p2 = 9 * cu * d
    pult = min(p1, p2) * pEleLength
    y50 = 2.5 * eps50 * d
    return pult, y50
########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
    ####################################################################################################################
    pyDepth, gamma, phiDegree, b, pEleLength, puSwitch, ASwitch, kSwitch, gwtSwitch = 30, 9.8098, 36.918, 1.8, 0.5, 1, 1, 2, 2
    py_sand(pyDepth, gamma, phiDegree, b, pEleLength, puSwitch, ASwitch, kSwitch, gwtSwitch)
