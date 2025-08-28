#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-08-27
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import math
import pandas as pd
import time
import openseespy.opensees as ops
import sys
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
############################----import auxiliary modules---#############################################################
###---responseSpectraCalculation is a class for acceleration,velocity and displacement response spectra calculation.
from auxiliaryModules.mainMod import responseSpectraCalculation
########################################################################################################################
########################################################################################################################
########################################################################################################################
for i in range(160):
    horizontalPath = os.path.join("inputGroundMotion/horizontal", str(i + 1) + ".txt")
    verticalPath = os.path.join("inputGroundMotion/vertical", str(i + 1) + ".txt")
    deltaT = np.loadtxt("deltaT.txt")
    dt=deltaT[i]
    horizontalAcc = np.loadtxt(horizontalPath)
    verticalAcc = np.loadtxt(verticalPath)

    logP = [0.05 * x for x in range(-40, 21, 1)]
    period = [10 ** y for y in logP]
    beta=0.05

    logSpectrumSa = []
    logSpectrumSv = []
    logSpectrumSd = []
    logT =[math.log(each,10) for each in period]
    # startTime = time.time()
    SaHo,SvHo,SdHo=responseSpectraCalculation(horizontalAcc, dt,period, beta)
    # endTime=time.time()
    # print("run time: ",endTime-startTime)
    SaVe, SvVe, SdVe = responseSpectraCalculation(verticalAcc, dt, period, beta)
    logSpectrumSa = [math.log(math.sqrt(each1*each2), 10) for each1,each2 in zip(SaHo,SaVe)]
    logSpectrumSv = [math.log(math.sqrt(each1 * each2), 10) for each1, each2 in zip(SvHo, SvVe)]
    logSpectrumSd = [math.log(math.sqrt(each1 * each2), 10) for each1, each2 in zip(SdHo, SdVe)]

    spectrumMatSa = np.asmatrix(logSpectrumSa).T
    spectrumMatSv = np.asmatrix(logSpectrumSv).T
    spectrumMatSd = np.asmatrix(logSpectrumSd).T

    logTMat = np.asmatrix(logT).T
    aa = np.hstack((logTMat, spectrumMatSa))
    bb = np.hstack((logTMat, spectrumMatSv))
    cc = np.hstack((logTMat, spectrumMatSd))

    savePath1 = os.path.join("5%geometricMeanSa/" + str(i + 1) + ".txt")
    np.savetxt(savePath1, aa, fmt="%f")

    savePath2 = os.path.join("5%geometricMeanSv/" + str(i + 1) + ".txt")
    np.savetxt(savePath2, bb, fmt="%f")

    savePath3 = os.path.join("5%geometricMeanSd/" + str(i + 1) + ".txt")
    np.savetxt(savePath3, cc, fmt="%f")
