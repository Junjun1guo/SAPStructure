#-*-coding: UTF-8-*-
#########################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 27/08/2025
#  Environemet: Successfully excucted in python 3.11
#########################################################################
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import numba
###########################################
import numpy as np
import numba
from numba import njit, prange

@njit(fastmath=True, parallel=True)
def SaSvSd(acc: np.ndarray, dt: float, T: np.ndarray, beta: float):
    m = 1.0
    g0 = 9.81
    r, b = 0.5, 0.25
    nT = len(T)
    n = len(acc)
    accConvert = acc * g0
    SaArray = np.zeros(nT)
    SvArray = np.zeros(nT)
    SdArray = np.zeros(nT)

    ### use prange and parallel can take advantage of multi CPU
    for i1 in prange(nT):
        ti = T[i1]
        k = (2.0 * np.pi / ti) ** 2 * m
        w0 = 2.0 * np.pi / ti
        c = 2.0 * m * w0 * beta

        velList = np.zeros(n)
        dispList = np.zeros(n)
        accList = np.zeros(n)

        accList[0] = (accConvert[0] - c * velList[0] - k * dispList[0]) / m

        a1 = m / (b * dt * dt) + r * c / (b * dt)
        a2 = m / (b * dt) + (r / b - 1.0) * c
        a3 = (1.0 / (2.0 * b) - 1.0) * m + dt * c * (r / (2.0 * b) - 1.0)
        k1 = k + a1

        for j1 in range(n - 1):
            ptemp = -accConvert[j1 + 1] + a1 * dispList[j1] + a2 * velList[j1] + a3 * accList[j1]
            disptemp = ptemp / k1
            dispList[j1 + 1] = disptemp

            velTemp = (r * (dispList[j1 + 1] - dispList[j1]) / (b * dt)
                       + (1.0 - r / b) * velList[j1]
                       + dt * accList[j1] * (1.0 - r / (2.0 * b)))
            velList[j1 + 1] = velTemp

            accTemp = ((dispList[j1 + 1] - dispList[j1]) / (b * dt * dt)
                       - velList[j1] / (b * dt)
                       - accList[j1] * (1.0 / (2.0 * b) - 1.0))
            accList[j1 + 1] = accTemp

        maxAcc = 0.0
        maxVel = 0.0
        maxDisp = 0.0

        for i2 in range(n):
            absVel = abs(velList[i2])
            absDisp = abs(dispList[i2])
            absAcc = abs(accList[i2] + accConvert[i2])
            if absVel > maxVel:
                maxVel = absVel
            if absDisp > maxDisp:
                maxDisp = absDisp
            if absAcc > maxAcc:
                maxAcc = absAcc

        SaArray[i1] = maxAcc / g0
        SvArray[i1] = maxVel * 100.0
        SdArray[i1] = maxDisp * 100.0

    return SaArray, SvArray, SdArray
if __name__ == "__main__":
	#####################################################################
	pass
