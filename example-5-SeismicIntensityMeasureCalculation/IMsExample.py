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
from auxiliaryModules.mainMod import CalculateGroundMotionIMs
########################################################################################################################
########################################################################################################################
########################################################################################################################
acc=np.loadtxt("1.txt")
dt=0.01
initiation=CalculateGroundMotionIMs(acc,dt)
imInstance=initiation.calIMs() ###---return the instance of the class IMs()
# print(help(imInstance))  ###---help document for each method
print("PGA:",imInstance.PGA())
print("PGV:", imInstance.PGV())
print("PGD:", imInstance.PGD())
print("vmax/amax:", imInstance.VmaxDivAmax())
print("aRMS:", imInstance.aRMS())
print("vRMS:", imInstance.vRMS())
print("dRMS:", imInstance.dRMS())
print("AI:", imInstance.AI())
print("Ic:", imInstance.Ic())
print("SED:", imInstance.SED())
print("CAV:", imInstance.CAV())
print("ASI:", imInstance.ASI())
print("VSI:", imInstance.VSI())
print("SMA:", imInstance.SMA())
print("SMV:", imInstance.SMV())
print("A95:", imInstance.A95())
print("Ia:", imInstance.Ia())
print("FI:", imInstance.FI())
print("Iv:", imInstance.Iv())
print("Id:", imInstance.Id())
