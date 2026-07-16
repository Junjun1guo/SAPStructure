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
numSimMotions=100 ###-total number of simulated ground motions
Vp=None ###-velocity pulse magnitude (cm/s), None means no specific value for Vp
Tp=4.5 ###-velocity pulse period (s), None means no specific value for Tp
########################################################################################################################
########################################################################################################################
# folderName=f"Vp={Vp}_Tp={Tp}" ###-folder name
# if os.path.exists(folderName):
#     shutil.rmtree(folderName)
# os.makedirs(folderName)
# for i1 in range(numSimMotions):
# 	print(f"Current simulation number {i1+1}!")
# 	while True:
# 		###-moment magnitude of an Earthquake, [5.5,7.9]
# 		moment = np.random.uniform(5.5, 7.9)
# 		###-unit,km,the closest distance from the site to the fault rupture in kilometers, [0.07,31]
# 		Rrup = np.random.uniform(5,5.1)
# 		###-unit,km,directivity parameters s or d (input the large of s and d), [1.2,135]
# 		s_or_d = np.random.uniform(30, 30)
# 		###-unit,degreens,directivity angle theta or phi (input corresponding to s_or_d) [0,90]
# 		theta_or_phi = np.random.uniform(45, 45)
# 		###-unit:m, site soil shear wave average velocity over the top 30 meters in meters per second,[139,2016]
# 		Vs30 = np.random.uniform(1000,1001)
# 		###-type of faulting,including ["strikeSlip","reverseAndReverseOblique"]
# 		faultType = "strikeSlip"
# 		###-unit:km, depth to the top of the rupture plane in kilometers, larger or equal to 0
# 		zTOR = 0
# 		numSimMotions_perLoop = 1
# 		nearFaultMotionSimInstance =NearFaultGroundMotionSim_twoComp_specificPulseParas(faultType,moment, zTOR, Rrup,
# 									Vs30, s_or_d,theta_or_phi,numSimMotions_perLoop)
# 		pulseComponnetsList, pulseParas = nearFaultMotionSimInstance.runSpecificPulseMotion(Vp=Vp,Tp=Tp)
# 		if pulseComponnetsList is not None:
# 			np.savetxt(folderName+"/"+str(i1+1)+"_pulseComp.txt",np.asmatrix(pulseComponnetsList[0][0]).T,fmt="%.6f")
# 			np.savetxt(folderName+"/"+str(i1+1)+"_orthogComp.txt",np.asmatrix(pulseComponnetsList[0][1]).T,fmt="%.6f")
# 			np.savetxt(folderName+"/"+str(i1+1)+"_pulseModel.txt",np.asmatrix(pulseComponnetsList[0][2]).T,fmt="%.6f")
# 			np.savetxt(folderName + "/" + str(i1 + 1) + "_pulseParas.txt", np.asmatrix(pulseParas).T,fmt="%.6f")
# 			break
########################################################################################################################
########################################################################################################################
###---plot time history
# folderName=os.path.abspath(f"Vp={Vp}_Tp={Tp}")
# motionNumber=6
# NearFaultGroundMotionSim_twoComp_specificPulseParas.timeHistoryPlot(folderName,motionNumber)
########################################################################################################################
########################################################################################################################
###---plot response spectra
# folderName=os.path.abspath(f"Vp={Vp}_Tp={Tp}")
# NearFaultGroundMotionSim_twoComp_specificPulseParas.responseSpectraPlot(folderName=folderName,maxPeriod=6,beta=0.05,
#                                     logCoord=True,SaLim=[0.01,4],SvLim=[0.01,500],SdLim=[0.001,500])

