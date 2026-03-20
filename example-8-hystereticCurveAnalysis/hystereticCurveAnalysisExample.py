#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2026-03-19
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import math
import pandas as pd
import time
import openseespy.opensees as ops
import sys
import pandas as pd
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
############################----import auxiliary modules---#############################################################
from auxiliaryModules.mainMod import HystereticCurveAnalysis
########################################################################################################################
df=pd.read_csv('cableBearingData.txt',sep='\t')
data = np.loadtxt('cableBearingData.txt')
xData,yData= data[:, 0], data[:, 1]
########################################################################################################################
hystereticCurveInstance=HystereticCurveAnalysis(xData,yData)
# hystereticCurveInstance.plotHystereticCurve(xlabel="disp.(mm)",ylabel="force(kN)",multiColors=True)
# hystereticCurveInstance.yAxisDataTranslation(startX=-25,endX=25)
# hystereticCurveInstance.plotHystereticCurve(xlabel="disp.(mm)",ylabel="force(kN)",multiColors=True)
hystereticCurveInstance.yAxisDataTranslation(startX=-10,endX=10)
hystereticCurveInstance.plotHystereticCurve()
###############################################
# hystereticCurveInstance.plotLoop_doubleDirection(loopNumber=1)
###############################################
posValue,negValue=hystereticCurveInstance.skeletonCurve_doubleDirection()
########################################################################################################################



