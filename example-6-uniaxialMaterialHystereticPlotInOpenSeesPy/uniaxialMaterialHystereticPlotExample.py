#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-08-28
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
###---Hysteretic curve plot for uniaxial materials in OpenSeesPy software.
from auxiliaryModules.mainMod import MaterialTest
from auxiliaryModules.mainMod import DisplacementHistory
########################################################################################################################
########################################################################################################################
########################################################################################################################
dispInstance=DisplacementHistory()
#######################---monotonic displacement history
# dispHistory=dispInstance.monotonicDisp(numDivide=100,maxDisp=1,scaleFactor=1,plotDispHistory=True)
#######################---cyclicOneSideConstPeak displacement history
# dispHistory=dispInstance.cyclicOneSideConstPeak(numCycle=10,numDivide=100,maxDisp=-0.3,scaleFactor=2,
#                            plotDispHistory=True)
#######################---cyclicOneSideConstPeak displacement history
# dispHistory = dispInstance.cyclicOneSideLinearIncrPeak(numCycle=10, numDivide=15, maxDisp=-0.5,
#                                                                 scaleFactor=2,plotDispHistory=True)
#######################---cyclicOneSideLinearIncrPeak_nConst displacement history
# dispHistory=dispInstance.cyclicOneSideLinearIncrPeak_nConst(numCycle= 10, numDivide= 100, maxDisp= -1,scaleFactor= 1,
#                                                 numLocalConst= 3, plotDispHistory= True)
#######################---cyclicTwoSideConstPeak displacement history
# dispHistory= dispInstance.cyclicTwoSideConstPeak(numCycle= 2, numDivide= 10, maxDisp= 1,scaleFactor= 1.5,plotDispHistory= True)
#######################---cyclicTwoSideConstPeak displacement history
dispHistory= dispInstance.cyclicTwoSideLinearIncrPeak(numCycle=10,numDivide=20,maxDisp=-0.1,scaleFactor=1,plotDispHistory=True)
######################---cyclicTwoSideConstPeak displacement history
# dispHistory= dispInstance.cyclicTwoSideLinearIncrPeak_nConst(numCycle=5,numDivide=10,maxDisp=1,scaleFactor=1,
#                                        numLocalConst=2,plotDispHistory=True)
########################################################################################################################
materialTestInstance=MaterialTest()
##################################
# dispList,forceList=materialTestInstance.uniaixalMaterialTest("uniaxialMaterial('PySimple1',3000,2,815,0.0142771,0.3)",
#                                                              dispHistory)
# dispList,forceList=materialTestInstance.uniaixalMaterialTest("uniaxialMaterial('QzSimple1',6000,2,1800,0.0142)",
#                                                              dispHistory)
dispList,forceList=materialTestInstance.uniaixalMaterialTest("uniaxialMaterial('TzSimple1',1,2,1.8*80,0.0014588)",
                                                             dispHistory)
