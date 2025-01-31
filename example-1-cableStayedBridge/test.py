#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2023-04-03
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import math
import pandas as pd
from time import time
import openseespy.opensees as ops
import sys
sys.path.append("..")
############################----import auxiliary modules---#######################################
###---CalculateGroundMotionIMs is a class for calculate ground motion intensity measure,please use the command
###---print(help(CalculateGroundMotionIMs)) to check the structure and the usage of the class
from auxiliaryModules.mainMod  import CalculateGroundMotionIMs
###---GroundMotionProcess is a class for ground motion baseline correction ,fltering, and conversion among acceleration,
###---velocity and displacement, please use the command print(help(GroundMotionProcess)) to check the structure and
###---the usage of the class
from auxiliaryModules.mainMod import GroundMotionProcess
###---OpenSeesPyX is a class for the visualization and quick construction of OpenSeesPy model. please use the command
###---print(help(OpenSeesPyX)) to check the structure and the usage of the class
from auxiliaryModules.mainMod import OpenSeesPyX
###---SectionPropertyCalculate is class for calculating the section properties. please use the command
###---print(help(SectionPropertyCalculate)) to check the structure and the usage of the class
from auxiliaryModules.mainMod import SectionPropertyCalculate
###---SectionFiberDivide is a class for section fiber divide, used for quickly construct fiber-based nonlinear elements.
###---please use the command print(help(SectionFiberDivide)) to check the structure and the usage of the class
from auxiliaryModules.mainMod import SectionFiberDivide
###---SectMCAnalysis is a class for section moment curvature analysis. please use the command
###---print(help(SectMCAnalysis)) to check the structure and the usage of the class
from auxiliaryModules.mainMod import SectMCAnalysis
###---ExciteAnyDirectionOpenSees is a class for horizontally rotate FE model,it is convinient to get the rotated node
###---coordinates use this class. please use the command print(help(ExciteAnyDirectionOpenSees)) to check the structure
###---and the usage of the class
from auxiliaryModules.mainMod import ExciteAnyDirectionOpenSees
###---PythonInteractSAP2000 is a class for python  interacting with the SAP2000 program.
from auxiliaryModules.mainMod import PythonInteractSAP2000
########################################################################################################################
########################################################################################################################
########################################################################################################################
##################---初始化辅助类OpenSeesPyX用于OPenSeesPy模型可视化及后处理---############
opsX = OpenSeesPyX(dataBaseName="cableStayedBrideModel")  ###初始化OpenSeesPyX类
########################################################################################################################
##########################################---定义模型维数及自由度数---######################################################


# opsX.auxiliary_writeModelInformationToDB() ###---将模型信息写入数据库，以便在SAPStructure中显示模型
########################################################################################################################
