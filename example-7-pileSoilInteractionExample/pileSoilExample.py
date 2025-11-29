#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Haofan Peng, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2025-11-29
########################################################################################################################
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import math
import pandas as pd
from time import time
import openseespy.opensees as ops
import sys
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
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
from auxiliaryModules.py_tz_qz_SoilParameters.pyModel import (py_sand,py_softClay_freeWater,
                            py_stiffClay_freeWater,py_stiffClay_withoutFreeWater)
from auxiliaryModules.py_tz_qz_SoilParameters.tzModel import tz_sand,tz_clay
from auxiliaryModules.py_tz_qz_SoilParameters.qzModel import qz_sand,qz_clay
########################################################################################################################
########################################################################################################################
##################---初始化辅助类OpenSeesPyX用于OPenSeesPy模型可视化及后处理---###############################################
opsX = OpenSeesPyX(dataBaseName="PileSoilExample")  ###初始化OpenSeesPyX类
########################################################################################################################
#########################################Scenario Parameters##########################################################
# N = 2  ## Num. of pile rows
# d = 1.6  ## Diameter of piles (m)
# rho_pile = 0.015  ## Longitudinal reinforcement ratio of piles
# alpha = 0.11  ## Axial ratio of piles
# p2pDis = 2.7  ## Pile-to-Pile distance (d)
# Dr = 0.52  ## Soil Relative Density
# ScourDepth = 2.0  ## Scour depth (m)
# Hp = 10.5  ## Pier height (m)
# DT = 4.8  ## Pier size (m) along the loading direction
# DL = 3.2  ## Pier size (m) perpendicular to the loading direction
# rho_pier = 0.03 ## Longitudinal reinforcement ratio of pier
# rhoTran_pile = 0.012  ## Transverse reinforcement ratio of pile
# fy_pile = 471	  ## steel yield strength (MPa)
# fc_pile = 42  ## concrete cover peak strength (MPa)
# LP = 19.5  ## Pile total length (m)
########################################################################################################################
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
###################################### define pier steel/concrete material##############################################
ops.uniaxialMaterial('Concrete04',9913,-42000.0,-0.002,-0.005,30459481.28251694)  ###cover
ops.uniaxialMaterial('Concrete04',9914,-50130.62184414807,-0.003935862343844779,-0.012420003272895231,33277401.28882108)  ###core
ops.uniaxialMaterial('MinMax',9903,9913,'-min',-0.005)  ###cover
ops.uniaxialMaterial('MinMax',9904,9914,'-min',-0.012420003272895231)  ###core
####steel
ops.uniaxialMaterial('Steel02',9905,471000,200000000,0.01,20,0.925,0.15,0,1,0,1)
###################################### define pile steel/concrete material##############################################
ops.uniaxialMaterial('Concrete04',9933,-42000.0,-0.002,-0.005,30459481.28251694)  ###cover
ops.uniaxialMaterial('Concrete04',9934,-53953.94126414787,-0.004846176491463779,-0.013388007402836017,34523072.9009604)  ###core
ops.uniaxialMaterial('MinMax',9923,9933,'-min',-0.005)  ###cover
ops.uniaxialMaterial('MinMax',9924,9934,'-min',-0.013388007402836017)  ###core
###################################### define pile ground effect ##############################################
p2pDis = 2.7 ## Pile-to-Pile distance (d)
P_leading = 0.100 * (p2pDis - 3.0) + 0.8
P_middle = 0.225 * (p2pDis - 3.0) + 0.4
P_trailing = 0.200 * (p2pDis - 3.0) + 0.3
###################################### define soil material ##############################################
gamma = 9.535830426995947  # effective density
phi = 32.8148 # friction angle DDS
Gsoil =  61653.18781754988 # shear modulus DDS
###################################### define pier fiber ##############################################
d =1.6    ## Diameter of piles (m)
DT = 4.8  ## Pier size (m) along the loading direction
DL = 3.2  ## Pier size (m) perpendicular to the loading direction
cover = 0.1  ## Concrete cover thickness
thickness = d/3.0
Cor11Y = DT / 2.0
Cor11Z = DL / 2.0
Cor12Y = -Cor11Y
Cor12Z = Cor11Z
Cor13Y = -Cor11Y
Cor13Z = -Cor11Z
Cor14Y = Cor11Y
Cor14Z = -Cor11Z

Cor21Y = Cor11Y - cover
Cor21Z = Cor11Z - cover
Cor22Y = -Cor21Y
Cor22Z = Cor21Z
Cor23Y = -Cor21Y
Cor23Z = -Cor21Z
Cor24Y = Cor21Y
Cor24Z = -Cor21Z

Cor31Y = Cor21Y - thickness
Cor31Z = Cor21Z - thickness
Cor32Y = -Cor31Y
Cor32Z = Cor31Z
Cor33Y = -Cor31Y
Cor33Z = -Cor31Z
Cor34Y = Cor31Y
Cor34Z = -Cor31Z

Cor41Y = Cor31Y - cover
Cor41Z = Cor31Z - cover
Cor42Y = -Cor41Y
Cor42Z = Cor41Z
Cor43Y = -Cor41Y
Cor43Z = -Cor41Z
Cor44Y = Cor41Y
Cor44Z = -Cor41Z
rho_pier = 0.03   ## Longitudinal reinforcement ratio of pier
areaBar = 4 * (Cor21Y * Cor21Z - Cor31Y * Cor31Z) * rho_pier / 182.0

opsX.section('Fiber', 1, '-GJ', 1.0e10)
# Outer cover
opsX.patch('rect', 9903, int(cover / 0.1), int(2 * Cor11Z / 0.1), Cor24Y, Cor14Z, Cor11Y, Cor11Z)  # right
opsX.patch('rect', 9903, int(cover / 0.1), int(2 * Cor11Z / 0.1), Cor13Y, Cor13Z, Cor22Y, Cor12Z)  # left
opsX.patch('rect', 9903, int(2 * Cor21Y / 0.1), int(cover / 0.1), Cor22Y, Cor22Z, Cor21Y, Cor11Z)  # top
opsX.patch('rect', 9903, int(2 * Cor21Y / 0.1), int(cover / 0.1), Cor23Y, Cor13Z, Cor24Y, Cor24Z)  # bottom
# Core concrete
opsX.patch('rect', 9904, int(thickness / 0.1), int(2 * Cor21Z / 0.1), Cor34Y, Cor24Z, Cor21Y, Cor21Z)  # right
opsX.patch('rect', 9904, int(thickness / 0.1), int(2 * Cor21Z / 0.1), Cor23Y, Cor23Z, Cor32Y, Cor22Z)  # left
opsX.patch('rect', 9904, int(2 * Cor31Y / 0.1), int(thickness / 0.1), Cor32Y, Cor32Z, Cor31Y, Cor21Z)  # top
opsX.patch('rect', 9904, int(2 * Cor31Y / 0.1), int(thickness / 0.1), Cor33Y, Cor23Z, Cor34Y, Cor34Z)  # bottom
# Inner cover
opsX.patch('rect', 9903, int(cover / 0.1), int(2 * Cor31Z / 0.1), Cor44Y, Cor34Z, Cor31Y, Cor31Z)  # right
opsX.patch('rect', 9903, int(cover / 0.1), int(2 * Cor31Z / 0.1), Cor33Y, Cor33Z, Cor42Y, Cor32Z)  # left
opsX.patch('rect', 9903, int(2 * Cor41Y / 0.1), int(cover / 0.1), Cor42Y, Cor42Z, Cor41Y, Cor31Z)  # top
opsX.patch('rect', 9903, int(2 * Cor41Y / 0.1), int(cover / 0.1), Cor43Y, Cor33Z, Cor44Y, Cor44Z)  # bottom
# Outer rebars
opsX.layer('straight', 9905, 36, areaBar, Cor22Y, Cor22Z, Cor21Y, Cor21Z)  # top
opsX.layer('straight', 9905, 36, areaBar, Cor23Y, Cor23Z, Cor24Y, Cor24Z)  # bottom
opsX.layer('straight', 9905, 17, areaBar, Cor22Y, Cor22Z - 0.15, Cor23Y, Cor23Z + 0.15)  # left
opsX.layer('straight', 9905, 17, areaBar, Cor21Y, Cor21Z - 0.15, Cor24Y, Cor24Z + 0.15)  # right
# Inner rebars
opsX.layer('straight', 9905, 27, areaBar, Cor32Y, Cor32Z, Cor31Y, Cor31Z)  # top
opsX.layer('straight', 9905, 27, areaBar, Cor33Y, Cor33Z, Cor34Y, Cor34Z)  # bottom
opsX.layer('straight', 9905, 11, areaBar, Cor32Y, Cor32Z - 0.15, Cor33Y, Cor33Z + 0.15)  # left
opsX.layer('straight', 9905, 11, areaBar, Cor31Y, Cor31Z - 0.15, Cor34Y, Cor34Z + 0.15)  # right

# pltHandle = opsX.auxiliary_fiberSectionPlot(sectTag=1)
# pltHandle.savefig("sect2.eps")
# pltHandle.show()
###################################### define pile fiber ##############################################
opsX.section('Fiber',2,'-GJ',1.0e10)
opsX.patch('circ',9923,10,2,0,0,0.7,0.8,0,360)
opsX.patch('circ',9924,10,8,0,0, 0 , 0.7 ,0,360)
opsX.layer('circ',9905,30,0.0007696902001294995,0,0,0.7)

# pltHandle = opsX.auxiliary_fiberSectionPlot(sectTag=2)
# pltHandle.savefig("sect2.eps")
# pltHandle.show()
##########################################---定义桩节点---######################################################
pileNodes = np.loadtxt('ModelInfo/pileNodes.txt')
for each in pileNodes:
    NodeTag = int(each[0])
    xCoord = float(each[1])
    yCoord = float(each[2])
    zCoord = float(each[3])
    nodeMass = float(each[4])
    opsX.node(NodeTag,xCoord,yCoord,zCoord,'-mass',nodeMass,nodeMass,nodeMass,0.0,0.0,0.0)
    ops.fix(NodeTag,0,0,1,1,1,0)

##########################################---定义承台节点---######################################################
capNodes = np.loadtxt('ModelInfo/capNodes.txt')
for each in capNodes:
    NodeTag = int(each[0])
    xCoord = float(each[1])
    yCoord = float(each[2])
    zCoord = float(each[3])
    nodeMass = float(each[4])
    opsX.node(NodeTag,xCoord,yCoord,zCoord,'-mass',nodeMass,nodeMass,nodeMass,0.0,0.0,0.0)
    ops.fix(NodeTag,0,0,1,1,1,0)

##########################################---定义桥墩节点---######################################################
pierNodes = np.loadtxt('ModelInfo/pierNodes.txt')
for each in pierNodes:
    NodeTag = int(each[0])
    xCoord = float(each[1])
    yCoord = float(each[2])
    zCoord = float(each[3])
    nodeMass = float(each[4])
    opsX.node(NodeTag,xCoord,yCoord,zCoord,'-mass',nodeMass,nodeMass,nodeMass,0.0,0.0,0.0)
    ops.fix(NodeTag,0,0,1,1,1,0)

##########################################---定义支座---######################################################
ops.uniaxialMaterial('Elastic', 9901, 1E10)
opsX.element('zeroLength', 9999, 122, 123, '-mat', 9901, 9901, 9901, 9901, 9901, 9901,
             '-dir', 1, 2, 3, 4, 5, 6)

##########################################---定义局部坐标转换---######################################################
opsX.geomTransf('Corotational', 1, 0, 0, -1)  # Considering the P-△ effect
opsX.geomTransf('Linear', 2, 0, 1, 0)  # For pile-cap connections

##########################################---定义桩单元---######################################################
integrationPoint=5
pileElement = np.loadtxt('ModelInfo/pileElement.txt')
for each in pileElement:
    EleTag = int(each[0])
    NodeI = int(each[1])
    NodeJ = int(each[2])
    opsX.element('nonlinearBeamColumn',EleTag,NodeI,NodeJ,integrationPoint,2,1)

##########################################---定义承台单元---######################################################
capElement = np.loadtxt('ModelInfo/capElement.txt')
for each in capElement:
    EleTag = int(each[0])
    NodeI = int(each[1])
    NodeJ = int(each[2])
    A = float(each[3])
    E = float(each[4])
    G = float(each[5])
    J = float(each[6])
    Iy = float(each[7])
    Iz = float(each[8])
    Transf = int(each[9])
    opsX.element('elasticBeamColumn', EleTag, NodeI, NodeJ, A, E, G, J, Iy, Iz, Transf)

##########################################---定义桥墩单元---######################################################
pierElement = np.loadtxt('ModelInfo/pierElement.txt')
integrationPoint=5
for each in pierElement:
    EleTag = int(each[0])
    EleNodeI = int(each[1])
    EleNodeJ = int(each[2])
    EleSection = int(each[3])
    EleGeomTransf = int(each[4])
    opsX.element('nonlinearBeamColumn', EleTag, EleNodeI, EleNodeJ, integrationPoint,EleSection, EleGeomTransf)

##########################################---定义土节点---######################################################
soilNodes = np.loadtxt('ModelInfo/soilNodes.txt')
for each in soilNodes:
    NodeTag = int(each[0])
    xCoord = float(each[1])
    yCoord = float(each[2])
    zCoord = float(each[3])
    nodeMass = float(each[4])
    opsX.node(NodeTag,xCoord,yCoord,zCoord,'-mass',nodeMass,nodeMass,nodeMass,0.0,0.0,0.0)
    ops.fix(NodeTag,1,1,1,1,1,1)

##########################################---定义py,qz弹簧---######################################################
pySpring = np.loadtxt('ModelInfo/py,qzUniaxialMat.txt')
for each in pySpring:
    pyDepth = float(each[0])
    pyTag1 = int(each[1])
    pyTag3 = int(each[2])
    qzTag = int(each[3])
    pult0,y50 = py_sand(pyDepth,gamma,phi,1.6,0.5,1,2,1,2)
    pult_leadingPile=pult0*P_trailing
    pult_trailingPile=pult0 *P_trailing
    ops.uniaxialMaterial('PySimple1',pyTag1,2,pult_leadingPile,y50,0.3)
    ops.uniaxialMaterial('PySimple1',pyTag3,2,pult_trailingPile,y50,0.3)
    pult_active=pult0*(P_leading-P_trailing)
    ops.uniaxialMaterial('QzSimple1',qzTag,2,pult_active,y50)

##########################################---定义土py,qz零长度单元---######################################################
soilZeroElement = np.loadtxt('ModelInfo/soilZeroElement1.txt')
for each in soilZeroElement:
    EleTagPy = int(each[0])
    EleTagQz = int(each[1])
    EleNodeI = int(each[2])
    EleNodeJ = int(each[3])
    EleMaterialTagPy = int(each[4])
    EleMaterialTagQz = int(each[5])
    opsX.element('zeroLength',EleTagPy,EleNodeI,EleNodeJ,'-mat',EleMaterialTagPy,'-dir',1)
    opsX.element('zeroLength',EleTagQz,EleNodeI,EleNodeJ,'-mat',EleMaterialTagQz,'-dir',1)

soilZeroElement = np.loadtxt('ModelInfo/soilZeroElement2.txt')
for each in soilZeroElement:
    EleTagPy = int(each[0])
    EleTagQz = int(each[1])
    EleNodeI = int(each[2])
    EleNodeJ = int(each[3])
    EleMaterialTagPy = int(each[4])
    EleMaterialTagQz = int(each[5])
    opsX.element('zeroLength',EleTagPy,EleNodeI,EleNodeJ,'-mat',EleMaterialTagPy,'-dir',1,'-orient',-1,0,0,0,1,0)
    opsX.element('zeroLength',EleTagQz,EleNodeI,EleNodeJ,'-mat',EleMaterialTagQz,'-dir',1,'-orient',-1,0,0,0,1,0)

##########################################---定义tz弹簧---######################################################
tzSpring = np.loadtxt('ModelInfo/tzUniaxialMat.txt')
for each in tzSpring:
    pyDepth = float(each[0])
    tzTag = int(each[1])
    tult, z50 = tz_sand(phi, 1.6, gamma, pyDepth, 0.5)
    ops.uniaxialMaterial('TzSimple1', tzTag, 2, 1.8 * tult, z50)
ops.uniaxialMaterial('ENT',10001,10000000)

##########################################---定义tz零长度单元---######################################################
tzZeroElement = np.loadtxt('ModelInfo/tzZeroElement.txt')
for each in tzZeroElement:
    EleTag = int(each[0])
    EleNodeI = int(each[1])
    EleNodeJ = int(each[2])
    EleMaterialTagtz = int(each[3])
    opsX.element('zeroLength',EleTag,EleNodeI,EleNodeJ,'-mat',EleMaterialTagtz,'-dir',2)

# ##########################################---施加重力荷载---############################################################
nodesTags=ops.getNodeTags()
nodesMass=[ops.nodeMass(each)[0] for each in nodesTags]
nodeList=[[node,mass] for node,mass in zip(nodesTags, nodesMass)]
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
for each in nodeList:
    ops.load(int(each[0]), 0.0, -each[1] * 9.81, 0.0, 0.0, 0.0, 0.0)

# ##########################################---施加重力荷载---############################################################
ops.test('NormDispIncr', 1.0e-4, 50)
ops.algorithm('NewtonLineSearch', 0.8)
ops.integrator('Newmark', 0.5, 0.25)
ops.system('BandGeneral')
ops.constraints('Transformation')
ops.numberer('RCM')
ops.analysis('Transient')
ok = ops.analyze(1, 1)
if ok == 0:
    print("--- Structure gravity analysis succeed!")
else:
    print("--- Structure gravity analysis disconvergence!!!")
ops.setTime(0.0)
ops.wipeAnalysis()
ops.remove("recorders")
############################################################################################################
opsX.auxiliary_writeModelInformationToDB()  ###---将模型信息写入数据库，以便在SAPStructure中显示模型
###########################################---重力分析---#################################################################
# opsX.integration_analysisGravity(totalStep=1,recordList=None) ###---采用集成式方法进行自重分析(包装了自重分析相关命令
# ################################################################################################################
# opsX.integration_analysisModal(numModes=50)  ###---采用集成式方法进行自重分析(包装了模态分析相关命令)
########################################
### Beginning Modal Analysis..."
### determine Natural Period, Frequency & damping parameters for SDOF
### set xDamp 0.05; # damping ratio (0.02-0.05-typical)
m_order = 10  # model order
lamiga = ops.eigen(m_order)
periodList = []
for i in range(m_order):
    omegaa = (lamiga[i]) ** 0.5
    Ti = 2 * 3.1415926 / omegaa
    periodList.append(Ti)
    print(f"{i + 1}th period is: {Ti}")

# ################################################################################################################
periodPath = os.path.join("Period", f"{1}")
os.makedirs(periodPath, exist_ok=True)
opsX.auxiliary_writeDataIntoTxtFile(savePath=periodPath, filename='period', listData=[periodList], decimals=[6])
################################################################################################################
xDamp = 0.05
MpropSwitch = 1.0
KcurrSwitch = 1.0
KcommSwitch = 0.0
KinitSwitch = 0.0
nEigenI = 1
nEigenJ = 2
lambdaI = lamiga[nEigenI - 1]
lambdaJ = lamiga[nEigenJ - 1]
omegaI = lambdaI ** 0.5
omegaJ = lambdaJ ** 0.5
alphaM = MpropSwitch * xDamp * (2 * omegaI * omegaJ) / (omegaI + omegaJ)  # M-prop. damping; D = alphaM*M
betaKcurr = KcurrSwitch * 2. * xDamp / (omegaI + omegaJ)  # current-K;      +beatKcurr*KCurrent
betaKcomm = KcommSwitch * 2. * xDamp / (omegaI + omegaJ)  # last-committed K;   +betaKcomm*KlastCommitt
betaKinit = KinitSwitch * 2. * xDamp / (omegaI + omegaJ)
print("alphaM,betaK=", alphaM, betaKcurr)
### "*** Beginning Time History Analysis..."
motionDT =0.005
tsTag=10
Scale = 1
ops.timeSeries('Path',tsTag,'-dt',motionDT,'-filePath',f'inputGroundMotion/waves/31.txt','-factor',9.81*Scale)
ops.pattern('UniformExcitation',400,1,'-accel',tsTag)
# Load Recorder
recDT=motionDT
#### pile and cap response
ops.recorder('Node','-file',f'Data/31-1-10000.out','-timeSeries',tsTag,'-time',
             '-dT',recDT,'-node',611,'-dof',1,'accel') ###cap acc. (-timeSeries: absolute = gm + relative)
ops.recorder('Node','-file',f'Data/31-1-11000.out','-time','-dT',recDT,
             '-node',611,'-dof',1,2,6,'disp') #cap disp. and rotation
ops.recorder('Node','-file',f'Data/31-1-11010.out','-time','-dT',recDT,
             '-node',610,'-dof',1,2,6,'disp') #cap disp. and rotation
for row in range(2):
    pilenum_ctr=(2 *row + 1)*1000
    pilenum_ctr_Ele=(2 *row + 1)*1100
    ops.recorder('Node','-file',f'Data/31-1-{12000+row*20000}.out','-time',
                 '-dT',recDT,'-nodeRange',pilenum_ctr+1,pilenum_ctr+40,'-dof',1,'disp')
    ops.recorder('Node','-file',f'Data/31-1-{15000 + row * 20000}.out',
                 '-time','-dT',recDT,'-nodeRange',pilenum_ctr+1,pilenum_ctr+40,'-dof',2,'disp')
    ops.recorder('Element','-file',f'Data/31-1-{13000 + row * 20000}.out',
                 '-time','-dT',recDT,'-eleRange',pilenum_ctr_Ele+1,pilenum_ctr_Ele+39,'localForce')
    ops.recorder('Element','-file',f'Data/31-1-{14000 +row * 20000}.out',
                 '-time','-dT',recDT,'-eleRange',pilenum_ctr_Ele+1,pilenum_ctr_Ele+39,'section',5,'deformation')
### pier and superstructure response
ops.recorder('Node','-file',f'Data/31-1-15000.out','-timeSeries',tsTag,
             '-time','-dT',recDT,'-node',123,'-dof',1,'accel') ### superstruc. acc.
ops.recorder('Node','-file',f'Data/31-1-16000.out','-time','-dT',recDT,
             '-node',123,'-dof',1,2,6,'disp')
ops.recorder('Node','-file',f'Data/31-1-17000.out','-time','-dT',recDT,
             '-nodeRange',101,122,'-dof',1,'disp') ### pier node lateral disp.
ops.recorder('Node','-file',f'Data/31-1-17500.out','-time','-dT',recDT,
             '-nodeRange',101,122,'-dof',2,'disp')
ops.recorder('Element','-file',f'Data/31-1-18000.out','-time','-dT',
             recDT,'-eleRange',101,121,'localForce')
ops.recorder('Element','-file',f'Data/31-1-19000.out','-time','-dT',recDT,
             '-eleRange',101,121,'section',1,'deformation')
### reaction force and displacment in bearing
ops.recorder('Element','-file',f'Data/31-1-19999.out','-time','-dT',recDT,
             '-ele',55,66,'localForce')
#######################################################################################################################
ops.constraints('Penalty',1.e20,1.e20)
ops.test('NormDispIncr',1.0e-4,35,0)
ops.algorithm('KrylovNewton')
ops.numberer('RCM')
ops.system('ProfileSPD')
ops.integrator('Newmark',0.5,0.25)
ops.rayleigh(alphaM,betaKcurr,betaKinit,betaKcomm)
ops.analysis('Transient')
motionSteps = 11007
tFinal=motionSteps *motionDT
tFinal=10
tCurrent=ops.getTime()
ok=0
step=1
while (ok==0)and(tCurrent<tFinal):
    ok=ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print("Trying Broyden...")
        ops.algorithm('Broyden',8)
        ok=ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print("Trying Newton ..")
        ops.algorithm('Newton')
        ok=ops.analyze(1,motionDT)
    if ok!=0:
        ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print("Trying Newton With LineSearch ..")
        ops.algorithm('NewtonLineSearch')
        ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print("Trying Broyden-Fletcher-Goldfarb-Shanno (BFGS) ..")
        ops.algorithm('BFGS')
        ok=ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print('Trying SecantNewton ..')
        ops.algorithm('SecantNewton')
        ok=ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok!=0:
        print("Trying ModifiedNewton ..")
        ops.algorithm('ModifiedNewton','-initial')
        ok=ops.analyze(1,motionDT)
    if ok!=0:
        ok=ops.analyze(1,motionDT/2.0)
    if ok!=0:
        ok=ops.analyze(1,motionDT/4.0)
    if ok==0:
        print(f"Succeed! Step {step}/{motionSteps} of Acc 31 with Scalefactor {Scale} in Case 1 ")
        ops.algorithm('KrylovNewton')
    if ok!=0:
        print(f"Disconvergence! Step {step}/{motionSteps} of Acc 31 with Scalefactor {Scale} in Case 1")
    step=step + 1
    tCurrent=ops.getTime()

