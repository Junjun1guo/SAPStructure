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
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
############################----import auxiliary modules---#############################################################
from auxiliaryModules.mainMod import SectionPropertyCalculate
########################################################################################################################
###-use simple section to check the scaleFactor is right
sectionPropertyInstance=SectionPropertyCalculate()
A, Iyy, Izz, J, Cy, Cz,outSideNode,outSideEle,inSideNode,inSideEle= sectionPropertyInstance.dxf_sectionproperties(
dxfFileName="sections.dxf",layerName="pier",scaleFactor=1000,meshSize=0.01,numCircleSeg=50,numArcSeg=10,numEllipseSeg=20,numSplineSeg=20)
print("A=", A, " Iyy=", Iyy, " Izz=", Izz, " J=", J, " Cy=", Cy, " Cz=",Cz)
########################################################################################################################



