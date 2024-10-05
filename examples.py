# -*-coding: UTF-8-*-
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully excucted in python 3.11
#  Date: 2021.04.27
# import necessary modules
from auxiliaryModules.mainMod import SectionPropertyCalculate
sectionPropertyInstance=SectionPropertyCalculate()
A, Iyy, Izz, J, Cy, Cz,outSideNode,outSideEle,inSideNode,inSideEle= sectionPropertyInstance.dxf_sectionproperties\
    ("00.dxf","0",scaleFactor=1000,meshSize=0.0005)
print("A=", A, " Iyy=", Iyy, " Izz=", Izz, " J=", J, " Cy=", Cy, " Cz=",Cz)

