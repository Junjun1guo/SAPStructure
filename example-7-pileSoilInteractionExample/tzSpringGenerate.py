#-*-coding: UTF-8-*-
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from tzSpringModel import sandTz,clayTz

#######################################
depth=[0.5*(i+1) for i in range(100)]
sandDepth=[k1 for k1 in depth if k1<=23]
clayDepth=[k2 for k2 in depth if k2>23]
frictionAngle=30
sandUnitWeight=19
D=2
clayUnitWeight=17
shearStrength=50
eps50=0.01
shearModulus=1e5
shearVelocity=200
J=0.5
pult=[]
y50=[]
soilType=[]
materialNumber=[]
c=[]

for i1 in range(len(sandDepth)):
	pult1,y1=sandTz (frictionAngle,D,sandUnitWeight,sandDepth[i1],0.5)
	pult.append(pult1)
	y50.append(y1)
	materialNumber.append(i1+80001)
	c.append(sandUnitWeight/10.0*D*shearVelocity)
	soilType.append(2)
for i2 in range(len(clayDepth)):
	pult2,y2=clayTz (shearStrength,D,clayUnitWeight,clayDepth[i1],0.5)
	pult.append(pult2)
	y50.append(y2)
	soilType.append(1)
	materialNumber.append(i2+len(sandDepth)+80001)
	c.append(clayUnitWeight/10.0*D*shearVelocity)

matTag=np.mat(materialNumber).T
pyType=np.mat(soilType).T
pultMat=np.mat(pult).T
y50Mat=np.mat(y50).T
cMat=np.mat(c).T
pySave=np.hstack((matTag,pyType,pultMat,y50Mat,cMat))
np.savetxt("tzMaterial.txt",pySave,fmt="%d %d %.8f %.8f %.8f")






