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
@numba.jit(nopython=True)
def SaSvSd(acc:list,dt:float,T:list,beta:float):
	"""
	acceleration,velocity and displacement response spectra calculation
	acc:acceleration time history(g)
	dt:time interval(s)
	T:periods list
	beta:damping ratio
	return:Sa(g),Sv(cm/s),Sd(cm)
	"""
	m = 1
	accConvert=[each*9.81 for each in acc]
	r,b=0.5,0.25
	SaArray=np.zeros(len(T))
	SvArray=np.zeros(len(T))
	SdArray=np.zeros(len(T))
	for i1 in range(len(T)):
		k=(2*np.pi/T[i1])**2*m
		w0=2.0*np.pi/T[i1]
		c=2.0*m*w0*beta
		accList=np.zeros((1,len(acc)))
		velList=np.zeros((1,len(acc)))
		dispList=np.zeros((1,len(acc)))
		velList[0,0]=0.0
		dispList[0,0]=0.0
		acc0=(accConvert[0]-c*velList[0,0]-k*dispList[0,0])/m
		accList[0,0]=acc0
		a1=m/(b*dt*dt)+r*c/(b*dt)
		a2=m/(b*dt)+(r/b-1.0)*c
		a3=(1.0/(2.0*b)-1.0)*m+dt*c*(r/(2.0*b)-1.0)
		k1=k+a1
		num1=len(accConvert)-1
		for j1 in range(num1):
			ptemp=-accConvert[j1+1]+a1*dispList[0,j1]+a2*velList[0,j1]+a3*accList[0,j1]
			disptemp=ptemp/k1
			dispList[0,j1+1]=disptemp
			velTemp=r*(dispList[0,j1+1]-dispList[0,j1])/(b*dt)+(1.0-r/b)*velList[0,j1]+dt*accList[0,j1]*(1.0-r/(2.0*b))
			velList[0,j1+1]=velTemp
			accTemp=(dispList[0,j1+1]-dispList[0,j1])/(b*dt*dt)-velList[0,j1]/(b*dt)-accList[0,j1]*(1.0/(2.0*b)-1.0)
			accList[0,j1+1]=accTemp
		maxAcc,maxVel,maxDisp=0.0,0.0,0.0
		for i2 in range(len(accConvert)):
			if np.fabs(velList[0,i2])>maxVel:
				maxVel=np.fabs(velList[0,i2])
			if np.fabs(dispList[0,i2])>maxDisp:
				maxDisp=np.fabs(dispList[0,i2])
			if np.fabs(accList[0,i2]+accConvert[i2])>maxAcc:
				maxAcc=np.fabs(accList[0,i2]+accConvert[i2])
		SaArray[i1]=maxAcc/9.81
		SvArray[i1]=maxVel*100
		SdArray[i1]=maxDisp*100
	return SaArray,SvArray,SdArray
if __name__ == "__main__":
	#####################################################################
	pass
