#-*-coding: UTF-8-*-
#########################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 25/02/2019
#  Environemet: Successfully excucted in python 3.11
#########################################################################
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import norm
from responseSpectMain import SaSvSd
###########################################
f=[np.exp(i0*0.035-3) for i0 in range(200)]##frequencies(Hz)0.05~50Hz
Tarray=[1.0/each1 for each1 in f]
Tarray.reverse()
beta=0.05
##################################################
dt1=np.loadtxt("nonPulseMotions/dt.txt")
length1=np.loadtxt("nonPulseMotions/length.txt")
saNon=[]
svNon=[]
for i1 in range(10):
	acc1=np.loadtxt("nonPulseMotions/FN/"+str(i1+1)+".txt")
	sa1, sv1, sd1 = SaSvSd(acc1, dt1[i1], Tarray, beta)
	saNon.append(sa1)
	svNon.append(sv1)

dt1=np.loadtxt("pulseMotions/dt.txt")
length1=np.loadtxt("pulseMotions/length.txt")
saPulse=[]
svPulse=[]
for i1 in range(10):
	acc1=np.loadtxt("pulseMotions/FN/"+str(i1+1)+".txt")
	sa1, sv1, sd1 = SaSvSd(acc1, dt1[i1], Tarray, beta)
	saPulse.append(sa1)
	svPulse.append(sv1)
#####################################################
fig,((ax1,ax2))=plt.subplots(1,2)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.grid(which="both")
ax1.set_xlim(0.01896,20)
ax1.set_ylim(10**(-4),10**1)
for i1 in range(10):
	ax1.plot(Tarray,saPulse[i1],"k")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.grid(which="both")
ax2.set_xlim(0.01896,20)
ax2.set_ylim(10**(-2),10**3)
for i2 in range(10):
	ax2.plot(Tarray,svPulse[i2],"k")
plt.savefig("spectComapre.eps")
plt.show()