#-*-coding: UTF-8-*-
#########################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-08-27
#########################################################################
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib
matplotlib.use('TkAgg')
###########################################
yvalue=[]
for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean1=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean1real=[10**x for x in mean1]
var1=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std1=[x**0.5  for x in var1]
var1real=[10**y for y in std1]

upper1=[]
lower1=[]
for i1 in range(len(mean1)):
	inorm=norm.ppf(0.975, loc=mean1[i1], scale=std1[i1])
	upper1.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean1[i1], scale=std1[i1])
	lower1.append(lnorm)

upperReal1=[10**x for x in upper1]
lowerReal1=[10**x for x in lower1]
###########################################
yvalue=[]
for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean2=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean2real=[10**x for x in mean2]
var2=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std2=[x**0.5  for x in var2]
var2real=[10**y for y in std2]

upper2=[]
lower2=[]
for i1 in range(len(mean2)):
	inorm=norm.ppf(0.975, loc=mean2[i1], scale=std2[i1])
	upper2.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean2[i1], scale=std2[i1])
	lower2.append(lnorm)

upperReal2=[10**x for x in upper2]
lowerReal2=[10**x for x in lower2]
###########################################
yvalue=[]
for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean3=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean3real=[10**x for x in mean3]
var3=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std3=[x**0.5  for x in var3]
var3real=[10**y for y in std3]

upper3=[]
lower3=[]
for i1 in range(len(mean3)):
	inorm=norm.ppf(0.975, loc=mean3[i1], scale=std3[i1])
	upper3.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean3[i1], scale=std3[i1])
	lower3.append(lnorm)

upperReal3=[10**x for x in upper3]
lowerReal3=[10**x for x in lower3]
###########################################
yvalue=[]
for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean4=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean4real=[10**x for x in mean4]
var4=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std4=[x**0.5  for x in var4]
var4real=[10**y for y in std4]

upper4=[]
lower4=[]
for i1 in range(len(mean4)):
	inorm=norm.ppf(0.975, loc=mean4[i1], scale=std4[i1])
	upper4.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean4[i1], scale=std4[i1])
	lower4.append(lnorm)

upperReal4=[10**x for x in upper4]
lowerReal4=[10**x for x in lower4]
###########################################
yvalue=[]
for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean5=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean5real=[10**x for x in mean5]
var5=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std5=[x**0.5  for x in var5]
var5real=[10**y for y in std5]

upper5=[]
lower5=[]
for i1 in range(len(mean5)):
	inorm=norm.ppf(0.975, loc=mean5[i1], scale=std5[i1])
	upper5.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean5[i1], scale=std5[i1])
	lower5.append(lnorm)

upperReal5=[10**x for x in upper5]
lowerReal5=[10**x for x in lower5]
###########################################
yvalue=[]
for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean6=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean6real=[10**x for x in mean6]
var6=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std6=[x**0.5  for x in var6]
var6real=[10**y for y in std6]

upper6=[]
lower6=[]
for i1 in range(len(mean6)):
	inorm=norm.ppf(0.975, loc=mean6[i1], scale=std6[i1])
	upper6.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean6[i1], scale=std6[i1])
	lower6.append(lnorm)

upperReal6=[10**x for x in upper6]
lowerReal6=[10**x for x in lower6]
###########################################
yvalue=[]
for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean7=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean7real=[10**x for x in mean7]
var7=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std7=[x**0.5  for x in var7]
var7real=[10**y for y in std7]

upper7=[]
lower7=[]
for i1 in range(len(mean7)):
	inorm=norm.ppf(0.975, loc=mean7[i1], scale=std7[i1])
	upper7.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean7[i1], scale=std7[i1])
	lower7.append(lnorm)

upperReal7=[10**x for x in upper7]
lowerReal7=[10**x for x in lower7]
###########################################
yvalue=[]
for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean8=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean8real=[10**x for x in mean8]
var8=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std8=[x**0.5  for x in var8]
var8real=[10**y for y in std8]

upper8=[]
lower8=[]
for i1 in range(len(mean8)):
	inorm=norm.ppf(0.975, loc=mean8[i1], scale=std8[i1])
	upper8.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean8[i1], scale=std8[i1])
	lower8.append(lnorm)

upperReal8=[10**x for x in upper8]
lowerReal8=[10**x for x in lower8]
###########################################
yvalue=[]
for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean9=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean9real=[10**x for x in mean9]
var9=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std9=[x**0.5  for x in var9]
var9real=[10**y for y in std9]

upper9=[]
lower9=[]
for i1 in range(len(mean9)):
	inorm=norm.ppf(0.975, loc=mean9[i1], scale=std9[i1])
	upper9.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean9[i1], scale=std9[i1])
	lower9.append(lnorm)

upperReal9=[10**x for x in upper9]
lowerReal9=[10**x for x in lower9]
###########################################
yvalue=[]
for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean10=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean10real=[10**x for x in mean10]
var10=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std10=[x**0.5  for x in var10]
var10real=[10**y for y in std10]

upper10=[]
lower10=[]
for i1 in range(len(mean10)):
	inorm=norm.ppf(0.975, loc=mean10[i1], scale=std10[i1])
	upper10.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean10[i1], scale=std10[i1])
	lower10.append(lnorm)

upperReal10=[10**x for x in upper10]
lowerReal10=[10**x for x in lower10]
###########################################
yvalue=[]
for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean11=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean11real=[10**x for x in mean11]
var11=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std11=[x**0.5  for x in var11]
var11real=[10**y for y in std11]

upper11=[]
lower11=[]
for i1 in range(len(mean11)):
	inorm=norm.ppf(0.975, loc=mean11[i1], scale=std11[i1])
	upper11.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean11[i1], scale=std11[i1])
	lower11.append(lnorm)

upperReal11=[10**x for x in upper11]
lowerReal11=[10**x for x in lower11]
###########################################
yvalue=[]
for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xlog=nploadtxt[:,0]
	ylog=nploadtxt[:,1]
	yvalue.append(ylog)
mean12=np.mean(np.asmatrix(yvalue),axis=0).tolist()[0]
mean12real=[10**x for x in mean12]
var12=np.var(np.asmatrix(yvalue),axis=0).tolist()[0]
std12=[x**0.5  for x in var12]
var12real=[10**y for y in std12]

upper12=[]
lower12=[]
for i1 in range(len(mean12)):
	inorm=norm.ppf(0.975, loc=mean12[i1], scale=std12[i1])
	upper12.append(inorm)
	lnorm=norm.ppf(0.025, loc=mean12[i1], scale=std12[i1])
	lower12.append(lnorm)

upperReal12=[10**x for x in upper12]
lowerReal12=[10**x for x in lower12]
########################################################################################################################
fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(ax10,ax11,ax12))=plt.subplots(4,3)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.grid(which="both")
ax1.set_ylim(10**(-4),10**1)

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.grid(which="both")
ax2.set_ylim(10**(-4),10**4)

ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.grid(which="both")
ax3.set_ylim(10**(-5),10**4)

ax4.set_xscale("log")
ax4.set_yscale("log")
ax4.grid(which="both")
ax4.set_ylim(10**(-4),10**1)

ax5.set_xscale("log")
ax5.set_yscale("log")
ax5.grid(which="both")
ax5.set_ylim(10**(-4),10**4)


ax6.set_xscale("log")
ax6.set_yscale("log")
ax6.grid(which="both")
ax6.set_ylim(10**(-5),10**4)

ax7.set_xscale("log")
ax7.set_yscale("log")
ax7.grid(which="both")
ax7.set_ylim(10**(-4),10**1)


ax8.set_xscale("log")
ax8.set_yscale("log")
ax8.grid(which="both")
ax8.set_ylim(10**(-4),10**4)

ax9.set_xscale("log")
ax9.set_yscale("log")
ax9.grid(which="both")
ax9.set_ylim(10**(-5),10**4)

ax10.set_xscale("log")
ax10.set_yscale("log")
ax10.grid(which="both")
ax10.set_ylim(10**(-4),10**1)

ax11.set_xscale("log")
ax11.set_yscale("log")
ax11.grid(which="both")
ax11.set_ylim(10**(-4),10**4)

ax12.set_xscale("log")
ax12.set_yscale("log")
ax12.grid(which="both")
ax12.set_ylim(10**(-5),10**4)

for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax1.plot(xreal,yreal,"k")
ax1.plot(xreal,mean1real,"r",linewidth=2)
ax1.plot(xreal,upperReal1,"b",linewidth=2)
ax1.plot(xreal,lowerReal1,"b",linewidth=2)


for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax2.plot(xreal,yreal,"k")
ax2.plot(xreal,mean2real,"r",linewidth=2)
ax2.plot(xreal,upperReal2,"b",linewidth=2)
ax2.plot(xreal,lowerReal2,"b",linewidth=2)

for i1 in range(40):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax3.plot(xreal,yreal,"k")
ax3.plot(xreal,mean3real,"r",linewidth=2)
ax3.plot(xreal,upperReal3,"b",linewidth=2)
ax3.plot(xreal,lowerReal3,"b",linewidth=2)

for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax4.plot(xreal,yreal,"k")
ax4.plot(xreal,mean4real,"r",linewidth=2)
ax4.plot(xreal,upperReal4,"b",linewidth=2)
ax4.plot(xreal,lowerReal4,"b",linewidth=2)

for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax5.plot(xreal,yreal,"k")
ax5.plot(xreal,mean5real,"r",linewidth=2)
ax5.plot(xreal,upperReal5,"b",linewidth=2)
ax5.plot(xreal,lowerReal5,"b",linewidth=2)

for i1 in range(40,80,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax6.plot(xreal,yreal,"k")
ax6.plot(xreal,mean6real,"r",linewidth=2)
ax6.plot(xreal,upperReal6,"b",linewidth=2)
ax6.plot(xreal,lowerReal6,"b",linewidth=2)

for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax7.plot(xreal,yreal,"k")
ax7.plot(xreal,mean7real,"r",linewidth=2)
ax7.plot(xreal,upperReal7,"b",linewidth=2)
ax7.plot(xreal,lowerReal7,"b",linewidth=2)

for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax8.plot(xreal,yreal,"k")
ax8.plot(xreal,mean8real,"r",linewidth=2)
ax8.plot(xreal,upperReal8,"b",linewidth=2)
ax8.plot(xreal,lowerReal8,"b",linewidth=2)


for i1 in range(80,120,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax9.plot(xreal,yreal,"k")
ax9.plot(xreal,mean9real,"r",linewidth=2)
ax9.plot(xreal,upperReal9,"b",linewidth=2)
ax9.plot(xreal,lowerReal9,"b",linewidth=2)

for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSa/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax10.plot(xreal,yreal,"k")
ax10.plot(xreal,mean10real,"r",linewidth=2)
ax10.plot(xreal,upperReal10,"b",linewidth=2)
ax10.plot(xreal,lowerReal10,"b",linewidth=2)

for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSv/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax11.plot(xreal,yreal,"k")
ax11.plot(xreal,mean11real,"r",linewidth=2)
ax11.plot(xreal,upperReal11,"b",linewidth=2)
ax11.plot(xreal,lowerReal11,"b",linewidth=2)


for i1 in range(120,160,1):
	loadPath=os.path.join("5%geometricMeanSd/"+str(i1+1)+".txt")
	nploadtxt=np.loadtxt(loadPath)
	xload=nploadtxt[:,0]
	yload=nploadtxt[:,1]
	xreal=[10**x for x in xload]
	yreal=[10**y for y in yload]
	ax12.plot(xreal,yreal,"k")
ax12.plot(xreal,mean12real,"r",linewidth=2)
ax12.plot(xreal,upperReal12,"b",linewidth=2)
ax12.plot(xreal,lowerReal12,"b",linewidth=2)

plt.savefig("groundMotionSpectrum.eps")
plt.show()

	
	
	
		
		
		
