# -*- coding: utf-8 -*-
import numpy as np
import math
from .responseSpectraCalculate import SaSvSd
import time
########################################################################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 27/08/2025
#  Environemet: Successfully excucted in python 3.13
########################################################################################################################
class IMs():
	"""
	A class for ground motion intensity measures calculation.
	"""
	def __init__(self,acc,dt):
		"""
		acc(g): acceleration time history
		dt(s): time interval
		"""
		self.acc=acc
		self.t=dt
		self.num=len(self.acc)
	####################################################################################################################
	def PGA(self):
		"""
		return-PGA(g): peak ground acceleration (PGA)
		"""
		pga=np.fabs(self.acc).max()
		return pga
	####################################################################################################################
	def AcctoVelocity (self):
		"""
		return-vel(cm/s): convert acceleration to velocity
		"""
		vel=[0]
		acc=self.acc
		# for i in range(self.num-1):
		# 	velocity=(acc[i]+acc[i+1])*self.t/2*981+vel[-1]
		# 	vel.append(velocity)
		[[velocity:=(acc[i]+acc[i+1])*self.t/2*981+vel[-1],vel.append(velocity)] for i in range(self.num-1)]
		return vel
	####################################################################################################################
	def PGV (self):
		"""
		return-PGV(cm/s): peack ground velocity (PGV)
		"""
		veloc=self.AcctoVelocity()
		pgv=np.fabs(veloc).max()
		return pgv
	####################################################################################################################
	def AcctoDisp(self):
		"""
		return-disp(cm): convert acceleration time history to displacement
		"""
		disp = [0]
		veld = self.veloc = self.AcctoVelocity()
		# for i in range(self.num - 1):
		# 	displacement = (veld[i] + veld[i + 1]) * self.t / 2 + disp[-1]
		# 	disp.append(displacement)
		[[displacement:= (veld[i] + veld[i + 1]) * self.t / 2 + disp[-1],disp.append(displacement)] for i in range(self.num - 1)]
		return disp
	####################################################################################################################
	def PGD(self):
		"""
		return-PGD(cm): Peak ground displacement
		"""
		dispc = self.AcctoDisp()
		pgd = np.fabs(dispc).max()
		return pgd
	####################################################################################################################
	def VmaxDivAmax(self):
		"""
		return-vmax/amax(s):ratio between PGV and PGA
		"""
		vmax = self.PGV()
		amax = self.PGA() * 981
		ratio = vmax / amax
		return ratio
	####################################################################################################################
	def aRMS(self):
		"""
		return-arms(g): root mean square acceleration
		"""
		a2 = np.array(self.acc ** 2)
		ttol = self.num * self.t
		atol = sum(a2) * self.t
		arms = (atol / ttol) ** 0.5
		return arms
	####################################################################################################################
	def vRMS(self):
		"""
		return-vrms(cm/s):root mean square velocity
		"""
		vel = np.array(self.AcctoVelocity())
		vel2 = vel ** 2
		ttol = self.num * self.t
		vtol = sum(vel2) * self.t
		vrms = (vtol / ttol) ** 0.5
		return vrms
	####################################################################################################################
	def dRMS(self):
		"""
		return-drms(cm): root mean square displacement
		"""
		disp = np.array(self.AcctoDisp())
		disp2 = disp ** 2
		ttol = self.num * self.t
		dtol = sum(disp2) * self.t
		drms = (dtol / ttol) ** 0.5
		return drms
	####################################################################################################################
	def AI(self):
		"""
		return-AI(m/s): Arias intensity measure
		"""
		a2 = self.acc * 9.81
		a3 = np.array(a2)
		a3tol = sum(a3 ** 2) * self.t
		ai = a3tol * math.pi / (2 * 9.81)
		return ai
	####################################################################################################################
	def Ic(self):
		"""
		return-Ic(g3/2*s1/2):Characteristic intensity
		"""
		ttol = self.num * self.t
		Ic = (self.aRMS()) ** (1.5) * (ttol) ** 0.5
		return Ic
	####################################################################################################################
	def SED(self):
		"""
		return-sed(cm2/s): Specific energy density
		"""
		vel = self.AcctoVelocity()
		velarr = np.array(vel)
		sed = sum(velarr ** 2) * self.t
		return sed
	####################################################################################################################
	def CAV(self):
		"""
		return-cav(cm/s):Cumulative absolute velocity
		"""
		accarr = np.array(self.acc)
		accabs = np.fabs(accarr) * 981
		cav = sum(accabs) * self.t
		return cav
	####################################################################################################################
	def DisptoVelocity(self, displacement, t):
		"""
		Convert displacement(cm) to velocity(cm/s)
		--------------------------------------------------
		displacement(cm): displacement time history
		t(s): time interval
		return-vel(cm/s): velocity time history
		"""
		n = len(displacement)
		vel = [0]
		disp = displacement
		# for i in range(n - 1):
		# 	vell = 2 * (disp[i + 1] - disp[i]) / t - vel[-1]
		# 	vel.append(vell)
		[[vell:= 2 * (disp[i + 1] - disp[i]) / t - vel[-1],vel.append(vell)] for i in range(n - 1)]
		return vel
	####################################################################################################################
	def VeltoAccele (self,vel,t):
		"""
		Convert velocity(cm/s) to acceleration(g)
		----------------------------------------------------
		vel(cm/s): velocity time history
		t(s): time interval
		return-acc(g): acceleration time history
		"""
		n = len(vel)
		acc = [0]
		# for i in range(n - 1):
		# 	accel = (vel[i + 1] - vel[i]) / (981 * float(dt))
		# 	acc.append(accel)
		[[accel:= (vel[i + 1] - vel[i]) / (981 * float(dt)),acc.append(accel)] for i in range(n - 1)]
		return acc
	####################################################################################################################
	def ASI(self):
		"""
		return-ASI(g*s): Acceleration spectral intensity
		"""
		T=np.arange(0.1,0.52,0.02)
		# startTime=time.time()
		Tr,Sar,Svr,Sdr= self.calResponseSpe(T, 0.05)
		# endTime=time.time()
		# print("run time=",endTime-startTime)
		asi = sum(Sar) * 0.02
		return asi
	####################################################################################################################
	def calResponseSpe (self,T,beta):
		"""
		Response spectral calculation
		------------------------------------------------
		SaSvSd(acc:list,dt:float,T:list,beta:float)
		T(s):Period list
		beta:Damping ratio
		return:Period Value(s)，absolute acceleratio spectrum（g)，relative velocity spectrum (cm/s)
		and relative displacement spectrum (cm)
		"""
		sa,sv,sd = SaSvSd(self.acc,self.t, T, beta)
		return T,sa,sv,sd
	####################################################################################################################
	def VSI(self):
		"""
		return-VSI(cm):Velocity spectral intensity
		"""
		T = np.arange(0.1, 2.52, 0.02)
		Tr, Sar, Svr, Sdr = self.calResponseSpe(T, 0.05)
		vsi = sum(Svr) * 0.02
		return vsi
	####################################################################################################################
	def __SMR (self,response):
		#this parameter gives the sustained maximum acceleration/velocity during three cycles,
		#and is defined as the third highest absolute value of acceleration/velocity in the time-history
		#(note: in order for an absolute value to be considered as a "maximum", it must be larger than values 20 steps before and 20 steps after).
		#Nuttli O.W. [1979] "The relation of sustained maximum ground acceleration and velocity to earthquake intensity and magnitude,
		#" Miscellaneous Paper S-71-1, Report 16, U.S. Army Corps of Engineers, Waterways Experiment Station, Vicksburg, Mississippi.
		a=response
		aBefore=a[0:19]
		aAfter=a[-20:-1]
		aMiddle=a[20:-21]
		b=[]
		c=[]
		# for i in range(len(aMiddle)-1):
		# 	if aMiddle[i]*aMiddle[i+1]<0:
		# 		b.append(i)
		[b.append(i) for i in range(len(aMiddle)-1) if aMiddle[i]*aMiddle[i+1]<0 ]

		# for i in np.arange(0,len(b)-1,1):
		# 	c.append(np.fabs(aMiddle[b[i]:b[i+1]]).max())
		[c.append(np.fabs(aMiddle[b[i]:b[i+1]]).max()) for i in np.arange(0,len(b)-1,1)]
		c.sort()
		return c[-3]
	####################################################################################################################
	def SMA(self):
		"""
		return-SMA(g): 持续最大加速度(g),加速度时程绝对值第三高值
		"""
		sma = self.__SMR(self.acc)
		return sma
	####################################################################################################################
	def SMV (self):
		"""
		:return: Sustained Maximum Acceleration, Third Highest Absolute Value in Acceleration Time History
		"""
		smv=self.__SMR(self.AcctoVelocity ())
		return smv
	####################################################################################################################
	def A95 (self):
		#The acceleration level below which 95% of the total Arias intensity is contained.
		#In other words, if the entire accelerogram yields a value of Ia equal to 100, the
		#A95 parameter is the threshold of acceleration such that integrating all the values
		#of the accelerogram below it, one gets an Ia=95. (g)
		ai=self.AI()
		threshold=0.95*ai
		acc=np.fabs(self.acc*9.81)
		t=self.t
		acc.sort()
		acclist=[0]
		# for i in range(len(acc)):
		# 	acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1])
		[acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) for i in range(len(acc))]
		# for i in range(len(acclist)):
		# 	if acclist[i]>=threshold:
		# 		a95=acc[i-1]/9.81
		# 		break
		a95 = next((acc[i - 1] / 9.81 for i in range(len(acclist)) if acclist[i] >= threshold), None)
		return a95
	####################################################################################################################
	def Ia (self):
		#Compound acc.-related IM  (g.s**1/3)
		#Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn30(12):1791–1816
		#Ia=PGA*(td**1/3) td=t2-t1   t1=t(5%AI)  t2=t(95%AI)
		ai=self.AI()
		thresholdt1=0.05*ai
		thresholdt2=0.95*ai
		t=self.t
		acc=self.acc*9.81
		acclist=[0]
		t1=[]
		t2=[]
		# for i in range(len(acc)):
		# 	acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1])
		[acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) for i in range(len(acc))]
		# for i1 in range(len(acclist)):
		# 	if acclist[i1]>=thresholdt1:
		# 		t1.append(i1*t)
		# 		break
		t1_value=next((i1*t for i1 in range(len(acclist)) if acclist[i1]>=thresholdt1 ), None)
		if t1_value is not None:
			t1.append(t1_value)
		# for i2 in range(len(acclist)):
		# 	if acclist[i2]>=thresholdt2:
		# 		t2.append(i2*t)
		# 		break
		t2_value = next((i2 * t for i2 in range(len(acclist)) if acclist[i2] >= thresholdt2), None)
		if t2_value is not None:
			t2.append(t2_value)
		td=t2[0]-t1[0]
		pga=self.PGA()
		Ia=pga*(td**(1.0/3.0))
		return Ia
	####################################################################################################################
	def FI (self):
		#Fajfar intensity  (cm/s**0.75)
		#Fajfar P, Vidic T, Fischinger M (1990) A measure of earthquake motion capacity to damage medium-period
		#structures. Soil Dyn Earthq Eng 9(5):236–242
		#FI=PGV*td**0.25   td=t2-t1   t1=t(5%AI)  t2=t(95%AI)
		ai=self.AI()
		thresholdt1=0.05*ai
		thresholdt2=0.95*ai
		t=self.t
		acc=self.acc*9.81
		acclist=[0]
		t1=[]
		t2=[]
		# for i in range(len(acc)):
		# 	acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1])
		[acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) for i in range(len(acc))]
		# for i1 in range(len(acclist)):
		# 	if acclist[i1]>=thresholdt1:
		# 		t1.append(i1*t)
		# 		break
		t1_value = next((i1 * t for i1 in range(len(acclist)) if acclist[i1] >= thresholdt1), None)
		if t1_value is not None:
			t1.append(t1_value)
		# for i2 in range(len(acclist)):
		# 	if acclist[i2]>=thresholdt2:
		# 		t2.append(i2*t)
		# 		break
		t2_value = next((i2 * t for i2 in range(len(acclist)) if acclist[i2] >= thresholdt2), None)
		if t2_value is not None:
			t2.append(t2_value)
		td=t2[0]-t1[0]
		pgv=self.PGV()
		fi=pgv*(td**(0.25))
		return fi
	####################################################################################################################
	def Iv (self):
		#Compound vel.-related IM (cm**2/3)/s**1/3
		#Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn
		#30(12):1791–1816
		#IvPGV**(2/3)*td(1/3)
		ai=self.AI()
		thresholdt1=0.05*ai
		thresholdt2=0.95*ai
		t=self.t
		acc=self.acc*9.81
		acclist=[0]
		t1=[]
		t2=[]
		# for i in range(len(acc)):
		# 	acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1])
		[acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) for i in range(len(acc))]
		# for i1 in range(len(acclist)):
		# 	if acclist[i1]>=thresholdt1:
		# 		t1.append(i1*t)
		# 		break
		t1_value = next((i1 * t for i1 in range(len(acclist)) if acclist[i1] >= thresholdt1), None)
		if t1_value is not None:
			t1.append(t1_value)
		# for i2 in range(len(acclist)):
		# 	if acclist[i2]>=thresholdt2:
		# 		t2.append(i2*t)
		# 		break
		t2_value = next((i2 * t for i2 in range(len(acclist)) if acclist[i2] >= thresholdt2), None)
		if t2_value is not None:
			t2.append(t2_value)
		td=t2[0]-t1[0]
		pgv=self.PGV()
		Iv=pgv**(2.0/3.0)*(td**(1.0/3.0))
		return Iv
	####################################################################################################################
	def Id (self):
		#Compound disp.-related IM (cm.s**1/3)
		##Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn
		#30(12):1791–1816
		#Id=PGD*td**1/3
		ai=self.AI()
		thresholdt1=0.05*ai
		thresholdt2=0.95*ai
		t=self.t
		acc=self.acc*9.81
		acclist=[0]
		t1=[]
		t2=[]
		# for i in range(len(acc)):
		# 	acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1])
		[acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) for i in range(len(acc))]
		# for i1 in range(len(acclist)):
		# 	if acclist[i1]>=thresholdt1:
		# 		t1.append(i1*t)
		# 		break
		t1_value = next((i1 * t for i1 in range(len(acclist)) if acclist[i1] >= thresholdt1), None)
		if t1_value is not None:
			t1.append(t1_value)
		# for i2 in range(len(acclist)):
		# 	if acclist[i2]>=thresholdt2:
		# 		t2.append(i2*t)
		# 		break
		t2_value = next((i2 * t for i2 in range(len(acclist)) if acclist[i2] >= thresholdt2), None)
		if t2_value is not None:
			t2.append(t2_value)
		td=t2[0]-t1[0]
		pgd=self.PGD()
		Id=pgd*(td**(1.0/3.0))
		return Id
########################################################################################################################
########################################################################################################################
if __name__=='__main__':
	acc = np.loadtxt("1.txt")
	imInstance = IMs(acc, 0.01)
	print("PGA:",imInstance.PGA())
	print("PGV:", imInstance.PGV())
	print("PGD:", imInstance.PGD())
	print("vmax/amax:", imInstance.VmaxDivAmax())
	print("aRMS:", imInstance.aRMS())
	print("vRMS:", imInstance.vRMS())
	print("dRMS:", imInstance.dRMS())
	print("AI:", imInstance.AI())
	print("Ic:", imInstance.Ic())
	print("SED:", imInstance.SED())
	print("CAV:", imInstance.CAV())
	print("ASI:", imInstance.ASI())
	print("VSI:", imInstance.VSI())
	print("SMA:", imInstance.SMA())
	print("SMV:", imInstance.SMV())
	print("A95:", imInstance.A95())
	print("Ia:", imInstance.Ia())
	print("FI:", imInstance.FI())
	print("Iv:", imInstance.Iv())
	print("Id:", imInstance.Id())





































