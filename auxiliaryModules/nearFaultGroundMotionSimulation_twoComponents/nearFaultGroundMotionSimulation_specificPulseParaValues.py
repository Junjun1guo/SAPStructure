# -*-coding: UTF-8-*-
########################################################################################################################
#  Author: Junjun Guo,Tongji University. https://github.com/Junjun1guo
#  E-mail: guojj@tongji.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.8
#  Date: 2022-10-02
########################################################################################################################
########################---导入必要模块--#################################
########################################################################
# import necessary modules
import os
from tqdm.std import trange
import cmath
from multiprocessing import Pool ##多进程
from sko.DE import DE
import time
import random
from numba import jit,prange
from scipy.stats import beta,norm,uniform
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from pathlib import Path
matplotlib.use('TkAgg')
#####---对于类中方法的加速可以在类中保留该方法的接口，然后将实现部分写成numba加速的函数进行调用，从而保持接口的完整性与性能的一致性
#########################################################################
########################################################################################################################
########################################################################################################################
def accToVelocity(dt, acc):
	"""
	加速度时程转换为速度时程
	from acceleration (g) to velocity (cm/s)
	dt:time interval (s)
	acc: acceleration time history (g/s2)
	vel: velocity time history (cm/s)
	output:vel-velocity time history (cm/s)
	"""
	vel = [0]
	num = len(acc)
	for i in range(num - 1):
		velocity = (acc[i] + acc[i + 1]) * dt / float(2) * 981 + vel[-1]
		vel.append(velocity)
	return vel
########################################################################################################################
########################################################################################################################
def velToDisplacement (dt,vel):
	"""
	from velocity (cm/s) to displacement (cm)
	input:dt-time interval(s)
	vel-velocity time history(cm/s)
	output:disp-displacement time history(cm)
	"""
	disp=[0]
	num=len(vel)
	for i in range(num-1):
		displacement=(vel[i]+vel[i+1])*dt/float(2)+disp[-1]
		disp.append(displacement)
	return disp
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
def motionFourierSpectra(acc,dt):
	"""
	---地震动傅里叶谱幅值及相位计算---
	注意：这里傅里叶幅值谱经FFT变换后除以采样频率得到（即乘以地震动时间间隔dt)
	输入:
		acc(float list)-加速度时程列表(g)
		dt(float)-时间间隔
	输出：
		freq(Hz),xabs，phase(弧度）-傅里叶幅值谱频率，幅值及相位角
	"""
	xw=np.fft.fft(acc)
	xabs=np.abs(xw)
	magnitudeSpectra=xabs*dt
	freq=np.fft.fftfreq(len(acc),dt).tolist()
	phase = [cmath.phase(each) for each in xw]
	return freq,magnitudeSpectra,phase
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
@jit(nopython=True)
def ResponseCal(beta,dt,Tarray,accArray):
	""""""
	accm = accArray * 9.81
	accLength = len(accm)
	r = 0.5
	b = 0.25
	saArray=np.zeros(len(Tarray))
	svArray = np.zeros(len(Tarray))
	sdArray = np.zeros(len(Tarray))
	for i1 in prange(len(Tarray)):
		k=(2.0*np.pi/Tarray[i1])**2*1
		w0=2.0*np.pi/Tarray[i1]
		c=2.0*1*w0*beta
		acc=np.zeros(accLength)
		vel=np.zeros(accLength)
		disp=np.zeros(accLength)
		acc0=(accm[0]-c*vel[0]-k*disp[0])/1.0
		acc[0]=acc0
		a1=1.0/(b*dt**2)+r*c/(b*dt)
		a2=1.0/(b*dt)+(r/b-1.0)*c
		a3=(1.0/(2.0*b)-1)*1.0+dt*c*(r/(2.0*b)-1.0)
		k1=k+a1
		num1=accLength-1
		for j1 in prange(num1):
			ptemp=-accm[j1+1]+a1*disp[j1]+a2*vel[j1]+a3*acc[j1]
			disptemp=ptemp/float(k1)
			disp[j1+1]=disptemp
			velTemp=r*(disp[j1+1]-disp[j1])/(b*dt)+(1.0-r/b)*vel[j1]+dt*acc[j1]*(1.0-r/(2.0*b))
			vel[j1+1]=velTemp
			accTemp=(disp[j1+1]-disp[j1])/(b*dt**2)-vel[j1]/(b*dt)-acc[j1]*(1.0/(2.0*b)-1.0)
			acc[j1+1]=accTemp
		accAbs=accm+acc
		maxAcc=np.max(np.abs(accAbs))
		maxVel=np.max(np.abs(vel))
		maxDisp=np.max(np.abs(disp))
		saArray[i1]=maxAcc/9.81
		svArray[i1]=maxVel*100.0
		sdArray[i1]=maxDisp*100.0
	return saArray,svArray,sdArray
########################################################################################################################
########################################################################################################################
class NearFaultMotionSim():
	"""
	---基于场地及震源特性的随机近场地震动模拟类---
	参考文献：
	1. Dabaghi and Der Kiureghian (2014 PEER report) "Stochastic Modeling and Simulation of Near-Fault Ground Motions for
	   Performance-Based Earthquake Engineering
	2. Dabaghi and Der Kiureghian (2017) "Stochastic model for simulation of NF GMs"
	3. Dabaghi and Der Kiureghian (2018) "Simulation of orthogonal horizontal components of near-fault ground motion
	   for specified EQ source and site characteristics"
	4. Dabaghi and Der Kiureghian (2022) “Simulation of near-fault ground motions for randomized directivity parameters”
	"""
	def __init__(self,faultType,simulationType,momentMag,zTOR,Rrup,Vs30,s_or_d,theta_or_phi,numSimMotions):
		"""
		faultType(str): type of faulting,including ["strikeSlip","reverseAndReverseOblique"]
		simulationType(str): type of simulatated grond motions, including ["pulseAndNoPulse","onlyPulse","onlyNoPulse"]
		momentMag(float): moment magnitude of an Earthquake
		zTOR(float,km): depth to the top of the rupture plane in kilometers
		Rrup(float,km): the closest distance from the site to the fault rupture in kilometers
		Vs30(float,m/s): site soil shear wave average velocity over the top 30 meters in meters per second
		s_or_d(float,km): directivity parameters s or d (input the large of s and d)
		theta_or_phi(float,degrees): directivity angle theta or phi (input corresponding to s_or_d)
		numSimMotions(int): the number of simulated ground motion time histories
		"""
		self.faultType=faultType
		if not (self.faultType=="strikeSlip" or self.faultType=="reverseAndReverseOblique"):
			raise ValueError("fault type should be 'strikeSlip' or 'reverseAndReverseOblique'")
		self.simulationType=simulationType
		self.momentMag=momentMag
		if not (5.5<=self.momentMag<=7.9):
			raise ValueError("moment magnitude should within (5.5,7.9)")
		self.zTOR=zTOR
		if self.zTOR<0:
			raise ValueError("zToR should larger or equal to 0")
		self.Rrup=Rrup
		if not (0.07<=self.Rrup<=31):
			raise ValueError("Rrup should within [0.07,31]")
		self.Vs30=Vs30
		if not (139<=self.Vs30<=2016):
			raise ValueError("Vs30 should within [139,2016]")
		self.s_or_d=s_or_d
		if not (1.2<=self.s_or_d<=135):
			raise ValueError("sord should within [1.2,135]")
		self.theta_or_phi=theta_or_phi
		if not (0.0<=self.theta_or_phi<=90):
			raise ValueError("theta_or_phi should within [0,90]")
		self.numSimMotions=numSimMotions
		self.timeStep=0.005
		self.numPulseMotion=0
		self.start_time=0.0
		if self.simulationType=="onlyPulse":
			self.numPulseMotion=self.numSimMotions
		elif self.simulationType=="onlyNoPulse":
			self.numPulseMotion=0
		elif self.simulationType=="pulseAndNoPulse":
			self.numPulseMotion=self.numPulseMotionCal(self.numSimMotions)
		else:
			raise ValueError("simulationType should be 'onlyPulse','onlyNoPulse' or 'pulseAndNoPulse'")
		self.numNoPulseMotion=int(self.numSimMotions-self.numPulseMotion)
		###Pulse parameters---Vp,Tp,gamma,niu,Do_maxp,Ia_res,D595_res,D05_res,D030_res,fmid_res,f'res,zeta_res,Ia_po,D595_po,
		###---D05_po,D030_po,fmid_po,f'po,zeta_po
		###---
		self.sigma_pulse=np.array([0.385316782551070,0.580604607504062,1.000000000000000,1.000000000000000,
		                      0.468600626648626,0.781462017439926,0.371951335342951,0.442124880737131,
		                      0.393548709857826,0.409988222892834,0.820134132413354,1.096436125406160,
		                      0.746552904533869,0.402440900041447,0.461424491954810,0.407724539607907,
		                      0.440166826670740,0.824632603897275,0.961997443697123])
		###---noPulse parameters:
		###---Ia_NP1,D595_NP1,D05_NP1,D030_NP1,fmid_NP1,f'NP1,zeta_NP1,Ia_NP2,D595_NP2,D05_NP2,D030_NP2,fmid_NP2,f'NP2,
		###---zeta_NP2
		self.sigma_noPulse=np.array([1.052723262620090,0.398427668080412,0.456828618083131,0.305727090125879,
		                             0.447517114210032,0.941288665677244,1.007680943597980,1.028226030318770,
		                             0.375877126809162,0.458413470215522,0.294118965466636,0.399966941388161,
		                             0.831550984874095,0.887870513796394])
		self.corrMatrix_pulse=np.array([[1, -0.175836500638571, -0.0203457508302324, 0.173589795302921, 0.191081284058678,
		0.446601858355920, 0.0422292189436895, 0.0268120615665615, 0.120003074985399, -0.382241273128019, 0.0573276270760638,
		0.155115849407343, 0.407222009805941, -0.0405527073664880, 0.0164573906036360, 0.0415103476380082, -0.282855958356553,
		0.106173613413053, 0.0486136094283240],[-0.175836500638571, 1, 0.183359484001117, 0.00131202858838653, 0.431552926960257,
		-0.0802776119560438, 0.0983744630339078, 0.310456055209018, 0.368157334142620, 0.0509849248659661, 0.0223164167958617,
		0.175239604674009, -0.0979434412969958, 0.143078236501468, 0.295214504980143, 0.370873286184142, 0.0409704140512328,
		-0.0681917932721656, 0.243299650091412],[-0.0203457508302324, 0.183359484001117, 1, -0.190019984820145, 0.242770958876810,
		0.178170099097870, 0.107202069371700, 0.150512140384424, 0.236433714526518, -0.114026379004051, 0.0535511100783777,
		0.0642193542887409, 0.0688913090123055, 0.127858468767928, 0.0893601270266779, 0.213748741406994, -0.0697188946975121,
		0.0212082250410815, 0.123098089235859],[0.173589795302921, 0.00131202858838653, -0.190019984820145, 1, 0.119159474433286,
		-0.0812496120969922, 0.0911626190533515, 0.0653720081018262, 0.0718151835177181, -0.133641807530098, -0.0920397447082436,
		0.0299407456292105, -0.0242282623085006, -0.0172855366304964, 0.0663348465393307, 0.0729217102589892, -0.146245197915904,
		0.0999563418090508, -0.0409433584122388],[0.191081284058678, 0.431552926960257, 0.242770958876810, 0.119159474433286,
		1, 0.0597439686298471, 0.163724290728144, 0.733127256394521, 0.788656583220928, -0.0312901152364670, -0.153794052040040,
		0.126072923696987, 0.0157717534215698, 0.186868291148981, 0.680765575675777, 0.749289248752773, 0.00752109582582447,
		-0.165776964800244, 0.192516732945902],[0.446601858355920, -0.0802776119560438, 0.178170099097870, -0.0812496120969922,
		0.0597439686298471, 1, -0.0244049816649646, 0.0584250580150467, 0.0856617566986387, 0.0783654782068241, 0.0990943110343599,
		0.0264806656713516, 0.837794567493788, 0.0235744551632278, 0.00722635162195151, 0.0801743456052000, 0.125480219427063,
		0.0813027074983168, 0.0454360407125637],[0.0422292189436895, 0.0983744630339078, 0.107202069371700, 0.0911626190533515,
		0.163724290728144, -0.0244049816649646, 1, 0.0612360822187756, 0.245411252666782, -0.0247015935697729, -0.224809121780326,
		0.0459354003374228, 0.0364365286784528, 0.760546075583923, 0.0370009049319755, 0.205658174145097, -0.0837872744722257,
		-0.0338684418879410, -0.0274880324550587],[0.0268120615665615, 0.310456055209018, 0.150512140384424, 0.0653720081018262,
		0.733127256394521, 0.0584250580150467, 0.0612360822187756, 1, 0.855219194729329, 0.0120435632525462, -0.0489647334211007,
		0.199955997824984, 0.0212627853872894, 0.146801029873847, 0.931624963947970, 0.850430406207709, 0.0232039506387295,
		-0.00394570150177181, 0.244390432559234],[0.120003074985399, 0.368157334142620, 0.236433714526518, 0.0718151835177181,
		0.788656583220928, 0.0856617566986387, 0.245411252666782, 0.855219194729329, 1, 0.0298843738756568, -0.0618147333807716,
		0.197760496777865, 0.104861038989822, 0.262606914312613, 0.795537444898935, 0.906602423826413, 0.0381266607245457,
		-0.0143594270096637, 0.226001111347866],[-0.382241273128019, 0.0509849248659661, -0.114026379004051, -0.133641807530098,
		-0.0312901152364670, 0.0783654782068241, -0.0247015935697729, 0.0120435632525462, 0.0298843738756568, 1, -0.241519621897699,
		0.112473855187119, 0.171776282660699, -0.0492035952139502, 0.0493755848461730, 0.112703656481410, 0.864715714588433,
		-0.286206615545551, 0.157882561174870],[0.0573276270760638, 0.0223164167958617, 0.0535511100783777, -0.0920397447082436,
		-0.153794052040040, 0.0990943110343599, -0.224809121780326, -0.0489647334211007, -0.0618147333807716, -0.241519621897699,
		1, 0.112365366315021, -0.0106706754632641, -0.0488287220881385, -0.0635241398373312, -0.0911374530847290, -0.0885687207002145,
		0.421522224993818, 0.239492035016932],[0.155115849407343, 0.175239604674009, 0.0642193542887409, 0.0299407456292105,
		0.126072923696987, 0.0264806656713516, 0.0459354003374228, 0.199955997824984, 0.197760496777865, 0.112473855187119,
		0.112365366315021, 1, -0.0274158393051499, 0.0871234642667006, 0.160919524797887, 0.256368904115107, 0.275537783357955,
		-0.175867301345906, 0.792491388890200],[0.407222009805941, -0.0979434412969958, 0.0688913090123055, -0.0242282623085006,
		0.0157717534215698, 0.837794567493788, 0.0364365286784528, 0.0212627853872894, 0.104861038989822, 0.171776282660699,
		-0.0106706754632641, -0.0274158393051499, 1, -0.168436208860883, 0.0139228196755138, 0.0606940022435021, 0.0950730856215810,
		0.166344473008140, -0.0377746475967110],[-0.0405527073664880, 0.143078236501468, 0.127858468767928, -0.0172855366304964,
		0.186868291148981, 0.0235744551632278, 0.760546075583923, 0.146801029873847, 0.262606914312613, -0.0492035952139502,
		-0.0488287220881385, 0.0871234642667006, -0.168436208860883, 1, 0.0584697041090115, 0.243068478704018, 0.0101378284587601,
		-0.134310788077535, 0.0858141211009980],[0.0164573906036360, 0.295214504980143, 0.0893601270266779, 0.0663348465393307,
		0.680765575675777, 0.00722635162195151, 0.0370009049319755, 0.931624963947970, 0.795537444898935, 0.0493755848461730,
		-0.0635241398373312, 0.160919524797887, 0.0139228196755138, 0.0584697041090115, 1, 0.841226866362572, 0.0319804533198723,
		0.0395868643324250, 0.183336741524880],[0.0415103476380082, 0.370873286184142, 0.213748741406994, 0.0729217102589892,
		0.749289248752773, 0.0801743456052000, 0.205658174145097, 0.850430406207709, 0.906602423826413, 0.112703656481410,
		-0.0911374530847290, 0.256368904115107, 0.0606940022435021, 0.243068478704018, 0.841226866362572, 1, 0.0999775444207582,
		-0.0848763381048173, 0.277349661980356],[-0.282855958356553, 0.0409704140512328, -0.0697188946975121, -0.146245197915904,
		0.00752109582582447, 0.125480219427063, -0.0837872744722257, 0.0232039506387295, 0.0381266607245457, 0.864715714588433,
		-0.0885687207002145, 0.275537783357955, 0.0950730856215810, 0.0101378284587601, 0.0319804533198723, 0.0999775444207582,
		1, -0.432951535761925, 0.262376572921422],[0.106173613413053, -0.0681917932721656, 0.0212082250410815, 0.0999563418090508,
		-0.165776964800244, 0.0813027074983168, -0.0338684418879410, -0.00394570150177181, -0.0143594270096637, -0.286206615545551,
		0.421522224993818, -0.175867301345906, 0.166344473008140, -0.134310788077535, 0.0395868643324250, -0.0848763381048173,
		-0.432951535761925, 1, -0.184861197592255],[0.0486136094283240, 0.243299650091412, 0.123098089235859, -0.0409433584122388,
		0.192516732945902, 0.0454360407125637, -0.0274880324550587, 0.244390432559234, 0.226001111347866, 0.157882561174870, 0.239492035016932,
		0.792491388890200, -0.0377746475967110, 0.0858141211009980, 0.183336741524880, 0.277349661980356, 0.262376572921422, -0.184861197592255, 1]])
		
		self.corrMatrix_noPulse=np.array([[1, -0.183620641202513, 0.0890171218487119, 0.104132896092390, 0.0143281984142704,
		0.202871723469377, -0.151909317725644, 0.945163870283100, -0.0778432911362303, 0.0495683691288216, 0.0966843496273208,
		0.0917721771113965, 0.103741261286614, -0.121195511065596],[-0.183620641202513, 1, 0.0854423655718587, 0.307373593686928,
		-0.0150490152853094, -0.149397471797000, 0.0888348816281089, -0.0794047082384602, 0.848433024094089, 0.0950830201555768,
		0.288741920713302, -0.0596971902001111, -0.0160895956989375, 0.113905908782625],[0.0890171218487119, 0.0854423655718587,
		1, 0.813213357087557, -0.225438606252670, 0.000735703328121357, -0.0840944291284059, 0.0562899783045387, 0.137192754102300,
		0.907605550316775, 0.788263163970441, -0.189167012249831, -0.0219422025252862, -0.0931560965857281],[0.104132896092390,
		0.307373593686928, 0.813213357087557, 1, -0.162877782549727, -0.0946212742252604, -0.0192432432422716, 0.131014027317413,
		0.289495591671556, 0.752836252137069, 0.908489822222397, -0.151488045845577, -0.0637984858287899, -0.0498615936932623],
		[0.0143281984142704, -0.0150490152853094, -0.225438606252670, -0.162877782549727, 1, -0.187527904866745, -0.163661457269968,
		0.0720174301613615, -0.0805658506819631, -0.173027594836660, -0.165135405413428, 0.897082085156820, -0.0778400029130002,
		 -0.00420375269007833],[0.202871723469377, -0.149397471797000, 0.000735703328121357, -0.0946212742252604, -0.187527904866745,
		1, -0.0853760806838177, 0.151981779909482, -0.0282514238500333, 0.00217240817323823, -0.0879254146125414, -0.0874961774514421,
		0.647468312996265, -0.157070775414507],[-0.151909317725644, 0.0888348816281089, -0.0840944291284059, -0.0192432432422716,
		-0.163661457269968, -0.0853760806838177, 1, -0.109821275189057, 0.0605607584025393, -0.0674419385042876, -0.0178734643092523,
		-0.0879585887923949, -0.105931812299965, 0.761324011431321],[0.945163870283100, -0.0794047082384602, 0.0562899783045387,
		0.131014027317413, 0.0720174301613615, 0.151981779909482, -0.109821275189057, 1, -0.0744735914486929, 0.0477307773362732,
		0.116976058986177, 0.104731056969540, 0.137507782979518, -0.107158677162180],[-0.0778432911362303, 0.848433024094089,
		0.137192754102300, 0.289495591671556, -0.0805658506819631, -0.0282514238500333, 0.0605607584025393, -0.0744735914486929,
		1, 0.0782297693189997, 0.294723991917604, -0.0895932216480470, -0.0501878780366773, 0.0967842822751738],[0.0495683691288216,
		0.0950830201555768, 0.907605550316775, 0.752836252137069, -0.173027594836660, 0.00217240817323823, -0.0674419385042876,
		0.0477307773362732, 0.0782297693189997, 1, 0.786122088469745, -0.177521125226899, 0.00868592321024064, -0.0671981184080449],
		[0.0966843496273208, 0.288741920713302, 0.788263163970441, 0.908489822222397, -0.165135405413428, -0.0879254146125414,
		-0.0178734643092523, 0.116976058986177, 0.294723991917604, 0.786122088469745, 1, -0.168356869639510, -0.0773913106180435,
		-0.0274804568206274],[0.0917721771113965, -0.0596971902001111, -0.189167012249831, -0.151488045845577, 0.897082085156820,
		-0.0874961774514421, -0.0879585887923949, 0.104731056969540, -0.0895932216480470, -0.177521125226899, -0.168356869639510,
		1, -0.183773515755380, 0.00693211686662153],[0.103741261286614, -0.0160895956989375, -0.0219422025252862, -0.0637984858287899,
		-0.0778400029130002, 0.647468312996265, -0.105931812299965, 0.137507782979518, -0.0501878780366773, 0.00868592321024064,
		-0.0773913106180435, -0.183773515755380, 1, -0.110874872058067],[-0.121195511065596, 0.113905908782625, -0.0931560965857281,
		-0.0498615936932623, -0.00420375269007833, -0.157070775414507, 0.761324011431321, -0.107158677162180, 0.0967842822751738,
		-0.0671981184080449, -0.0274804568206274, 0.00693211686662153, -0.110874872058067, 1]])
		
		###---const,M,(M-6.5)I(M>6.5),lnsqrt(Rrup2+h2),Mlnsqrt(Rrup2+h2),F.ffilt_z,lnVs30,s_or_d
		self.regressionCoeff_pulse=np.array([[1.69862554416145, 0.608190177030033, -0.608190177030033, -0.576217471720854,
		0, 0.183071013159864, -0.0939319189357984, 0.00657091132855174],[-2.47924338395758, 0.670395625327187, 0, 0, 0,
		-0.263957799575519, -0.232548659903406, 0.00791988432530859],[0, 0, 0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0, 0, 0],
		[-4.24873260256537, 0.852185710838635, 0, 0.389602053633798, 0, -0.380323064996537, -0.0880126751081291, 0],
		[-2.11599237942931, 1.47405211856417, -1.37810504103118, -1.07311968166742, 0, 0.336513504829882, 0, 0],
		[-0.381092000896248, 0.732824259094057, 0, 0.216502277780155, 0, -0.162653108888503, -0.426573570884097, 0],
		[-5.56310544580757, 0.905239391642438, 0, 0.385150957140062, 0, -0.282428134685156, 0, 0],[-4.77682417817789,
		0.879981585446539, 0, 0.310609998745732, 0, -0.339225914155431, 0, 0],[0.966712608298821, -0.110938996751065,
		0, 0, 0, 0, 0.183289269601829, 0],[-2.16587686889173, 0.321501356495567, 0, 0, 0, 0, 0, 0],[-1.70734873800232,
		0.433032709755610, 0, -0.412648447269525, 0, 0, 0, 0],[-0.263198077701693, 1.13060571952101, -1.16957372638359,
		-1.65164476300789, 0.104746885300915, 0.404058507352901, 0, 0],[-0.515969600508314, 0.754135122993588, 0,
		0.191575083960373, 0, -0.121665118341219, -0.423844926774291, 0],[-5.77208004012831, 0.923144661954273, 0, 0.402940883581593,
		0, -0.238200773846646, 0, 0],[-5.01588271867143, 0.905027154000921, 0, 0.326841161038211, 0, -0.328283913521590, 0, 0],
		[0.434339606308037, -0.125225415996048, 0, 0, 0, 0, 0.301631247865572, 0],[-2.87544520181323, 0.415682222485807,
		0, 0, 0, 0, 0, 0],[-1.86755738290362, 0.457448335779201, 0, -0.501103981295545, 0, 0, 0, 0]])
		
		self.regressionCoeff_noPulse=np.array([[8.09695881287823, 1.00609515629221, -1.39347614723327, -4.85869770683701,
		0.472644100309933, 0.434550762616159, -0.862562872197509, 0],[-1.03473761679032, 0.769091178587874, 0, 0.412237308297152,
		0, -0.377739650769220, -0.424234099315427, 0],[-4.72728279119446, 0.709717476708319, 0, 0.470974168011549, 0,
		-0.123518047425648, 0, 0],[-4.44400222195478, 0.798093247074753, 0, 0.345405210060350, 0, -0.230823340141895, 0, 0],
		[0.247133528936450, -0.149209862203390, 0, 0, 0, 0, 0.377202902904920, 0],[-1.44302935447839, 0.223053706671624, 0, 0,
		0, 0, 0, 0],[-0.380413278316438, 0.159342468070527, 0, -0.298208438215333, 0, 0, 0, 0],[7.30682757526241, 0.999256668956432,
		-1.33082594407524, -4.95306361630276, 0.490554994733579, 0.442502068793772, -0.835310070621911, 0],[-0.403711730133755,
		0.672375321924977, 0, 0.335372498461681, 0, -0.330322239630250, -0.366700025738387, 0],[-4.79820204505010, 0.709160958437296,
		0, 0.472560804537015, 0, -0.0755764830052928, 0, 0],[-4.35041760661412, 0.785290791385159, 0, 0.325462132630085, 0,
		-0.221525656800750, 0, 0],[0.424849811725595, -0.181204207590470, 0, 0, 0, 0, 0.401549107903204, 0],[-2.97911606394595,
		0.420016455603546, 0, 0, 0, 0, 0, 0],[-0.703694160589291, 0.160571013696218, 0, -0.145792047865653, 0, 0, 0, 0]])
		
		self.parameters_lower_bound=np.array([0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -3.5, -4.7, 0, 0, 0, 0, 0, -3.5, -4.7])
		self.parameters_upper_bound=np.array([0, 0, 3.2, 2, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 0, 1.5, 0])
		
		self.fitted_params_1=np.array([0, 0, 1.30326178289206, 0, 0, 0, 0, 0, 0, 0, 14.2935537214223, 5.33551936215137, 0, 0,
		0, 0, 0, 14.2935537214223, 5.33551936215137])
		self.fitted_params_2=np.array([0, 0, 3.96858083951547, 0, 0, 0, 0, 0, 0, 0, 6.40242376475815, 3.82954843573707, 0, 0,
		0, 0, 0, 6.40242376475815, 3.82954843573707])
		self.fitted_params_3=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.42179588354923, 0, 0, 0, 0, 0, 0, 4.42179588354923, 0])
			
	def numPulseMotionCal(self,numMotions):
		"""
		---速度脉冲地震动模拟数量计算方法---
		输入:
			numMotions(int): the total number of the near fault ground motions need to be simulated
		返回:
			numPulseMotion(int): the number of pulse ground motions
		"""
		###基于Shahi & Baker(2014)建立的模型计算平面双向地震动是否是脉冲地震动的概率
		pulseProb=0.0
		if self.faultType=="strikeSlip":
			pulseProb=1.0/(1.0+np.exp(0.457+0.126*self.Rrup-0.244*self.s_or_d**0.5+0.013*self.theta_or_phi))
		elif self.faultType=="reverseAndReverseOblique":
			pulseProb=1.0/(1.0+np.exp(0.304+0.072*self.Rrup-0.208*self.s_or_d**0.5+0.021*theta_or_phi))
		else:
			raise ValueError("faultType should be 'strikeSlip' or 'reverseAndReverseOblique'")
		uniformRandoms=np.random.uniform(0,1,numMotions)
		numPulseMotion=np.sum([1 if pulseProb>eachItem else 0 for eachItem in uniformRandoms ])
		return numPulseMotion
	
	def inv_double_exp(self,probability,param_a,param_b,param_c,lower_bound):
		"""双指数函数反函数值求取方法"""
		if probability<0 or probability>1:
			raise ValueError("probability argument less than 0 or greater than 1")
		location_inv=(1.0/param_b)*np.log((param_b/param_c)*probability+np.exp(param_b*lower_bound))
		if location_inv<lower_bound or location_inv>0.0:
			location_inv=-(1.0/param_a)*np.log((param_a/param_b)*(1.0-np.exp(param_b*lower_bound))-
			                                   (param_a/param_c)*probability+1.0)
		return location_inv
		
	def transformParametersToOriginalSpace(self,pulse,modelParams):
		"""将模拟生成的参数变换回初始空间方法"""
		backModelParams=np.zeros((np.size(modelParams)))
		if pulse:
			indices=[0,1,4,5,6,7,8,9,12,13,14,15,16]
			for i1 in range(len(indices)):
				indexValue=indices[i1]
				backModelParams[indexValue]=np.exp(modelParams[0,indexValue])
			####--计算gamma
			gammaP=norm.cdf(modelParams[0,2])
			gamaDist=beta.ppf(gammaP,self.fitted_params_1[2],self.fitted_params_2[2])##反函数
			backGamma=gamaDist*(self.parameters_upper_bound[2]-self.parameters_lower_bound[2])+self.parameters_lower_bound[2]
			backModelParams[2]=backGamma
			###--计算niu
			niuP = norm.cdf(modelParams[0, 3])
			uniformDist = uniform.ppf(niuP,self.parameters_lower_bound[3],self.parameters_upper_bound[3])
			backModelParams[3]=uniformDist
			##--计算f',res
			doubleExpP=norm.cdf(modelParams[0, 10])
			invValue=self.inv_double_exp(doubleExpP,self.fitted_params_1[10],self.fitted_params_2[10],
			                             self.fitted_params_3[10],self.parameters_lower_bound[10])
			backModelParams[10]=invValue
			##--计算zeta,res
			zetaP = norm.cdf(modelParams[0, 11])##cdf累计概率函数
			zetaDist = beta.ppf(zetaP, self.fitted_params_1[11], self.fitted_params_2[11])##ppf已知概率求临界值
			backZeta = zetaDist * (self.parameters_upper_bound[11] - self.parameters_lower_bound[11]) + \
			            self.parameters_lower_bound[11]
			backModelParams[11] = np.exp(backZeta)
			##--计算f'Ia_po
			Ia_poP = norm.cdf(modelParams[0, 17])
			invValueIa_po = self.inv_double_exp(Ia_poP, self.fitted_params_1[17], self.fitted_params_2[17],
			                               self.fitted_params_3[17], self.parameters_lower_bound[17])
			backModelParams[17] = invValueIa_po
			##--计算zeta_po
			zeta_poP = norm.cdf(modelParams[0, 18])
			zeta_poDist = beta.ppf(zeta_poP, self.fitted_params_1[18], self.fitted_params_2[18])
			backZeta_po = zeta_poDist * (self.parameters_upper_bound[18] - self.parameters_lower_bound[18]) + \
			           self.parameters_lower_bound[18]
			backModelParams[18] = np.exp(backZeta_po)
		else:
			indices = [0, 1, 2,3,4,7,8,9,10,11]
			for i1 in range(len(indices)):
				indexValue=indices[i1]
				backModelParams[indexValue]=np.exp(modelParams[0,indexValue])
			##--计算f'_component 1
			doubleExpP = norm.cdf(modelParams[0, 5])
			invValue = self.inv_double_exp(doubleExpP, self.fitted_params_1[10], self.fitted_params_2[10],
			                               self.fitted_params_3[10], self.parameters_lower_bound[10])
			backModelParams[5] = invValue
			##--计算zeta_componnet 1
			zetaP = norm.cdf(modelParams[0, 6])
			zetaDist = beta.ppf(zetaP, self.fitted_params_1[11], self.fitted_params_2[11])
			backZeta = zetaDist * (self.parameters_upper_bound[11] - self.parameters_lower_bound[11]) + \
			           self.parameters_lower_bound[11]
			backModelParams[6] = np.exp(backZeta)
			##--计算f'_component 2
			Ia_poP = norm.cdf(modelParams[0, 12])
			invValueIa_po = self.inv_double_exp(Ia_poP, self.fitted_params_1[17], self.fitted_params_2[17],
			                                    self.fitted_params_3[17], self.parameters_lower_bound[17])
			backModelParams[12] = invValueIa_po
			##--计算zeta_component 2
			zeta_poP = norm.cdf(modelParams[0, 13])
			zeta_poDist = beta.ppf(zeta_poP, self.fitted_params_1[18], self.fitted_params_2[18])##类似于标准分布，x取值在[0,1]之间
			backZeta_po = zeta_poDist * (self.parameters_upper_bound[18] - self.parameters_lower_bound[18]) + \
			              self.parameters_lower_bound[18]
			backModelParams[13] = np.exp(backZeta_po)
		return backModelParams
		
	def modelParametersSimulate(self,numMotion,pulse=True):
		"""模拟生成参数方法"""
		if pulse:
			rows,cols=np.shape(self.corrMatrix_pulse)
			covarianceMatrix=np.zeros((rows,cols))
			for i1 in range(rows):
				for j1 in range(cols):
					covarianceMatrix[i1,j1]=self.corrMatrix_pulse[i1,j1]*self.sigma_pulse[i1]*self.sigma_pulse[j1]##方差矩阵
		else:
			rows, cols = np.shape(self.corrMatrix_noPulse)
			covarianceMatrix = np.zeros((rows, cols))
			for i1 in range(rows):
				for j1 in range(cols):
					covarianceMatrix[i1, j1] = self.corrMatrix_noPulse[i1, j1] * self.sigma_noPulse[i1] * self.sigma_noPulse[j1]
		if pulse:
			simulatedParas=np.zeros((numMotion,19))
		else:
			simulatedParas=np.zeros((numMotion,14))
		depth_parameter=self.zTOR if self.zTOR<1.0 else 1.0
		site_parameter=self.Vs30 if self.Vs30<=1100.0 else 1100.0
		lnVs30=np.log(site_parameter)
		fault_parameter=0 if self.faultType=="strikeSlip" else 1
		M65=self.momentMag-6.5 if self.momentMag>6.5 else 0
		Fflt=fault_parameter*depth_parameter
		lnR=np.log((self.Rrup**2+6**2)**0.5)
		MlnR=self.momentMag*np.log((self.Rrup**2+6**2)**0.5)
		itemParams=np.array([1.0,self.momentMag,M65,lnR,MlnR,Fflt,lnVs30,self.s_or_d])
		if pulse:
			expectedParamsValues=np.dot(self.regressionCoeff_pulse,itemParams)##各个强度指标与输入参数回归关系
		else:
			expectedParamsValues = np.dot(self.regressionCoeff_noPulse, itemParams)
		error_mean=np.zeros(19) if pulse else np.zeros(14)
		for i1 in range(numMotion):
			test=-1
			while (test<0.0):##确保脉冲模型起始位置不能小于0
				multiNormRandom=np.random.multivariate_normal(error_mean,covarianceMatrix,size=1)
				epsilon=multiNormRandom/self.sigma_pulse if pulse else multiNormRandom/self.sigma_noPulse
				max_epsilon=np.max(np.abs(epsilon))##抽取的样本到均值距离与标准差的比值
				while max_epsilon>2.0:##避免随机生成的参数距离目标分布均值过远
					multiNormRandom = np.random.multivariate_normal(error_mean, covarianceMatrix, size=1)
					epsilon = multiNormRandom / self.sigma_pulse if pulse else multiNormRandom / self.sigma_noPulse
					max_epsilon = np.max(np.abs(epsilon))
				model_params=expectedParamsValues+multiNormRandom##均值加上随机误差项
				backTransfParams=self.transformParametersToOriginalSpace(pulse,model_params)
				if pulse:
					test=backTransfParams[4]-0.5*backTransfParams[1]*backTransfParams[2]
				else:
					test=1.0
			simulatedParas[i1,:]=backTransfParams
		return simulatedParas
	
	def sumFunc(self,p):
		""""""
		alpha,beta,t_max_q,d05,d030,d095=p
		t0=0
		t5_fit =t0 + (((5.0 / 100.0) * (t_max_q - t0)**(2.0 * alpha)) *((t_max_q - t0) + (2.0 * alpha + 1.0) / (2.0 * beta)))\
		        **(1.0 / (2.0 * alpha + 1.0))
		t30_fit =t0 + (((30.0 / 100.0) * (t_max_q - t0)**(2.0 * alpha)) *((t_max_q - t0) + (2.0 * alpha + 1.0) / (2.0 * beta)))\
		        **(1.0 / (2.0 * alpha + 1.0))
		t95_fit = t0 + (((95.0 / 100.0) * (t_max_q - t0) ** (2.0 * alpha)) * ((t_max_q - t0) + (2.0 * alpha + 1.0) / (2.0 * beta))) \
		          ** (1.0 / (2.0 * alpha + 1.0))
		if t5_fit>t_max_q:
			t5_fit = t_max_q -(1.0 / (2.0 * beta)) *np.log(((100.0 - 5.0) / 100.0) *((t_max_q - t0) * (2.0 * beta) / (2.0 * alpha + 1.0) +1.0))
		if t30_fit > t_max_q:
			t30_fit = t_max_q - (1.0 / (2.0 * beta)) * np.log(((100.0 - 30.0) / 100.0) * ((t_max_q - t0) * (2.0 * beta) / (2.0 * alpha + 1.0) + 1.0))
		if t95_fit > t_max_q:
			t95_fit = t_max_q - (1.0 / (2.0 * beta)) * np.log(((100.0 - 95.0) / 100.0) * ((t_max_q - t0) * (2.0 * beta) / (2.0 * alpha + 1.0) + 1.0))
		d05_fit = t5_fit - t0
		d030_fit = t30_fit - t0
		d095_fit = t95_fit - t0
		errorSum=(d05-d05_fit)**2 +(d030-d030_fit)**2 +(d095-d095_fit)**2
		return errorSum
		
	def backcalculate_modulating_params(self,q_params,t0):
		"""基于差分进化算法的调制函数参数优化方法"""
		arias_intensity = q_params[0] / 981
		d595 = q_params[1]
		d05 = q_params[2]
		d030 = q_params[3]
		t30 = t0 + d030
		t95=d595+d05+t0
		d095=d595+d05
		########alpha,beta,tmax,q,c
		alpha = [0.01, 10]
		beta = [0.01, 2]
		tmaxq = [t0, t95*2]
		
		de = DE(func=self.sumFunc, n_dim=6, size_pop=50, max_iter=150, lb=[alpha[0], beta[0], tmaxq[0],d05,d030,d095],
		        ub=[alpha[1], beta[1], tmaxq[1],d05,d030,d095])
		best_x, best_y = de.run()
		# plt.plot(de.generation_best_Y)
		# plt.show()
		alpha_opt,beta_opt,tmaxq_opt,_,_,_=best_x
		c_opt= (arias_intensity /((np.pi / 2.0) * ((tmaxq_opt- t0) / (2.0 * alpha_opt + 1.0) +1.0 / (2.0 * beta_opt))))**0.5
		return alpha_opt,beta_opt,tmaxq_opt,c_opt

	@staticmethod
	@jit(nopython=True)
	def calc_time_to_intensity(acceleration, percentage, time_step):  ##percentage是乘以100后的数据
		"""不同Arias强度对应的时间函数"""
		t01_sum = [each ** 2 for each in acceleration]
		t01_sumCum = [0.0]
		[t01_sumCum.append(t01_sumCum[-1] + t01_sum[i1]) for i1 in range(len(t01_sum))]
		t01_sum = t01_sumCum[1:]
		t01_sum = [each * 100.0 / t01_sum[-1] for each in t01_sum]  ##百分数
		tIndex = 0
		for i2 in prange(len(t01_sum)):
			if t01_sum[i2] >= percentage:
				tIndex = i2 + 1
				break
		return tIndex * time_step


	@staticmethod
	@jit(nopython=True)
	def calc_linear_filter(num_steps, filter_params, t01, tmid, t99, time_step_):
		""""""
		filter_func = np.zeros(num_steps)
		##--Mininum frequency in Hz
		min_freq = 0.3
		##--Frequency at tmid, in Hz
		mid_freq = filter_params[0]
		##--Slope of frequency assumed constant
		freq_slope = filter_params[1]
		for i in prange(num_steps):
			current_time = i * time_step_
			if (current_time < t01):
				filter_func[i] = min_freq * 2.0 * np.pi if min_freq > (mid_freq + freq_slope * (t01 - tmid)) \
					else (mid_freq + freq_slope * (t01 - tmid)) * 2.0 * np.pi
			elif (current_time <= t99):
				filter_func[i] = min_freq * 2.0 * np.pi if min_freq > mid_freq + freq_slope * (current_time - tmid) \
					else (mid_freq + freq_slope * (current_time - tmid)) * 2.0 * np.pi
			else:
				filter_func[i] = min_freq * 2.0 * np.pi if min_freq > mid_freq + freq_slope * (t99 - tmid) \
					else (mid_freq + freq_slope * (t99 - tmid)) * 2.0 * np.pi
		return filter_func

	@staticmethod
	@jit(nopython=True)
	def calc_impulse_response_filter(num_steps, input_filter, zeta, time_step_):
		"""计算脉冲响应函数"""
		impulse_response = np.zeros((num_steps, num_steps))  ##矩阵
		for i in prange(num_steps):
			omega = input_filter[i]
			times = np.array([j * time_step_ for j in range(num_steps - i)])

			res = ((omega / np.sqrt(1.0 - zeta * zeta)) * np.exp(-zeta * omega * times) *
				   np.sin((omega * np.sqrt(1.0 - zeta * zeta) * times)))
			for i1 in prange(len(times)):
				impulse_response[i, i + i1] = res[i1]
		denominator1 = np.power(impulse_response, 2).sum(axis=0)

		denominator = np.power(denominator1, 0.5)
		denominator[0] = 0.1
		for i2 in prange(num_steps):
			for j2 in range(num_steps):
				impulse_response[i2, j2] = impulse_response[i2, j2] / denominator[j2]
		return impulse_response


	def simulate_white_noise(self,modulating_params,filter_params,num_steps,time_step_,t0):
		"""滤过白噪声模拟方法"""
		###---计算调幅函数
		mod_func_vals=np.zeros(num_steps)
		for i in range(num_steps):
			time=i*time_step_
			if time<t0:
				mod_func_vals[i]=0.0
			elif time<modulating_params[2]:
				mod_func_vals[i]=modulating_params[3]*((time - t0) / (modulating_params[2] - t0))**modulating_params[0]
			else:
				mod_func_vals[i] =modulating_params[3]*np.exp(-modulating_params[1]*(time - modulating_params[2]))
		###--计算频率函数
		t01 = NearFaultMotionSim.calc_time_to_intensity(mod_func_vals, 1.0,time_step_)##1%对应的Arias强度对应的时间
		tmid = NearFaultMotionSim.calc_time_to_intensity(mod_func_vals, 30.0,time_step_)##30%对应的Arias强度对应的时间
		t99 = NearFaultMotionSim.calc_time_to_intensity(mod_func_vals, 99.0,time_step_)
		###--定义滤过频率及带宽
		frequency_filter =NearFaultMotionSim.calc_linear_filter(num_steps, filter_params, t01, tmid, t99,time_step_)
		###---生成白噪声
		white_noise= np.random.normal(loc=0, scale=1, size=num_steps)##随机正态分布数
		
		###---计算单位脉冲响应
		impulse_response = NearFaultMotionSim.calc_impulse_response_filter(num_steps, frequency_filter, filter_params[2],time_step_)##脉冲响应
		freq_func=np.dot(white_noise,impulse_response)##第一个点只有第一个脉冲有贡献，最后一个点前面所有单位脉冲都有贡献
		####---生成滤过白噪声
		filtered_white_noise=freq_func*mod_func_vals##经过时、频域调节的地震动时程
		return filtered_white_noise
	
	def filter_acceleration(self,accel_history, freq_corner, filter_order):
		"""基于快速傅里叶变换的高通滤波函数"""
		###---计算归一化的阶段频率
		accel_fft= np.fft.fft(accel_history)
		freq = np.fft.fftfreq(len(accel_history), self.timeStep).tolist()
		accelFftAbs=abs(accel_fft)
		num_samples=len(accel_fft)
		
		filter=np.zeros(int((num_samples / 2) +1))
		for i2 in range(len(filter)):
			tempFreq=(freq[i2]/freq_corner)**(2.0*filter_order)
			filter[i2]=np.sqrt(tempFreq/(1.0+tempFreq))
		highpass_filter = np.zeros(2 * len(filter) - 2)
		highpass_filter[:len(filter)] = filter
		highpass_filter[len(filter):] = filter[1:(len(filter) - 1)][::-1]
		accel_fft_filt=accel_fft*highpass_filter
		filtered_acc=np.real(np.fft.ifft(accel_fft_filt))
		return filtered_acc

	@staticmethod
	@jit(nopython=True)
	def calc_pulse_acceleration(num_steps, parameters, start_time_, time_step_):
		"""低频脉冲分量时程计算函数"""
		pulse_velocity = parameters[0]
		pulse_frequency = 1.0 / parameters[1]  ##脉冲频率
		oscillation_param = parameters[2]  ##震荡参数
		phase_angle = parameters[3] * np.pi  ##相角
		peak_time = start_time_ + parameters[4]  ##脉冲峰值对应的时刻
		sinplus = np.sin(phase_angle + oscillation_param * np.pi)
		sinMinus = np.sin(phase_angle - oscillation_param * np.pi)
		gamma2 = 1 - oscillation_param ** 2
		resp_disp = pulse_velocity * (sinplus - sinMinus) / (4.0 * np.pi * gamma2 * pulse_frequency)
		##--计算速度脉冲时程
		velocity_history = np.zeros(num_steps)
		for i in prange(num_steps):
			time = i * time_step_
			time1 = peak_time - 0.5 * oscillation_param / pulse_frequency
			time2 = peak_time + 0.5 * oscillation_param / pulse_frequency
			if (time > time1) and (time <= time2):
				velocity_history[i] = (0.5 * pulse_velocity * np.cos(
					2.0 * np.pi * pulse_frequency * (time - peak_time) +
					phase_angle) - resp_disp * pulse_frequency / oscillation_param) * (
												  1.0 + np.cos(2.0 * np.pi * pulse_frequency *
															   (time - peak_time) / oscillation_param))
		acc_pulse = derivative(velocity_history, 1.0 / (981.0 * time_step_), True)
		return acc_pulse
		
	def simulate_near_fault_ground_motion(self,pulse,parameters):
		"""模拟近断层地震动方法"""
		alpha_1=np.zeros(7)##高频部分1
		alpha_2=np.zeros(7)##高频部分2
		if pulse:
			alpha_1=parameters[5:12]
			alpha_2=parameters[12:]
		else:
			alpha_1=parameters[0:7]
			alpha_2=parameters[7:14]
		###---设置调幅及滤过参数-基于DE算法
		modulating_params_1 =self.backcalculate_modulating_params(alpha_1[0:4],self.start_time)
		modulating_params_2 =self.backcalculate_modulating_params(alpha_2[0:4],self.start_time)
		filter_params_1 = alpha_1[4:7]
		filter_params_2 = alpha_2[4:7]
		###--确定模拟地震动长度
		time1=self.start_time+alpha_1[1]+alpha_1[2]
		time2=self.start_time+alpha_2[1]+alpha_2[2]
		t95=time1 if time1>time2 else time2
		
		num_steps=int(np.ceil(2.5 * t95 /self.timeStep))
		num_steps=num_steps+1 if num_steps%2==1 else num_steps ##偶数个数
		###---生成调幅滤过白噪声
		white_noise_1 =self.simulate_white_noise(modulating_params_1, filter_params_1, num_steps,
		                                                       self.timeStep,self.start_time)
		white_noise_2 = self.simulate_white_noise(modulating_params_2, filter_params_2, num_steps,
		                                                        self.timeStep,self.start_time)
		###---高通滤波及加0
		freq_corner =10.0**(1.4071 - 0.3452 * self.momentMag)  ##In Hz0
		filter_order = 4
		padding_duration = 0.5 * 1.5 * filter_order / freq_corner
		num_pads =int(np.ceil(padding_duration/self.timeStep))
		accel_padded_1=np.zeros(num_pads+num_steps+num_pads)
		accel_padded_2 = np.zeros(num_pads + num_steps + num_pads)
		accel_padded_1[num_pads:(num_steps+num_pads)]=white_noise_1
		accel_padded_2[num_pads:(num_steps + num_pads)] = white_noise_2
		###---加0后加速度时程滤波
		accel_comp_1=self.filter_acceleration(accel_padded_1, freq_corner, filter_order)##FFT高通滤波
		accel_comp_2=self.filter_acceleration(accel_padded_2, freq_corner, filter_order)
		###---基于AI强度一致缩放时程
		target_ai_1 = alpha_1[0] / 981.0
		target_ai_2 = alpha_2[0] / 981.0
		##--计算Arias强度
		arias_intensity_1=[each**2*self.timeStep*np.pi/2.0 for each in accel_comp_1]
		arias_intensity_2 = [each ** 2 * self.timeStep * np.pi / 2.0 for each in accel_comp_2]
		arias_intensity_1Sum = [0.0]
		[arias_intensity_1Sum.append(arias_intensity_1Sum[-1] + arias_intensity_1[i1]) for i1 in range(len(arias_intensity_1))]
		Arias_comp1=arias_intensity_1Sum[-1]
		arias_intensity_2Sum = [0.0]
		[arias_intensity_2Sum.append(arias_intensity_2Sum[-1] + arias_intensity_2[i1]) for i1 in range(len(arias_intensity_2))]
		Arias_comp2 = arias_intensity_2Sum[-1]
		scale_factor_1=(target_ai_1/Arias_comp1)**0.5
		scale_factor_2 = (target_ai_2 / Arias_comp2) ** 0.5
		accel_comp_1_scale=[each*scale_factor_1 for each in accel_comp_1]
		accel_comp_2_scale = [each * scale_factor_2 for each in accel_comp_2]
		###--如果是脉冲地震动，将脉冲加速度时程叠加到component_1上
		#################################################
		pulseZeroList=[]
		if pulse:
			pulse_accel = NearFaultMotionSim.calc_pulse_acceleration(num_steps, parameters,self.start_time,self.timeStep)##低频脉冲分量时程
			##--将加速度脉冲分量叠加到component 1
			pulseZeroList=[0.0 for ij in range(len(accel_comp_1_scale))]

			for j in range(len(pulse_accel)):
				accel_comp_1_scale[j + num_pads - 1] =accel_comp_1_scale[j + num_pads - 1] + pulse_accel[j]
				pulseZeroList[j + num_pads - 1]=pulseZeroList[j + num_pads - 1]+pulse_accel[j]
		###--地震动截断处理
		vel1=accToVelocity(self.timeStep,accel_comp_1_scale)##---加速度转换为速度
		disp1=velToDisplacement(self.timeStep,vel1)##---速度转换为位移
		maxdisp=0.01 ##单位cm
		index1=next((n1 for n1, v1 in enumerate(disp1) if abs(v1) >= maxdisp))##---生成器式查找
		######
		vel2 = accToVelocity(self.timeStep, accel_comp_2_scale)
		disp2 = velToDisplacement(self.timeStep, vel2)
		index2 = next((n2 for n2, v2 in enumerate(disp2) if abs(v2) >= maxdisp))
		indexValue=index1 if index1<index2 else index2
		#################################################
		returnPulseModelAcc=None
		if len(pulseZeroList)>0:
			returnPulseModelAcc=pulseZeroList[indexValue:-num_pads]
		# np.savetxt("pulse_acc.txt",pulseZeroList[indexValue:-num_pads],fmt="%.6f")
		# np.savetxt("pulseMotion_acc.txt",accel_comp_1_scale[indexValue:-num_pads],fmt="%.6f")
		##################################################
		return accel_comp_1_scale[indexValue:-num_pads],accel_comp_2_scale[indexValue:-num_pads],returnPulseModelAcc
	
	def runSimulation(self):
		"""
		---模拟生成地震动时程方法---
		"""
		#####---1.模拟脉冲与非脉冲参数
		#########################################################
		###---参数修改部分
		parameters_pulse=self.modelParametersSimulate(self.numPulseMotion,pulse=True)
		# if np.abs(parameters_pulse[0,1]-5)<1000:
		# 	pass
		# else:
		# 	return None,None
		print("parameters_pulse",parameters_pulse)
		#########################################################
		#########################################################
		parameters_noPulse = self.modelParametersSimulate(self.numNoPulseMotion, pulse=False)
		#####---2.模拟脉冲地震动
		##########################
		pulseComponnetsList=[]
		noPulseComponnetsList=[]
		for i0 in trange(self.numSimMotions):
			if i0<self.numPulseMotion:##脉冲地震动模拟
				pulse_comp1, pulse_comp2 = self.simulate_near_fault_ground_motion(True, parameters_pulse[i0, :])
				pulseComponnetsList.append([pulse_comp1,pulse_comp2])
			else:##非脉冲地震动模拟
				noPulse_comp1, noPulse_comp2 = self.simulate_near_fault_ground_motion(False, parameters_noPulse[int(i0-self.numPulseMotion), :])
				noPulseComponnetsList.append([noPulse_comp1,noPulse_comp2])
		return pulseComponnetsList,noPulseComponnetsList##---加速度时程列表[[comp1_acc],...]

	def runSpecificPulseMotion(self,Vp=None,Tp=None):
		"""
		---模拟特定脉冲幅值，周期地震动---
		Vp(cm/s)-速度脉冲幅值
		Tp(s)-脉冲周期
		"""
		#####---1.模拟脉冲与非脉冲参数
		#########################################################
		parameters_pulse = self.modelParametersSimulate(self.numPulseMotion, pulse=True)
		if Vp is not None:
			if np.abs(parameters_pulse[0,0]-Vp)<0.01*Vp:
				pass
			else:
				return None,None
		if Tp is not None:
			if np.abs(parameters_pulse[0,1]-Tp)<0.01*Tp:
				pass
			else:
				return None,None
		pulseParas=parameters_pulse[0,:5]
		pulseComponnetsList = []
		for i0 in trange(self.numSimMotions):
			pulse_comp1, pulse_comp2,returnPulseModelAcc= self.simulate_near_fault_ground_motion(True, parameters_pulse[i0, :])
			pulseComponnetsList.append([pulse_comp1,pulse_comp2,returnPulseModelAcc])
		return pulseComponnetsList,pulseParas
########################################################################################################################
########################################################################################################################
@jit(nopython=True)
def derivative(coefficients,constant_factor,add_zero):
	"""求导计算函数"""
	if add_zero:
		derivative=np.zeros(len(coefficients))
		derivative[0] = coefficients[0] * constant_factor
		for i in prange(1,len(derivative)):
			derivative[i] = (coefficients[i] - coefficients[i - 1]) * constant_factor
		return derivative
	else:
		derivative=np.zeros(len(coefficients) - 1)
		for i in prange(len(derivative)):
			derivative[i] = (coefficients[i + 1] - coefficients[i]) * constant_factor
		return derivative
########################################################################################################################
########################################################################################################################
def responseSpectraPlot(folderName,maxPeriod=10,beta=0.05,logCoord=True,SaLim=[0.01,10],SvLim=[0.01,500],SdLim=[0.001,500]):
	"""
	plot acceleration, velocity and displacement response spectra. time interval is 0.005 s
	----------------------------------------------------------------------------------------------------------------
	folderName(str)-the folder name that sotore generated motions
	maxPeriod(float)-the maximum period for response spectra plotting
	beta(float)-damping ratio, default value is 0.05
	lodCoord(bool)-Base-10 logarithmic coordinates, the default Value is True
	SaLim,SvLim,SdLim(list)-the minimum and maximum value for y coordinate
	"""
	f = [np.exp(i0 * 0.035 - 3) for i0 in range(200)]  ##frequencies(Hz)0.05~50Hz
	TList = [1.0 / each1 for each1 in f]
	TList.reverse()
	relaFolderName = folderName.split('\\')[-1]
	##############
	##############
	files = [f.name for f in Path(folderName).iterdir() if f.is_file()]
	motionNumbers=list(set([each.split('_')[0] for each in files]))
	totalMotion = int(len(files) / 4)
	dt = 0.005  ###---time interval 0.005 s
	pulseComp_Sa=[]
	pulseComp_Sv=[]
	pulseComp_Sd=[]
	orthogComp_Sa=[]
	orthogComp_Sv=[]
	orthogComp_Sd=[]
	pulseModel_Sa=[]
	pulseMode_Sv=[]
	pulseModel_Sd=[]
	for i1 in range(totalMotion):
		###---pulseComp acc(g)
		pulseCompAcc = np.loadtxt(folderName + f"/{motionNumbers[i1]}_pulseComp.txt")
		saArray, svArray, sdArray=ResponseCal(beta=beta,dt=dt,Tarray=TList,accArray=pulseCompAcc)
		pulseComp_Sa.append(saArray)
		pulseComp_Sv.append(svArray)
		pulseComp_Sd.append(sdArray)
		###---orthogComp acc(g)
		orthogCompAcc = np.loadtxt(folderName + f"/{motionNumbers[i1]}_orthogComp.txt")
		saArray, svArray, sdArray = ResponseCal(beta=beta,dt=dt,Tarray=TList,accArray=pulseCompAcc)
		orthogComp_Sa.append(saArray)
		orthogComp_Sv.append(svArray)
		orthogComp_Sd.append(sdArray)
		###---pulseModel acc(g)
		pulseModelAcc = np.loadtxt(folderName + f"/{motionNumbers[i1]}_pulseModel.txt")
		saArray, svArray, sdArray = ResponseCal(beta=beta,dt=dt,Tarray=TList,accArray=pulseCompAcc)
		pulseModel_Sa.append(saArray)
		pulseMode_Sv.append(svArray)
		pulseModel_Sd.append(sdArray)
	############################################
	fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3,figsize=(20,15))
	if logCoord:###-log-log coordinate
		####################################################################################
		###-ax1 set
		ax1.set_title('PulseComp', fontsize=11)
		ax1.set_ylabel('Sa (g)', fontsize=10)
		ax1.grid(which="both")
		ax1.set_xlim(TList[0],maxPeriod)
		ax1.set_ylim(SaLim[0],SaLim[1])
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		for i2 in range(totalMotion):
			ax1.plot(TList, pulseComp_Sa[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax1_mean = np.mean(pulseComp_Sa, axis=0)
		ax1_std = np.std(pulseComp_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax1_upper = ax1_mean + ax1_std
		ax1_lower = ax1_mean - ax1_std
		ax1.plot(TList, ax1_mean, "r", linewidth=1.5, label='mean')
		ax1.plot(TList, ax1_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax1.plot(TList, ax1_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax2 set
		ax2.set_title('OrthogComp', fontsize=11)
		ax2.set_ylabel('Sa (g)', fontsize=10)
		ax2.grid(which="both")
		ax2.set_xlim(TList[0], maxPeriod)
		ax2.set_ylim(SaLim[0], SaLim[1])
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		for i2 in range(totalMotion):
			ax2.plot(TList, orthogComp_Sa[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax2_mean = np.mean(orthogComp_Sa, axis=0)
		ax2_std = np.std(orthogComp_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax2_upper = ax2_mean + ax2_std
		ax2_lower = ax2_mean - ax2_std
		ax2.plot(TList, ax2_mean, "r", linewidth=1.5, label='mean')
		ax2.plot(TList, ax2_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax2.plot(TList, ax2_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax3 set
		ax3.set_title('pulseModel', fontsize=11)
		ax3.set_ylabel('Sa (g)', fontsize=10)
		ax3.grid(which="both")
		ax3.set_xlim(TList[0], maxPeriod)
		ax3.set_ylim(SaLim[0], SaLim[1])
		ax3.set_xscale("log")
		ax3.set_yscale("log")
		for i2 in range(totalMotion):
			ax3.plot(TList, pulseModel_Sa[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax3_mean = np.mean(pulseModel_Sa, axis=0)
		ax3_std = np.std(pulseModel_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax3_upper = ax3_mean + ax3_std
		ax3_lower = ax3_mean - ax3_std
		ax3.plot(TList, ax3_mean, "r", linewidth=1.5, label='mean')
		ax3.plot(TList, ax3_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax3.plot(TList, ax3_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax4 set
		# ax4.set_title('pulseModel', fontsize=11)
		ax4.set_ylabel('Sv (cm/s)', fontsize=10)
		ax4.grid(which="both")
		ax4.set_xlim(TList[0], maxPeriod)
		ax4.set_ylim(SvLim[0], SvLim[1])
		ax4.set_xscale("log")
		ax4.set_yscale("log")
		for i2 in range(totalMotion):
			ax4.plot(TList, pulseComp_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax4_mean = np.mean(pulseComp_Sv, axis=0)
		ax4_std = np.std(pulseComp_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax4_upper = ax4_mean + ax4_std
		ax4_lower = ax4_mean - ax4_std
		ax4.plot(TList, ax4_mean, "r", linewidth=1.5, label='mean')
		ax4.plot(TList, ax4_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax4.plot(TList, ax4_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax5 set
		# ax5.set_title('pulseModel', fontsize=11)
		ax5.set_ylabel('Sv (cm/s)', fontsize=10)
		ax5.grid(which="both")
		ax5.set_xlim(TList[0], maxPeriod)
		ax5.set_ylim(SvLim[0], SvLim[1])
		ax5.set_xscale("log")
		ax5.set_yscale("log")
		for i2 in range(totalMotion):
			ax5.plot(TList, orthogComp_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax5_mean = np.mean(orthogComp_Sv, axis=0)
		ax5_std = np.std(orthogComp_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax5_upper = ax5_mean + ax5_std
		ax5_lower = ax5_mean - ax5_std
		ax5.plot(TList, ax5_mean, "r", linewidth=1.5, label='mean')
		ax5.plot(TList, ax5_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax5.plot(TList, ax5_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax6 set
		# ax6.set_title('pulseModel', fontsize=11)
		ax6.set_ylabel('Sv (cm/s)', fontsize=10)
		ax6.grid(which="both")
		ax6.set_xlim(TList[0], maxPeriod)
		ax6.set_ylim(SvLim[0], SvLim[1])
		ax6.set_xscale("log")
		ax6.set_yscale("log")
		for i2 in range(totalMotion):
			ax6.plot(TList, pulseMode_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax6_mean = np.mean(pulseMode_Sv, axis=0)
		ax6_std = np.std(pulseMode_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax6_upper = ax6_mean + ax6_std
		ax6_lower = ax6_mean - ax6_std
		ax6.plot(TList, ax6_mean, "r", linewidth=1.5, label='mean')
		ax6.plot(TList, ax6_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax6.plot(TList, ax6_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax7 set
		# ax6.set_title('pulseModel', fontsize=11)
		ax7.set_ylabel('Sd (cm)', fontsize=10)
		ax7.set_xlabel('Period (s)', fontsize=10)
		ax7.grid(which="both")
		ax7.set_xlim(TList[0], maxPeriod)
		ax7.set_ylim(SdLim[0], SdLim[1])
		ax7.set_xscale("log")
		ax7.set_yscale("log")
		for i2 in range(totalMotion):
			ax7.plot(TList, pulseComp_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax7_mean = np.mean(pulseComp_Sd, axis=0)
		ax7_std = np.std(pulseComp_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax7_upper = ax7_mean + ax7_std
		ax7_lower = ax7_mean - ax7_std
		ax7.plot(TList, ax7_mean, "r", linewidth=1.5, label='mean')
		ax7.plot(TList, ax7_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax7.plot(TList, ax7_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax8 set
		# 8.set_title('pulseModel', fontsize=11)
		ax8.set_ylabel('Sd (cm)', fontsize=10)
		ax8.set_xlabel('Period (s)', fontsize=10)
		ax8.grid(which="both")
		ax8.set_xlim(TList[0], maxPeriod)
		ax8.set_ylim(SdLim[0], SdLim[1])
		ax8.set_xscale("log")
		ax8.set_yscale("log")
		for i2 in range(totalMotion):
			ax8.plot(TList, orthogComp_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax8_mean = np.mean(orthogComp_Sd, axis=0)
		ax8_std = np.std(orthogComp_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax8_upper = ax8_mean + ax8_std
		ax8_lower = ax8_mean - ax8_std
		ax8.plot(TList, ax8_mean, "r", linewidth=1.5, label='mean')
		ax8.plot(TList, ax8_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax8.plot(TList, ax8_lower, "b--", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax9 set
		# 8.set_title('pulseModel', fontsize=11)
		ax9.set_ylabel('Sd (cm)', fontsize=10)
		ax9.set_xlabel('Period (s)', fontsize=10)
		ax9.grid(which="both")
		ax9.set_xlim(TList[0], maxPeriod)
		ax9.set_ylim(SdLim[0], SdLim[1])
		ax9.set_xscale("log")
		ax9.set_yscale("log")
		for i2 in range(totalMotion):
			ax9.plot(TList, pulseModel_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax9_mean = np.mean(pulseModel_Sd, axis=0)
		ax9_std = np.std(pulseModel_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax9_upper = ax9_mean + ax9_std
		ax9_lower = ax9_mean - ax9_std
		ax9.plot(TList, ax9_mean, "r", linewidth=1.5, label='mean')
		ax9.plot(TList, ax9_upper, "b--", linewidth=1.5, label=f'mean+std')
		ax9.plot(TList, ax9_lower, "b--", linewidth=1.5, label=f'mean-std')


	else: ###-natural coordinate
		####################################################################################
		###-ax1 set
		ax1.set_title('PulseComp', fontsize=11)
		ax1.set_ylabel('Sa (g)', fontsize=10)
		ax1.grid(which="both")
		ax1.set_xlim(TList[0], maxPeriod)
		ax1.set_ylim(SaLim[0], SaLim[1])
		for i2 in range(totalMotion):
			ax1.plot(TList,pulseComp_Sa[i2], "grey",linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2) ###-pulseComp Sa(g)
		ax1_mean = np.mean(pulseComp_Sa, axis=0)
		ax1_std  = np.std(pulseComp_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax1_upper=ax1_mean+ax1_std
		ax1_lower=ax1_mean-ax1_std
		ax1.plot(TList, ax1_mean, "r", linewidth=1.5, label='mean')
		ax1.plot(TList,ax1_upper, "b", linewidth=1.5,label=f'mean+std')
		ax1.plot(TList, ax1_lower, "b", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax2 set
		ax2.set_title('OrthogComp', fontsize=11)
		ax2.set_ylabel('Sa (g)', fontsize=10)
		ax2.grid(which="both")
		ax2.set_xlim(TList[0], maxPeriod)
		ax2.set_ylim(SaLim[0], SaLim[1])
		for i2 in range(totalMotion):
			ax2.plot(TList, orthogComp_Sa[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax2_mean = np.mean(orthogComp_Sa, axis=0)
		ax2_std = np.std(orthogComp_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax2_upper = ax2_mean + ax2_std
		ax2_lower = ax2_mean - ax2_std
		ax2.plot(TList, ax2_mean, "r", linewidth=1.5, label='mean')
		ax2.plot(TList, ax2_upper, "b", linewidth=1.5, label=f'mean+std')
		ax2.plot(TList, ax2_lower, "b", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax3 set
		ax3.set_title('pulseModel', fontsize=11)
		ax3.set_ylabel('Sa (g)', fontsize=10)
		ax3.grid(which="both")
		ax3.set_xlim(TList[0], maxPeriod)
		ax3.set_ylim(SaLim[0], SaLim[1])
		for i2 in range(totalMotion):
			ax3.plot(TList, pulseModel_Sa[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax3_mean = np.mean(pulseModel_Sa, axis=0)
		ax3_std = np.std(pulseModel_Sa, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax3_upper = ax3_mean + ax3_std
		ax3_lower = ax3_mean - ax3_std
		ax3.plot(TList, ax3_mean, "r", linewidth=1.5, label='mean')
		ax3.plot(TList, ax3_upper, "b", linewidth=1.5, label=f'mean+std')
		ax3.plot(TList, ax3_lower, "b", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax4 set
		# ax4.set_title('pulseModel', fontsize=11)
		ax4.set_ylabel('Sv (cm/s)', fontsize=10)
		ax4.grid(which="both")
		ax4.set_xlim(TList[0], maxPeriod)
		ax4.set_ylim(SvLim[0], SvLim[1])
		for i2 in range(totalMotion):
			ax4.plot(TList, pulseComp_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax4_mean = np.mean(pulseComp_Sv, axis=0)
		ax4_std = np.std(pulseComp_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax4_upper = ax4_mean + ax4_std
		ax4_lower = ax4_mean - ax4_std
		ax4.plot(TList, ax4_mean, "r", linewidth=1.5, label='mean')
		ax4.plot(TList, ax4_upper, "b", linewidth=1.5, label=f'mean+std')
		ax4.plot(TList, ax4_lower, "b", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax5 set
		# ax5.set_title('pulseModel', fontsize=11)
		ax5.set_ylabel('Sv (cm/s)', fontsize=10)
		ax5.grid(which="both")
		ax5.set_xlim(TList[0], maxPeriod)
		ax5.set_ylim(SvLim[0], SvLim[1])
		for i2 in range(totalMotion):
			ax5.plot(TList, orthogComp_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax5_mean = np.mean(orthogComp_Sv, axis=0)
		ax5_std = np.std(orthogComp_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax5_upper = ax5_mean + ax5_std
		ax5_lower = ax5_mean - ax5_std
		ax5.plot(TList, ax5_mean, "r", linewidth=1.5, label='mean')
		ax5.plot(TList, ax5_upper, "b", linewidth=1.5, label=f'mean+std')
		ax5.plot(TList, ax5_lower, "b", linewidth=1.5, label=f'mean-std')
		####################################################################################
		###-ax6 set
		# ax6.set_title('pulseModel', fontsize=11)
		ax6.set_ylabel('Sv (cm/s)', fontsize=10)
		ax6.grid(which="both")
		ax6.set_xlim(TList[0], maxPeriod)
		ax6.set_ylim(SvLim[0], SvLim[1])
		for i2 in range(totalMotion):
			ax6.plot(TList, pulseMode_Sv[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax6_mean = np.mean(pulseMode_Sv, axis=0)
		ax6_std = np.std(pulseMode_Sv, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax6_upper = ax6_mean + ax6_std
		ax6_lower = ax6_mean - ax6_std
		ax6.plot(TList, ax6_mean, "r", linewidth=1.5, label='mean')
		ax6.plot(TList, ax6_upper, "b", linewidth=1.5, label=f'mean+std')
		ax6.plot(TList, ax6_lower, "b", linewidth=1.5, label=f'mean-std')
		###-ax7 set
		# ax7.set_title('pulseModel', fontsize=11)
		ax7.set_ylabel('Sd (cm)', fontsize=10)
		ax7.grid(which="both")
		ax7.set_xlim(TList[0], maxPeriod)
		ax7.set_ylim(SdLim[0], SdLim[1])
		for i2 in range(totalMotion):
			ax7.plot(TList, pulseComp_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax7_mean = np.mean(pulseComp_Sd, axis=0)
		ax7_std = np.std(pulseComp_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax7_upper = ax7_mean + ax7_std
		ax7_lower = ax7_mean - ax7_std
		ax7.plot(TList, ax7_mean, "r", linewidth=1.5, label='mean')
		ax7.plot(TList, ax7_upper, "b", linewidth=1.5, label=f'mean+std')
		ax7.plot(TList, ax7_lower, "b", linewidth=1.5, label=f'mean-std')
		###-ax8 set
		# ax8.set_title('pulseModel', fontsize=11)
		ax8.set_ylabel('Sd (cm)', fontsize=10)
		ax8.grid(which="both")
		ax8.set_xlim(TList[0], maxPeriod)
		ax8.set_ylim(SdLim[0], SdLim[1])
		for i2 in range(totalMotion):
			ax8.plot(TList, orthogComp_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax8_mean = np.mean(orthogComp_Sd, axis=0)
		ax8_std = np.std(orthogComp_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax8_upper = ax8_mean + ax8_std
		ax8_lower = ax8_mean - ax8_std
		ax8.plot(TList, ax8_mean, "r", linewidth=1.5, label='mean')
		ax8.plot(TList, ax8_upper, "b", linewidth=1.5, label=f'mean+std')
		ax8.plot(TList, ax8_lower, "b", linewidth=1.5, label=f'mean-std')
		###-ax9 set
		# ax8.set_title('pulseModel', fontsize=11)
		ax9.set_ylabel('Sd (cm)', fontsize=10)
		ax9.grid(which="both")
		ax9.set_xlim(TList[0], maxPeriod)
		ax9.set_ylim(SdLim[0], SdLim[1])
		for i2 in range(totalMotion):
			ax9.plot(TList, pulseModel_Sd[i2], "grey", linewidth=1,label=f'motion {motionNumbers[i2]}',picker=2)  ###-pulseComp Sa(g)
		ax9_mean = np.mean(pulseModel_Sd, axis=0)
		ax9_std = np.std(pulseModel_Sd, axis=0, ddof=0)  # ddof=1 为样本标准差；若要总体则 ddof=0
		ax9_upper = ax9_mean + ax9_std
		ax9_lower = ax9_mean - ax9_std
		ax9.plot(TList, ax9_mean, "r", linewidth=1.5, label='mean')
		ax9.plot(TList, ax9_upper, "b", linewidth=1.5, label=f'mean+std')
		ax9.plot(TList, ax9_lower, "b", linewidth=1.5, label=f'mean-std')


	############################
	# ax1.legend()
	# ax2.legend()
	# ax3.legend()
	# ax4.legend()
	# ax5.legend()
	# ax6.legend()
	# ax7.legend()
	# ax8.legend()
	# ax9.legend()
	#######################################################
	plt.savefig(f'{relaFolderName}_response spectra.jpg', dpi=600, bbox_inches='tight')
	plt.savefig(f'{relaFolderName}_response spectra.eps', dpi=600, bbox_inches='tight')
	##########################################################
	def on_pick(event):
		line=event.artist
		print("clicked motion number is: ",line.get_label())
	fig.canvas.mpl_connect('pick_event', on_pick)

	plt.show()
########################################################################################################################
########################################################################################################################
def timeHistoryPlot(folderName,motionNumber=1):
	"""
	plot acceleration, velocity and displacement time histories for the pulse component, orthogonal component and
	pulse model. time interval is 0.005 s
	----------------------------------------------------------------------------------------------------------------
	folderName(str)-the folder name that sotore generated motions
	motionNumber(int)-the number of motion to be plotted
	"""
	motionNumberList=[]
	files = [f.name for f in Path(folderName).iterdir() if f.is_file()]
	relaFolderName= folderName.split('\\')[-1]
	totalMotion=int(len(files)/4)
	dt=0.005 ###---time interval 0.005 s
	if motionNumber is not None: ###-绘制特定编号地震动
		fig, axes = plt.subplots(3, 3, figsize=(20, 5))
		fig.suptitle(f'Ground motion {motionNumber} time history', fontsize=16)
		###---first column, pulse component
		###---acc(g)
		pulseCompAcc=np.loadtxt(folderName+f"/{motionNumber}_pulseComp.txt").tolist()
		timeHistory=[i1*0.005 for i1 in range(len(pulseCompAcc))]
		axes[0, 0].plot(timeHistory, pulseCompAcc, 'b', linewidth=1)
		axes[0, 0].set_ylabel('Acc(g)', fontsize=10)
		axes[0, 0].set_title('PulseComp', fontsize=11)
		axes[0, 0].grid(True,which='both',linestyle='--',color='grey')
		axes[0,0].set_xlim([0,timeHistory[-1]])
		absMax=max(np.abs(pulseCompAcc))*1.1
		axes[0,0].set_ylim([-absMax,absMax])
		###---vel(cm/s)
		pulseCompVel=accToVelocity(dt,pulseCompAcc)
		timeHistory = [i1 * 0.005 for i1 in range(len(pulseCompVel))]
		axes[1, 0].plot(timeHistory, pulseCompVel, 'b', linewidth=1)
		axes[1, 0].set_ylabel('Vel(cm/s)', fontsize=10)
		# axes[1, 0].set_title('PulseComp Acc', fontsize=11)
		axes[1, 0].grid(True, which='both', linestyle='--', color='grey')
		axes[1, 0].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(pulseCompVel)) * 1.1
		axes[1, 0].set_ylim([-absMax, absMax])
		###---disp(cm)
		pulseCompDisp =velToDisplacement(dt,pulseCompVel)
		timeHistory = [i1 * 0.005 for i1 in range(len(pulseCompDisp))]
		axes[2, 0].plot(timeHistory, pulseCompDisp, 'b', linewidth=1)
		axes[2, 0].set_ylabel('Disp(cm)', fontsize=10)
		axes[2, 0].set_xlabel('Time (s)', fontsize=10)
		# axes[2, 0].set_title('PulseComp Acc', fontsize=11)
		axes[2, 0].grid(True, which='both', linestyle='--', color='grey')
		axes[2, 0].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(pulseCompDisp)) * 1.1
		axes[2, 0].set_ylim([-absMax, absMax])

		###---second column, orthogonal component
		###---acc(g)
		orthogCompAcc = np.loadtxt(folderName + f"/{motionNumber}_orthogComp.txt").tolist()
		timeHistory = [i1 * 0.005 for i1 in range(len(orthogCompAcc))]
		axes[0, 1].plot(timeHistory, orthogCompAcc, 'b', linewidth=1)
		axes[0, 1].set_ylabel('Acc (g)', fontsize=10)
		axes[0, 1].set_title('OrthogComp', fontsize=11)
		axes[0, 1].grid(True, which='both', linestyle='--', color='grey')
		axes[0, 1].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(orthogCompAcc)) * 1.1
		axes[0, 1].set_ylim([-absMax, absMax])
		###---vel(cm/s)
		orthogCompVel=accToVelocity(dt,orthogCompAcc)
		timeHistory = [i1 * 0.005 for i1 in range(len(orthogCompVel))]
		axes[1, 1].plot(timeHistory, orthogCompVel, 'b', linewidth=1)
		axes[1, 1].set_ylabel('Vel (cm/s)', fontsize=10)
		# axes[0, 1].set_title('orthogComp Acc', fontsize=11)
		axes[1, 1].grid(True, which='both', linestyle='--', color='grey')
		axes[1, 1].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(orthogCompVel)) * 1.1
		axes[1, 1].set_ylim([-absMax, absMax])
		###---disp(cm)
		orthogCompDisp = velToDisplacement(dt, orthogCompVel)
		timeHistory = [i1 * 0.005 for i1 in range(len(orthogCompDisp))]
		axes[2, 1].plot(timeHistory, orthogCompDisp, 'b', linewidth=1)
		axes[2, 1].set_ylabel('Disp (cm)', fontsize=10)
		axes[2, 1].set_xlabel('Time (s)', fontsize=10)
		# axes[0, 1].set_title('orthogComp Acc', fontsize=11)
		axes[2, 1].grid(True, which='both', linestyle='--', color='grey')
		axes[2, 1].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(orthogCompDisp)) * 1.1
		axes[2, 1].set_ylim([-absMax, absMax])


		###---third column, pulse model component
		###---acc(g)
		pulseModelAcc = np.loadtxt(folderName + f"/{motionNumber}_pulseModel.txt").tolist()
		timeHistory = [i1 * 0.005 for i1 in range(len(pulseModelAcc))]
		axes[0, 2].plot(timeHistory, pulseModelAcc, 'b', linewidth=1)
		axes[0, 2].set_ylabel('Acc (g)', fontsize=10)
		axes[0, 2].set_title('PulseModel', fontsize=11)
		axes[0, 2].grid(True, which='both', linestyle='--', color='grey')
		axes[0, 2].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(pulseModelAcc)) * 1.1
		axes[0, 2].set_ylim([-absMax, absMax])
		###---vel(cm/s)
		pulseModelVel =accToVelocity(dt,pulseModelAcc)
		timeHistory = [i1 * 0.005 for i1 in range(len(pulseModelVel))]
		axes[1, 2].plot(timeHistory, pulseModelVel, 'b', linewidth=1)
		axes[1, 2].set_ylabel('Vel (cm/s)', fontsize=10)
		# axes[0, 2].set_title('pulseModel Acc', fontsize=11)
		axes[1, 2].grid(True, which='both', linestyle='--', color='grey')
		axes[1, 2].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(pulseModelVel)) * 1.1
		axes[1, 2].set_ylim([-absMax, absMax])
		###---disp(cm)
		pulseModelDisp = velToDisplacement(dt, pulseModelVel)
		timeHistory = [i1 * 0.005 for i1 in range(len(pulseModelDisp))]
		axes[2, 2].plot(timeHistory, pulseModelDisp, 'b', linewidth=1)
		axes[2, 2].set_ylabel('Disp (cm)', fontsize=10)
		axes[2, 2].set_xlabel('Time (s)', fontsize=10)
		# axes[0, 2].set_title('pulseModel Acc', fontsize=11)
		axes[2, 2].grid(True, which='both', linestyle='--', color='grey')
		axes[2, 2].set_xlim([0, timeHistory[-1]])
		absMax = max(np.abs(pulseModelDisp)) * 1.1
		axes[2, 2].set_ylim([-absMax, absMax])
	plt.savefig(f'{relaFolderName}_ground motion {motionNumber} time history.jpg',dpi=600,bbox_inches='tight')
	plt.savefig(f'{relaFolderName}_ground motion {motionNumber} time history.eps', dpi=600, bbox_inches='tight')
	plt.show()
########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
	#################################################################################
	#######################------################################
	i11=0
	j11=1
	while True:
		print(f"第{i11+1}次模拟!")
		moment=np.random.uniform(6,7.9)
		Rrup=np.random.uniform(0.1,30)
		s_or_d=np.random.uniform(1.2,135)
		theta_or_phi=np.random.uniform(0,90)
		Vs30=np.random.uniform(700,2000)


		faultType="strikeSlip" ##Including: ["strikeSlip","reverseAndReverseOblique"]
		simulationType="onlyPulse"  ##Including: ["pulseAndNoPulse","onlyPulse","onlyNoPulse"]
		momentMag=moment
		zTOR=0
		Rrup=Rrup
		Vs30=Vs30
		s_or_d=s_or_d
		theta_or_phi=theta_or_phi
		numSimMotions=1
		nearFaultMotionSimInstance=NearFaultMotionSim(faultType,simulationType,momentMag,zTOR,Rrup,Vs30,s_or_d,
		theta_or_phi,numSimMotions)
		pulseComponnetsList,noPulseComponnetsList=nearFaultMotionSimInstance.runSimulation()
		if len(pulseComponnetsList)>0:
			np.savetxt("Vp=Tp=5s/"+str(j11)+".txt",np.asmatrix(pulseComponnetsList[0][0]).T,fmt="%.6f")
			j11+=1
		i11+=1
	###############################################
	################---反应谱绘制
	
	###################
	f = [np.exp(i0 * 0.035 - 3) for i0 in range(200)]  ##frequencies(Hz)0.05~50Hz
	TList = [1.0 / each1 for each1 in f]
	TList.reverse()
	beta = 0.05
	dt = 0.005
	# acc=np.loadtxt("maxVel_acceleration_interpolation/9.txt")
	# targetSa, targetSv, targetSd = ResponseCal(beta, dt, np.array(TList), np.array(acc))
	#################
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax2.set_xscale("log")
	ax2.set_yscale("log")
	ax3.set_xscale("log")
	ax3.set_yscale("log")
	ax4.set_xscale("log")
	ax4.set_yscale("log")
	ax1.grid(which="both")
	ax2.grid(which="both")
	ax3.grid(which="both")
	ax4.grid(which="both")
	ax1.set_ylim(10 ** (-2), 5*10 ** 0)
	ax2.set_ylim(10 ** (-2), 5* 10 ** 2)
	ax3.set_ylim(5*10 ** (-4), 5 * 10 ** 2)
	ax1.set_xlim(0.01896, 20)
	ax2.set_xlim(0.01896, 20)
	ax3.set_xlim(0.01896, 20)
	############
	for i1 in range(len(pulseComponnetsList)):
		SaArray1, SvArray1, SdArray1 = ResponseCal(beta, dt, np.array(TList), np.array(pulseComponnetsList[i1][0]))
		SaArray2, SvArray2, SdArray2 = ResponseCal(beta, dt, np.array(TList), np.array(pulseComponnetsList[i1][1]))
		SaGeomAve=[(each1*each2)**0.5 for each1,each2 in zip(SaArray1,SaArray2)]
		SvGeomAve = [(each1*each2)**0.5 for each1, each2 in zip(SvArray1, SvArray2)]
		SdGeomAve = [(each1*each2)**0.5 for each1, each2 in zip(SdArray1, SdArray2)]
		ax1.plot(TList, SaGeomAve, "k")
		ax2.plot(TList, SvGeomAve, "k")
		ax3.plot(TList,SdGeomAve, "k")
	for i1 in range(len(noPulseComponnetsList)):
		SaArray1, SvArray1, SdArray1 = ResponseCal(beta, dt, np.array(TList), np.array(noPulseComponnetsList[i1][0]))
		SaArray2, SvArray2, SdArray2 = ResponseCal(beta, dt, np.array(TList), np.array(noPulseComponnetsList[i1][1]))
		SaGeomAve=[(each1*each2)**0.5 for each1,each2 in zip(SaArray1,SaArray2)]
		SvGeomAve = [(each1*each2)**0.5 for each1, each2 in zip(SvArray1, SvArray2)]
		SdGeomAve = [(each1*each2)**0.5 for each1, each2 in zip(SdArray1, SdArray2)]
		ax1.plot(TList, SaGeomAve, "k")
		ax2.plot(TList, SvGeomAve, "k")
		ax3.plot(TList,SdGeomAve, "k")
	# ax1.plot(TList, targetSa, "r")
	# ax2.plot(TList, targetSv, "r")
	# ax3.plot(TList, targetSd, "r")
	plt.show()
	################################################################################################################
	# ################################################---时程绘制
	# pulseComp=np.loadtxt("pulse_acc.txt")
	# totalMotion=np.loadtxt("pulseMotion_acc.txt")
	#
	# fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
	# ax1.grid(which="both")
	# ax2.grid(which="both")
	# ax3.grid(which="both")
	# ax4.grid(which="both")
	# ax5.grid(which="both")
	# ax6.grid(which="both")
	# dt = 0.005
	# times=[j1*dt for j1 in range(len(pulseComp))]
	# acc1=pulseComp
	# acc2 = totalMotion
	# vel1=accToVelocity(dt,acc1)
	# vel2=accToVelocity(dt,acc2)
	# disp1=velToDisplacement(dt,vel1)
	# disp2=velToDisplacement(dt,vel2)
	# ax1.plot(times,acc1)
	# ax2.plot(times,acc2)
	# ax3.plot(times, vel1)
	# ax4.plot(times, vel2)
	# ax5.plot(times, disp1)
	# ax6.plot(times, disp2)
	# plt.show()
	#