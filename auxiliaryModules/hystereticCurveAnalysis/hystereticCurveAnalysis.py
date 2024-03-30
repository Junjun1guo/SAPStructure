#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
#####Units: Length-mm, Force-N, mass-ton, Stress-Mpa, g=9810mm/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: guojj@tongji.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2024-03-22
########################################################################################################################
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import openseespy.opensees as ops
########################################################################################################################
########################################################################################################################
class HystereticCurveAnalysis():
	"""
	--------------------------------------------------------------------------------------------------------------------
	A class for hysteretic curve analyses (version:0.6.0)
	Environemet: Successfully executed in python 3.11
	Date: 2024-03-22
	--------------------------------------------------------------------------------------------------------------------
	    ** ********************************************************************************  **
	    ** (C) Copyright 2024, School of Civil Engineering, Beijing Jiaotong University      **
	    ** All Rights Reserved.                                                              **
	    **                                                                                   **
	    ** Commercial use of this program is strictly prohibited.                            **
	    **                                                                                   **
	    ** Developed by:                                                                     **
	    **   Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo           **
	    **   jjguo2@bjtu.edu.cn/guojj_ce@163.com                                             **
	    ** ********************************************************************************  **
	"""
	def __init__(self,xDataList,yDataList):
		"""
        ----------------------------------------------------------------------------------------------------------------
        Initialize the class
        ------------------------------------------
        Inputs:
        xDataList,yDataList(list)-The x and y values of the hysteretic curve data
        ----------------------------------------------------------------------------------------------------------------
		"""
		self.xValueList,self.yValueList=xDataList,yDataList
		self.loopInfoDict={} ###---loopNum:[startIndex,endIndex]
		self.tensionOnly=True

	def plotHystereticCurve(self,saveFig=False,xlabel="x",ylabel='y',title='Hysteretic curve',multiColors=False):
		"""
        Plot the original hysteretic curve
        ------------------------------------------
        Inputs:
        saveFig(bool)-Save the plot to eps figure, default not saved
        xDataList,yDataList(list)-The x and y values of the hysteretic curve data
        multiColors(bool)-Whether distinguish different loops with different colors, default value is False
		"""
		ax=plt.subplot(1,1,1)
		if not multiColors:
			plt.plot(self.xValueList,self.yValueList, ls='-', color='b', lw=1.5)
		else:
			self._hystereticLoopDecomp()
			legendList=list(self.loopInfoDict.keys())
			legendList=[str(each) for each in legendList]
			for key,value in self.loopInfoDict.items():
				start,end=value
				plt.plot(self.xValueList[start:end],self.yValueList[start:end],ls='-',lw=1.5)
			ax.legend(labels=legendList)
		plt.xlabel(f'{xlabel}', fontsize=15)
		plt.ylabel(f'{ylabel}', fontsize=15)
		plt.title(f'{title}', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		if saveFig:
			plt.savefig(title+".eps")
			plt.savefig(title + ".jpg")
		plt.show()

	def skeletonCurve(self,saveData=False,saveFig=False):
		"""
		Get and plot the skeleton curve of the hysteretic data
        ------------------------------------------
        Inputs:
        saveData(bool)-Whether save the skeleton curve data, default value is False
        saveFig(bool)-Whether save the skeleton curve figure, default value is False
		"""
		posIndex,negIndex=[],[]
		self._hystereticLoopDecomp()
		for eachI in self.loopInfoDict.values():
			posValue=max(self.yValueList[eachI[0]:eachI[1]])
			indexValue=[index for index,value in enumerate(self.yValueList[eachI[0]:eachI[1]])
						if value==posValue][0]+eachI[0]
			posIndex.append(indexValue)
		if not self.tensionOnly:
			for eachI in self.loopInfoDict.values():
				negValue = min(self.yValueList[eachI[0]:eachI[1]])
				indexValue = [index for index, value in enumerate(self.yValueList[eachI[0]:eachI[1]])
							  if value == posValue][0] + eachI[0]
				negIndex.append(indexValue)
		if len(negIndex)>0:
			negIndex.reverse()
		plt.subplot(1, 1, 1)
		plt.plot(self.xValueList, self.yValueList, ls='-', color='b', lw=1)
		plt.xlabel(f'x', fontsize=15)
		plt.ylabel(f'y', fontsize=15)
		plt.title(f'Hysteretic curve', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		negX, negY = [], []
		if len(negIndex)>0:
			[[negX.append(self.xValueList[each]),negY.append(self.yValueList[each])] for each in negIndex]
			for item in range(len(self.xValueList)):
				indexI=negIndex[-1]+item
				if self.yValueList[indexI]>0:
					break
				else:
					negX.append(self.xValueList[indexI])
					negY.append(self.yValueList[indexI])
			plt.plot(negX, negY, ls='-', color='r', lw=1.5)
		posX,posY=[],[]
		for i1 in range(posIndex[0]):
			posX.append(self.xValueList[i1])
			posY.append(self.yValueList[i1])
		for i2 in range(len(posIndex)-1):
			posX.append(self.xValueList[posIndex[i2+1]])
			posY.append(self.yValueList[posIndex[i2+1]])
		plt.plot(posX, posY, ls='-', color='r', lw=1.5)

		if saveData:
			xData=negX+posX
			yData=negY+posY
			saveDataValue=[]
			for eachx,eachy in zip(xData,yData):
				saveDataValue.append([eachx,eachy])
			np.savetxt("skeletonData.txt",saveDataValue, fmt="%.6f %.6f")
		if saveFig:
			plt.savefig("skeleton.eps")
			plt.savefig("skeleton.jpg")
		plt.show()

	def _hystereticLoopDecomp(self):
		"""
		Decompositon of the hysteretic curves
        ------------------------------------------
        Inputs:
		"""
		loopInfoList=[]
		startIndex=0
		xValue,yValue=self.xValueList,self.yValueList
		xValue[0]=0.0
		if min(xValue)<-0.01*max(xValue): ###---for tensile and compresive loops
			self.tensionOnly=False
			for i1 in range(startIndex, len(xValue)-1):
				if ((xValue[i1+1]-xValue[i1])>0)and(yValue[i1]*yValue[i1+1]<=0):
					loopInfoList.append(i1)
		else:###---for tension only loops
			for i1 in range(1,len(xValue)-1):
				if (xValue[i1]-xValue[i1-1]<0)and(xValue[i1+1]-xValue[i1]>0):
					loopInfoList.append(i1)
		loopInfoList=[0]+loopInfoList+[len(xValue)-1]
		for num in range(len(loopInfoList)-1):
			self.loopInfoDict[num+1]=[loopInfoList[num],loopInfoList[num+1]]

	def plotLoop(self,loopNumber,saveData=False,saveFig=False,dottedLine=True):
		"""
		Plot each hysteretic loop
        ------------------------------------------
        Inputs:
        loopNumber(int)-The hysteretic loop number
        saveData(bool)-Whether save the hysteretic loop data, default value is False
        saveFig(bool)-Whether save the hysteretic loop figure, default value is False
		"""
		self._hystereticLoopDecomp()
		start,end=self.loopInfoDict[loopNumber]
		plt.subplot(1, 1, 1)
		if dottedLine:
			plt.plot(self.xValueList[start:end], self.yValueList[start:end], '.', color='b',markersize=5)
		else:
			plt.plot(self.xValueList[start:end], self.yValueList[start:end], ls='-', color='b', lw=1)
		plt.xlabel(f'x', fontsize=15)
		plt.ylabel(f'y', fontsize=15)
		plt.title(f'loop-{loopNumber}', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		if saveFig:
			plt.savefig(f"loop-{loopNumber}.eps")
			plt.savefig(f"loop-{loopNumber}.jpg")
		if saveData:
			saveDataValue=[]
			for eachx,eachy in zip(self.xValueList[start:end], self.yValueList[start:end]):
				saveDataValue.append([eachx,eachy])
			np.savetxt(f"loop-{loopNumber}.txt",saveDataValue, fmt="%.6f %.6f")
		plt.show()

########################################################################################################################
########################################################################################################################
if __name__=="__main__":
	testData = np.loadtxt('testSMACableData.txt')
	testDisp = list(testData[:, 0])
	testForce = list(testData[:, 1] * 1000)
	# plt.plot(testDisp, testForce)
	# plt.show()
	instance=HystereticCurveAnalysis(xDataList=testDisp,yDataList=testForce)
	instance.plotHystereticCurve(saveFig=False,multiColors=False)
	# instance.skeletonCurve(saveData=True,saveFig=False)
	# instance.plotLoop(loopNumber=2,saveData=False,saveFig=False)












