#-*-coding: UTF-8-*-
########################################################################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 28/08/2025
#  Environemet: Successfully excucted in python 3.13
########################################################################################################################
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import numpy as np
import numba
from matplotlib.lines import lineStyles
from numba import njit, prange
import time
########################################################################################################################
class DisplacementHistory():
    """
    A class for generating several types of displacement history
    """
    def __init__(self):
        """
        """
        pass

    def _plotDispHistory(self,timeList:list,dispList:list):
        """
        ----------------------------------------------------------------------------------------------------------------
        This method plots the displacement history versus pseudo time history
        ----------------------------------------------------------------------------------------------------------------
        :param timeList(list): pseudo time history list
        :param dispList(list): displacement history list
        :return:
        """
        fig,ax=plt.subplots(figsize=(6,1),dpi=600,constrained_layout=True)
        ax.plot(timeList,dispList,linewidth=0.6,color='b')
        ax.set_xlabel('Pseudo time',fontsize=8)
        ax.set_ylabel('Displacement',fontsize=8)
        ax.minorticks_on()
        ax.grid(which='major',linestyle='--',color='gray',linewidth=0.15)
        ax.grid(which='minor', linestyle='--', color='gray', linewidth=0.15)
        ax.tick_params(labelsize=8)
        ax.set_xlim(0,1)
        ax.set_ylim(1.1*min(dispList),1.1*max(dispList))
        fig.savefig('displacementHistory.png')
        fig.savefig('displacementHistory.eps')
        plt.show()


    def monotonicDispHistory(self,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for monotonic loading
        ----------------------------------------------------------------------------------------------------------------
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        peakValue=maxDisp*scaleFactor
        incrValue=peakValue/numDivide
        dispHistoryList=[i1*incrValue for i1 in range(numDivide+1)]
        timeHistoryList=[i2/numDivide for i2 in range(numDivide+1)]
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        return dispHistoryList

    def cyclicOneSideConstPeak(self,numCycle:int=1,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,
                               plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for constant peak values on positive or negtive side
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        # startTime=time.time()
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        dispHistoryList=[incrValue*i1 for i1 in range(numDivide)]+[incrValue*j1 for j1 in range(numDivide,0,-1)]
        dispHistoryList=dispHistoryList*numCycle
        dispHistoryList.append(0)
        numPts=len(dispHistoryList)
        timeHistoryList=[j2/(numPts-1) for j2 in range(numPts)]
        # endTime = time.time()
        # print("total time:", endTime - startTime)
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        # endTime=time.time()
        # print("total time:",endTime-startTime)
        return dispHistoryList

    def cyclicOneSideLinearIncrPeak(self,numCycle:int=1,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,
                               plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for linear incremental peak on positive or negtive side
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        peakValuesList=[(i1+1)*peakValue/numCycle for i1 in range(numCycle)]
        dispHistoryList=[]
        for i2 in range(numCycle):
            numPts=int(peakValuesList[i2]/incrValue)
            incr=peakValuesList[i2]/numPts
            [dispHistoryList.append(i3*incr) for i3 in range(numPts)]
            [dispHistoryList.append(i4*incr) for i4 in range(numPts,0,-1)]
        dispHistoryList.append(0)
        numPts = len(dispHistoryList)
        timeHistoryList = [j2 / (numPts - 1) for j2 in range(numPts)]
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        return dispHistoryList

    def cyclicOneSideLinearIncrPeak_nConst(self,numCycle:int=1,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,
                                           numLocalConst:int=2,plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for linear incremental peak on positive or negtive side, and has
        n local const peak
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param numLocalConst(int): number of local constant peaks, default=2
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        peakValuesList=[(i1+1)*peakValue/numCycle for i1 in range(numCycle)]
        dispHistoryList=[]
        for i2 in range(numCycle):
            localDispList=[]
            numPts=int(peakValuesList[i2]/incrValue)
            incr=peakValuesList[i2]/numPts
            [localDispList.append(i3*incr) for i3 in range(numPts)]
            [localDispList.append(i4*incr) for i4 in range(numPts,0,-1)]
            dispHistoryList+=localDispList*numLocalConst
        dispHistoryList.append(0)
        numPts = len(dispHistoryList)
        timeHistoryList = [j2 / (numPts - 1) for j2 in range(numPts)]
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        return dispHistoryList

    def cyclicTwoSideConstPeak(self, numCycle: int = 1, numDivide: int = 100, maxDisp: float = 1,
                               scaleFactor: float = 1,plotDispHistory: bool = False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for constant peak values on both positive and negtive sides
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        # startTime=time.time()
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        # print("peakValue:",peakValue)
        # print("incrValue:",incrValue)
        dispHistoryList = ([incrValue * i1 for i1 in range(numDivide)] +
        [incrValue * j1 for j1 in range(numDivide, -numDivide, -1)]+[incrValue * i3 for i3 in range(-numDivide,0,1)])
        dispHistoryList = dispHistoryList * numCycle
        dispHistoryList.append(0)
        numPts = len(dispHistoryList)
        timeHistoryList = [j2 / (numPts - 1) for j2 in range(numPts)]
        # endTime = time.time()
        # print("total time:", endTime - startTime)
        if plotDispHistory == True:
            self._plotDispHistory(timeHistoryList, dispHistoryList)
        # endTime=time.time()
        # print("total time:",endTime-startTime)
        return dispHistoryList

    def cyclicTwoSideLinearIncrPeak(self,numCycle:int=1,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,
                               plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for linear incremental peak on both positive and negative sides
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        peakValuesList=[(i1+1)*peakValue/numCycle for i1 in range(numCycle)]
        dispHistoryList=[]
        for i2 in range(numCycle):
            numPts=int(peakValuesList[i2]/incrValue)
            incr=peakValuesList[i2]/numPts
            [dispHistoryList.append(i3*incr) for i3 in range(numPts)]
            [dispHistoryList.append(i4*incr) for i4 in range(numPts,-numPts,-1)]
            [dispHistoryList.append(i5*incr) for i5 in range(-numPts,0,1)]
        dispHistoryList.append(0)
        numPts = len(dispHistoryList)
        timeHistoryList = [j2 / (numPts - 1) for j2 in range(numPts)]
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        return dispHistoryList

    def cyclicTwoSideLinearIncrPeak_nConst(self,numCycle:int=1,numDivide:int=100,maxDisp:float=1,scaleFactor:float=1,
                                           numLocalConst:int=2,plotDispHistory:bool=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        This class generates a displacement history for linear incremental peak on both positive and negative sides,
        and has n local const peak
        ----------------------------------------------------------------------------------------------------------------
        :param numCycle(int): number of cycles for the disp history, default=1
        :param numDivide(int):number of divisions for the maximum displacement, default=100
        :param maxDisp(float):maximum displacement, positive means positive tangent for staring line,default=1
        :param scaleFactor(float):scale factor for the maximum displacement,default=1
        :param numLocalConst(int): number of local constant peaks, default=2
        :param plotDispHistory(bool): whether to plot displacement history, default=False
        :return:displacement history
        """
        peakValue = maxDisp * scaleFactor
        incrValue = peakValue / numDivide
        peakValuesList=[(i1+1)*peakValue/numCycle for i1 in range(numCycle)]
        dispHistoryList=[]
        for i2 in range(numCycle):
            localDispList=[]
            numPts=int(peakValuesList[i2]/incrValue)
            incr=peakValuesList[i2]/numPts
            [localDispList.append(i3 * incr) for i3 in range(numPts)]
            [localDispList.append(i4 * incr) for i4 in range(numPts, -numPts, -1)]
            [localDispList.append(i5 * incr) for i5 in range(-numPts, 0, 1)]
            dispHistoryList+=localDispList*numLocalConst
        dispHistoryList.append(0)
        numPts = len(dispHistoryList)
        timeHistoryList = [j2 / (numPts - 1) for j2 in range(numPts)]
        if plotDispHistory==True:
            self._plotDispHistory(timeHistoryList,dispHistoryList)
        return dispHistoryList

########################################################################################################################
########################################################################################################################
if __name__ == "__main__":
    dispInstance=DisplacementHistory()
    #######################---monotonic displacement history
    # monotonicDispHistory=dispInstance.monotonicDispHistory(numDivide=100,maxDisp=1,scaleFactor=1,plotDispHistory=True)
    #######################---cyclicOneSideConstPeak displacement history
    # cyclicOneSideConstPeakIns=dispInstance.cyclicOneSideConstPeak(numCycle=10,numDivide=100,maxDisp=-0.3,scaleFactor=2,
    #                            plotDispHistory=True)
    #######################---cyclicOneSideConstPeak displacement history
    # cyclicOneSideConstPeakIns = dispInstance.cyclicOneSideLinearIncrPeak(numCycle=10, numDivide=15, maxDisp=-0.5,
    #                                                                 scaleFactor=2,plotDispHistory=True)
    #######################---cyclicOneSideLinearIncrPeak_nConst displacement history
    # dispInstance.cyclicOneSideLinearIncrPeak_nConst(numCycle= 10, numDivide= 100, maxDisp= -1,scaleFactor= 1,
    #                                                 numLocalConst= 3, plotDispHistory= True)
    #######################---cyclicTwoSideConstPeak displacement history
    # dispInstance.cyclicTwoSideConstPeak(numCycle= 2, numDivide= 10, maxDisp= 1,scaleFactor= 1.5,plotDispHistory= True)
    #######################---cyclicTwoSideConstPeak displacement history
    # dispInstance.cyclicTwoSideLinearIncrPeak(numCycle=10,numDivide=10,maxDisp=-1,scaleFactor=1,plotDispHistory=True)
    #######################---cyclicTwoSideConstPeak displacement history
    dispInstance.cyclicTwoSideLinearIncrPeak_nConst(numCycle=5,numDivide=10,maxDisp=1,scaleFactor=1,
                                           numLocalConst=2,plotDispHistory=True)

