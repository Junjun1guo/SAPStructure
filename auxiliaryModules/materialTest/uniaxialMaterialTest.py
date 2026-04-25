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
import openseespy.opensees as ops
import re
########################################################################################################################
def uniaxialMatTest(matList,disp,saveFig,figureName):
    """
    A function for performing uniaxial material test.
    """
    for each in matList:
        match = re.search(r'\b\d+\b',each)
        first_number=None
        if match:
            first_number = match.group()
        else:
            print("mateiral tag should be integer number!")
    ######################################################################################
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)  ###for 3D model
    ops.node(1, 0)
    ops.node(2, 0)
    ops.fix(1, 1)
    for each in matList:
        eval(f"ops.{each}")
    ops.uniaxialMaterial('Elastic', 11111099,1e-10)
    ops.uniaxialMaterial('Parallel',21111099,eval(first_number),11111099)
    ###############################################################
    # ops.equalDOF(1, 2, 2, 3, 4, 5, 6)
    ops.element('zeroLength', 1, 1, 2, '-mat',21111099, '-dir', 1)
    ###############################################################
    ops.loadConst('-time', 0.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1)
    ###############################################################
    ops.system('UmfPack')
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.test('NormDispIncr', 1.0e-1, 2000)
    ops.algorithm('KrylovNewton')
    ops.analysis('Static')
    ###############################################################
    resultList = []
    dispIncrList=[0]+[disp[i0+1]-disp[i0] for i0 in range(len(disp)-1)]
    ################################################################
    for i1 in range(len(dispIncrList)):
        ops.integrator('DisplacementControl', 2, 1, dispIncrList[i1])
        ops.analyze(1)
        dispValues = round(ops.nodeDisp(2, 1), 8)
        # print(dispValues)
        forceValues = round(ops.eleForce(1, 1) * -1, 8)
        resultList.append([dispValues, forceValues])
    dispList = [each[0] for each in resultList]
    forceList = [each[1] for each in resultList]
    ################################################################
    fig, ax = plt.subplots(figsize=(4, 2.5), dpi=600, constrained_layout=True)
    ax.plot(dispList,forceList, linewidth=0.4, color='b')
    ax.set_xlabel('displacement', fontsize=8)
    ax.set_ylabel('force', fontsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    # ax.minorticks_on()
    majorGrid=ax.grid(which='major', linestyle='--', color='gray', linewidth=0.15)
    # ax.grid(which='minor', linestyle='--', color='gray', linewidth=0.15)
    ax.tick_params(labelsize=8)
    ax.set_xlim(1.1*min(dispList), 1.1*max(dispList))
    ax.set_ylim(1.1 * min(forceList), 1.1 * max(forceList))
    if saveFig:
        fig.savefig(f'{figureName}.png')
        fig.savefig(f'{figureName}.eps')
    plt.show()
    ##########################################################################
    return dispList,forceList
########################################################################################################################
########################################################################################################################
if __name__ == "__main__":
    from displacementHistory import DisplacementHistory
    dispInstance=DisplacementHistory()
    # dispHistory=dispInstance.monotonicDispHistory(numDivide=100,maxDisp=-0.003,scaleFactor=1,plotDispHistory=True)
    dispHistory=dispInstance.cyclicTwoSideLinearIncrPeak(numCycle=10, numDivide=20, maxDisp=-0.004, scaleFactor=1, plotDispHistory=False)
    frictionCoeff = 0.05  #####-支座摩擦副摩擦系取为0.05
    shearKeyLongK = 1650 / 0.002  ####-取位移2mm发生剪断
    ops.uniaxialMaterial('Elastic', 101, shearKeyLongK)
    ops.uniaxialMaterial('MinMax', 102, 101, '-min', -0.002, '-max', 0.002)  ###-位移超过2mm让剪力键失效
    materialList=[f"uniaxialMaterial('Elastic', 101, {shearKeyLongK})",
                  f"uniaxialMaterial('MinMax', 102, 101, '-min', -0.002, '-max', 0.002) "]
    uniaxialMatTest(matList=materialList,disp=dispHistory,saveFig=False,figureName="hystereticCurve")