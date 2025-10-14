#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
#####Units: Length-mm, Force-N, mass-ton, Stress-Mpa, g=9810mm/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2022-01-29
########################################################################################################################
########################---import modules---#################################
import h5py
import os
import numpy as np
import fnmatch
import datetime
########################################################################################################################
class h5pyDB(object):
    """Save the data to h5py database"""
    def __init__(self,dbPath):
        """
        Initialize the class
        Inputs:
            dbPath(str)-the path of the database
        """
        self._dbPath = dbPath
        self.values=0
        self.resultFileName=dbPath

    @classmethod
    def initDB(self,dbPath):
        """Initialize the database"""
        if os.path.exists(dbPath):
            os.remove(dbPath)
        with h5py.File(dbPath, 'w') as f:
            pass  ###---create an empty h5 database

    def saveResult(self,dataName,resultsList,headNameList,operationIndexStr='replace'):
        """
        ---A general save template for hdf5 database ---
        Save results to database, resultsList=[[result0_0,result0_1,],[],...]
        headNameList=[headName_0,headName_1,...]
        dataName(str)
        operationIndexStr='replace' or 'append' 'replace' means replace the data, 'append' means append data after the last row
        """
        if len(resultsList)>0:
            list0 = resultsList[0]
            saveTypeList = []
            typeDict = {"int": "np.int32", "float": "np.float32", "str": "h5py.string_dtype(encoding='utf-8')"}
            saveTypeList=[typeDict["int"] if isinstance(eachValue, (int,np.int64,np.uint64)) else
                typeDict["float"] if isinstance(eachValue, (float, np.float64)) else
                typeDict["str"] if isinstance(eachValue, str) else None  for eachValue in list0]

            dtypeStr = "np.dtype(["
            dtypeStr += ''.join([f"('{headNameList[i1]}',{saveTypeList[i1]})," for i1 in range(len(list0) - 1)])
            dtypeStr += f"('{headNameList[-1]}',{saveTypeList[-1]})])"
            dtype = eval(dtypeStr)
            structured_data = np.zeros(len(resultsList), dtype=dtype)
            for i2 in range(len(headNameList)):
                structured_data[headNameList[i2]] = [each[i2] for each in resultsList]

            with h5py.File(self.resultFileName, 'a',libver='latest') as f:
                dataset = f.get(dataName)
                if dataset is None:
                    dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype,chunks=True)
                    new_size = len(resultsList)
                    dataset.resize(new_size, axis=0)
                    dataset[0:new_size] = structured_data
                else:
                    if operationIndexStr == 'replace':
                        del f[dataName]
                        dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype,chunks=True)
                        new_size = len(resultsList)
                        dataset.resize(new_size, axis=0)
                        dataset[0:new_size] = structured_data
                    elif operationIndexStr == 'append':
                        original_size = dataset.shape[0]
                        new_size = original_size + len(resultsList)
                        dataset.resize(new_size, axis=0)
                        dataset[original_size:new_size] = structured_data

########################################################################################################################
########################################################################################################################
# if __name__ == '__main__':
#     dbPath = "resultsDB_1_114.h5"
#     dbInstance=h5pyDB(dbPath)
#     # dbInstance.initDB(dbPath)
#     # nodeList=[[4,1.1,1.2,3],[6,3,4.3,4],[7,3,4.3,4]]
#     # dbInstance.saveNodes("nodes",nodeList)
#     # dbInstance.getNodes('nodes')
#     # dbInstance.saveEles('zeroEle',nodeList)
#     # dbInstance.getEles('zeroEle')
#     modeList=[[1,1,2,3,4,5.5,6],[2,1,2,3,4,5.5,6],[3,1,2,3,4,5.5,6]]
#     # dbInstance.saveModes('mdoe-2',modeList)
#     # dbInstance.getModes('mdoe-3')
#     periodList=[[1,0.3],[2,0.6]]
#     # dbInstance.savePeriod(periodList)
#     # dbInstance.getPeriod()
#     geomfList=[[1,1.0,0.2,0.3],[2,0.7,0.3,0.5]]
#     # dbInstance.saveGeomTransf('ele1Geomf',geomfList)
#     # dbInstance.getGeomTransf("ele1Geomff")
#     localEleCoord=[[1,2,4],[3,4,5]]
#     localEleCoordSpecial=[[1,2,3.4,2.3,4.5,3.4,5.6,7.6],[2,2,3.4,2.3,4.5,3.4,5.6,7.6]]
#     # dbInstance.saveEleLocalCoordSys('1_realEle',localEleCoord)
#     # dbInstance.saveEleLocalCoordSys('1_specialEle', localEleCoordSpecial)
#     # dbInstance.getEleLocalCoordSys("1_realEle")
#     nodeTimeHisList=[[0.0,1.2,3.3,4.2],[0.1,1.2,3.3,4.2]]
#     # dbInstance.saveNodeTimeHistory('node_disp_1',nodeTimeHisList)
#     # dbInstance.getNodeTimeHistory(1,'disp',3)
#     trussTimeList=[[0.1,2.3],[0.3,0.6]]
#     # dbInstance.saveTrussEleResponseTimeHistory('trussEle_axialForce_1',trussTimeList)
#     # dbInstance.getTrussEleResponseTimeHistory(1,'axialForce')
#     zeroTimeList=[[0.1,3.2,3.3,4.0],[0.2,3.3,7.3,5.0]]
#     dirList=[1,4,3]
#     # dbInstance.saveZeroEleResponseTimeHistory('zeroEle_localForce_11',zeroTimeList,dirList)
#     # dbInstance.getZeroEleResponseTimeHistory(11,'localForce',4)
#     nonEleSectionList=[[1,0.2,0.4,0.6,0.8],[2,0.6,0.8,0.9,1.0]]
#     # dbInstance.saveNonEleSectResponseTimeHistory('nonEle_sectionForce_13',nonEleSectionList)
#     # dbInstance.getNonEleSectResponseTimeHistory(13,'sectionForce',4)
#     nonZeroEleList=[[1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.10,0.11,0.12],
#                     [2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.10,0.11,0.12]]
#     # dbInstance.saveNonZeroEleResponseTimeHistory('nonZeroEle_localForce_11',nonZeroEleList)
#     # dbInstance.getNonZeroEleResponseTimeHistory(11,'localForce','J',1)
#
#
