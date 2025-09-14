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
########################################################################################################################
class h5pyDB(object):
    """Save the data to sqlite database"""
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

    def saveNodes(self,nodesSaveName,nodeList):
        """Save nodes to database, [[nodeTag,xCoord,yCoord,zCoord...],[],...]"""
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('modelInfo/'+nodesSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('modelInfo/'+nodesSaveName, [1, 4], maxshape=[None,4],
                                           chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['nodeTag:int','xCoord:float','yCoord:float','zCoord:float']
                newShape=(len(nodeList),4)
                dataSet.resize(newShape)
                dataSet[:]=nodeList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(nodeList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=nodeList

    def getNodes(self,saveNodeName):
        """
        return nodes from database
        saveNodeName(str)-the table name of saved nodes
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('modelInfo/'+saveNodeName)
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet
            else:
                print(f'''table {saveNodeName} doesn't exitst!''')
                return None

    def saveEles(self,elesSaveName,elesList):
        """Save nodes to database, [[eleTag,nodeI,nodeJ,...,'1D'],[],...]"""
        eleDimStr=elesList[0][-1]
        elesList=[each[:-1] for each in elesList]
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('modelInfo/'+elesSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('modelInfo/'+elesSaveName, [1, len(elesList[0])],
                                           maxshape=[None,len(elesList[0])],
                                           chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['eleTag:int','nodeI:int','nodeJ:int','...']
                dataSet.attrs['eleDim']=eleDimStr
                newShape=(len(elesList),len(elesList[0]))
                dataSet.resize(newShape)
                dataSet[:]=elesList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(elesList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=elesList

    def getEles(self,saveElesName):
        """
        return elements from database
        saveElesName(str)-the table name of saved eles
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('modelInfo/'+saveElesName)
            eleDimStr=dataSet.attrs['eleDim']
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet,eleDimStr
            else:
                print(f'''table {saveElesName} doesn't exitst!''')
                return None,None

    def saveModes(self,modesName,modesList):
        """Save modes to database, [[nodeTag,[dof1value,dof2value],...],[],...]"""
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('modalInfo/'+modesName)
            if dataSet is None:
                dataSet = f.create_dataset('modalInfo/'+modesName, [1, len(modesList[0])],
                                           maxshape=[None,len(modesList[0])],
                                           chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['nodeTag:int','dof1Value:float','dof2Value:float','...']
                newShape=(len(modesList),len(modesList[0]))
                dataSet.resize(newShape)
                dataSet[:]=modesList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(modesList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=modesList


    def getModes(self,saveModesName):
        """
        return modes from database
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('modalInfo/'+saveModesName)
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet
            else:
                print(f'''table {saveModesName} doesn't exitst!''')
                return None

    def savePeriod(self,periodList):
        """Save periods to database, [[periodNum,value],[],...]"""
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('periodInfo/period')
            if dataSet is None:
                dataSet = f.create_dataset('periodInfo/period', [1,2],maxshape=[None,2],
                                           chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['periodNum:int','periodValue:float']
                newShape=(len(periodList),len(periodList[0]))
                dataSet.resize(newShape)
                dataSet[:]=periodList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(periodList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=periodList

    def getPeriod(self):
        """
        return periods from database
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('periodInfo/period')
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet
            else:
                print(f'''table period doesn't exitst!''')
                return None

    def saveGeomTransf(self,geomTransfSaveName,geomfList):
        """
        Save the geomTransf to Database
        geomTransfSaveName(str)-the name of the saved table
        geomfList(list)-[[geomfTag1,localZX_1,localZY_1,localZZ_1],[geomfTag2,localZX_2,localZY_2,localZZ_2],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('geomTransf/'+geomTransfSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('geomTransf/'+geomTransfSaveName, [1,len(geomfList[0])],
                maxshape=[None,len(geomfList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['geomfTag:int','localZX_1','localZY_1','...']
                newShape=(len(geomfList),len(geomfList[0]))
                dataSet.resize(newShape)
                dataSet[:]=geomfList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(geomfList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=geomfList

    def getGeomTransf(self, saveGeomTransfName):
        """
        return geomTransf from database
        saveGeomTransfName(str)-the table name of saved geomTransf
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('geomTransf/'+saveGeomTransfName)
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet
            else:
                print(f'''table {saveGeomTransfName} doesn't exitst!''')
                return None

    def saveEleLocalCoordSys(self,SaveName,EleLocalCoordSys):
        """
        Save element local coordinate systems to database
        localZSaveName(str)-the name of the saved table
        localZList(list)-for real length element, [['realEle',nodeI,nodeJ,localZTag]]
                        -for zeroLength ele or node, [['specialEle',nodeI,nodeJ,(localX_x,localX_y,localX_z),
                        (localY_x,localY_y,localY_z)]]
        """
        eleTypeName=EleLocalCoordSys[0][0]
        EleLocalCoordSys=[each[1:] for each in EleLocalCoordSys]
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('eleLocalCoordSys/'+SaveName)
            if dataSet is None:
                dataSet = f.create_dataset('eleLocalCoordSys/'+SaveName, [1,len(EleLocalCoordSys[0])],
                maxshape=[None,len(EleLocalCoordSys[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['nodeI:int','nodeJ:int','localZTag/localX,localY']
                eleCol=len(EleLocalCoordSys[0])
                dataSet.attrs['eleType']=eleTypeName
                newShape=(len(EleLocalCoordSys),len(EleLocalCoordSys[0]))
                dataSet.resize(newShape)
                dataSet[:]=EleLocalCoordSys
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(EleLocalCoordSys),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=EleLocalCoordSys

    def getEleLocalCoordSys(self, saveEleLocalCoordSysName):
        """
        return EleLocalCoordSys from database
        saveEleLocalCoordSysName(str)-the table name of saved EleLocalCoordSysName
        """
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('eleLocalCoordSys/'+saveEleLocalCoordSysName)
            eleTypeStr=dataSet.attrs['eleType']
            if dataSet is not None:
                returnDataSet=[list(each) for each in dataSet]
                return returnDataSet,eleTypeStr
            else:
                print(f'''table {saveEleLocalCoordSysName} doesn't exitst!''')
                return None,None

    def saveNodeTimeHistory(self,nodeSaveName,nodeHistoryList):
        """
        ---Save node time history responses to database---
        Inputs:
            nodeSaveName(str)-a table name for the saved responses, e.g,'node_disp_1'
            nodeHistoryList(list)-e.g.,[[time0,U1_0,U2_0,U3_0],[time1,U1_1,U2_1,U3_1],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('timeHistoryResponse/'+nodeSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('timeHistoryResponse/'+nodeSaveName, [1,len(nodeHistoryList[0])],
                maxshape=[None,len(nodeHistoryList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['time:float','dof1Value:float','dof2Value:float','dof3Value:float']
                newShape=(len(nodeHistoryList),len(nodeHistoryList[0]))
                dataSet.resize(newShape)
                dataSet[:]=nodeHistoryList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(nodeHistoryList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=nodeHistoryList

    def getNodeTimeHistory(self,nodeTag,resType,dof):
        """
        ---return node time history response---
        Inputs:
            nodeTag(int)-the tag of the inquried node
            resType(str)-node response type, resType='disp','vel','accel' or 'reaction'
            dof(int)-return the corresponding dof response,dof=1,2 or 3
        Outputs:
            historyResList(list),[times0,times1,times2,...],[response0,response1,response2,...]
        """
        tableName = 'node_' + resType + '_' + str(nodeTag)
        with h5py.File(self._dbPath, 'r') as f:
            dataSet=f.get('timeHistoryResponse/'+tableName)
            if dataSet is not None:
                timesList=[]
                responseList=[]
                [[timesList.append(each[0]),responseList.append(each[dof])] for each in dataSet]
                return timesList,responseList
            else:
                print(f'''table {tableName} doesn't exitst!''')
                return None,None

    def saveTrussEleResponseTimeHistory(self,eleSaveName,eleHistoryList):
        """
        ---Save truss element axial force and deformation database---
        elesSaveName(str)-a table name for saved truss element responses, e.g. 'element_axialForce_1','element_axialDeform_1'
        eleHistoryList(list)-e.g.,[[time0,resValue0],[time1,resValue1],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('timeHistoryResponse/'+eleSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('timeHistoryResponse/'+eleSaveName, [1,len(eleHistoryList[0])],
                maxshape=[None,len(eleHistoryList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['time:float','resValue:float']
                newShape=(len(eleHistoryList),len(eleHistoryList[0]))
                dataSet.resize(newShape)
                dataSet[:]=eleHistoryList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(eleHistoryList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=eleHistoryList

    def getTrussEleResponseTimeHistory(self,eleTag,resType):
        """
        ---Return the truss element time history response---
        Inputs:
            eleTag(int)-the tag of the inquried element
            resType(str)-element response type, resType='axialForce','axialDeform'
        Outputs:
            timesList,responsesList

        """
        tableName='trussEle_'+resType+'_'+str(eleTag)
        with h5py.File(self._dbPath, 'r') as f:
            dataSet = f.get('timeHistoryResponse/' + tableName)
            if dataSet is not None:
                timesList = []
                responseList = []
                [[timesList.append(each[0]), responseList.append(each[1])] for each in dataSet]
                return timesList, responseList
            else:
                print(f'''table {tableName} doesn't exitst!''')
                return None,None

    def saveZeroEleResponseTimeHistory(self,eleSaveName,eleHistoryList,dirList):
        """
        ---Save zeroLength element responses database---
        elesSaveName(str)-a table name for saved truss element responses, e.g. 'zeroEle_deformation_1','zeroEle_localForce_1'
        eleHistoryList(list)-e.g.,[[time0,resValue0_1,resValue0_2,resValue0_3],[],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('timeHistoryResponse/'+eleSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('timeHistoryResponse/'+eleSaveName, [1,len(eleHistoryList[0])],
                maxshape=[None,len(eleHistoryList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['[time,resValue0_1,resValue0_2,resValue0_3],[],...]']
                dataSet.attrs['dirList']=dirList
                newShape=(len(eleHistoryList),len(eleHistoryList[0]))
                dataSet.resize(newShape)
                dataSet[:]=eleHistoryList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(eleHistoryList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=eleHistoryList

    def getZeroEleResponseTimeHistory(self,eleTag,resType,dof):
        """
        ---Return the zeroLength element time history response---
        Inputs:
            eleTag(int)-the tag of the inquried element
            resType(str)-element response type, resType='localForce','deformation'
            dof(int)-return the corresponding dof response,dof=1,2,or 3
        Outputs:
            timesList,responsesList
        """
        tableName = 'zeroEle_' + resType + '_' + str(eleTag)
        with h5py.File(self._dbPath, 'r') as f:
            dataSet = f.get('timeHistoryResponse/' + tableName)
            if dataSet is not None:
                timesList = []
                responseList = []
                dirList=[each for each in dataSet.attrs['dirList']]
                if dof not in dirList:
                    [[timesList.append(each[0]),responseList.append(0)] for each in dataSet]
                    return timesList,responseList
                dofIndex=dirList.index(dof)+1
                [[timesList.append(each[0]), responseList.append(each[dofIndex])] for each in dataSet]
                return timesList, responseList
            else:
                print(f'''table {tableName} doesn't exitst!''')
                return None,None

    def saveNonEleSectResponseTimeHistory(self,eleSaveName,eleHistoryList):
        """
        ---Save nonlinear element section responses database---
        elesSaveName(str)-a table name for saved truss element responses, e.g. 'nonEle_sectionForce_1','nonEle_sectionDeformation_1'
        eleHistoryList(list)-e.g.,[[time0,response0_1,response0_2,response0_3,responses0_4],[],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('timeHistoryResponse/'+eleSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('timeHistoryResponse/'+eleSaveName, [1,len(eleHistoryList[0])],
                maxshape=[None,len(eleHistoryList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['[time,resValue0_1,resValue0_2,resValue0_3,response0_4],[],...]']
                newShape=(len(eleHistoryList),len(eleHistoryList[0]))
                dataSet.resize(newShape)
                dataSet[:]=eleHistoryList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(eleHistoryList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=eleHistoryList

    def getNonEleSectResponseTimeHistory(self,eleTag,resType,dof):
        """
        ---Return the nonlinear element section time history response---
        Inputs:
            eleTag(int)-the tag of the inquried element
            resType(str)-element response type, resType='sectionForce','sectionDeformation'
            dof(int)-return the corresponding dof response,dof=1,2,3,or 4
                   1-axial direction, 2-rotate about local z,3-rotate about local y, 4-rotate about local x
        Outputs:
            timesList,responsesList
        """
        tableName = 'nonEle_' + resType + '_' + str(eleTag)
        with h5py.File(self._dbPath, 'r') as f:
            dataSet = f.get('timeHistoryResponse/' + tableName)
            if dataSet is not None:
                timesList = []
                responseList = []
                [[timesList.append(each[0]), responseList.append(each[dof])] for each in dataSet]
                return timesList, responseList
            else:
                print(f'''table {tableName} doesn't exitst!''')
                return None,None

    def saveNonZeroEleResponseTimeHistory(self,eleSaveName,eleHistoryList):
        """
        ---Save non zerolength element localForce responses database---
        elesSaveName(str)-a table name for saved truss element responses, e.g. 'nonZeroEle_localForce_1'
        eleHistoryList(list)-e.g.,[[time0,response0_I1,response0_I2,response0_I3,response0_I4,response0_I5,
        response0_I6,response0_J1,response0_J2,response0_J3,response0_J4,response0_J5,response0_J6],[],...]
        """
        with h5py.File(self._dbPath,'a') as f:
            dataSet=f.get('timeHistoryResponse/'+eleSaveName)
            if dataSet is None:
                dataSet = f.create_dataset('timeHistoryResponse/'+eleSaveName, [1,len(eleHistoryList[0])],
                maxshape=[None,len(eleHistoryList[0])],chunks=True, compression='gzip',compression_opts=7)
                dataSet.attrs['column_names']=['[time,resValue0_I1,resValue0_I2,resValue0_I3,resValue0_I4,'
                'resValue0_I5,resValue0_I6,resValue0_J1,resValue0_J2,resValue0_J3,resValue0_J4,resValue0_J5,'
                                               'resValue0_J6],[],...]']
                newShape=(len(eleHistoryList),len(eleHistoryList[0]))
                dataSet.resize(newShape)
                dataSet[:]=eleHistoryList
            else:
                currentShape=dataSet.shape
                newShape=(currentShape[0]+len(eleHistoryList),currentShape[1])
                dataSet.resize(newShape)
                dataSet[currentShape[0]:]=eleHistoryList

    def getNonZeroEleResponseTimeHistory(self,eleTag,resType,eleEnd,dof):
        """
        ---Return the non zerolength element time history response---
        Inputs:
            eleTag(int)-the tag of the inquried element
            resType(str)-element response type, resType='localForce'
            eleEnd(str)-the end of the element, 'I' or 'J'
            dof(int)-return the corresponding dof response,dof=1,2,3,4,5 or 6
        Outputs:
            timesList,responsesList
        """
        tableName = 'nonZeroEle_' + resType + '_' + str(eleTag)
        with h5py.File(self._dbPath, 'r') as f:
            dataSet = f.get('timeHistoryResponse/' + tableName)
            if dataSet is not None:
                timesList = []
                responseList = []
                if eleEnd=='I':
                    indexValue=dof
                elif eleEnd=="J":
                    indexValue=6+dof
                [[timesList.append(each[0]), responseList.append(each[indexValue])] for each in dataSet]
                return timesList, responseList
            else:
                print(f'''table {tableName} doesn't exitst!''')
                return None,None
    ####################################################################################################################
    ####################################################################################################################
    @classmethod
    def dataSetsInGroup(cls, dbPath, groupName):
        """
        """
        with h5py.File(dbPath, 'r') as f:
            if groupName in f:
                group = f[groupName]
                datasets = list(group.keys())
                return datasets
            else:
                print(f"Group {groupName} does not exist.")
                return None

    @classmethod
    def partialMatchGroups(cls, dbPath, partialGroupName):
        """"""
        with h5py.File(dbPath, 'r') as f:
            matching_groups = []
            for group_name in f.keys():
                partialName = f"{partialGroupName}*"
                if fnmatch.fnmatch(group_name, partialName):
                    matching_groups.append(group_name)
            return matching_groups

    @classmethod
    def partialMatchDataSets(cls, dbPath, groupName, partialDataSetName):
        """"""
        with h5py.File(dbPath, 'r') as f:
            group = f[groupName]
            partial_name = f'{partialDataSetName}*'
            matching_datasets = []
            for name, item in group.items():
                if isinstance(item, h5py.Dataset):
                    if fnmatch.fnmatch(name, partial_name):
                        matching_datasets.append(name)
            return matching_datasets

    @classmethod
    def deleteData(cls, dbPath, dataSetName):
        """"""
        with h5py.File(dbPath, 'a') as f:
            if dataSetName in f:
                del f[dataSetName]

    def saveResult(self,dataName,resultsList,headNameList,operationIndexStr='replace'):
        """
        ---A general save template for hdf5 database ---
        Save results to database, resultsList=[[result0_0,result0_1,],[],...]
        headNameList=[headName_0,headName_1,...]
        dataName(str)
        operationIndexStr='replace' or 'append'
                          'replace' means replace the data, 'append' means append data after the last row
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

            with h5py.File(self.resultFileName, 'a') as f:
                dataset = f.get(dataName)
                if dataset is None:
                    dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype,chunks=True,
                                               compression='gzip',compression_opts=9,shuffle=True)
                    new_size = len(resultsList)
                    dataset.resize(new_size, axis=0)
                    dataset[0:new_size] = structured_data
                else:
                    if operationIndexStr == 'replace':
                        del f[dataName]
                        dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype,chunks=True,
                                               compression='gzip',compression_opts=9,shuffle=True)
                        new_size = len(resultsList)
                        dataset.resize(new_size, axis=0)
                        dataset[0:new_size] = structured_data
                    elif operationIndexStr == 'append':
                        original_size = dataset.shape[0]
                        new_size = original_size + len(resultsList)
                        dataset.resize(new_size, axis=0)
                        dataset[original_size:new_size] = structured_data

    @classmethod
    def getResult(cls, dbPath, dataName):
        """
        dataName(str)---multi-level data table name,eg. 'group/dataSetName'
        """
        with h5py.File(dbPath, 'r') as f:
            dataSet = f.get(dataName)
            if dataSet is not None:
                returnList = []
                [[tempList := [], [tempList.append(each.decode("utf-8")) if isinstance(each, bytes) else
                                   tempList.append(float(each)) if isinstance(each, (np.float32, np.float64)) else
                                   tempList.append(int(each)) if isinstance(each, np.int32) else None
                                   for each in eachRow], returnList.append(tempList)] for eachRow in dataSet[:]]
                return returnList
            else:
                print(f'''table {dataName} doesn't exitst!''')
                return None
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
