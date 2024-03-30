########################################################################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 03/30/2024
########################################################################################################################
import numpy as np
import os
import math
import pyvista as pv
import matplotlib.pyplot as plt
import sys
########################################################################################################################
class EleMeshPlotAndSelect():
    """
    --------------------------------------------------------------------------------------------------------------------
    A class for visualizing finite elements meshed with pygmsh module and selecting mesh components (version:0.1.0)
    Environemet: Successfully executed in python 3.11
    Date: 2024-01-15
    --------------------------------------------------------------------------------------------------------------------
    Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
    --------------------------------------------------------------------------------------------------------------------
        ** **************************************************************************** **
        ** (C) Copyright 2024, School of Civil Engineering,Beijing Jiaotong University  **
        ** All Rights Reserved.                                                         **
        **                                                                              **
        ** Commercial use of this program is strictly prohibited.                       **
        **                                                                              **
        ** Developed by:                                                                **
        **   Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo      **
        **   jjguo2@bjtu.edu.cn/guojj_ce@163.com                                        **
        ** **************************************************************************** **
    """
    def __init__(self):
        """
        ----------------------------------------------------------------------------------------------------------------
        Initialize the class
        ------------------------------------------
        ----------------------------------------------------------------------------------------------------------------
        """
        self.cellTypeDict_1D = {'line': pv.CellType.LINE,'line3':pv.CellType.CUBIC_LINE}
        self.cellTypeDict_2D = {'triangle': pv.CellType.TRIANGLE,'quad':pv.CellType.QUAD,
                                'triangle6':pv.CellType.POLYGON}
        self.cellTypeDict_3D = {'tetra': pv.CellType.TETRA,'hexahedron':pv.CellType.HEXAHEDRON,
                                'wedge':pv.CellType.WEDGE}

        self.meshList = []  ###---mesh对象list
        self.nodeStartNumberList = []  ###---每个mesh部分节点编号开始数
        self.elementStartNumberList = []  ###---每个mesh部分单元编号开始数
        self.pvPlotGridList = []
        self.pointsTagDict = {}
        self.elemenTagDict = {}

        self.plotter = pv.Plotter(title='Element Mesh Visualization')
        self.plotter.show_axes()
        # self.plotter.background_color="gray"
        self.plotter.title

    def addMesh(self,mesh,nodeStartNumber,elementStartNumber):
        """
        ----------------------------------------------------------------------------------------------------------------
        Add pygmsh generated mesh
        -----------------------------
        Inputs:
            nodeStartNumber(int)-The start number of the meshed nodes
            elementStartNumber(int)-The start number of the meshed elements
        ----------------------------------------------------------------------------------------------------------------
        """
        self.meshList.append(mesh)
        self.nodeStartNumberList.append(nodeStartNumber)
        self.elementStartNumberList.append(elementStartNumber)

    def meshProcess(self,showPointTag=False):
        """
        ----------------------------------------------------------------------------------------------------------------
        Add pygmsh generated mesh
        -----------------------------
        Inputs:
        ----------------------------------------------------------------------------------------------------------------
        """
        for i0 in range(len(self.meshList)):
            eachMesh=self.meshList[i0]
            points = eachMesh.points
            nodesLabelList = [self.nodeStartNumberList[i0] + i1 for i1 in range(len(points))]
            pointsList=points.tolist()
            for eachTag,eachValue in zip(nodesLabelList,pointsList):
                self.pointsTagDict[eachTag]=eachValue
            cellDict=eachMesh.cells_dict
            keys=cellDict.keys()
            eleStartTag=self.elementStartNumberList[i0]
            eleIncrTag=eleStartTag
            for eachKey in keys:
                if eachKey!='vertex':
                    eachTypeEle=cellDict[eachKey]
                    for eachEle in eachTypeEle:
                        eleNodes=(eachEle+self.nodeStartNumberList[i0]).tolist()
                        self.elemenTagDict[eleIncrTag]=[eachKey,eleNodes]
                        eleIncrTag+=1

            def pointCallback(pickedNode):
                """"""
                pickedNodeList=pickedNode.tolist()
                for key,value in self.pointsTagDict.items():
                    nodeDist=np.sqrt((pickedNodeList[0]-value[0])**2+(pickedNodeList[1]-value[1])**2+
                                     (pickedNodeList[2]-value[2])**2)
                    if nodeDist<=1.0e-16:
                        self.plotter.add_title(f"Node tag: {key}     coordinates:{(round(value[0],3),round(value[1],3),round(value[2],3))}",
                        font='times', color='red', font_size=10)
                        break

            for each3D in self.cellTypeDict_3D:
                if each3D in keys:
                    ele_3D = cellDict[each3D].tolist()
                    elesMesh = [[len(each)] + each for each in ele_3D]
                    eleTypesList= [self.cellTypeDict_3D[each3D] for each in elesMesh]
                    grid = pv.UnstructuredGrid(np.array(elesMesh), np.array(eleTypesList), points)
                    self.pvPlotGridList.append(grid)
                    self.plotter.add_mesh(grid, show_edges=True, color='lightBlue', style='surface', point_size=2,
                                          line_width=1,opacity=1)

                    self.plotter.add_mesh(points, color="k", point_size=5, render_points_as_spheres=True)
                    self.plotter.enable_point_picking(left_clicking=True,callback=pointCallback,show_message=False)
            if showPointTag:
                nodesLabelList=[self.nodeStartNumberList[i0]+i1 for i1 in range(len(points))]
                self.plotter.add_point_labels(points, nodesLabelList, bold=False, point_size=0.001,
                                              font_size=16, text_color='blue', font_family="times",show_points=False,
                                              shape=None,tolerance=1.0)
            if len(self.pvPlotGridList)>0:
                continue

            for each2D in self.cellTypeDict_2D:
                if each2D in keys:
                    if each2D=='triangle6':
                        ele_2D = cellDict[each2D].tolist()
                        elesMesh = [[len(each[:3])] + each[:3] for each in ele_2D]
                        eleTypesList = [pv.CellType.TRIANGLE for each in elesMesh]
                        grid = pv.UnstructuredGrid(np.array(elesMesh), np.array(eleTypesList), points)
                        self.pvPlotGridList.append(grid)
                        self.plotter.add_mesh(grid, show_edges=True, color='lightBlue', style='surface', point_size=2,
                                              line_width=1, opacity=1)

                        self.plotter.add_mesh(points, color="k", point_size=5, render_points_as_spheres=True)
                        self.plotter.enable_point_picking(left_clicking=True, callback=pointCallback,
                                          show_message=False)
                    else:
                        ele_2D = cellDict[each2D].tolist()
                        elesMesh = [[len(each)] + each for each in ele_2D]
                        eleTypesList= [self.cellTypeDict_2D[each2D] for each in elesMesh]
                        grid = pv.UnstructuredGrid(np.array(elesMesh), np.array(eleTypesList), points)
                        self.pvPlotGridList.append(grid)
                        self.plotter.add_mesh(grid, show_edges=True, color='lightBlue', style='surface', point_size=2,
                                              line_width=1,opacity=1)

                        self.plotter.add_mesh(points, color="k", point_size=5, render_points_as_spheres=True)
                        self.plotter.enable_point_picking(left_clicking=True,callback=pointCallback,show_message=False)
            if showPointTag:
                nodesLabelList=[self.nodeStartNumberList[i0]+i1 for i1 in range(len(points))]
                self.plotter.add_point_labels(points, nodesLabelList, bold=False, point_size=0.001,
                                              font_size=16, text_color='blue', font_family="times",show_points=False,
                                              shape=None,tolerance=1.0)
            if len(self.pvPlotGridList)>0:
                continue

            for each1D in self.cellTypeDict_1D:
                if each1D in keys:
                    ele_1D = cellDict[each1D].tolist()
                    elesMesh = [[len(each)] + each for each in ele_1D]
                    eleTypesList= [self.cellTypeDict_1D[each1D] for each in elesMesh]
                    grid = pv.UnstructuredGrid(np.array(elesMesh), np.array(eleTypesList), points)
                    self.pvPlotGridList.append(grid)
                    self.plotter.add_mesh(grid, show_edges=True, color='lightBlue', style='surface', point_size=2,
                                          line_width=1,opacity=1)

                    self.plotter.add_mesh(points, color="k", point_size=5, render_points_as_spheres=True)
                    self.plotter.enable_point_picking(left_clicking=True,callback=pointCallback,show_message=False)
            if showPointTag:
                nodesLabelList=[self.nodeStartNumberList[i0]+i1 for i1 in range(len(points))]
                self.plotter.add_point_labels(points, nodesLabelList, bold=False, point_size=0.001,
                                              font_size=16, text_color='blue', font_family="times",show_points=False,
                                              shape=None,tolerance=1.0)
            if len(self.pvPlotGridList)>0:
                continue

    def meshPlot(self):
        """
        ----------------------------------------------------------------------------------------------------------------
        Show the element mesh
        -----------------------------
        Inputs:
        ----------------------------------------------------------------------------------------------------------------
        """
        self.plotter.show()

    def saveNodes(self,saveName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Save the points as txt format, [nodeTag,X,Y,Z]
        -----------------------------
        Inputs:
            saveName(str)-the Name of the file, eg."nodes"
        ----------------------------------------------------------------------------------------------------------------
        """
        nodesDict=self.pointsTagDict
        saveNodesList=[[key]+values for key,values in nodesDict.items()]
        np.savetxt(f"{saveName}.txt",saveNodesList,fmt="%d    %.6f    %.6f    %.6f")

    def saveElements(self,saveName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Save the elements as txt format, different element type will be identified with different postfix
        -----------------------------
        Inputs:
            saveName(str)-the Name of the file, eg."savedElements"
        ----------------------------------------------------------------------------------------------------------------
        """
        elesDict=self.elemenTagDict
        eleTypes=set([each[0] for each in elesDict.values()])
        for eachType in eleTypes:
            saveList=[]
            for key,value in elesDict.items():
                if value[0]==eachType:
                    saveList.append([key]+value[1:][0])
            np.savetxt(f"{saveName}_{eachType}.txt",saveList,fmt="%d")

    def selectNodes_inLine(self,startNodeTag,endNodeTag,saveSetName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Select nodes in a line segment after providing start and end node Tags
        -----------------------------
        Inputs:
            startNodeTag(int)-the tag of the start node
            endNodeTag(int)-the tag of the end node
            saveSetName(str)-the name of the saved file
        ----------------------------------------------------------------------------------------------------------------
        """
        starNodeCoord=self.pointsTagDict[startNodeTag]
        endNodeCoord=self.pointsTagDict[endNodeTag]
        x1,y1,z1=starNodeCoord
        x2,y2,z2=endNodeCoord
        minX,maxX=min(x1,x2),max(x1,x2)
        minY,maxY=min(y1,y2),max(y1,y2)
        minZ,maxZ=min(z1,z2),max(z1,z2)

        m=x2-x1
        n=y2-y1
        p=z2-z1
        def twoNodeEq(xi,yi,zi):
            """"""
            MINVALUE=1.0e-8
            if (np.abs(m)>MINVALUE)and(np.abs(n)>MINVALUE)and(np.abs(p)>MINVALUE):
                value1=(xi-x1)/float(m)
                value2=(yi-y1)/float(n)
                value3=(zi-z1)/float(p)
                if (np.abs(value1-value2)<MINVALUE)and(np.abs(value1-value3)<MINVALUE)and(np.abs(value2-value3)<MINVALUE):
                    return True
            if (np.abs(m)<MINVALUE)and(np.abs(n)>MINVALUE)and(np.abs(p)>MINVALUE):
                value2 = (yi - y1) / float(n)
                value3 = (zi - z1) / float(p)
                if (np.abs(value2 - value3) < MINVALUE):
                    return True
            if (np.abs(m)>MINVALUE)and(np.abs(n)<MINVALUE)and(np.abs(p)>MINVALUE):
                value1 = (xi - x1) / float(m)
                value3 = (zi - z1) / float(p)
                if (np.abs(value1 - value3) < MINVALUE):
                    return True
            if (np.abs(m)>MINVALUE)and(np.abs(n)>MINVALUE)and(np.abs(p)<MINVALUE):
                value1 = (xi - x1) / float(m)
                value2 = (yi - y1) / float(n)
                if (np.abs(value1 - value2) < MINVALUE):
                    return True
            if (np.abs(m)<MINVALUE)and(np.abs(n)<MINVALUE)and(np.abs(p)>MINVALUE):
                if (minZ<=zi<=maxZ):
                    return True
            if (np.abs(m)<MINVALUE)and(np.abs(n)>MINVALUE)and(np.abs(p)<MINVALUE):
                if (minY<=yi<=maxY):
                    return True
            if (np.abs(m)>MINVALUE)and(np.abs(n)<MINVALUE)and(np.abs(p)<MINVALUE):
                if (minX<=xi<=maxX):
                    return True
            return False
        selectedNodesList=[]
        savedNodeTag=[]
        for eachNodeTag,nodeCoord in self.pointsTagDict.items():
            xi,yi,zi=nodeCoord
            if(minX<=xi<=maxX)and(minY<=yi<=maxY)and(minZ<=zi<=maxZ)and(twoNodeEq(xi,yi,zi)):
                selectedNodesList.append([eachNodeTag]+nodeCoord)
                savedNodeTag.append(eachNodeTag)
        np.savetxt(f"{saveSetName}_nodeSet.txt",savedNodeTag,fmt="%d")
        nodesCoordList = [each[1:] for each in selectedNodesList]
        vertices = np.array(nodesCoordList)
        surf = pv.PolyData(vertices)  # vtk format data
        self.plotter.add_mesh(surf, color="red", point_size=8, render_points_as_spheres=True)

    def selectNodes_inPlane(self,planeNode1Tag,planeNode2Tag,planeNode3Tag,saveSetName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Select all nodes in a plane that defined with three different nodes in the plane
        -----------------------------
        Inputs:
            planeNode1Tag(int)-the tag of the first node in the plane
            planeNode2Tag(int)-the tag of the second node in the plane
            planeNode3Tag(int)-the tag of the third node in the plane
            saveSetName(str)-the name of the saved file
        ----------------------------------------------------------------------------------------------------------------
        """
        MINVALUE=1.0E-8
        nodeCoord_1 = self.pointsTagDict[planeNode1Tag]
        nodeCoord_2= self.pointsTagDict[planeNode2Tag]
        nodeCoord_3 = self.pointsTagDict[planeNode3Tag]
        x1, y1, z1 = nodeCoord_1
        x2, y2, z2 = nodeCoord_2
        x3, y3, z3 = nodeCoord_3
        vector1 = [x2 - x1, y2 - y1, z2 - z1]
        vector2 = [x3 - x1, y3 - y1, z3 - z1]
        normal_vector = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                         vector1[2] * vector2[0] - vector1[0] * vector2[2],
                         vector1[0] * vector2[1] - vector1[1] * vector2[0]]
        A, B, C = normal_vector
        D = - (A * x1 + B * y1 + C * z1)
        def nodeInPlane(xi,yi,zi):
            """"""
            calValue=A*xi+B*yi+C*zi+D
            if np.abs(calValue)<MINVALUE:
                return True
            else:
                return False

        selectedNodesList = []
        savedNodeTag = []
        for eachNodeTag, nodeCoord in self.pointsTagDict.items():
            xi, yi, zi = nodeCoord
            if (nodeInPlane(xi, yi, zi)):
                selectedNodesList.append([eachNodeTag] + nodeCoord)
                savedNodeTag.append(eachNodeTag)
        np.savetxt(f"{saveSetName}_nodeSet.txt", savedNodeTag, fmt="%d")
        nodesCoordList = [each[1:] for each in selectedNodesList]
        vertices = np.array(nodesCoordList)
        surf = pv.PolyData(vertices)  # vtk format data
        self.plotter.add_mesh(surf, color="red", point_size=8, render_points_as_spheres=True)

    def selectNodes_XYZRanges(self,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,saveSetName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Select all nodes within a box with length in [xmin,xmax], widht in [ymin,ymax], height in [zmim,zmax]
        -----------------------------
        Inputs:
            Xmin,Xmax(float)-the lower and upper bounds for the box length
            Ymin,Ymax(float)-the lower and upper bounds for the box width
            Zmin,Zmax(float)-the lower and upper bounds for the box height
            saveSetName(str)-the name of the saved file
        ----------------------------------------------------------------------------------------------------------------
        """
        selectedNodesList = []
        savedNodeTag = []
        for eachNodeTag, nodeCoord in self.pointsTagDict.items():
            xi, yi, zi = nodeCoord
            if (Xmin<=xi<=Xmax)and(Ymin<=yi<=Ymax)and(Zmin<=zi<=Zmax):
                selectedNodesList.append([eachNodeTag] + nodeCoord)
                savedNodeTag.append(eachNodeTag)
        np.savetxt(f"{saveSetName}_nodeSet.txt", savedNodeTag, fmt="%d")
        nodesCoordList = [each[1:] for each in selectedNodesList]
        vertices = np.array(nodesCoordList)
        surf = pv.PolyData(vertices)  # vtk format data
        self.plotter.add_mesh(surf, color="red", point_size=8, render_points_as_spheres=True)

    def selectNodes_betweenTwoConcentricCylinders(self,circleCenter,radius1,radius2,Zmin,Zmax,saveSetName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Select all nodes within two concentric cylinders
        -----------------------------
        Inputs:
            circleCenter(float,float)-the X and Y coordinates of the center of the cylinder
            radius1(float)-the radius of the internal circle
            radius2(float)-the radius of the external circle
            Zmin(float)-the bottom height of the cylinder
            Mmax(float)-the top height of the cylinder
            saveSetName(str)-the name of the saved file
        ----------------------------------------------------------------------------------------------------------------
        """
        selectedNodesList = []
        savedNodeTag = []
        x0,y0=circleCenter
        minR,maxR=min(radius1,radius2),max(radius1,radius2)
        for eachNodeTag, nodeCoord in self.pointsTagDict.items():
            xi, yi, zi = nodeCoord
            distL=np.sqrt((xi-x0)**2+(yi-y0)**2)
            if (minR<=distL<=maxR)and(Zmin<=zi<=Zmax):
                selectedNodesList.append([eachNodeTag] + nodeCoord)
                savedNodeTag.append(eachNodeTag)
        np.savetxt(f"{saveSetName}_nodeSet.txt", savedNodeTag, fmt="%d")
        nodesCoordList = [each[1:] for each in selectedNodesList]
        vertices = np.array(nodesCoordList)
        surf = pv.PolyData(vertices)  # vtk format data
        self.plotter.add_mesh(surf, color="red", point_size=8, render_points_as_spheres=True)


    def selectFaces_inPlane(self,planeNode1Tag,planeNode2Tag,planeNode3Tag,saveSetName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Select all faces within a plane that defined with three different nodes
        -----------------------------
        Inputs:
            planeNode1Tag(int)-the tag of the first node in the plane
            planeNode2Tag(int)-the tag of the second node in the plane
            planeNode3Tag(int)-the tag of the third node in the plane
            saveSetName(str)-the name of the saved file
        ----------------------------------------------------------------------------------------------------------------
        """
        MINVALUE = 1.0E-8
        nodeCoord_1 = self.pointsTagDict[planeNode1Tag]
        nodeCoord_2 = self.pointsTagDict[planeNode2Tag]
        nodeCoord_3 = self.pointsTagDict[planeNode3Tag]
        x1, y1, z1 = nodeCoord_1
        x2, y2, z2 = nodeCoord_2
        x3, y3, z3 = nodeCoord_3
        vector1 = [x2 - x1, y2 - y1, z2 - z1]
        vector2 = [x3 - x1, y3 - y1, z3 - z1]
        normal_vector = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                         vector1[2] * vector2[0] - vector1[0] * vector2[2],
                         vector1[0] * vector2[1] - vector1[1] * vector2[0]]
        A, B, C = normal_vector
        D = - (A * x1 + B * y1 + C * z1)

        def nodeInPlane(xi, yi, zi):
            """"""
            calValue = A * xi + B * yi + C * zi + D
            if np.abs(calValue) < MINVALUE:
                return True
            else:
                return False
        elesMesh=[]
        for eleTag,eleValues in self.elemenTagDict.items():
            eleNodes=eleValues[1]
            eachEleList=[]
            if len(eleNodes)>=3:
                for eachNode in eleNodes:
                    eachNodeCoord=self.pointsTagDict[eachNode]
                    xi,yi,zi=eachNodeCoord
                    if nodeInPlane(xi,yi,zi):
                        eachEleList.append(eachNode)
            if len(eachEleList)>=3:
                elesMesh.append([len(eachEleList)]+eachEleList)
        eleTypesList = [pv.CellType.POLYGON for each in elesMesh]
        seletedNodesList=[]
        for eachEle in elesMesh:
            for eachNode in eachEle[1:]:
                if eachNode not in seletedNodesList:
                    seletedNodesList.append(eachNode)
        selectedNodesDict={seletedNodesList[i1]:i1 for i1 in range(len(seletedNodesList))}
        points=np.array([self.pointsTagDict[each] for each in seletedNodesList])
        meshes=[]
        for eacE in elesMesh:
            eTag=eacE[0]
            eachList=[]
            for eachN in eacE[1:]:
                eachList.append(selectedNodesDict[eachN])
            meshes.append([eTag]+eachList)
        grid = pv.UnstructuredGrid(np.array(meshes), np.array(eleTypesList), points)
        self.pvPlotGridList.append(grid)
        self.plotter.add_mesh(grid, show_edges=True, color='red', style='surface', point_size=2, line_width=1,
                              opacity=1)
########################################################################################################################
if __name__ == "__main__":
    pass







