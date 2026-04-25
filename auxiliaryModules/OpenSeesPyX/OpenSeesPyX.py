########################################################################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 05/02/2020
########################################################################################################################
import h5py
import os
import numpy as np
import openseespy.opensees as ops
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Wedge
from .modelResultsDB import h5pyDB  ##h5py database to store opensees results
###############################
import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QToolTip, QFileDialog, QPushButton, QVBoxLayout, \
    QHBoxLayout
from PySide6.QtCore import Qt, QPointF, QSize, QSizeF, QMarginsF
from PySide6.QtGui import QPainter, QPen, QBrush, QColor, QFont, QImage, QPainter as QPainterExport, QPageSize, \
    QPageLayout
from PySide6.QtPrintSupport import QPrinter
import math
########################################################################################################################
class OpenSeesPyX():
    """
    --------------------------------------------------------------------------------------------------------------------
    A class for openSeesPy visualization and structural model analysis (version:0.1.0)
    Environemet: Successfully executed in python 3.11
    Date: 2023-04-07
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
    def __init__(self,dataBaseName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Initialize the class
        ------------------------------------------
        Inputs:
            dataBaseName(str)-the name of the database
        ----------------------------------------------------------------------------------------------------------------
        """
        self.nodeSetNameList=[]
        self.eleSetNameList=[]
        self.modalNameList=[]
        self.specialEleTransfList = []
        self.EleLocalCoordSysSetNameList=[]
        self.localTransfNameList=[]
        self.materialNumberDict={}
        self.nodeLinkDict={}
        self.eleLinkDict={}
        self.nodeMassDict={}
        self.zeroEleDirDict={}
        self.dbPath = dataBaseName+".h5"
        self.saveInstance = h5pyDB(self.dbPath)
        h5pyDB.initDB(self.dbPath)
        #########----存储模型节点，单元等列表
        self.saveNodeList = []
        self.saveGeomfList=[]
        self.saveEleList=[]
        self.EleLocalCoordSys=[]
        #######----存储纤维截面信息
        self.currentSectTag=None
        self.fiberSectDict={}

    def auxiliary_vectXY(self,a):
        """
        ----------------------------------------------------------------------------------------------------------------
        obtain the normalized xvector and yvector by providing the axial vector of the element
        a-(elemental aixal vector),eg.(x2-x1,y2-y1,z2-z1)
        ------------------------------------------------------------
        return
        xVector-normalized axial vector
        yvector-normalized y vector that perpendicular to xVector
        ----------------------------------------------------------------------------------------------------------------
        """
        ax, ay, az = a

        # 选取绝对值最小分量对应的基向量
        if abs(ax) <= abs(ay) and abs(ax) <= abs(az):
            b = (1, 0, 0)
        elif abs(ay) <= abs(ax) and abs(ay) <= abs(az):
            b = (0, 1, 0)
        else:
            b = (0, 0, 1)

        def cross_product(a, b):
            return (
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
            )
        ######################################################################################
        v = cross_product(a, b)
        vectorNorma = np.linalg.norm(a)
        vectorNormv = np.linalg.norm(v)
        xVect = [a[0] / float(vectorNorma), a[1] / float(vectorNorma), a[2] / float(vectorNorma)]
        yvect = [v[0] / float(vectorNormv), v[1] / float(vectorNormv), v[2] / float(vectorNormv)]
        return xVect, yvect

    def auxiliary_localZvector(self,eleVector,refVector):
        """
        ----------------------------------------------------------------------------------------------------------------
        get local Z vector that perperticular to two vetors in a plane,one is a vector based on element(from I to J),
        and the other is easily specified,eg. 0,0,1
        -----------------------------
        Inputs:
            eleVector(list)-A element vector made up of I and J nodes, eg.(xj-xi,yj-yi,zj-zi)
            refVector(list)-A reference vector that in the same plane with eleVector, and perperticular to localZvector
        ----------------------------------------------------------------------------------------------------------------
        """
        a = np.array(eleVector)
        b = np.array(refVector)
        c = np.cross(a, b)
        vectorNorm = np.linalg.norm(c)
        localzVector = (c[0] / float(vectorNorm), c[1] / float(vectorNorm), c[2] / float(vectorNorm))
        return localzVector
    def auxiliary_materialReNumber(self,materialName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Register a material name in the domain and return a unique number of the mateiral
        -----------------------------
        Inputs:
            materialName(str)-the name tag of the material
        ----------------------------------------------------------------------------------------------------------------
        """
        if materialName not in self.materialNumberDict.keys():
            self.materialNumberDict[materialName]=len(self.materialNumberDict.keys())+1000000
        return self.materialNumberDict[materialName]

    def auxiliary_nodeReNumber(self,nodeName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Register a node name in the domain and return a unique number of the node
        -----------------------------
        Inputs:
            nodeName(str)-the name tag of the node
        ----------------------------------------------------------------------------------------------------------------
        """
        if nodeName not in self.nodeLinkDict.keys():
            self.nodeLinkDict[nodeName] = len(self.nodeLinkDict.keys()) + 1000000
        return self.nodeLinkDict[nodeName]

    def auxiliary_eleReNumber(self,eleName):
        """
        ----------------------------------------------------------------------------------------------------------------
        Register a element name in the domain and return a unique number of the element
        -----------------------------
        Inputs:
            eleName(str)-the name tag of the element
        ----------------------------------------------------------------------------------------------------------------
        """
        if eleName not in self.eleLinkDict.keys():
            self.eleLinkDict[eleName]=len(self.eleLinkDict.keys())+1000000
        return self.eleLinkDict[eleName]
    def auxiliary_nodeMassGenerate(self,nodeI,nodeJ,eleLength,massUnitLength):
        """
        ----------------------------------------------------------------------------------------------------------------
        node mass generate
        -----------------------------
        Inputs:
        nodeI,nodeJ(int)-girder end nodes
        eleLength(float)-girder element length
        massUnitLength(float)-mass per unit length
        Outputs:
        finally return the lumped node masses with call nodeMassDict
        ----------------------------------------------------------------------------------------------------------------
        """
        totalMass = eleLength * massUnitLength
        if nodeI not in self.nodeMassDict.keys():
            self.nodeMassDict[nodeI] = round(0.5 * totalMass, 6)
        else:
            self.nodeMassDict[nodeI] = round(self.nodeMassDict[nodeI] + 0.5 * totalMass, 6)
        if nodeJ not in self.nodeMassDict.keys():
            self.nodeMassDict[nodeJ] = round(0.5 * totalMass, 6)
        else:
            self.nodeMassDict[nodeJ] = round(self.nodeMassDict[nodeJ] + 0.5 * totalMass, 6)

    def auxiliary_concreteMaterialProperty_JTG2004(self,concreteTag="C40"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the concrete material properties based on Chinese JTG standard (JTG D62-2004 C15)---
        concreteMatProDict=["concreteTag":[unitWeight(kN/m3),E,poission,coefficient of thermal expansion,G]]
        ----------------------------------------------------------------------------------------------------------------
        """
        concreteMatProDict = {"C15": [25, 22000000, 0.2, 1.000E-05, 9166667],
                              "C20": [25, 25500000, 0.2, 1.000E-05, 10625000],
                              "C25": [25, 28000000, 0.2, 1.000E-05, 11666667],
                              "C30": [25, 30000000, 0.2, 1.000E-05, 12500000],
                              "C35": [25, 31500000, 0.2, 1.000E-05, 13125000],
                              "C40": [25, 32500000, 0.2, 1.000E-05, 13541667],
                              "C45": [25, 33500000, 0.2, 1.000E-05, 13958333],
                              "C50": [25, 34500000, 0.2, 1.000E-05, 14375000],
                              "C55": [25, 35500000, 0.2, 1.000E-05, 14791667],
                              "C60": [25, 36000000, 0.2, 1.000E-05, 15000000],
                              "C65": [25, 36500000, 0.2, 1.000E-05, 15208333],
                              "C70": [25, 37000000, 0.2, 1.000E-05, 15416667],
                              "C75": [25, 37500000, 0.2, 1.000E-05, 15625000],
                              "C80": [25, 38000000, 0.2, 1.000E-05, 15833333]}
        return concreteMatProDict[concreteTag]

    def auxiliary_rebarMaterialProperty_GB(self,rebarTag="HRB400"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the rebar material properties based on Chinese GB standard (GB50010)---
        reBarMatProDict={"rebarTag":[unitWeight,E,poission,coefficient of thermal expansion,G,Fy,Fu]}
        ----------------------------------------------------------------------------------------------------------------
        """
        reBarMatProDict = {"HPB300": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 300000, 420000],
                           "HRB335": [77, 2.000E+08, 0.3, 1.170E-05, 76923077, 335000, 455000],
                           "HRB400": [77, 2.000E+08, 0.3, 1.170E-05, 76923077, 400000, 540000],
                           "HRB500": [77, 2.000E+08, 0.3, 1.170E-05, 76923077, 500000, 630000]}
        return reBarMatProDict[rebarTag]

    def auxiliary_rebarMaterialProperty_JTG2004(self,rebarTag="R235"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the rebar material properties based on Chinese JTG standard (JTG D62-2004)---
        reBarMatProDict={"rebarTag":[unitWeight,E,poission,coefficient of thermal expansion,G,Fy,Fu]}
        ----------------------------------------------------------------------------------------------------------------
        """
        reBarMatProDict = {"R235": [77, 2.100E+08, 0.3, 1.170E-05, 80769231,235000,290000],
                           "HRB335": [77, 2.000E+08, 0.3, 1.170E-05,76923077,335000,420000],
                           "HRB400": [77, 2.000E+08, 0.3, 1.170E-05, 76923077,400000,630000],
                           "KL400": [77, 2.000E+08, 0.3, 1.170E-05, 76923077, 400000, 630000]}
        return reBarMatProDict[rebarTag]

    def auxiliary_steelMaterialProperty_JTG2004(self,steelTag="Q235q"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the steel material properties based on Chinese JTG standard (JTG D62-2004)---
        steelMatProDict={"rebarTag":[unitWeight,E,poission,coefficient of thermal expansion,G,Fy,Fu]}
        ----------------------------------------------------------------------------------------------------------------
        """
        steelMatProDict = {"Q235q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 235000,400000],
                           "Q345q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 345000,510000],
                           "Q370q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 370000, 510000],
                           "Q420q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 420000, 540000],
                           "Q460q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 460000, 570000],
                           "Q500q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 500000, 600000],
                           "Q550q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 550000, 660000],
                           "Q620q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 620000, 720000],
                           "Q690q": [77, 2.100E+08, 0.3, 1.170E-05, 80769231, 690000, 770000]}
        return steelMatProDict[steelTag]

    def auxiliary_tendonMaterialProperty_JTG2004(self,tendonTag="fpk1470"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the tendon material properties based on Chinese JTG standard (JTG D62-2004)---
        tendonMatProDict={"tendonTag":[unitWeight,E,poission,coefficient of thermal expansion,G,Fy,Fu]}
        ----------------------------------------------------------------------------------------------------------------
        """
        tendonMatProDict = {"fpk1470": [77, 1.95E+08, 0.3, 1.170E-05,75000000,1470000,1690000],
                            "fpk1570": [77, 1.95E+08, 0.3, 1.170E-05, 75000000,1570000,1800000],
                            "fpk1720": [77, 1.95E+08, 0.3, 1.170E-05, 75000000,1720000,1980000],
                            "fpk1860": [77, 1.95E+08, 0.3, 1.170E-05, 75000000,1860000,2140000]}
        return tendonMatProDict[tendonTag]
    ####################################################################################################################
    def auxiliary_concreteUltimateStrength_TB10092_2017(self,concreteTag="C25"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the concrete ultimate strength based on Code for Design of Concrete Structures of Railway Bridge and Culvert
        (TB 10092-2017)---
        Input:
            concrete Tag(str), one of the grade "C25","C30","C35","C40","C45","C50","C55" or "C60"
        OutPut:
            ultimate strength (axial compressive, axial tensile) (KPa)
        ----------------------------------------------------------------------------------------------------------------
        """
        concreteUltimateStrengthDict = {"C25": [17.0*1000,2.0*1000],"C30": [20.0*1000,2.2*1000],
                                        "C35": [23.5*1000,2.5*1000],"C40": [27.0*1000,2.7*1000],
                                        "C45": [30.0*1000,2.9*1000],"C50": [33.5*1000,3.1*1000],
                                        "C55": [37.0*1000,3.3*1000],"C60": [40.0*1000,3.5*1000]}
        return concreteUltimateStrengthDict[concreteTag]
    ####################################################################################################################
    def auxiliary_concreteEGPoissonRation_TB10092_2017(self,concreteTag="C25"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---the concrete elastic modulus, shear modulus and Possion ration based on Code for Design of Concrete Structures
        of Railway Bridge and Culvert(TB 10092-2017)---
        Input:
            concrete Tag(str), one of the grade "C25","C30","C35","C40","C45","C50","C55" or "C60"
        OutPut:
            [elasticModulus(kPa), shearModulus(kPa), PossionRation]
        ----------------------------------------------------------------------------------------------------------------
        """
        concreteEGPDict = {"C25": [3.0*10**7,0.43*3.0*10**7,0.2],"C30": [3.2*10**7,0.43*3.2*10**7,0.2],
                           "C35": [3.3*10**7,0.43*3.3*10**7,0.2],"C40": [3.4*10**7,0.43*3.4*10**7,0.2],
                           "C45": [3.45*10**7,0.43*3.45*10**7,0.2],"C50": [3.55*10**7,0.43*3.55*10**7,0.2],
                           "C55": [3.6*10**7,0.43*3.6*10**7,0.2],"C60": [3.65*10**7,0.43*3.65*10**7,0.2]}
        return concreteEGPDict[concreteTag]
    ####################################################################################################################
    def auxiliary_tendonStandardTensileStrength_TB10092_2017(self,tendonCode="HPB300"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---tendon standard tensile strength based on Code for Design of Concrete Structures of Railway Bridge and
        Culvert(TB 10092-2017)---
        Input:
            tendon code(str), one of the grade "HPB300","HRB400","HRB500","C40","PSB830","PSB980"
        OutPut:
            [standardTensileStrength(KPa)]
        ----------------------------------------------------------------------------------------------------------------
        """
        retrunDict = {"HPB300": 300*10**3,"HRB400": 400*10**3,"HRB500": 500*10**3,"PSB830": 830*10**3,"PSB980": 980*10**3}
        return retrunDict[tendonCode]
    ####################################################################################################################
    def auxiliary_steelElasticModulus_TB10092_2017(self,steelCode="wire"):
        """
        ----------------------------------------------------------------------------------------------------------------
        ---steel elastic modulus based on Code for Design of Concrete Structures of Railway Bridge andCulvert(TB 10092-2017)---
        Input:
            steel code(str), one of the grade "wire","wireStrand","prestressedBar","HPB300","HRB400","HRB500"
        OutPut:
            [elasticModulus(KPa)]
        ----------------------------------------------------------------------------------------------------------------
        """
        retrunDict = {"wire":2.05*10**5,"wireStrand":1.95*10**5,"prestressedBar":2.0*10**5,"HPB300":2.1*10**5,
                      "HRB400":2.0*10**5,"HRB500":2.0*10**5}
        return retrunDict[steelCode]
    ####################################################################################################################
    def auxiliary_fiberSectionPlot(self,sectTag, fillFlag=1):
        """
        ----------------------------------------------------------------------------------------------------------------
        Used for plotting fiber sections,the commands are not added in the OpenSeesPy
        ----------------------------------------------------------------------------------------------------------------
        1. Inspired by plotSection matlab function written by D. Vamvatsikos available at
           http://users.ntua.gr/divamva/software.html (plotSection.zip)
        2. Inspired by the module opsvis developed by sewkokot available at
           https://github.com/sewkokot/opsvis
        ----------------------------------------------------------------------------------------------------------------
        Inputs:
            sectTag(int)-the tag of the fiber section to be plotted
            fillFlag (int): 1 - filled fibers with color specified in matcolor
                            0 - no color, only the outline of fibers
        Output:
            plt-plot handle
        Notes:
            To use this auxiliary function, the section fiber should be constructed with the OpenSeesPyX command
        ----------------------------------------------------------------------------------------------------------------
        ###---Example:
                pltHandle=opsX.auxiliary_fiberSectionPlot(sectTag=110001)
                pltHandle.show()
        ----------------------------------------------------------------------------------------------------------------
        """
        fiberSectList=self.fiberSectDict[sectTag]
        fig, ax = plt.subplots()
        ax.set_xlabel('y')
        ax.set_ylabel('z')
        ax.grid(False)
        matTag=list(set([item[2] for item in fiberSectList][1:]))
        color_names = ["silver","green","orange","brown", "gray", "red", "yellow", "cyan", "magenta", "blue",
                        "purple", "turquoise",  "gold", "navy", "olive", "maroon", "beige"]
        matColorDict={matTag[i1]:color_names[i1] for i1 in range(len(matTag))}
        for item in fiberSectList:
            if item[0] == 'layer':
                matTag = item[2]
                if item[1] == 'straight':
                    n_bars = item[3]
                    As = item[4]
                    Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                    r = np.sqrt(As / np.pi)
                    Y = np.linspace(Iy, Jy, n_bars)
                    Z = np.linspace(Iz, Jz, n_bars)
                    [[bar := Circle((yi,zi), r, ec='k', fc='k', zorder=10),ax.add_patch(bar)] for zi, yi in zip(Z, Y)]

                if item[1] == 'circ':
                    n_bars, As = item[3], item[4]
                    yC, zC, arc_radius = item[5], item[6], item[7]
                    if len(item) > 8:
                        a0_deg, a1_deg = item[8], item[9]
                        if ((a1_deg - a0_deg) >= 360. and n_bars > 0):
                            a1_deg = a0_deg + 360. - 360. / n_bars
                    else:
                        a0_deg, a1_deg = 0., 360. - 360. / n_bars
                    a0_rad, a1_rad = np.pi * a0_deg / 180., np.pi * a1_deg / 180.
                    r_bar = np.sqrt(As / np.pi)
                    thetas = np.linspace(a0_rad, a1_rad, n_bars)
                    Y = yC + arc_radius * np.cos(thetas)
                    Z = zC + arc_radius * np.sin(thetas)
                    [[bar := Circle((yi,zi), r_bar, ec='k', fc='k', zorder=10),ax.add_patch(bar)] for zi, yi in zip(Z, Y)]

            if (item[0] == 'patch' and (item[1] == 'quad' or item[1] == 'quadr' or
                                        item[1] == 'rect')):
                matTag, nIJ, nJK = item[2], item[3], item[4]

                if item[1] == 'quad' or item[1] == 'quadr':
                    Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                    Ky, Kz, Ly, Lz = item[9], item[10], item[11], item[12]

                if item[1] == 'rect':
                    Iy, Iz, Ky, Kz = item[5], item[6], item[7], item[8]
                    Jy, Jz, Ly, Lz = Ky, Iz, Iy, Kz
                # check for convexity (vector products)
                outIJxIK = (Jy - Iy) * (Kz - Iz) - (Ky - Iy) * (Jz - Iz)
                outIKxIL = (Ky - Iy) * (Lz - Iz) - (Ly - Iy) * (Kz - Iz)
                # check if I, J, L points are colinear
                outIJxIL = (Jy - Iy) * (Lz - Iz) - (Ly - Iy) * (Jz - Iz)
                # outJKxJL = (Ky-Jy)*(Lz-Jz) - (Ly-Jy)*(Kz-Jz)
                if outIJxIK <= 0 or outIKxIL <= 0 or outIJxIL <= 0:
                    print('Warning! Patch quad is non-convex or counter-clockwise defined or has at least 3'
                          'colinear points in line')  # noqa: E501

                IJz, IJy = np.linspace(Iz, Jz, nIJ + 1), np.linspace(Iy, Jy, nIJ + 1)
                JKz, JKy = np.linspace(Jz, Kz, nJK + 1), np.linspace(Jy, Ky, nJK + 1)
                LKz, LKy = np.linspace(Lz, Kz, nIJ + 1), np.linspace(Ly, Ky, nIJ + 1)
                ILz, ILy = np.linspace(Iz, Lz, nJK + 1), np.linspace(Iy, Ly, nJK + 1)

                if fillFlag:
                    Z = np.zeros((nIJ + 1, nJK + 1))
                    Y = np.zeros((nIJ + 1, nJK + 1))
                    for j in range(nIJ + 1):
                        Z[j, :] = np.linspace(IJz[j], LKz[j], nJK + 1)
                        Y[j, :] = np.linspace(IJy[j], LKy[j], nJK + 1)
                    [[[zy := np.array([[Y[j, k],Z[j, k]],
                                           [Y[j, k + 1],Z[j, k + 1]],
                                           [Y[j + 1, k + 1],Z[j + 1, k + 1]],
                                           [Y[j + 1, k],Z[j + 1, k]]]),poly := Polygon(zy, closed=True, ec='k',
                                           fc=matColorDict[matTag]),
                       ax.add_patch(poly)] for k in range(nJK)] for j in range(nIJ)]
                else:
                    # horizontal lines
                    [plt.plot([ay, by],[az, bz], 'b-', zorder=1) for az, bz, ay, by in zip(IJz, LKz, IJy, LKy)]
                    # vertical lines
                    [plt.plot( [ay, by],[az, bz], 'b-', zorder=1) for az, bz, ay, by in zip(JKz, ILz, JKy, ILy)]

            if item[0] == 'patch' and item[1] == 'circ':
                matTag, nc, nr = item[2], item[3], item[4]
                yC, zC, ri, re = item[5], item[6], item[7], item[8]
                a0, a1 = item[9], item[10]
                dr = (re - ri) / nr
                dth = (a1 - a0) / nc
                [[rj := ri + j * dr,rj1 := rj + dr,[[thi:= a0 + i * dth,thi1:= thi + dth,wedge:= Wedge((yC,zC),
                    rj1, thi, thi1, width=dr, ec='k',lw=1, fc=matColorDict[matTag]),ax.add_patch(wedge)]
                                                    for i in range(nc)]] for j in range(nr)]
                ax.axis('equal')
        plt.axis('equal')
        plt.title(f"Section fibers-{sectTag}")
        return plt
    ####################################################################################################################
    def auxiliary_dicreteFiberPointsPlot(self,data_groups,colors,radius=2):
        """
        plot discrete fiber points (ylocal,zlocal,area)
        data_groups-fiber data, eg. [[[bartag1,y1,z1,A1],[bartag2,y2,z2,A2],...],[[core1,y1,z1,A1],[core2,y2,z2,A2],...],...]
        colors-color list, eg. ["red","green","blue"]
        "red", "green", "blue", "yellow", "orange","purple", "pink", "cyan", "magenta", "brown","lime", "navy", "teal",
        "olive", "maroon","gray", "grey", "black", "white", "gold"
        """
        visualizer = DataVisualizer()
        visualizer.show_scatter_plot(data_groups=data_groups,colors=colors,radius=radius,window_title="Data Visualization")
    ####################################################################################################################
    def auxiliary_writeModelInformationToDB(self):
        """
        When call this method all nodes,elements..., that need to be visualized by SAPStructure will be written to the
        result database
        """
        nodeDict={}
        for eachItem in self.saveNodeList:
            if eachItem[1] not in nodeDict.keys():
                nodeDict[eachItem[1]]=[]
                nodeDict[eachItem[1]].append(eachItem[0])
            else:
                nodeDict[eachItem[1]].append(eachItem[0])
        # for eachType in nodeDict.keys():
        #     saveName = f"modelInfo/{eachType}"
        #     saveValueList=nodeDict[eachType]
        #     headNameList=["nodeTag","xCoord","yCoord","zCoord"]
        #     operationIndexStr = 'replace'
        #     self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)
        [[saveName:= f"modelInfo/{eachType}",saveValueList:=nodeDict[eachType],
              headNameList:=["nodeTag","xCoord","yCoord","zCoord"],operationIndexStr:= 'append',
              self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)]
             for eachType in nodeDict.keys()]
        self.saveNodeList=[]###---empty the list after writing data to database
        ###############################################################################################################
        ###############################################################################################################
        geomDict={}
        for eachItem in self.saveGeomfList:
            if eachItem[1] not in geomDict.keys():
                geomDict[eachItem[1]] = []
                geomDict[eachItem[1]].append(eachItem[0])
            else:
                geomDict[eachItem[1]].append(eachItem[0])
        # for eachGeomType in geomDict.keys():
        #     saveName = f"geomTransf/{eachGeomType}"
        #     saveValueList=geomDict[eachGeomType]
        #     headNameList=["geomTransfTag","localZx","localZy","localZz"]
        #     operationIndexStr = 'replace'
        #     self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)
        [[saveName:= f"geomTransf/{eachGeomType}",saveValueList:=geomDict[eachGeomType],
          headNameList:=["geomTransfTag","localZx","localZy","localZz"],operationIndexStr:= 'append',
          self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)] for eachGeomType in geomDict.keys()]
        self.saveGeomfList=[]
        ###############################################################################################################
        ###############################################################################################################
        eleDict={}
        for eachItem in self.saveEleList:
            if eachItem[1] not in eleDict.keys():
                eleDict[eachItem[1]] = []
                eleDict[eachItem[1]].append(eachItem[0])
            else:
                eleDict[eachItem[1]].append(eachItem[0])
        # for eachEleType in eleDict.keys():
        #     saveName = f"modelInfo/{eachEleType}"
        #     saveValueList=eleDict[eachEleType]
        #     totalItems=len(eleDict[eachEleType][0])-2
        #     headNameList=["eleTag"]
        #     for i11 in range(totalItems):
        #         headNameList.append(f"node_{i11+1}")
        #     headNameList.append("eleDim")
        #     operationIndexStr = 'replace'
        #     self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)
        [[saveName:= f"modelInfo/{eachEleType}",saveValueList:=eleDict[eachEleType],
        totalItems:=len(eleDict[eachEleType][0])-2,headNameList:=["eleTag"],[headNameList.append(f"node_{i11+1}")
        for i11 in range(totalItems)],headNameList.append("eleDim"),operationIndexStr:= 'append',
            self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)] for eachEleType in eleDict.keys()]
        self.saveEleList=[]
        ###############################################################################################################
        ###############################################################################################################
        eleLocCoordDict={}
        for eachItem in self.EleLocalCoordSys:
            if eachItem[1] not in eleLocCoordDict.keys():
                eleLocCoordDict[eachItem[1]] = []
                eleLocCoordDict[eachItem[1]].append(eachItem[0])
            else:
                eleLocCoordDict[eachItem[1]].append(eachItem[0])
        # for eachLocCoord in eleLocCoordDict.keys():
        #     saveName = f"eleLocalCoordSys/{eachLocCoord}"
        #     saveValueList = eleLocCoordDict[eachLocCoord]
        #     headNameList=None
        #     if len(eleLocCoordDict[eachLocCoord][0])==4:
        #         headNameList = ["eleType", "nodeI", "nodeJ", "geomTransfTag"]
        #     elif len(eleLocCoordDict[eachLocCoord][0])==9:
        #         headNameList = ["eleType", "nodeI", "nodeJ", "localX_x", "localX_y", "localX_z","localY_x","localY_y","localY_z"]
        #     operationIndexStr = 'replace'
        #     self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)
        [[saveName:= f"eleLocalCoordSys/{eachLocCoord}",saveValueList:= eleLocCoordDict[eachLocCoord],headNameList:=None,
          length:= len(eleLocCoordDict[eachLocCoord][0]),headNameList:= {4: ["eleType", "nodeI", "nodeJ", "geomTransfTag"],
            9: ["eleType", "nodeI", "nodeJ","localX_x", "localX_y", "localX_z","localY_x", "localY_y", "localY_z"],}.get(length, []),
          operationIndexStr:= 'append',self.saveInstance.saveResult(saveName, saveValueList, headNameList, operationIndexStr)]
         for eachLocCoord in eleLocCoordDict.keys()]
        self.EleLocalCoordSys=[]
        ###############################################################################################################
        ###############################################################################################################
        #####---for node informaiton display
        nodeTags=ops.getNodeTags()
        # for eachNode in nodeTags:
        #     saveName=f"nodeInformationDisplay/{eachNode}"
        #     savedValueList = [[eachNode,str(ops.nodeCoord(eachNode)),str(ops.nodeMass(eachNode))]]
        #     headNameList = ["nodeTag","nodeCoord.","nodeMass"]
        #     operationIndexStr = 'replace'
        #     self.saveInstance.saveResult(saveName, savedValueList, headNameList, operationIndexStr)
        [[saveName:=f"nodeInformationDisplay/{eachNode}",savedValueList:= [[eachNode,str(ops.nodeCoord(eachNode)),
        str(ops.nodeMass(eachNode))]],headNameList:=["nodeTag","nodeCoord.","nodeMass"],
          operationIndexStr:= 'replace',self.saveInstance.saveResult(saveName, savedValueList, headNameList, operationIndexStr)
          ] for eachNode in nodeTags]
        #####-for element information display
        ############################################
        #############################################
        def key_from_list(lst):
            return tuple(sorted(lst))
        #############################################
        # geomTagInfoDict={}
        # for eachType in geomDict.keys():
        #     for eachValue in geomDict[eachType]:
        #         geomTagInfoDict[eachValue[0]]=[eachType]+eachValue
        geomTagInfoDict = {eachValue[0]: [eachType] + list(eachValue) for eachType, values in geomDict.items()
            for eachValue in values}
        #####################
        # eleTagInfoDict={}
        # for eachEleType in eleDict.keys():
        #     for eachEle in eleDict[eachEleType]:
        #         eleTagInfoDict[eachEle[0]]=[eachEleType]+eachEle[1:]
        eleTagInfoDict = {eachEle[0]: [eachEleType] + list(eachEle[1:]) for eachEleType, eles in eleDict.items()
            for eachEle in eles}
        ###################
        # nodesTagLocalCoordDict={}
        # for eachCordType in eleLocCoordDict.keys():
        #     for eachNodes in eleLocCoordDict[eachCordType]:
        #         nodesTagLocalCoordDict[str(eachNodes[1:3])]=[eachCordType]+eachNodes
        nodesTagLocalCoordDict = {key_from_list(nodes[1:3]): [cordType] + list(nodes) for cordType, nodes_list in eleLocCoordDict.items()
            for nodes in nodes_list}
        ###############################################################################################################
        ###############################################################################################################
        eleTags=ops.getEleTags()
        for eachEle in eleTags:
            try:
                saveName = f"elementInformationDisplay/{eachEle}"
                sortedNodes = key_from_list(ops.eleNodes(eachEle))
                eleType=eleTagInfoDict[eachEle][0]
                eleDim=eleTagInfoDict[eachEle][-1]
                try:
                    specialEleIndex=nodesTagLocalCoordDict[sortedNodes][1]
                except:
                    specialEleIndex=None
                localXYVector=None
                geomTransform=None
                localStr=None
                dirExpress=None
                if eleDim=='1D':
                    if specialEleIndex=='specialEle':
                        localXYVector0=nodesTagLocalCoordDict[sortedNodes][-6:]
                        localXYVector=[float(each) for each in localXYVector0]
                    elif specialEleIndex=='realEle':
                        geomTag=nodesTagLocalCoordDict[sortedNodes][-1]
                        geomTransform=geomTagInfoDict[geomTag]
                if localXYVector is not None:
                    localStr='eleLocalXYvector:'
                    dirExpress=localXYVector
                elif geomTransform is not None:
                    localStr='eleGeomTransf:'
                    dirExpress=geomTransform
                else:
                    localStr=None
                    dirExpress=None
                localDirString=f"{localStr}{dirExpress}"
                savedValueList = [[eachEle,str(sortedNodes),eleType,localDirString]]
                headNameList = ["eleTag","endNodes","eleType","localDirection"]
                operationIndexStr = 'replace'
                self.saveInstance.saveResult(saveName, savedValueList, headNameList, operationIndexStr)
            except:
                pass
    ####################################################################################################################
    def auxiliary_writeDataIntoTxtFile(self,savePath,filename,listData,decimals):
        """
        ----------------------------------------------------------------------------------------------------------------
        Write list data into txt file
        ----------------------------------------------------------------------------------------------------------------
        Inputs:
            savePath(str)-the path of the saved file, e.g. 'Period\1'
            filename(str)-the name of the file, e.g. 'period'
            listData(list[list])-the nested list to store the data, e.g. [[1,2,3,4],[0.12,0.22,0.333,0.412]]
                                the first list ([1,2,3,4]) will be saved as the first column in the file
            decimals(list[int])-to specify the decimal for each column data in the file, e.g. [0,2]
        ----------------------------------------------------------------------------------------------------------------
        """
        if not listData:
            raise ValueError("listData cannot be empty!")
        ncols = len(listData)
        nrows = len(listData[0])
        if any(len(c) != nrows for c in listData):
            raise ValueError("All columns should be the same size!")
        if isinstance(decimals, int):
            decimals = [decimals] * ncols
        if len(decimals) != ncols:
            raise ValueError("decimals must be a single integer or a list of integers the same as the number of columns!")
        fmts = [f"{{:.{d}f}}" for d in decimals]
        with open(f"{savePath}/{filename}.txt", 'w+', encoding='utf-8') as f:
            # for row in zip(*listData):
            #     line = ' '.join(fmt.format(val) for fmt, val in zip(fmts, row))
            #     f.write(line + '\n')
            [[line:= ' '.join(fmt.format(val) for fmt, val in zip(fmts, row)),f.write(line + '\n')]
             for row in zip(*listData)]
    ####################################################################################################################
    def node(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        node(nodeTag, *crds, '-ndf', ndf, '-mass', *mass, '-disp', *disp, '-vel', *vel, '-accel', *accel)
        Create a OpenSees node.
        nodeTag (int)	node tag.
        crds (list (float))	nodal coordinates.
        ndf (float)	nodal ndf. (optional)
        mass (list (float))	nodal mass. (optional)
        vel (list (float))	nodal velocities. (optional)
        accel (list (float))	nodal accelerations. (optional)
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.node{args}") ###---create node in OpenSeesPy
        nodeArgsList=list(args)
        nodeTagValue=nodeArgsList[0]
        coords = ops.nodeCoord(nodeTagValue)
        if len(coords)==1:
            coords=coords+[0.0,0.0]
        elif len(coords)==2:
            coords=coords+[0.0]
        else:
            pass
        tipValue = 'node_node'
        self.saveNodeList.append([[nodeTagValue]+coords,tipValue])
        self.nodeSetNameList.append(tipValue)
        # print("nodeNameList=",self.nodeSetNameList)


    def geomTransf(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        geomTransf(transfType, transfTag, *transfArgs)
        ----------------------------------------------
        The geometric-transformation command is used to construct a coordinate-transformation (CrdTransf) object,
        which transforms beam element stiffness and resisting force from the basic system to the global-coordinate
        system. The command has at least one argument, the transformation type.
        ----------------------------------------------
        transfType (str)	geomTransf type
        transfTag (int)	geomTransf tag.
        transfArgs (list)	a list of geomTransf arguments, must be preceded with *.
        ----------------------------------------------
        The following contain information about available transfType:
        Linear Transformation
        PDelta Transformation
        Corotational Transformation
        ----------------------------------------------
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.geomTransf{args}")  ###---create node in OpenSeesPy
        geomTransfArgsList = list(args)
        geomTransfType = geomTransfArgsList[0]
        geomTransfTag = geomTransfArgsList[1]
        tipsString =geomTransfType
        geomVect=[0.0,0.0,0.0]
        if geomTransfType in ["Linear","linear"]:
            if "-jntOffset" in geomTransfArgsList:
                indexValue=geomTransfArgsList.index("-jntOffset")
                if indexValue-1<=1:
                    geomVect=[0,0,1]
                else:
                    for i1 in range(indexValue-2):
                        geomVect[i1]=geomTransfArgsList[2+i1]
            else:
                totalNum=len(geomTransfArgsList)
                if totalNum==2:
                    geomVect=[0,0,1]
                else:
                    for i1 in range(totalNum-2):
                        geomVect[i1]=geomTransfArgsList[2+i1]
        elif geomTransfType in ["PDelta","pDelta"]:
            if "-jntOffset" in geomTransfArgsList:
                indexValue=geomTransfArgsList.index("-jntOffset")
                if indexValue-1<=1:
                    geomVect=[0,0,1]
                else:
                    for i1 in range(indexValue-2):
                        geomVect[i1]=geomTransfArgsList[2+i1]
            else:
                totalNum=len(geomTransfArgsList)
                if totalNum==2:
                    geomVect=[0,0,1]
                else:
                    for i1 in range(totalNum-2):
                        geomVect[i1]=geomTransfArgsList[2+i1]
        else: ###---Corotational
            if "-jntOffset" in geomTransfArgsList:
                indexValue = geomTransfArgsList.index("-jntOffset")
                if indexValue - 1 <= 1:
                    geomVect = [0, 0, 1]
            else:
                totalNum = len(geomTransfArgsList)
                if totalNum == 2:
                    geomVect = [0, 0, 1]
                else:
                    for i1 in range(totalNum - 2):
                        geomVect[i1] = geomTransfArgsList[2 + i1]
        saveGeomfList = []
        self.saveGeomfList.append([[geomTransfTag, geomVect[0],geomVect[1],geomVect[2]],tipsString + "_geomTransf"])
        self.localTransfNameList.append(tipsString + "_geomTransf")


    def equalDOF(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        equalDOF(rNodeTag, cNodeTag, *dofs)
        ----------------------------------------------
        Create a multi-point constraint between nodes.
        ----------------------------------------------
        rNodeTag (int)	integer tag identifying the retained, or master node.
        cNodeTag (int)	integer tag identifying the constrained, or slave node.
        dofs (list (int))	nodal degrees-of-freedom that are constrained at the cNode to be the same as those at
        the rNode Valid range is from 1 through ndf, the number of nodal degrees-of-freedom.
        ----------------------------------------------
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.equalDOF{args}")
        equalDOFArgsList = list(args)
        nodeI,nodeJ=equalDOFArgsList[0],equalDOFArgsList[1]
        nodeICoords = ops.nodeCoord(nodeI)
        nodeJCoords = ops.nodeCoord(nodeJ)
        equalLength=np.sqrt((nodeICoords[0]-nodeJCoords[0])**2+(nodeICoords[1]-nodeJCoords[1])**2+
                            (nodeICoords[2]-nodeJCoords[2])**2)
        tipValue = 'equalDOF'
        if equalLength<=1.0e-10:
            pass
        else:
            self.eleSetNameList.append(tipValue + "_ele")
            self.saveEleList.append([[0,nodeI,nodeJ],tipValue + "_ele"])



    def element(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        element(eleType, eleTag, *eleNodes, *eleArgs)
        Create a OpenSees element.
        eleType (str)	element type
        eleTag (int)	element tag.
        eleNodes (list (int))	a list of element nodes, must be preceded with *.
        eleArgs (list)	a list of element arguments, must be preceded with *.
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.element{args}")  ###---create node in OpenSeesPy
        argsList=list(args)
        eleType=argsList[0]
        tipValue =eleType
        eleTag=argsList[1]
        if eleType in ['elasticBeamColumn','ElasticBeamColumn']:
            returnValue=next(((index, item) for index, item in enumerate(argsList[1:]) if isinstance(item, str)), None)
            if returnValue==None:
                transTag=argsList[-1]
            else:
                transTag=argsList[returnValue[0]]
            eleNodes=ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle',eleNodes[0],eleNodes[1],transTag],tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['Truss','truss','TrussSection','trussSection','corotTruss','CorotTruss','corotTrussSection',
                         'CorotTrussSection']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['twoNodeLink','TwoNodeLink']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                    localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                    localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1]]+localXVector+localYVector,
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['Tri31','tri31']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellNLDKGQ','shellNLDKGQ']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['quad','Quad']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellMITC4','shellMITC4']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellDKGQ', 'shellDKGQ']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellDKGT', 'shellDKGT']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellNLDKGQ', 'sellNLDKGQ']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellNLDKGT', 'shellNLDKGT']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ShellNL', 'shellNL']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['bbarQuad', 'BbarQuad']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['enhancedQuad', 'EnhancedQuad']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SSPquad','sSPquad']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['MVLEM_3D','mVLEM_3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SFI_MVLEM_3D', 'sFI_MVLEM_3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['quadUP','QuadUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['bbarQuadUP','BbarQuadUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        # elif eleType in ['9_4_QuadUP']:
        #     tipValue="number_"+tipValue
        #     eleNodes = ops.eleNodes(eleTag)
        #     node1,node2,node3,node4=eleNodes[:4]
        #     eleNodes=[node1,node2,node3,node4]
        #     saveList = []
        #     saveList.append([eleTag] + eleNodes + ["2D"]) ####---去除最中心节点
        #     self.eleSetNameList.append(tipValue + "_ele")
        #     self.saveInstance.saveEles(elesSaveName=tipValue + "_ele", elesList=saveList)
        ################################################################################################################
        elif eleType in ['SSPquadUP','sSPquadUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SimpleContact2D','simpleContact2D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SimpleContact3D','simpleContact3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['BeamContact2D','beamContact2D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['BeamContact3D','beamContact3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['BeamEndContact3D','beamEndContact3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['VS3D4','vS3D4']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['AV3D4','aV3D4']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SurfaceLoad','surfaceLoad']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["2D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################



        ################################################################################################################
        elif eleType in ['FourNodeTetrahedron','fourNodeTetrahedron']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Tetrahedron"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['stdBrick','StdBrick']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['bbarBrick','BbarBrick']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['bbarBrick','BbarBrick']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SSPbrick','sSPbrick']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['brickUP','BrickUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['bbarBrickUP','BbarBrickUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SSPbrickUP','sSPbrickUP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['AC3D8','aC3D8']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ASI3D8','aSI3D8']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes + ["Brick"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################


        ################################################################################################################
        elif eleType in ['zeroLength', 'ZeroLength']:
            dirIndex = argsList.index('-dir')
            tempValue=[]
            for eachItem in argsList[(dirIndex+1):]:
                if isinstance(eachItem, int):
                    tempValue.append(eachItem)
                else:
                    break
            if eleTag not in self.zeroEleDirDict.keys():
                self.zeroEleDirDict[eleTag] = tempValue
            ###################################################
            eleNodes = ops.eleNodes(eleTag)
            localXVector = [1, 0, 0]
            localYVector = [0, 1, 0]
            localZVector = [0, 0, 1]
            if "-orient" in argsList:
                orientIndex = argsList.index('-orient')
                localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"] ,tipValue+ "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1]]+localXVector+localYVector,
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['zeroLengthND', 'ZeroLengthND']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector = [1, 0, 0]
            localYVector = [0, 1, 0]
            localZVector = [0, 0, 1]
            if "-orient" in argsList:
                orientIndex = argsList.index('-orient')
                localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['zeroLengthSection', 'ZeroLengthSection']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector = [1, 0, 0]
            localYVector = [0, 1, 0]
            localZVector = [0, 0, 1]
            if "-orient" in argsList:
                orientIndex = argsList.index('-orient')
                localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['CoupledZeroLength', 'coupledZeroLength']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['zeroLengthContact2D', 'ZeroLengthContact2D','zeroLengthContact3D','ZeroLengthContact3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['zeroLengthInterface2D', 'ZeroLengthInterface2D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['zeroLengthImpact3D', 'ZeroLengthImpact3D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        ################################################################################################################
        elif eleType in ['ModElasticBeam2d','modElasticBeam2d']:
            returnValue = next(((index, item) for index, item in enumerate(argsList[1:]) if isinstance(item, str)),None)
            if returnValue == None:
                transTag = argsList[-1]
            else:
                transTag = argsList[returnValue[0]]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['ElasticTimoshenkoBeam', 'elasticTimoshenkoBeam']:
            returnValue = next(((index, item) for index, item in enumerate(argsList[1:]) if isinstance(item, str)),None)
            if returnValue == None:
                transTag = argsList[-1]
            else:
                transTag = argsList[returnValue[0]]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],
                                          tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['dispBeamColumn', 'DispBeamColumn']:
            transTag=argsList[4]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],
                                          tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['forceBeamColumn', 'ForceBeamColumn']:
            transTag=argsList[4]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],
                                          tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['nonlinearBeamColumn', 'NonlinearBeamColumn']:
            returnValue = next(((index, item) for index, item in enumerate(argsList[1:]) if isinstance(item, str)),None)
            if returnValue == None:
                transTag = argsList[-1]
            else:
                transTag = argsList[returnValue[0]]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],
                                          tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['dispBeamColumnInt', 'DispBeamColumnInt']:
            transTag=argsList[6]
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['realEle', eleNodes[0], eleNodes[1], transTag],
                                          tipValue + "_eleLocCordSys"])
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
        ################################################################################################################
        elif eleType in ['MVLEM', 'mVLEM']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['SFI_MVLEM', 'sFI_MVLEM']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['beamColumnJoint','BeamColumnJoint']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['ElasticTubularJoint','elasticTubularJoint']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['Joint2D','joint2D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['Joint2D','joint2D']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['CatenaryCable','catenaryCable']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################



        ################################################################################################################
        elif eleType in ['elastomericBearingPlasticity','ElastomericBearingPlasticity']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]

            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['ElastomericBearingBoucWen','elastomericBearingBoucWen']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['flatSliderBearing','FlatSliderBearing']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['singleFPBearing','SingleFPBearing']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['TFP','tFP']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['TripleFrictionPendulum','tripleFrictionPendulum']:
            eleNodes = ops.eleNodes(eleTag)
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.eleSetNameList.append(tipValue + "_ele")
        ################################################################################################################
        elif eleType in ['multipleShearSpring','MultipleShearSpring']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['KikuchiBearing','kikuchiBearing']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['YamamotoBiaxialHDR','yamamotoBiaxialHDR']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['ElastomericX','elastomericX']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['LeadRubberX','leadRubberX']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['HDR','hDR']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['RJWatsonEqsBearing','rJWatsonEqsBearing']:
            eleNodes = ops.eleNodes(eleTag)
            localXVector=[1,0,0]
            localYVector=[0,1,0]
            localZVector=[0,0,1]
            nodeICoords = ops.nodeCoord(eleNodes[0])
            nodeJCoords = ops.nodeCoord(eleNodes[1])
            equalLength = np.sqrt((nodeICoords[0] - nodeJCoords[0]) ** 2 + (nodeICoords[1] - nodeJCoords[1]) ** 2 +
                                  (nodeICoords[2] - nodeJCoords[2]) ** 2)
            if equalLength<1.0e-10:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            else:
                if "-orient" in argsList:
                    orientIndex=argsList.index('-orient')
                    nextIndex=next((index for index, item in enumerate(argsList[int(orientIndex+1):])
                                    if isinstance(item, str)), None)
                    if nextIndex==None:
                        vectNum=len(argsList)-orientIndex-1
                        if vectNum==3:
                            localYVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localXVector=[nodeJCoords[0]-nodeICoords[0],nodeJCoords[1]-nodeICoords[1],
                            nodeJCoords[2]-nodeICoords[2]]
                        else:
                            localXVector=[argsList[orientIndex+1],argsList[orientIndex+2],argsList[orientIndex+3]]
                            localYVector=[argsList[orientIndex+4],argsList[orientIndex+5],argsList[orientIndex+6]]
                    else:
                        if nextIndex==3:
                            localYVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localXVector = [nodeJCoords[0] - nodeICoords[0], nodeJCoords[1] - nodeICoords[1],
                                            nodeJCoords[2] - nodeICoords[2]]
                        else:
                            localXVector = [argsList[orientIndex + 1], argsList[orientIndex + 2], argsList[orientIndex + 3]]
                            localYVector = [argsList[orientIndex + 4], argsList[orientIndex + 5], argsList[orientIndex + 6]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        elif eleType in ['FPBearingPTV','fPBearingPTV']:
            eleNodes = ops.eleNodes(eleTag)
            orientIndex=19
            localXVector = [argsList[orientIndex], argsList[orientIndex + 1], argsList[orientIndex + 2]]
            localYVector = [argsList[orientIndex + 3], argsList[orientIndex + 4], argsList[orientIndex + 5]]
            self.eleSetNameList.append(tipValue + "_ele")
            self.EleLocalCoordSysSetNameList.append(tipValue + "_eleLocCordSys")
            self.saveEleList.append([[eleTag] + eleNodes+["1D"],tipValue + "_ele"])
            self.EleLocalCoordSys.append([['specialEle', eleNodes[0], eleNodes[1], localXVector, localYVector],
                                          tipValue + "_eleLocCordSys"])
        ################################################################################################################
        else:
            print(123)
        ################################################################################################################

    def rigidLink(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        rigidLink(type, rNodeTag, cNodeTag)
        ----------------------------------------------
        Create a multi-point constraint between nodes.
        ----------------------------------------------
        type (str)
        string-based argument for rigid-link type:
        'bar': only the translational degree-of-freedom will be constrained to be exactly the same
        as those at the master node
        'beam': both the translational and rotational degrees of freedom are constrained.
        rNodeTag (int)	integer tag identifying the master node
        cNodeTag (int)	integar tag identifying the slave node
        ----------------------------------------------
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.rigidLink{args}")
        argsList = list(args)
        rigidType=argsList[0]
        rNode=argsList[1]
        cNode=argsList[2]
        rNodeCoords = ops.nodeCoord(rNode)
        cNodeCoords = ops.nodeCoord(cNode)
        equalLength = np.sqrt((rNodeCoords[0] -cNodeCoords[0]) ** 2 + (rNodeCoords[1] -cNodeCoords[1]) ** 2 +
                              (rNodeCoords[2] -cNodeCoords[2]) ** 2)
        tipValue = 'rigidLink'
        if equalLength <= 1.0e-10:
            pass
        else:
            self.eleSetNameList.append(tipValue + "_ele")
            self.saveEleList.append([[0, rNode, cNode]+["1D"],tipValue + "_ele"])
    ####################################################################################################################
    def section(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        section(secType, secTag, *secArgs)
        ----------------------------------------------
        Reference the useage in OpenSeesPy
        ----------------------------------------------
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.section{args}")
        self.currentSectTag=args[1]
        self.fiberSectDict[self.currentSectTag]=[]
        self.fiberSectDict[self.currentSectTag].append(['section']+list(args))
    ####################################################################################################################
    def patch(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        Reference the useage in OpenSeesPy
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.patch{args}")
        self.fiberSectDict[self.currentSectTag].append(['patch']+list(args))
    ####################################################################################################################
    def layer(self,*args):
        """
        ----------------------------------------------------------------------------------------------------------------
        Reference the useage in OpenSeesPy
        ----------------------------------------------------------------------------------------------------------------
        """
        eval(f"ops.layer{args}")
        self.fiberSectDict[self.currentSectTag].append(['layer']+list(args))
    ####################################################################################################################
    def integration_fiberSection(self,eleTag,fiberList,matTagList,GJValue=1.0e10):
        """
        ----------------------------------------------------------------------------------------------------------------
        Construct fiber seciton
        ------------------------------------------
        Inputs:
            eleTag(int)-the number of fiber section
            fiberList(list)-the fiber information list, eg. [[[yloc0_0,zloc0_0,A0_0],[yloc1_0,zloc01_0,A1_0],...],
                [[yloc0_1,zloc0_1,A0_1],[yloc1_1,zloc01_1,A1_1],...],...]
            matTagList(list)-material number list corresponding to each element in fiberList,eg.[1,2,...]
            GJValue-(float)-linear-elastic torsional stiffness assigned to the section (default value takes 1.0e10)
            tipsString(str)-print information on console
        ----------------------------------------------------------------------------------------------------------------
        """
        ops.section('Fiber', int(eleTag), '-GJ', GJValue)
        [[ops.fiber(eachItem[0], eachItem[1], eachItem[2], matTagList[i1]) for eachItem in fiberList[i1]]
         for i1 in range(len(fiberList))]
    ####################################################################################################################
    def _makeDirs(self,savePath):
        """
        ----------------------------------------------------------------------------------------------------------------
        Make directory if not exists,used only internally
        ------------------------------------------
        savePath(str)-the path of the directory
        ----------------------------------------------------------------------------------------------------------------
        """
        if os.path.exists(savePath):
            pass
        else:
            os.makedirs(savePath)

    ####################################################################################################################

    ####################################################################################################################
    def integration_recorderNode(self,savePath,filename,nodeLists,dofLists,responseType):
        """
        ----------------------------------------------------------------------------------------------------------------
        Records the response of a number of nodes at every converged step
        ------------------------------------------
        savePath(str)-the path of the directory,eg.'nodeDisp'
        fileName(str)-the name of the txt file, eg. 'case1'
        nodeLists(list)-nodes that need record responses, eg. [1,2,3,4]
        dofLists(list)-the specified dof at the nodes whose response is requested.eg. [1,2,3]
        responseType(str)-a string indicating response required
            including:
            'disp' displacement
            'vel' velocity
            'accel' acceleration
            'incrDisp' incremental displacement
            'reaction' nodal reaction
            'eigen i' eigenvector for mode i
            'rayleighForces' damping forces
        ----------------------------------------------------------------------------------------------------------------
        """
        self._makeDirs(savePath)
        fileName = savePath + '/' +filename+ '.txt'
        linkstr = f"ops.recorder('Node', '-file','{fileName}', '-time', '-node',"
        for each in nodeLists:
            linkstr+=f"{each}"+f","
        linkstr+=f"'-dof',"
        for each in dofLists:
            linkstr+=f"{each}"+f","
        linkstr+=f"'{responseType}')"
        eval(linkstr)
    ####################################################################################################################
    def integration_recorderElement(self,savePath,filename,eleList,responseTypeList):
        """
        ----------------------------------------------------------------------------------------------------------------
        Records the response of a number of elements at every converged step
        ------------------------------------------
        savePath(str)-the path of the directory,eg.'eleForce'
        fileName(str)-the name of the txt file, eg. 'case1'
        eleLists(list)-elements that need record responses, eg. [1,2,3,4]
        responseTypeList(list)-arguments which are passed to the setResponse()
            include:
            ['axialForce']-for truss element,1 column for each element
            ['section','1','force']-for nonlinear element force at integrationPoint 1, 4column for each element
            ['section', '1', 'deformation']-for nonlinear element deformation at integrationPoint 1,4column for each element
            ['localForce']-for elestic beamcolumn element and zerolength element force
            ['deformation']--for elestic beamcolumn element and zerolength element deformation
        ----------------------------------------------------------------------------------------------------------------
        """
        self._makeDirs(savePath)
        fileName = savePath + '/' +filename+ '.txt'
        linkstr = f"ops.recorder('Element', '-file','{fileName}', '-time', '-ele',"
        for each in eleList:
            linkstr+=f"{each}"+f","
        for i1 in range(len(responseTypeList)-1):
            linkstr += f"'{responseTypeList[i1]}'" + f","
        linkstr+=f"'{responseTypeList[-1]}'"+f")"
        eval(linkstr)
    ####################################################################################################################
    def integration_gravityLoad(self,nodesList):
        """
        ----------------------------------------------------------------------------------------------------------------
        Apply gravity load to associated nodes
        ------------------------------------------
        nodesList(list)-eg.[[[node1Tag,node1Mass],[],...],[],...]
        ----------------------------------------------------------------------------------------------------------------
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        for each in nodesList:
            ops.load(int(each[0]), 0.0, 0.0, -each[1] * 9.81, 0.0, 0.0, 0.0)
        ----------------------------------------------------------------------------------------------------------------
        """
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        [ops.load(int(each[0]), 0.0, 0.0, -each[1] * 9.81, 0.0, 0.0, 0.0) for each in nodesList]
    ####################################################################################################################
    def integration_analysisGravity(self,totalStep=1,recordList=None):
        """
        ----------------------------------------------------------------------------------------------------------------
        Static gravity load analysis
        ------------------------------------------
        totalStep(int)-Number of analysis steps to perform
        recordList(list)-responses need to be recordedd
        ----------------------------------------------------------------------------------------------------------------
        """
        ################################################################################################################
        def saveDataIntoOpenedDB(f,dataName,resultsList,headNameList,dirList=None):
            """
            Save results to database, resultsList=[[result0_0,result0_1,],[],...]
            headNameList=[headName_0,headName_1,...]
            dataName(str)
            operationIndexStr='replace' or 'append' 'replace'
            """
            if len(resultsList) > 0:
                # t1=time.time()
                list0 = resultsList[0]
                saveTypeList = []
                typeDict = {"int": "np.int32", "float": "np.float32", "str": "h5py.string_dtype(encoding='utf-8')"}
                saveTypeList = [typeDict["int"] if isinstance(eachValue, (int, np.int64, np.uint64)) else
                                typeDict["float"] if isinstance(eachValue, (float, np.float64)) else
                                typeDict["str"] if isinstance(eachValue, str) else None for eachValue in list0]

                dtypeStr = "np.dtype(["
                dtypeStr += ''.join([f"('{headNameList[i1]}',{saveTypeList[i1]})," for i1 in range(len(list0) - 1)])
                dtypeStr += f"('{headNameList[-1]}',{saveTypeList[-1]})])"
                dtype = eval(dtypeStr)
                structured_data = np.zeros(len(resultsList), dtype=dtype)
                for i2 in range(len(headNameList)):
                    structured_data[headNameList[i2]] = [each[i2] for each in resultsList]
                ########################################################################################################
                # t2=time.time()
                # print("time1=",t2-t1)
                dataset = f.get(dataName)
                # t3 = time.time()
                if dataset is None:
                    dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype, chunks=True,compression=None)
                    new_size = len(resultsList)
                    dataset.resize(new_size, axis=0)
                    dataset[0:new_size] = structured_data
                    if dirList is not None:
                        dataset.attrs['dirList'] =dirList
                else:
                    original_size = dataset.shape[0]
                    new_size = original_size + len(resultsList)
                    dataset.resize(new_size, axis=0)
                    dataset[original_size:new_size] = structured_data
                # t4 = time.time()
                # print("time2=", t4 - t3)
        ################################################################################################################
        ################################################################################################################
        f=None
        if recordList is None:
            totalStep = totalStep
            ops.system('UmfPack')
            ops.constraints('Transformation')
            ops.numberer('RCM')
            ops.test('NormDispIncr', 1.0e-3, 200)
            ops.algorithm('KrylovNewton')
            ops.integrator('LoadControl', 1.0 / float(totalStep))
            ops.analysis('Static')
            ops.analyze(int(totalStep))
            ops.loadConst('-time', 0.0)  ##-在后续分析步考虑重力作用(将重力视常荷载),并将拟时间置0
        else:
            totalStep = totalStep
            ops.system('UmfPack')
            ops.constraints('Transformation')
            ops.numberer('RCM')
            ops.test('NormDispIncr', 1.0e-3, 200)
            ops.algorithm('KrylovNewton')
            ###########################################################################################################
            dbPath = self.saveInstance._dbPath
            f = h5py.File(dbPath, 'a', libver='latest')
            ##############################################
            nodeDict = {}
            trussEleResponseDict = {}
            zeroEleResponseDict = {}
            zeroEleDirectionDict = {}
            nonEleSectResponsesDict = {}
            nonEleSectNumberDict = {}
            nonZeroEleResponsesDict = {}
            for each in recordList:
                if each[0] == 'node':
                    nodeIdenty, resType, nodeTags = each[0], each[1], each[2]
                    nodeItemDict = {(nodeIdenty + '_' + resType + '_' + str(eachNode)): [] for eachNode in nodeTags}
                    nodeDict = {**nodeDict, **nodeItemDict}  ##Merge two dicts
                elif each[0] == 'trussEle':
                    responseType, eleTags = each[1], each[2]
                    eleItemDict = {('trussEle_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    trussEleResponseDict = {**trussEleResponseDict, **eleItemDict}  ##Merge two dicts
                elif each[0] == 'zeroEle':
                    responseType, eleTags = each[1], each[2]
                    eleItemDict = {('zeroEle_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    zeroEleResponseDict = {**zeroEleResponseDict, **eleItemDict}  ##Merge two dicts
                    eleDirectDict = {('zeroEle_' + responseType + '_' + str(eachEle)): self.zeroEleDirDict[eachEle] for
                                     eachEle in eleTags}
                    zeroEleDirectionDict = {**zeroEleDirectionDict, **eleDirectDict}  ##Merge two dicts
                elif each[0] == 'nonEleSection':
                    responseType, sectNum, eleTags = each[1], each[2], each[3]
                    eleItemDict = {('nonEleSection_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    nonEleSectResponsesDict = {**nonEleSectResponsesDict, **eleItemDict}  ##Merge two dicts
                    sectNumDict = {('nonEleSection_' + responseType + '_' + str(eachEle)): sectNum for eachEle in
                                   eleTags}
                    nonEleSectNumberDict = {**nonEleSectNumberDict, **sectNumDict}  ##Merge two dicts
                elif each[0] == 'nonZeroEle':
                    responseType, eleTags = each[1], each[2]
                    eleItemDict = {('nonZeroEle_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    nonZeroEleResponsesDict = {**nonZeroEleResponsesDict, **eleItemDict}  ##Merge two dicts
        ####################---recorderProcess---###########
        if recordList != None:
            for i1 in range(totalStep):
                ops.integrator('LoadControl', 1.0 / float(totalStep))
                ops.analysis('Static')
                ops.analyze(1)
                tCurrent = ops.getTime()
                ######################################################
                if nodeDict:
                    nodeKeys = nodeDict.keys()
                    nodeResNameDict = {'disp': 'nodeDisp', 'vel': 'nodeVel', 'accel': 'nodeAccel',
                                       'reaction': 'nodeReaction'}
                    # for eachkey in nodeKeys:
                    #     resType= eachkey.split('_')[1]
                    #     nodeTag= eachkey.split('_')[2]
                    #     tempValue1=[round(tCurrent, 4)]
                    #     tempValue2= eval(f"ops.{nodeResNameDict[resType]}({nodeTag})")
                    #     tempValue= tempValue1 + tempValue2
                    #     nodeDict['node_' + resType + '_' + str(nodeTag)].append(tempValue)
                    [[resType:= eachkey.split('_')[1],nodeTag:= eachkey.split('_')[2],tempValue1:=[round(tCurrent, 4)],
                      tempValue2:= eval(f"ops.{nodeResNameDict[resType]}({nodeTag})"),tempValue:= tempValue1 + tempValue2,
                      nodeDict['node_' + resType + '_' + str(nodeTag)].append(tempValue)] for eachkey in nodeKeys]
                    ########################################
                    for eachkey in nodeKeys:
                        resType= eachkey.split("_")[1]
                        nodeTag= eachkey.split("_")[2]
                        saveValueList= nodeDict['node_' + resType + '_' + str(nodeTag)]
                        saveName = f"staticResponse_node/node_{resType}_{nodeTag}"
                        numArgs=len(saveValueList[0])-1
                        headNameList = ["pseudotime"]
                        for i2 in range(numArgs):
                            headNameList.append(f"value_{i2 + 1}")
                        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                        nodeDict['node_' + resType + '_' + str(nodeTag)]=[]
                #######################################################################################################
                #######################################################################################################
                if trussEleResponseDict:
                    eleKeys = trussEleResponseDict.keys()
                    eleResNameDict = {'axialForce': 'basicForce', 'axialDeform': 'basicDeformation'}
                    # for eachkey in eleKeys:
                    #     resType= eachkey.split("_")[1]
                    #     eleTag= eachkey.split("_")[2]
                    #     tempValue1 = [round(tCurrent, 4)]
                    #     tempValue2 = [eval(f"ops.{eleResNameDict[resType]}({eleTag})[0]")]
                    #     tempValue= tempValue1 + tempValue2
                    #     trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)].append(tempValue)
                    [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue1:= [round(tCurrent, 4)],
                      tempValue2:= [eval(f"ops.{eleResNameDict[resType]}({eleTag})[0]")],tempValue:= tempValue1 + tempValue2,
                      trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)].append(tempValue)] for eachkey in eleKeys]
                    ########################################################################################################
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        saveValueList= trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)]
                        saveName = f"staticResponse_element/trussEle_{resType}_{eleTag}"
                        numArgs = len(saveValueList[0]) - 1
                        headNameList = ["pseudotime"]
                        for i2 in range(numArgs):
                            headNameList.append(f"value_{i2 + 1}")
                        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                        trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)] = []
                ########################################################################################################
                ########################################################################################################
                if zeroEleResponseDict:
                    eleKeys = zeroEleResponseDict.keys()
                    # for eachkey in eleKeys:
                    #     resType= eachkey.split("_")[1]
                    #     eleTag= eachkey.split("_")[2]
                    #     tempValue1 = [round(tCurrent, 4)]
                    #     tempValue2= eval(f"ops.eleResponse({eleTag},'{resType}')")
                    #     tempValue= [tempValue1 + tempValue2] + [zeroEleDirectionDict[eachkey]]
                    #     zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)].append(tempValue)
                    [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue1:= [round(tCurrent, 4)],
                      tempValue2:= eval(f"ops.eleResponse({eleTag},'{resType}')"),
                      tempValue:= [tempValue1 + tempValue2] + [zeroEleDirectionDict[eachkey]],
                      zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)].append(tempValue)
                      ] for eachkey in eleKeys]
                    ################################################
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        saveValueList= [zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)][0][0]]
                        dirList= zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)][0][1]
                        saveName = f"staticResponse_element/zeroEle_{resType}_{eleTag}"
                        numArgs=len(saveValueList[0])-1
                        headNameList = ["pseudotime"]
                        for i2 in range(numArgs):
                            headNameList.append(f"value_{i2 + 1}")
                        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList,dirList)
                        zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)]=[]
                ########################################################################################################
                ########################################################################################################
                if nonEleSectResponsesDict:
                    eleKeys = nonEleSectResponsesDict.keys()
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        tempValue1 = [round(tCurrent, 4)]
                        tempValue = tempValue1 + [eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},1)"),
                            eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},2)"),
                            eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},3)"),
                            eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},4)")]
                        nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)].append(tempValue)
                    ##################################################
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        saveValueList= nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)]
                        saveName = f"staticResponse_element/nonEleSection_{resType}_{eleTag}"
                        numArgs = len(saveValueList[0]) - 1
                        headNameList = ["pseudotime"]
                        for i2 in range(numArgs):
                            headNameList.append(f"value_{i2 + 1}")
                        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                        nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)] = []
                ############################################################################################################
                ############################################################################################################
                if nonZeroEleResponsesDict:
                    eleKeys = nonZeroEleResponsesDict.keys()
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        tempValue1 = [round(tCurrent, 4)]
                        tempValue2= eval(f"ops.eleResponse({eleTag},'{resType}')")
                        tempValue= tempValue1 + tempValue2
                        nonZeroEleResponsesDict[eachkey].append(tempValue)
                    #########################################
                    for eachkey in eleKeys:
                        resType= eachkey.split("_")[1]
                        eleTag= eachkey.split("_")[2]
                        saveValueList= nonZeroEleResponsesDict[eachkey]
                        saveName = f"staticResponse_element/nonZeroEle_{resType}_{eleTag}"
                        numArgs=len(saveValueList[0])-1
                        headNameList = ["pseudotime"]
                        for i2 in range(numArgs):
                            headNameList.append(f"value_{i2 + 1}")
                        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                        nonZeroEleResponsesDict[eachkey]=[]
            ############################################################################################################
            ############################################################################################################
            ops.loadConst('-time', 0.0)  ##-在后续分析步考虑重力作用(将重力视常荷载),并将拟时间置0
            ############################################################################################################
        if f is not None:
            f.close()

    ####################################################################################################################
    def integration_analysisModal(self,numModes=10):
        """
        ----------------------------------------------------------------------------------------------------------------
        Modal analysis
        ------------------------------------------
        numModes(int)-number of eigenvalues required
        ----------------------------------------------------------------------------------------------------------------
        """
        ################################################################################################################
        def saveDataIntoOpenedDB(f,dataName,resultsList,headNameList,operationIndexStr='replace'):
            """
            Save results to database, resultsList=[[result0_0,result0_1,],[],...]
            headNameList=[headName_0,headName_1,...]
            dataName(str)
            operationIndexStr='replace' or 'append' 'replace'
            """
            if len(resultsList) > 0:
                list0 = resultsList[0]
                saveTypeList = []
                typeDict = {"int": "np.int32", "float": "np.float32", "str": "h5py.string_dtype(encoding='utf-8')"}
                saveTypeList = [typeDict["int"] if isinstance(eachValue, (int, np.int64, np.uint64)) else
                                typeDict["float"] if isinstance(eachValue, (float, np.float64)) else
                                typeDict["str"] if isinstance(eachValue, str) else None for eachValue in list0]

                dtypeStr = "np.dtype(["
                dtypeStr += ''.join([f"('{headNameList[i1]}',{saveTypeList[i1]})," for i1 in range(len(list0) - 1)])
                dtypeStr += f"('{headNameList[-1]}',{saveTypeList[-1]})])"
                dtype = eval(dtypeStr)
                structured_data = np.zeros(len(resultsList), dtype=dtype)
                for i2 in range(len(headNameList)):
                    structured_data[headNameList[i2]] = [each[i2] for each in resultsList]
                #################################################################################
                dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype, chunks=True,
                                           compression='gzip', compression_opts=9, shuffle=True)
                new_size = len(resultsList)
                dataset.resize(new_size, axis=0)
                dataset[0:new_size] = structured_data
            f.flush()
        ################################################################################################################
        dbPath=self.saveInstance._dbPath
        f=h5py.File(dbPath, 'a')
        allNodesTag = ops.getNodeTags()
        eigenValues = ops.eigen(numModes)
        ################################################################################################################
        for eachMode in range(numModes):
            saveList= []
            for eachNode in allNodesTag:
                nodeEigenValue= ops.nodeEigenvector(eachNode, int(eachMode + 1))
                saveList.append([eachNode] + nodeEigenValue)
            ########################
            dofNum=len(saveList[0])-1
            saveName = f"modalInfo/modal_mode_{eachMode+1}"
            saveValueList=saveList
            headNameList=["nodeTag"]
            for i2 in range(dofNum):
                headNameList.append(f"dof_{i2+1}")
            operationIndexStr = 'replace'
            saveDataIntoOpenedDB(f,saveName,saveValueList,headNameList,operationIndexStr)
        ########################
        savePeridList = []
        for i3 in range(numModes):
            periodT= 2.0 * 3.1415926 / float(eigenValues[i3] ** 0.5)
            savePeridList.append([i3 + 1, periodT])
        saveName = f"periodInfo/period"
        saveValueList =savePeridList
        headNameList = ["modeNum","period"]
        saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList, operationIndexStr)
        ###########################################################################
        f.close()
        ################################################################################################################
        if numModes>=10:
            for i2 in range(10):
                print(str(i2 + 1) + ' th period is: ' + str(savePeridList[i2]) + ' second')
        else:
            for i2 in range(numModes):
                print(str(i2 + 1) + ' th period is: ' + str(savePeridList[i2]) + ' second')
    ####################################################################################################################
    def integration_analysisModalProperties(self,numEigen,pflag=1,outname=None):
        """
        ----------------------------------------------------------------------------------------------------------------
            Details
            -------
                This script will return the modal properties of an OpenSeespy model.

            Information
            -----------
                Author: Volkan Ozsarac, Earthquake Engineering PhD Candidate
                Affiliation: University School for Advanced Studies IUSS Pavia
                e-mail: volkanozsarac@iusspavia.it

            References
            ----------
                Chopra, A.K. 2012. Dynamics of Structures: Theory and
                Applications to Earthquake Engineering, Prentice Hall.

            Notes
            -----
                Total (activated) mass is obtained by summing the masses assigned to the
                unrestrained degrees of freedoms (DOFs). Thus, it should not be confused
                with total mass assigned to all DOFs. Influence vectors for rotational
                excitation are not correct at the moment, this addition remains as future work.
                Which reference point to use is not clear for rotational excitations.
                SAP2000 and Seismostruct use different reference points.

            Parameters
            ----------
            numEigen : int
                Number of eigenvalues to calculate.
            pflag    : int (1 or 0)
                flag to print output information on screen
            outname  : str, optional (The default is None)
                if not None and pFlag==1, the modal properties for the
                first numEigen modes will be writtend into outname.csv.

            Returns
            -------
            T        : numpy.ndarray
                Period array for the first numEigen modes.
            Mratios  : dictionary
                Effective modal mass participation ratios for the first numEigen modes.
            Mfactors : dictionary
                Modal particpation factors for the first numEigen modes.
            Mtots    : dictionary
                Total activated masses.
        ----------------------------------------------------------------------------------------------------------------
            """
        import numpy as np
        import openseespy.opensees as op
        import sys

        op.wipeAnalysis()
        op.numberer("Plain")
        op.system('FullGeneral')
        op.algorithm('Linear')
        op.analysis('Transient')

        # Extract the Mass Matrix
        # Note that this is not the global mass matrix, but unrestrained part (Muu)
        op.integrator('GimmeMCK', 1.0, 0.0, 0.0)
        op.analyze(1, 0.0)
        # Number of equations in the model
        N = op.systemSize()  # Has to be done after analyze
        Mmatrix = op.printA('-ret')  # Or use op.printA('-file','M.out')
        Mmatrix = np.array(Mmatrix)  # Convert the list to an array
        Mmatrix.shape = (N, N)  # Make the array an NxN matrix
        print('\n************************************************************', \
              '\nExtracting the mass matrix, ignore the warnings...')

        # Determine maximum number of DOFs/node used in the system
        NDF = 0
        for node in op.getNodeTags():
            temp = len(op.nodeDOFs(node))
            if temp > NDF: NDF = temp

        DOFs = []  # List containing indices of unrestrained DOFs
        used = {}  # Dictionary with nodes and associated unrestrained DOFs
        ldict = {}  # Dictionary containing influence vectors
        Mratios = {}  # Dictionary containing effective modal masses ratios
        Mfactors = {}  # Dictionary containing modal participation factors
        for i in range(1, NDF + 1):
            ldict[i] = np.zeros([N, 1])
            Mratios[i] = np.zeros(numEigen)
            Mfactors[i] = np.zeros(numEigen)

        # Create the influence vectors, and get the unrestrained DOFs assigned to the nodes
        # TODO -1: The influence vectors are not correct in case of rotational excitations
        # One typical approach is to use center of mass on plane
        idx = 0  # Counter for unrestrained DOFs
        for node in op.getNodeTags():  # Start iterating over each node
            used[node] = []  # Unrestrain local DOF ids
            ndof = len(op.nodeDOFs(node))  # Total number of DOFs assigned
            for j in range(ndof):  # Iterate over each DOF
                temp = op.nodeDOFs(node)[j]  # Get the global DOF id (-1 if restrained)
                if temp not in DOFs and temp >= 0:  # Check if this DOF is unrestrained and is not known before
                    DOFs.append(temp)  # Save the global id of DOF
                    used[node].append(j + 1)  # Save the local id of DOF
                    ldict[j + 1][idx, 0] = 1  # Influence vectors for horizontal and vertical excitations
                    idx += 1  # Increase the counter

        # This does not seem necessary when numberer is "Plain"
        # But lets reorganize the mass matrix anyway
        Mmatrix = Mmatrix[DOFs, :][:, DOFs]

        # Calculate the total masses assigned to the unrestrained DOFs
        Mtots = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        for i in range(1, NDF + 1):
            Mtots[i] = (ldict[i].T @ Mmatrix @ ldict[i])[0, 0]

        # Perform eigenvalue analysis
        op.wipeAnalysis()
        listSolvers = ['-genBandArpack', '-fullGenLapack', '-symmBandLapack']
        ok = 1
        for s in listSolvers:
            print("Using %s as solver..." % s[1:])
            try:
                eigenValues = op.eigen(s, numEigen)
                catchOK = 0
                ok = 0
            except:
                catchOK = 1

            if catchOK == 0:
                for i in range(numEigen):
                    if eigenValues[i] < 0:
                        ok = 1
                if ok == 0:
                    print('Eigenvalue analysis is completed.')
                    break
        if ok != 0:
            print("Error on eigenvalue something is wrong...")
            sys.exit()
        else:
            Lambda = np.asarray(eigenValues)
            Omega = Lambda ** 0.5
            T = 2 * np.pi / Omega
            frq = 1 / T

        # Note: influence factors for rotational excitation is wrong!
        # Obtain modal properties
        for mode in range(1, numEigen + 1):
            idx = 0
            phi = np.zeros([N, 1])  # Eigen vector
            for node in used:
                for dof in used[node]:
                    phi[idx, 0] = op.nodeEigenvector(node, mode, dof)
                    idx += 1

            phi = phi / (phi.T @ Mmatrix @ phi) ** 0.5  # Normalize the eigen vector by modal mass
            Mn = phi.T @ Mmatrix @ phi  # Modal mass (should always be equal to 1)

            for j in range(1, NDF + 1):
                if Mtots[j] != 0:  # Check if any mass is assigned
                    Ln = phi.T @ Mmatrix @ ldict[j]  # Modal excitation factor
                    Mnstar = (Ln ** 2 / Mn)[0, 0]  # Effective modal mass
                    Mfactors[j][mode - 1] = Ln / Mn  # Modal participation factor
                    Mratios[j][mode - 1] = (Mnstar / Mtots[j] * 100)  # Effective modal mass participation ratio [%]

        for j in range(1, 7):
            try:
                Mratios[j]
            except:
                Mratios[j] = np.zeros(numEigen)
                Mfactors[j] = np.zeros(numEigen)

        # TODO-1: Results are not correct for rotational excitation cases, for now ignore those.
        del Mratios[6], Mratios[5], Mratios[4]
        del Mfactors[6], Mfactors[5], Mfactors[4]

        # Calculate cumulative modal mass participation ratio
        sM1 = np.cumsum(Mratios[1]);
        sM2 = np.cumsum(Mratios[2]);
        sM3 = np.cumsum(Mratios[3])

        # Print modal analysis results
        if pflag == 1:
            arguments = []
            arguments.append('Modal Periods and Frequencies')
            arguments.append('%4s|%8s|%10s|%12s|%12s' \
                             % ('Mode', 'T [sec]', 'f [Hz]', '\u03C9 [rad/sec]', '\u03BB [rad\u00b2/sec\u00b2]'))
            for mode in range(numEigen):
                arguments.append('%4s|%8s|%10s|%12s|%12s' \
                                 % ("{:.0f}".format(mode + 1), "{:.4f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                                    "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))
            arguments.append('Total Activated Masses')
            arguments.append('%8s|%8s|%8s' \
                             % ('M\u2081', 'M\u2082', 'M\u2083'))
            arguments.append('%8s|%8s|%8s' \
                             % ("{:.2f}".format(Mtots[1]), "{:.2f}".format(Mtots[2]), "{:.2f}".format(Mtots[3])))
            arguments.append('Modal Mass Participation Factors')
            arguments.append('%4s|%7s|%7s|%7s' \
                             % ('Mode', '\u0393\u2081', '\u0393\u2082', '\u0393\u2083'))
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode + 1), \
                                                      "{:.3f}".format(Mfactors[1][mode]),
                                                      "{:.3f}".format(Mfactors[2][mode]),
                                                      "{:.3f}".format(Mfactors[3][mode])))
            arguments.append('Effective Modal Mass Participation Ratios [%]')
            arguments.append('%4s|%7s|%7s|%7s' \
                             % ('Mode', 'U\u2081', 'U\u2082', 'U\u2083'))
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode + 1), \
                                                      "{:.3f}".format(Mratios[1][mode]),
                                                      "{:.3f}".format(Mratios[2][mode]),
                                                      "{:.3f}".format(Mratios[3][mode])))
            arguments.append('Cumulative Effective Modal Mass Participation Ratios [%]')
            arguments.append('%4s|%7s|%7s|%7s' \
                             % ('Mode', '\u2211U\u2081', '\u2211U\u2082', '\u2211U\u2083'))
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode + 1), \
                                                      "{:.3f}".format(sM1[mode]), "{:.3f}".format(sM2[mode]),
                                                      "{:.3f}".format(sM3[mode])))

                # To the screen
            arguments = '\n'.join(arguments);
            print(arguments)

            # To the .csv file
            if outname != None:
                with open(outname + '.csv', 'w', encoding='utf-32') as f:
                    f.write(arguments)

        return T, Mratios, Mfactors, Mtots
    ####################################################################################################################
    def integration_earthquakeExcite(self,RayleighDamping,waveLenth,dt,dirList,motionList,factor=9.81,recordList=None,waveNumber=None):
        """
        ----------------------------------------------------------------------------------------------------------------
        Apply a uniform excitation to a model acting in a certain direction
        ------------------------------------------
        RayleighDamping(list)-set the values of Rayleigh damping, with two options,
            option 1: RayleighDamping=['mode-1',dampingRatio,Tstart,Tend]
                    dampRatio(float)-the damping ratio for the structure,eg.0.05
                    Tstart,Tend(float)-the start and end periods for calculating rayleigh damping
            option 2: RayleighDamping=['mode-2',α,β1,β2,β3]
                    D=α×M＋β1×Kcurrent＋β2×Kinit＋β3×KlastCommit
        waveLenth(int)-the length of the ground motion
        dt(float)-the time interval of the motion
        dirList(list)-direction in which ground motion acts,eg. [1,3]
            1 corresponds to translation along the global X axis
            2 corresponds to translation along the global Y axis
            3 corresponds to translation along the global Z axis
            4 corresponds to rotation about the global X axis
            5 corresponds to rotation about the global Y axis
            6 corresponds to rotation about the global Z axis
        motionList(list)-grond motion paths corresponding to the dirList,eg.[path_acc_X,path_acc_Z]
        factor(float)-a value used to scale the ground motion time history
        waveNumber-ground motion number, used for print informaiton on screen
        ----------------------------------------------------------------------------------------------------------------
        """
        if RayleighDamping[0]=='mode-1':
            dampRatio,Tstart,Tend=RayleighDamping[1:]
            w1=2.0*np.pi/float(Tstart)
            w2=2.0*np.pi/float(Tend)
            a = dampRatio * 2.0 * w1 * w2 / float(w1 + w2)
            b = dampRatio * 2 / float(w1 +w2)
            ### D=α×M＋β1×Kcurrent＋β2×Kinit＋β3×KlastCommit Longitudinal direction
            ops.rayleigh(a, 0.0, 0.0,b)
            print('Rayleigh damping: ',a,0.0,0.0,b)
        elif RayleighDamping[0]=='mode-2':
            a,b1,b2,b3=RayleighDamping[1:]
            ops.rayleigh(a,b1,b2,b3)
            print('Rayleigh damping: ',a,b1,b2,b3)
        else:
            pass
        ops.loadConst('-time', 0.0)
        currentLength = waveLenth
        currentDt = dt
        dir_L, dir_T, dir_V = 1, 2, 3
        gmFact =factor
        for i1 in range(len(dirList)):
            ops.timeSeries('Path', int(i1+100), '-dt', currentDt, '-filePath',motionList[i1], '-factor', gmFact)
            ops.pattern('UniformExcitation', int(i1+1000), int(dirList[i1]), '-accel', int(i1+100))
        ######################################################
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        # ops.test('NormDispIncr', tol,maxNumIter)
        # ops.algorithm('KrylovNewton')
        # ops.integrator('Newmark', 0.5, 0.25)
        # ops.analysis('Transient')
        ################################################################################################################
        ################################################################################################################
        def saveDataIntoOpenedDB(f,dataName,resultsList,headNameList,dirList=None):
            """
            Save results to database, resultsList=[[result0_0,result0_1,],[],...]
            headNameList=[headName_0,headName_1,...]
            dataName(str)
            operationIndexStr='replace' or 'append' 'replace'
            """
            # if len(resultsList) > 0:
            #     # t1=time.time()
            #     list0 = resultsList[0]
            #     saveTypeList = []
            #     typeDict = {"int": "np.int32", "float": "np.float32", "str": "h5py.string_dtype(encoding='utf-8')"}
            #     saveTypeList = [typeDict["int"] if isinstance(eachValue, (int, np.int64, np.uint64)) else
            #                     typeDict["float"] if isinstance(eachValue, (float, np.float64)) else
            #                     typeDict["str"] if isinstance(eachValue, str) else None for eachValue in list0]
            #
            #     dtypeStr = "np.dtype(["
            #     dtypeStr += ''.join([f"('{headNameList[i1]}',{saveTypeList[i1]})," for i1 in range(len(list0) - 1)])
            #     dtypeStr += f"('{headNameList[-1]}',{saveTypeList[-1]})])"
            #     dtype = eval(dtypeStr)
            #     structured_data = np.zeros(len(resultsList), dtype=dtype)
            #     for i2 in range(len(headNameList)):
            #         structured_data[headNameList[i2]] = [each[i2] for each in resultsList]
            #     ########################################################################################################
            #     # t2=time.time()
            #     # print("time1=",t2-t1)
            #     dataset = f.get(dataName)
            #     # t3 = time.time()
            #     if dataset is None:
            #         dataset = f.create_dataset(dataName, shape=(0,), maxshape=(None,), dtype=dtype, chunks=True,compression=None)
            #         new_size = len(resultsList)
            #         dataset.resize(new_size, axis=0)
            #         dataset[0:new_size] = structured_data
            #         if dirList is not None:
            #             dataset.attrs['dirList'] =dirList
            #     else:
            #         original_size = dataset.shape[0]
            #         new_size = original_size + len(resultsList)
            #         dataset.resize(new_size, axis=0)
            #         dataset[original_size:new_size] = structured_data
            #     # t4 = time.time()
            #     # print("time2=", t4 - t3)
            ##################################################
            if len(resultsList) > 0:
                arr=np.array(resultsList).astype(np.float32,order='C')
                dataset = f.get(dataName)
                n_cols=len(resultsList[0])
                if dataset is None:
                    dataset = f.create_dataset(dataName,shape=(0, n_cols),maxshape=(None, n_cols),dtype=np.float32)
                    dataset.resize(len(resultsList), axis=0)
                    dataset[0:dataset.shape[0], :] = arr
                    if dirList is not None:
                        dataset.attrs['dirList'] =dirList
                else:
                    original_size = dataset.shape[0]
                    new_size = original_size + len(resultsList)
                    dataset.resize(new_size, axis=0)
                    dataset[original_size:new_size] =arr
        ################################################################################################################
        ################################################################################################################
        writeInterNum=200 ###---每200时间步将结果写入数据库一次
        f=None
        if recordList!=None:
            ######################################################################
            dbPath = self.saveInstance._dbPath
            f = h5py.File(dbPath, 'a', libver='latest')
            ######################################################################
            nodeDict={}
            trussEleResponseDict={}
            zeroEleResponseDict={}
            zeroEleDirectionDict={}
            nonEleSectResponsesDict={}
            nonEleSectNumberDict={}
            nonZeroEleResponsesDict={}
            for each in recordList:
                if each[0]=='node':
                    nodeIdenty, resType,nodeTags= each[0], each[1], each[2]
                    nodeItemDict={(nodeIdenty+'_'+resType+'_'+str(eachNode)):[] for eachNode in nodeTags}
                    nodeDict={**nodeDict,**nodeItemDict}##Merge two dicts
                elif each[0]=='trussEle':
                    responseType,eleTags =each[1],each[2]
                    eleItemDict={('trussEle_'+responseType+'_'+str(eachEle)):[] for eachEle in eleTags}
                    trussEleResponseDict = {**trussEleResponseDict, **eleItemDict}  ##Merge two dicts
                elif each[0]=='zeroEle':
                    responseType, eleTags = each[1], each[2]
                    eleItemDict = {('zeroEle_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    zeroEleResponseDict = {**zeroEleResponseDict, **eleItemDict}  ##Merge two dicts
                    eleDirectDict = {('zeroEle_' + responseType + '_' + str(eachEle)):self.zeroEleDirDict[eachEle] for eachEle in eleTags}
                    zeroEleDirectionDict = {**zeroEleDirectionDict, **eleDirectDict}  ##Merge two dicts
                elif each[0]=='nonEleSection':
                    responseType,sectNum,eleTags=each[1],each[2],each[3]
                    eleItemDict = {('nonEleSection_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    nonEleSectResponsesDict = {**nonEleSectResponsesDict, **eleItemDict}  ##Merge two dicts
                    sectNumDict = {('nonEleSection_' + responseType + '_' + str(eachEle)): sectNum for eachEle in eleTags}
                    nonEleSectNumberDict = {**nonEleSectNumberDict, **sectNumDict}  ##Merge two dicts
                elif each[0]=='nonZeroEle':
                    responseType,eleTags = each[1], each[2]
                    eleItemDict = {('nonZeroEle_' + responseType + '_' + str(eachEle)): [] for eachEle in eleTags}
                    nonZeroEleResponsesDict = {**nonZeroEleResponsesDict, **eleItemDict}  ##Merge two dicts
        ######################################################
        ######################################################
        startTime = time.perf_counter()
        tCurrent = ops.getTime()
        tFinal = currentLength * currentDt
        timeList = []
        maxNumIter=1000
        tol=1.0e-4
        deltaTList=[]
        deltaT=currentDt
        while (tCurrent < tFinal):
            deltaTList.append(deltaT)
            if deltaT<currentDt:
                if len(set(deltaTList[-400:]))==1:
                    deltaT=deltaT*2
                    deltaTList.append(deltaT)
                    print("The initial dt is:", currentDt, ", and the increase current dt is:", deltaT)
                if deltaT>=currentDt:
                    deltaT=currentDt
            ####['NormDispIncr', 'RelativeEnergyIncr', 'EnergyIncr', 'RelativeNormUnbalance',
            ##### 'RelativeNormDispIncr', 'NormUnbalance']
            ops.test('NormDispIncr', tol,maxNumIter)
            #####['KrylovNewton',, ['SecantNewton','-initial'], ['ModifiedNewton','-initial'],
            ##### ['RaphsonNewton','-initial'], 'PeriodicNewton','BFGS', 'Broyden', 'NewtonLineSearch'],前四个zengjia,'-initial'
            ops.algorithm('ModifiedNewton','-initial')  ##收敛性好于前两个
            # ops.algorithm('KrylovNewton')
            NewmarkGamma = 0.5
            NewmarkBeta = 0.25
            ops.integrator('Newmark', NewmarkGamma, NewmarkBeta)
            ops.analysis('Transient')
            ok = ops.analyze(1, deltaT)
            if (ok == 0):
                tCurrent = ops.getTime()
                timeList.append(tCurrent)
                endTime = time.perf_counter()
                realTime = endTime - startTime
                ########################################################################################################
                ####################---recorderProcess---###############################################################
                if recordList!=None:
                    if nodeDict:
                        nodeKeys=nodeDict.keys()
                        nodeResNameDict={'disp':'nodeDisp','vel':'nodeVel','accel':'nodeAccel','reaction':'nodeReaction'}
                        if (len(nodeDict[list(nodeKeys)[0]])>=writeInterNum) or (tCurrent>=tFinal):
                            # for eachkey in nodeKeys:
                            #     resType= eachkey.split("_")[1]
                            #     nodeTag= eachkey.split("_")[2]
                            #     saveValueList= nodeDict['node_' + resType + '_' + str(nodeTag)]
                            #     saveName = f"timeHistoryResponse_node/node_{resType}_{nodeTag}"
                            #     # saveName = f"modalInfo/modal_mode_{eachMode + 1}"
                            #     numArgs=len(saveValueList[0])-1
                            #     headNameList = ["time"]
                            #     for i2 in range(numArgs):
                            #         headNameList.append(f"value_{i2 + 1}")
                            #     saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                            [[resType:= eachkey.split("_")[1],nodeTag:= eachkey.split("_")[2],
                              saveValueList:= nodeDict['node_' + resType + '_' + str(nodeTag)],
                              saveName:= f"timeHistoryResponse_node/node_{resType}_{nodeTag}",
                              numArgs:=len(saveValueList[0])-1,headNameList:= ["time"],
                              [headNameList.append(f"value_{i2 + 1}") for i2 in range(numArgs)],operationIndexStr:= 'append',
                            saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)] for eachkey in nodeKeys]
                            ############################################################################################
                            f.flush()
                            for eachkey in nodeKeys:
                                nodeDict[eachkey] = [] #### return to null list
                        ################################################################################################
                        # for eachkey in nodeKeys:
                        #     resType= eachkey.split('_')[1]
                        #     nodeTag= eachkey.split('_')[2]
                        #     tempValue1= [round(tCurrent, 4)]
                        #     tempValue2= eval(f"ops.{nodeResNameDict[resType]}({nodeTag})")
                        #     tempValue3= [tempValue2[i1] for i1 in range(len(tempValue2))]
                        #     tempValue= tempValue1 + tempValue3
                        #     nodeDict['node_' + resType + '_' + str(nodeTag)].append(tempValue)
                        [[resType:= eachkey.split('_')[1],nodeTag:= eachkey.split('_')[2],tempValue1:= [round(tCurrent, 4)],
                          tempValue2:= eval(f"ops.{nodeResNameDict[resType]}({nodeTag})"),
                          tempValue3:= [tempValue2[i1] for i1 in range(len(tempValue2))],tempValue:= tempValue1 + tempValue3,
                          nodeDict['node_' + resType + '_' + str(nodeTag)].append(tempValue)] for eachkey in nodeKeys]
                    ###################################################################################################
                    ###################################################################################################
                    ###################################################################################################
                    if trussEleResponseDict:
                        eleKeys = trussEleResponseDict.keys()
                        eleResNameDict = {'axialForce': 'basicForce','axialDeform':'basicDeformation'}
                        if (len(trussEleResponseDict[list(eleKeys)[0]])>=writeInterNum) or (tCurrent>=tFinal):
                            # for eachkey in eleKeys:
                            #     resType= eachkey.split("_")[1]
                            #     eleTag= eachkey.split("_")[2]
                            #     saveValueList= trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)]
                            #     saveName = f"timeHistoryResponse_element/trussEle_{resType}_{eleTag}"
                            #     numArgs=len(saveValueList[0])-1
                            #     headNameList = ["time"]
                            #     for i2 in range(numArgs):
                            #         headNameList.append(f"value_{i2 + 1}")
                            #     saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                            [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],
                              saveValueList:= trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)],
                              saveName:= f"timeHistoryResponse_element/trussEle_{resType}_{eleTag}",
                              numArgs:=len(saveValueList[0])-1,headNameList:= ["time"],
                              [headNameList.append(f"value_{i2 + 1}") for i2 in range(numArgs)],
                              operationIndexStr:= 'append',
                              saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)]
                             for eachkey in eleKeys]
                            #############################################
                            f.flush()
                            for eachkey in eleKeys:
                                trussEleResponseDict[eachkey] = []
                            #############################################
                        # for eachkey in eleKeys:
                        #     resType= eachkey.split("_")[1]
                        #     eleTag= eachkey.split("_")[2]
                        #     tempValue1= [round(tCurrent, 4)]
                        #     tempValue2= [eval(f"ops.{eleResNameDict[resType]}({eleTag})[0]")]
                        #     tempValue= tempValue1 + tempValue2
                        #     trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)].append(tempValue)
                        [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue1:= [round(tCurrent, 4)],
                          tempValue2:= [eval(f"ops.{eleResNameDict[resType]}({eleTag})[0]")],tempValue:= tempValue1 + tempValue2,
                        trussEleResponseDict['trussEle_' + resType + '_' + str(eleTag)].append(tempValue)] for eachkey in eleKeys]
                    #####################################################################################################
                    ###################################################################################################
                    ###################################################################################################
                    if zeroEleResponseDict:
                        eleKeys = zeroEleResponseDict.keys()
                        if (len(zeroEleResponseDict[list(eleKeys)[0]])>=writeInterNum) or (tCurrent>=tFinal):
                            # for eachkey in eleKeys:
                            #     resType= eachkey.split("_")[1]
                            #     eleTag= eachkey.split("_")[2]
                            #     saveValueListSS= zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)]
                            #     dirValueList= saveValueListSS[0][1]
                            #     saveName = f"timeHistoryResponse_element/zeroEle_{resType}_{eleTag}"
                            #     saveValueList=[each[0] for each in saveValueListSS]
                            #     numArgs=len(saveValueList[0])-1
                            #     headNameList = ["time"]
                            #     for i2 in range(numArgs):
                            #         headNameList.append(f"value_{i2 + 1}")
                            #     saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList,dirValueList)
                            [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],
                              saveValueListSS:= zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)],
                              dirValueList:= saveValueListSS[0][1],
                              saveName:= f"timeHistoryResponse_element/zeroEle_{resType}_{eleTag}",
                              saveValueList:=[each[0] for each in saveValueListSS],
                              numArgs:=len(saveValueList[0])-1,headNameList:= ["time"],
                              [headNameList.append(f"value_{i2 + 1}") for i2 in range(numArgs)],
                              saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList,dirValueList)] for eachkey in eleKeys]
                            #################################################
                            f.flush()
                            for eachkey in eleKeys:
                                zeroEleResponseDict[eachkey] = []
                            ################################################
                        # for eachkey in eleKeys:
                        #     resType= eachkey.split("_")[1]
                        #     eleTag= eachkey.split("_")[2]
                        #     tempValue1= [round(tCurrent, 4)]
                        #     tempValue2= eval(f"ops.eleResponse({eleTag},'{resType}')")
                        #     tempValue3= tempValue2
                        #     tempValue = [tempValue1 + tempValue3] + [zeroEleDirectionDict[eachkey]]
                        #     zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)].append(tempValue)
                        [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue1:= [round(tCurrent, 4)],
                          tempValue2:= eval(f"ops.eleResponse({eleTag},'{resType}')"),tempValue3:= tempValue2,
                          tempValue:= [tempValue1 + tempValue3] + [zeroEleDirectionDict[eachkey]],
                          zeroEleResponseDict['zeroEle_' + resType + '_' + str(eleTag)].append(tempValue)] for eachkey in eleKeys]
                    ###################################################################################################
                    ###################################################################################################
                    ###################################################################################################
                    if nonEleSectResponsesDict:
                        eleKeys = nonEleSectResponsesDict.keys()
                        digitNumDict = {'sectionForce':6,'sectionDeformation':10}
                        if (len(nonEleSectResponsesDict[list(eleKeys)[0]])>=writeInterNum) or (tCurrent>=tFinal):
                            # for eachkey in eleKeys:
                            #     resType= eachkey.split("_")[1]
                            #     eleTag= eachkey.split("_")[2]
                            #     saveValueList= nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)]
                            #     saveName = f"timeHistoryResponse_element/nonEleSection_{resType}_{eleTag}"
                            #     numArgs=len(saveValueList[0])-1
                            #     headNameList = ["time"]
                            #     for i2 in range(numArgs):
                            #         headNameList.append(f"value_{i2 + 1}")
                            #     saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                            [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],
                              saveValueList:= nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)],
                              saveName:= f"timeHistoryResponse_element/nonEleSection_{resType}_{eleTag}",
                              numArgs:=len(saveValueList[0])-1,headNameList:= ["time"],
                              [headNameList.append(f"value_{i2 + 1}") for i2 in range(numArgs)],
                              saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)] for eachkey in eleKeys]

                            ###########################################################################
                            f.flush()
                            for eachkey in eleKeys:
                                nonEleSectResponsesDict[eachkey] = []
                            ########################################################################
                        # for eachkey in eleKeys:
                        #     resType= eachkey.split("_")[1]
                        #     eleTag= eachkey.split("_")[2]
                        #     tempValue= [round(tCurrent, 4)] + [
                        #         eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},1)"),
                        #         eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},2)"),
                        #         eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},3)"),
                        #         eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},4)")]
                        #     nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)].append(tempValue)
                        [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue:= [round(tCurrent, 4)]
                            + [eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},1)"),
                                eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},2)"),
                                eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},3)"),
                                eval(f"ops.{resType}({eleTag},{nonEleSectNumberDict[eachkey]},4)")],
                          nonEleSectResponsesDict['nonEleSection_' + resType + '_' + str(eleTag)].append(tempValue)] for eachkey in eleKeys]
                    ###################################################################################################
                    ###################################################################################################
                    ###################################################################################################
                    if nonZeroEleResponsesDict:
                        eleKeys = nonZeroEleResponsesDict.keys()
                        if (len(nonZeroEleResponsesDict[list(eleKeys)[0]])>=writeInterNum) or (tCurrent>=tFinal):
                            # for eachkey in eleKeys:
                            #     resType= eachkey.split("_")[1]
                            #     eleTag= eachkey.split("_")[2]
                            #     saveValueList= nonZeroEleResponsesDict[eachkey]
                            #     saveName = f"timeHistoryResponse_element/nonZeroEle_{resType}_{eleTag}"
                            #     numArgs=len(saveValueList[0])-1
                            #     headNameList = ["time"]
                            #     for i2 in range(numArgs):
                            #         headNameList.append(f"value_{i2 + 1}")
                            #     saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)
                            [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],
                              saveValueList:= nonZeroEleResponsesDict[eachkey],
                              saveName:= f"timeHistoryResponse_element/nonZeroEle_{resType}_{eleTag}",
                              numArgs:=len(saveValueList[0])-1,headNameList:= ["time"],
                              [headNameList.append(f"value_{i2 + 1}") for i2 in range(numArgs)],
                              [saveDataIntoOpenedDB(f, saveName, saveValueList, headNameList)]] for eachkey in eleKeys]
                            #####################################################
                            f.flush()
                            for eachkey in eleKeys:
                                nonZeroEleResponsesDict[eachkey] = []
                            ######################################################
                        # for eachkey in eleKeys:
                        #     resType= eachkey.split("_")[1]
                        #     eleTag= eachkey.split("_")[2]
                        #     tempValue1= [round(tCurrent, 4)]
                        #     tempValue2= eval(f"ops.eleResponse({eleTag},'{resType}')")
                        #     tempValue3= tempValue2
                        #     tempValue= tempValue1 + tempValue3
                        #     nonZeroEleResponsesDict[eachkey].append(tempValue)
                        [[resType:= eachkey.split("_")[1],eleTag:= eachkey.split("_")[2],tempValue1:= [round(tCurrent, 4)],
                          tempValue2:= eval(f"ops.eleResponse({eleTag},'{resType}')"),tempValue3:= tempValue2,
                          tempValue:= tempValue1 + tempValue3,nonZeroEleResponsesDict[eachkey].append(tempValue)] for eachkey in eleKeys]
                ########################################################################################################
                ########################################################################################################
                ########################################################################################################
                print('KrylovNewton',f'ground motion={waveNumber}','tol=',tol,'maxNumIter=',maxNumIter, 'totalTime=',
                      tFinal, 'tCurrent=',"{:.6f}".format(tCurrent),'time cost=', "{:.1f}".format(realTime), 'second')
            else:
                if deltaT>=0.000001*currentDt:
                    deltaT=deltaT*0.5
                    print("The initial dt is:",currentDt,", and the decrease current dt is:",deltaT)
                else:
                    print(f"ground moition {waveNumber}"," failed!")
                    break
            ############################################################################################################
        if f is not None:
            f.close()
########################################################################################################################
########################################################################################################################
class ScatterPlotWidget(QWidget):
    def __init__(self, data_groups, colors=None, radius=10, parent=None):
        """
        data_groups: list of lists, each sublist contains [id, y, z, area]
        colors: list of colors for each group (optional)
        radius: circle radius in pixels (all points same size)
        """
        super().__init__(parent)
        self.data_groups = data_groups
        self.radius = radius
        self.selected_point = None
        self.setMinimumSize(800, 600)
        self.setMouseTracking(True)

        # Flatten data for range calculation
        self.all_data = []
        self.group_indices = []
        for group_idx, group in enumerate(data_groups):
            for item in group:
                self.all_data.append(item)
                self.group_indices.append(group_idx)

        # Extract coordinates and areas
        self.y_coords = [item[1] for item in self.all_data]
        self.z_coords = [item[2] for item in self.all_data]
        self.areas = [item[3] for item in self.all_data]

        # Calculate data range
        self.y_min = min(self.y_coords)
        self.y_max = max(self.y_coords)
        self.z_min = min(self.z_coords)
        self.z_max = max(self.z_coords)

        # Calculate aspect ratio
        self.data_width = self.y_max - self.y_min if self.y_max != self.y_min else 1
        self.data_height = self.z_max - self.z_min if self.z_max != self.z_min else 1
        self.data_aspect_ratio = self.data_width / self.data_height

        # Store original color names
        self.color_names = []

        # Generate colors if not provided
        if colors is None or len(colors) < len(data_groups):
            self.colors = self.generate_colors(len(data_groups))
            self.color_names = [f"Auto Color {i + 1}" for i in range(len(data_groups))]
        else:
            self.colors = []
            for c in colors:
                parsed_color = self.parse_color(c)
                self.colors.append(parsed_color)
                # Store color name
                if isinstance(c, str):
                    self.color_names.append(c)
                elif isinstance(c, tuple):
                    self.color_names.append(f"RGB{c}")
                else:
                    self.color_names.append("Custom Color")

        # Margins
        self.margin_left = 80
        self.margin_right = 40
        self.margin_top = 40
        self.margin_bottom = 80

    def generate_colors(self, n):
        """Generate n distinct colors"""
        colors = []
        for i in range(n):
            hue = (i * 360 / n) % 360
            color = QColor.fromHsv(int(hue), 200, 220, 255)
            colors.append(color)
        return colors

    def parse_color(self, color):
        """Parse color from various formats"""
        # Predefined color mapping
        color_map = {
            "red": (255, 0, 0),
            "green": (0, 255, 0),
            "blue": (0, 0, 255),
            "yellow": (255, 255, 0),
            "orange": (255, 165, 0),
            "purple": (128, 0, 128),
            "pink": (255, 192, 203),
            "cyan": (0, 255, 255),
            "magenta": (255, 0, 255),
            "brown": (165, 42, 42),
            "lime": (0, 255, 0),
            "navy": (0, 0, 128),
            "teal": (0, 128, 128),
            "olive": (128, 128, 0),
            "maroon": (128, 0, 0),
            "gray": (128, 128, 128),
            "grey": (128, 128, 128),
            "black": (0, 0, 0),
            "white": (255, 255, 255),
            "gold": (255, 215, 0)
        }

        if isinstance(color, QColor):
            return color
        elif isinstance(color, str):
            # Check if it's a predefined color name
            color_lower = color.lower()
            if color_lower in color_map:
                rgb = color_map[color_lower]
                return QColor(rgb[0], rgb[1], rgb[2], 255)
            else:
                # Try to parse as hex or Qt color name
                c = QColor(color)
                if c.isValid():
                    c.setAlpha(255)
                    return c
                else:
                    # Default color if invalid
                    return QColor(100, 100, 200, 255)
        elif isinstance(color, tuple) and len(color) >= 3:
            if len(color) == 3:
                return QColor(color[0], color[1], color[2], 255)
            else:
                return QColor(color[0], color[1], color[2], color[3])
        else:
            return QColor(100, 100, 200, 255)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        # Background
        painter.fillRect(self.rect(), QColor(255, 255, 255))

        # Draw axes
        self.draw_axes(painter)

        # Draw grid
        self.draw_grid(painter)

        # Draw points
        self.draw_points(painter)

    def get_plot_rect(self):
        """Calculate plot area based on data aspect ratio"""
        available_width = self.width() - self.margin_left - self.margin_right
        available_height = self.height() - self.margin_top - self.margin_bottom
        available_aspect_ratio = available_width / available_height

        if self.data_aspect_ratio > available_aspect_ratio:
            plot_width = available_width
            plot_height = available_width / self.data_aspect_ratio
        else:
            plot_height = available_height
            plot_width = available_height * self.data_aspect_ratio

        x_offset = self.margin_left + (available_width - plot_width) / 2
        y_offset = self.margin_top + (available_height - plot_height) / 2

        return x_offset, y_offset, plot_width, plot_height

    def draw_grid(self, painter):
        """Draw grid lines"""
        x_offset, y_offset, plot_width, plot_height = self.get_plot_rect()

        pen = QPen(QColor(230, 230, 230), 1, Qt.DashLine)
        painter.setPen(pen)

        # Vertical grid lines
        for i in range(1, 10):
            x = x_offset + plot_width * i / 10
            painter.drawLine(int(x), int(y_offset),
                             int(x), int(y_offset + plot_height))

        # Horizontal grid lines
        for i in range(1, 10):
            y = y_offset + plot_height * i / 10
            painter.drawLine(int(x_offset), int(y),
                             int(x_offset + plot_width), int(y))

    def draw_axes(self, painter):
        """Draw coordinate axes"""
        x_offset, y_offset, plot_width, plot_height = self.get_plot_rect()

        pen = QPen(QColor(0, 0, 0), 2)
        painter.setPen(pen)

        # Y axis (left)
        painter.drawLine(int(x_offset), int(y_offset),
                         int(x_offset), int(y_offset + plot_height))

        # Z axis (bottom)
        painter.drawLine(int(x_offset), int(y_offset + plot_height),
                         int(x_offset + plot_width), int(y_offset + plot_height))

        # Draw ticks and labels
        font = QFont("Arial", 9)
        painter.setFont(font)

        # Z axis label (vertical) - closer to axis
        painter.save()
        painter.translate(25, int(y_offset + plot_height / 2))
        painter.rotate(-90)
        painter.drawText(0, 0, "Z")
        painter.restore()

        # Y axis label (horizontal) - closer to axis
        painter.drawText(int(x_offset + plot_width / 2 - 5),
                         self.height() - 25, "Y")

        # Z axis ticks (left side) - 2 decimal places
        num_ticks = 6
        for i in range(num_ticks):
            z_val = self.z_min + (self.z_max - self.z_min) * i / (num_ticks - 1)
            y_pos = y_offset + plot_height - (plot_height * i / (num_ticks - 1))
            painter.drawLine(int(x_offset - 5), int(y_pos),
                             int(x_offset), int(y_pos))
            painter.drawText(int(x_offset - 70), int(y_pos + 5), f"{z_val:.2f}")

        # Y axis ticks (bottom) - 2 decimal places
        for i in range(num_ticks):
            y_val = self.y_min + (self.y_max - self.y_min) * i / (num_ticks - 1)
            x_pos = x_offset + (plot_width * i / (num_ticks - 1))
            painter.drawLine(int(x_pos), int(y_offset + plot_height),
                             int(x_pos), int(y_offset + plot_height + 5))
            painter.drawText(int(x_pos - 20), int(y_offset + plot_height + 20),
                             f"{y_val:.2f}")

    def draw_points(self, painter):
        """Draw data points"""
        x_offset, y_offset, plot_width, plot_height = self.get_plot_rect()

        for i, item in enumerate(self.all_data):
            idx, y, z, area = item
            group_idx = self.group_indices[i]

            # Convert to screen coordinates
            x_screen = x_offset + ((y - self.y_min) / self.data_width) * plot_width
            y_screen = y_offset + plot_height - ((z - self.z_min) / self.data_height) * plot_height

            # Get color for this group
            color = self.colors[group_idx]

            # Draw circle with same color for fill and border (all same size)
            brush = QBrush(color)
            painter.setBrush(brush)
            pen = QPen(color, 1.5)
            painter.setPen(pen)
            painter.drawEllipse(QPointF(x_screen, y_screen), self.radius, self.radius)

    def mousePressEvent(self, event):
        """Mouse click event"""
        if event.button() == Qt.LeftButton:
            self.check_point_click(event.pos())

    def mouseMoveEvent(self, event):
        """Mouse move event"""
        self.check_point_hover(event.pos())

    def check_point_click(self, pos):
        """Check if a point was clicked"""
        x_offset, y_offset, plot_width, plot_height = self.get_plot_rect()

        for i, item in enumerate(self.all_data):
            idx, y, z, area = item
            group_idx = self.group_indices[i]

            # Convert to screen coordinates
            x_screen = x_offset + ((y - self.y_min) / self.data_width) * plot_width
            y_screen = y_offset + plot_height - ((z - self.z_min) / self.data_height) * plot_height

            # Calculate distance
            distance = math.sqrt((pos.x() - x_screen) ** 2 + (pos.y() - y_screen) ** 2)

            if distance <= self.radius:
                self.selected_point = i
                self.show_tooltip(pos, item, group_idx)
                return

    def check_point_hover(self, pos):
        """Check mouse hover"""
        x_offset, y_offset, plot_width, plot_height = self.get_plot_rect()

        for i, item in enumerate(self.all_data):
            idx, y, z, area = item

            # Convert to screen coordinates
            x_screen = x_offset + ((y - self.y_min) / self.data_width) * plot_width
            y_screen = y_offset + plot_height - ((z - self.z_min) / self.data_height) * plot_height

            # Calculate distance
            distance = math.sqrt((pos.x() - x_screen) ** 2 + (pos.y() - y_screen) ** 2)

            if distance <= self.radius:
                self.setCursor(Qt.PointingHandCursor)
                return

        self.setCursor(Qt.ArrowCursor)

    def show_tooltip(self, pos, item, group_idx):
        """Show tooltip - coordinates with 6 decimals, area with 9 decimals, and color name"""
        idx, y, z, area = item
        color_name = self.color_names[group_idx] if group_idx < len(self.color_names) else "Unknown"
        tooltip_text = f"ID: {idx}\nY: {y:.6f}\nZ: {z:.6f}\nArea: {area:.9f}\nColor: {color_name}"
        QToolTip.showText(self.mapToGlobal(pos), tooltip_text, self)

    def render_to_painter(self, painter, width, height):
        """Render the plot to a given painter with specified dimensions"""
        # Temporarily adjust widget size for rendering
        old_size = self.size()
        self.resize(width, height)

        # Draw on painter
        painter.fillRect(0, 0, width, height, QColor(255, 255, 255))
        self.draw_axes(painter)
        self.draw_grid(painter)
        self.draw_points(painter)

        # Restore original size
        self.resize(old_size)

    def export_to_png(self, filename, width=1920, height=1440, dpi=300):
        """Export plot to PNG file"""
        # Create image with specified size
        image = QImage(width, height, QImage.Format_RGB32)
        image.fill(QColor(255, 255, 255))

        # Set DPI
        dpi_to_dpm = dpi / 0.0254  # Convert DPI to dots per meter
        image.setDotsPerMeterX(int(dpi_to_dpm))
        image.setDotsPerMeterY(int(dpi_to_dpm))

        # Create painter for the image
        painter = QPainterExport(image)
        painter.setRenderHint(QPainter.Antialiasing)

        # Render plot
        self.render_to_painter(painter, width, height)

        painter.end()

        # Save image
        return image.save(filename, "PNG", 100)

    def export_to_jpg(self, filename, width=1920, height=1440, dpi=300, quality=95):
        """Export plot to JPG file"""
        # Create image with specified size
        image = QImage(width, height, QImage.Format_RGB32)
        image.fill(QColor(255, 255, 255))

        # Set DPI
        dpi_to_dpm = dpi / 0.0254  # Convert DPI to dots per meter
        image.setDotsPerMeterX(int(dpi_to_dpm))
        image.setDotsPerMeterY(int(dpi_to_dpm))

        # Create painter for the image
        painter = QPainterExport(image)
        painter.setRenderHint(QPainter.Antialiasing)

        # Render plot
        self.render_to_painter(painter, width, height)

        painter.end()

        # Save image with quality setting
        return image.save(filename, "JPG", quality)

    def export_to_eps(self, filename, width=1920, height=1440, dpi=300):
        """Export plot to EPS file (using PostScript printer)"""
        printer = QPrinter(QPrinter.HighResolution)
        printer.setOutputFileName(filename)
        printer.setOutputFormat(QPrinter.PdfFormat)  # Qt6 uses PDF format
        printer.setResolution(dpi)

        # Set page size (convert pixels to mm)
        width_mm = width * 25.4 / dpi
        height_mm = height * 25.4 / dpi
        page_size = QPageSize(QSizeF(width_mm, height_mm), QPageSize.Millimeter)
        printer.setPageSize(page_size)
        printer.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Millimeter)

        # Create painter for the printer
        painter = QPainterExport(printer)
        painter.setRenderHint(QPainter.Antialiasing)

        # Calculate scaling factor
        scale_x = printer.width() / width
        scale_y = printer.height() / height
        scale = min(scale_x, scale_y)

        # Scale and render
        painter.scale(scale, scale)
        self.render_to_painter(painter, width, height)

        painter.end()

        return True
########################################################################################################################
########################################################################################################################
class MainWindow(QMainWindow):
    def __init__(self, data_groups, colors=None, radius=10, window_title="Scatter Plot"):
        super().__init__()
        self.setWindowTitle(window_title)
        self.setGeometry(100, 100, 1000, 850)

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create scatter plot widget
        self.scatter_plot = ScatterPlotWidget(data_groups, colors, radius)
        layout.addWidget(self.scatter_plot)

        # Create button layout
        button_layout = QHBoxLayout()

        # Export buttons
        self.btn_export_png = QPushButton("Export PNG")
        self.btn_export_jpg = QPushButton("Export JPG")
        self.btn_export_eps = QPushButton("Export EPS")

        self.btn_export_png.clicked.connect(self.export_png)
        self.btn_export_jpg.clicked.connect(self.export_jpg)
        self.btn_export_eps.clicked.connect(self.export_eps)

        button_layout.addWidget(self.btn_export_png)
        button_layout.addWidget(self.btn_export_jpg)
        button_layout.addWidget(self.btn_export_eps)
        button_layout.addStretch()

        layout.addLayout(button_layout)

    def export_png(self):
        """Export to PNG"""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export PNG", "", "PNG Files (*.png)"
        )
        if filename:
            if not filename.lower().endswith('.png'):
                filename += '.png'
            try:
                success = self.scatter_plot.export_to_png(filename, width=1920, height=1440, dpi=300)
                if success:
                    print(f"Successfully exported to {filename}")
                else:
                    print(f"Failed to export to {filename}")
            except Exception as e:
                print(f"Error exporting PNG: {e}")

    def export_jpg(self):
        """Export to JPG"""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export JPG", "", "JPG Files (*.jpg *.jpeg)"
        )
        if filename:
            if not (filename.lower().endswith('.jpg') or filename.lower().endswith('.jpeg')):
                filename += '.jpg'
            try:
                success = self.scatter_plot.export_to_jpg(filename, width=1920, height=1440, dpi=300, quality=95)
                if success:
                    print(f"Successfully exported to {filename}")
                else:
                    print(f"Failed to export to {filename}")
            except Exception as e:
                print(f"Error exporting JPG: {e}")

    def export_eps(self):
        """Export to EPS (saved as PDF in Qt6)"""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export EPS", "", "EPS Files (*.eps);;PDF Files (*.pdf)"
        )
        if filename:
            # Accept both .eps and .pdf extensions
            if not (filename.lower().endswith('.eps') or filename.lower().endswith('.pdf')):
                filename += '.eps'
            try:
                success = self.scatter_plot.export_to_eps(filename, width=1920, height=1440, dpi=300)
                if success:
                    print(f"Successfully exported to {filename}")
                    if filename.lower().endswith('.eps'):
                        print("Note: File is saved in PDF format (Qt6 limitation)")
                else:
                    print(f"Failed to export to {filename}")
            except Exception as e:
                print(f"Error exporting EPS: {e}")
########################################################################################################################
########################################################################################################################
class DataVisualizer:
    """Class that contains the scatter plot visualization method"""

    def __init__(self):
        self.app = None
        self.window = None

    def show_scatter_plot(self, data_groups, colors=None, radius=10, window_title="Scatter Plot"):
        """
        Display scatter plot with given data

        Parameters:
        -----------
        data_groups : list
            List of groups, each group is a list of [id, y, z, area]
            Example: [
                [[1, 10, 20, 100], [2, 15, 25, 200]],  # Group 1
                [[3, 30, 35, 150], [4, 35, 40, 180]],  # Group 2
            ]
        colors : list, optional
            List of colors for each group. Supported formats:
            - Color names (string): "red", "green", "blue", "yellow", "orange",
                                   "purple", "pink", "cyan", "magenta", "brown",
                                   "lime", "navy", "teal", "olive", "maroon",
                                   "gray", "grey", "black", "white", "gold"
            - RGB tuple: (255, 0, 0)
            - RGBA tuple: (255, 0, 0, 255)
            - Hex string: "#FF0000"
            - QColor object
        radius : int, optional
            Circle radius in pixels (default: 10, all points same size)
        window_title : str, optional
            Window title (default: "Scatter Plot")
        """
        # Create QApplication if not exists
        if self.app is None:
            self.app = QApplication.instance()
            if self.app is None:
                self.app = QApplication(sys.argv)

        # Create main window with export buttons
        self.window = MainWindow(data_groups, colors, radius, window_title)

        # Show window
        self.window.show()

        # Execute event loop
        sys.exit(self.app.exec())
########################################################################################################################
########################################################################################################################
if __name__ == "__main__":
    pass







