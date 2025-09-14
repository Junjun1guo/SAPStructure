#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo, Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.13
#  Date: 2025-08-30
########################################################################################################################
########################---import modules---############################################################################
import numpy as np
import os
import math
import pandas as pd
from time import time
import openseespy.opensees as ops
import sys
########################################################################################################################
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
############################----import auxiliary modules---#############################################################
from auxiliaryModules.mainMod import OpenSeesPyX
from auxiliaryModules.py_tz_qz_SoilParameters.pyModel import (py_sand,py_softClay_freeWater,
                            py_stiffClay_freeWater,py_stiffClay_withoutFreeWater)
from auxiliaryModules.py_tz_qz_SoilParameters.tzModel import tz_sand,tz_clay
from auxiliaryModules.py_tz_qz_SoilParameters.qzModel import qz_sand,qz_clay
########################################################################################################################
########################################################################################################################
##################---初始化辅助类OpenSeesPyX用于OPenSeesPy模型可视化及后处理---###############################################
opsX = OpenSeesPyX(dataBaseName="pileSoilModel")  ###初始化OpenSeesPyX类
########################################################################################################################
ops.wipe()
for Scale in range(1,2,1):
    for caseNum in range(1,2,1):
        ops.wipe()
        with open("AccNum2.txt", "r+") as file_Acc:
            Acc_array = file_Acc.read().splitlines()
            AccNumber = Acc_array[caseNum-1]
        output = f"CaseNum({caseNum})"
        base_path = os.path.join("Data", output)
        pile_path = os.path.join(base_path, "Pile")
        strain_path = os.path.join(base_path, "Strain")
        log_path = os.path.join(base_path, "Log")
        os.makedirs(base_path, exist_ok=True)
        os.makedirs(pile_path, exist_ok=True)
        os.makedirs(strain_path, exist_ok=True)
        os.makedirs(log_path, exist_ok=True)
        ################################################################################################################
        ######################### Set Scenario Parameters
        with open("Cases.txt","r+") as casesParas:
            parasArray=casesParas.read().splitlines()[caseNum-1]
            parasSplit=parasArray.strip().split()
            N=int(parasSplit[0]) ## Num. of pile rows
            d=float(parasSplit[1]) ## Diameter of piles (m)
            rho_pile=float(parasSplit[2]) ## Longitudinal reinforcement ratio of piles
            alpha=float(parasSplit[3]) ## Axial ratio of piles
            p2pDis=float(parasSplit[4]) ## Pile-to-Pile distance (d)
            Dr=float(parasSplit[5]) ## Soil Relative Density
            ScourDepth=float(parasSplit[6]) ## Scour depth (m)
            Hp=float(parasSplit[7]) ## Pier height (m)
            DDT=float(parasSplit[8]) ## Pier size (d) along the loading direction
            rho_pier=float(parasSplit[9]) ## Longitudinal reinforcement ratio of pier
            rhoTran_pile=float(parasSplit[10]) ## Transverse reinforcement ratio of pile
            fy_pile=float(parasSplit[11]) ## steel yield strength (MPa)
            fc_pile=float(parasSplit[12]) ## concrete cover peak strength (MPa)
            LP=float(parasSplit[13]) ## Pile total length (m)
        p2sDis=1 ## Pile-to-side distance (d)
        Hc=2 ## Cap height (d)
        DDL=2 ## Pier size (d) perpendicular to the loading direction
        ########## steel/concrete material
        ### pile
        e_fc_pile=0.002  # concrete cover peak strength corresponding strain
        e_fcu_pile=0.005  # concrete cover ultimate strength corresponding strain
        E_y_pile=200.0  # steel young's modulus (GPa)
        e_fy_pile=fy_pile/E_y_pile/1000.0 #Yield strain of steel
        bs_pile=0.01  # hardening coefficient
        ke_pile=0.95
        fyh_pile=335.0
        fl_pile=0.5*ke_pile*rhoTran_pile*fyh_pile
        fcc_pile=(-1.254 + 2.254 * math.sqrt(1 + 7.94 * fl_pile / fc_pile) - 2 * fl_pile / fc_pile) * fc_pile
        e_fcc_pile=e_fc_pile * (1 + 5 * (fcc_pile / fc_pile - 1))
        e_fccu_pile=0.004 + 1.4 * rhoTran_pile * fyh_pile * 0.09 / fcc_pile
        ### pier
        fc_pier=fc_pile  # concrete cover peak strength (MPa)
        e_fc_pier=0.002  # concrete cover peak strength corresponding strain
        e_fcu_pier=0.005 # concrete cover ultimate strength corresponding strain
        fy_pier=fy_pile  # steel yield strength (MPa)
        E_y_pier=200.0   # steel young's modulus (GPa)
        bs_pier=0.01     # hardening coefficient
        ke_pier=0.75
        fyh_pier=335.0
        rhoTran_pier=0.01
        fl_pier=0.5*ke_pier*rhoTran_pier *fyh_pier
        fcc_pier=(-1.254 + 2.254 * math.sqrt(1 + 7.94 * fl_pier / fc_pier) - 2 * fl_pier / fc_pier) * fc_pier
        e_fcc_pier=e_fc_pier*(1+5*(fcc_pier/fc_pier-1))
        e_fccu_pier=0.004+1.4*rhoTran_pier*fyh_pier*0.09/fcc_pier
        ### Sand material
        Gs=2.650
        rw=1.000
        emax=0.894
        emin=0.516
        e=emax -Dr * (emax-emin)
        massDen=rw * (Gs+e) / (1.000 +e) ###---土的三相指标转换
        phi=16.000 *Dr *Dr + 0.170 *Dr + 28.400
        N60=60 *Dr *Dr
        Vs=85 * (2.5 + N60) ** 0.25
        Gsoil=massDen *Vs *Vs
        ### Pile group effect
        P_leading=0.100 * (p2pDis-3.0) + 0.8
        P_middle=0.225 * (p2pDis-3.0) + 0.4
        P_trailing=0.200 * (p2pDis-3.0) + 0.3
        ### Soil thickness
        fullDepth=LP -ScourDepth # full depth of soil layers
        ### cap mass
        Mcap=((N - 1) * p2pDis * d + 2 * p2sDis * d) * (p2pDis * d + 2 * p2sDis * d) * (Hc * d) * 2.6
        ### Record  case parameters
        casePara_path = os.path.join("Data",f"CaseNum({caseNum})","CaseParaFile")
        os.makedirs(casePara_path, exist_ok=True)
        writeList=[N,d,rho_pile,alpha,p2pDis,Dr,ScourDepth,Hp,DDT,rho_pier,rhoTran_pile,fy_pile,fc_pile,LP,
                   rhoTran_pier,p2sDis,Hc,DDL,massDen,phi,Gsoil,fcc_pile,e_fcc_pile,e_fccu_pile,fcc_pier,
                   e_fcc_pier,e_fccu_pier,e_fy_pile]
        with open(f"{casePara_path}/{AccNumber}-{Scale}.txt", "w+") as caseParasFile:
            #####
            for para in writeList:
                caseParasFile.write(str(para)+'\n')
        ################################################################################################################
        ######################### Define input ground motion
        with open("Waves/dt.txt", "r+") as file_deltaT:
            deltaT_array = file_deltaT.read().splitlines()[int(AccNumber)-1]
            motionDT=float(deltaT_array.strip().split()[0])

        with open("Waves/motionStep.txt", "r+") as file_gmlength:
            gmlength_array = file_gmlength.read().splitlines()[int(AccNumber)-1]
            motionSteps = int(gmlength_array.strip().split()[0])
        ################################################################################################################
        ######################### structure modeling
        ops.model('basic', '-ndm', 3, '-ndf', 6)
        ySize=0.5  # y (V) direction element size for piles and pier
        numYele=round(fullDepth / ySize) # number of embeded pile elements in y (V) direction, round down
        numYnode=numYele + 1
        numPele=round(LP/ySize) # number of total pile elements in y (V) direction
        numPnode=numPele + 1
        pi=3.1415926
        ######################### fiber pier/pile section
        DT=DDT *d
        DL=DDL *d
        thickness=d / 3.0
        cover=0.10
        ### Define Concrete-Steel Material Tag
        IDcover=9903  # material ID tag -- unconfined cover concrete
        IDcore=9904  # material ID tag -- confined core concrete
        IDreinf=9905  # material ID tag -- reinforcement
        # Concrete04
        # Cover Concrete 保护层混凝土材性参数
        coverfy_pier=-fc_pier * 1000
        coverdy_pier=-e_fc_pier
        coverdu_pier= -e_fcu_pier
        Ec_pier=4700 * math.sqrt(fc_pier) * 1000 # Concrete Elastic Modulus [kPa]
        # Core Concrete 核心混凝土材性参数
        corefy_pier= -fcc_pier * 1000 # concrete core strength as uncertainty
        coredy_pier= -e_fcc_pier
        coredu_pier= -e_fccu_pier
        Ecc_pier=4700 * math.sqrt(fcc_pier) * 1000
        # Establish concrete material 混凝土材料定义
        ops.uniaxialMaterial('Concrete04',9913,coverfy_pier,coverdy_pier,coverdu_pier,Ec_pier)
        ops.uniaxialMaterial('Concrete04', 9914, corefy_pier,coredy_pier,coredu_pier,Ecc_pier)
        ops.uniaxialMaterial('MinMax',IDcover,9913,'-min',coverdu_pier)
        ops.uniaxialMaterial('MinMax', IDcore, 9914, '-min', coredu_pier)
        # Longitudinal rebar 纵向钢筋材性参数(Steel02)
        Fy=fy_pier * 1000  # STEEL yield stress  read unit MPa
        Es=E_y_pier * 1000000  # modulus of steel read unit GPa
        Bs=bs_pier  # strain-hardening ratio
        R0=20  # control the transition from elastic to plastic branches
        cR1=0.925 # control the transition from elastic to plastic branches
        cR2=0.15
        a1=0
        a2=1
        a3=0
        a4=1
        # Establish Longitudinal rebar material 纵向钢筋定义
        ops.uniaxialMaterial('Steel02',IDreinf,Fy,Es,Bs,R0,cR1,cR2,a1,a2,a3,a4)
        # FIBER SECTION MESH
        # Corner coordinates
        Cor11Y=DT / 2.0
        Cor11Z=DL / 2.0
        Cor12Y=-Cor11Y
        Cor12Z=Cor11Z
        Cor13Y=-Cor11Y
        Cor13Z=-Cor11Z
        Cor14Y=Cor11Y
        Cor14Z=-Cor11Z

        Cor21Y=Cor11Y-cover
        Cor21Z=Cor11Z -cover
        Cor22Y=-Cor21Y
        Cor22Z=Cor21Z
        Cor23Y=-Cor21Y
        Cor23Z=-Cor21Z
        Cor24Y=Cor21Y
        Cor24Z=-Cor21Z

        Cor31Y=Cor21Y-thickness
        Cor31Z=Cor21Z -thickness
        Cor32Y=-Cor31Y
        Cor32Z=Cor31Z
        Cor33Y=-Cor31Y
        Cor33Z=-Cor31Z
        Cor34Y=Cor31Y
        Cor34Z=-Cor31Z

        Cor41Y=Cor31Y-cover
        Cor41Z=Cor31Z -cover
        Cor42Y=-Cor41Y
        Cor42Z=Cor41Z
        Cor43Y=-Cor41Y
        Cor43Z=-Cor41Z
        Cor44Y=Cor41Y
        Cor44Z=-Cor41Z

        areaBar=4 * (Cor21Y * Cor21Z - Cor31Y * Cor31Z) * rho_pier / 182.0

        opsX.section('Fiber',1,'-GJ', 1.0e10)
        # Outer cover
        opsX.patch('rect',IDcover,int(cover/0.1),int(2*Cor11Z/0.1),Cor24Y,Cor14Z,Cor11Y,Cor11Z) #right
        opsX.patch('rect',IDcover,int(cover/0.1),int(2*Cor11Z/0.1),Cor13Y,Cor13Z,Cor22Y,Cor12Z) #left
        opsX.patch('rect',IDcover,int(2*Cor21Y/0.1),int(cover/0.1),Cor22Y,Cor22Z,Cor21Y,Cor11Z) #top
        opsX.patch('rect', IDcover,int(2*Cor21Y/0.1),int(cover/0.1),Cor23Y,Cor13Z,Cor24Y,Cor24Z) #bottom
        # Core concrete
        opsX.patch('rect',IDcore,int(thickness/0.1),int(2*Cor21Z/0.1),Cor34Y,Cor24Z,Cor21Y,Cor21Z) #right
        opsX.patch('rect',IDcore,int(thickness/0.1),int(2*Cor21Z/0.1),Cor23Y,Cor23Z,Cor32Y,Cor22Z)  #left
        opsX.patch('rect',IDcore,int(2*Cor31Y/0.1),int(thickness/0.1),Cor32Y,Cor32Z,Cor31Y,Cor21Z )  #top
        opsX.patch('rect', IDcore,int(2*Cor31Y/0.1),int(thickness/0.1),Cor33Y,Cor23Z,Cor34Y,Cor34Z)  #bottom
        # Inner cover
        opsX.patch('rect', IDcover,int(cover/0.1),int(2*Cor31Z/0.1),Cor44Y,Cor34Z,Cor31Y,Cor31Z) # right
        opsX.patch('rect', IDcover,int(cover/0.1),int(2*Cor31Z/0.1),Cor33Y,Cor33Z,Cor42Y,Cor32Z)  #left
        opsX.patch('rect', IDcover,int(2*Cor41Y/0.1),int(cover/0.1),Cor42Y,Cor42Z,Cor41Y,Cor31Z) #top
        opsX.patch('rect', IDcover,int(2*Cor41Y/0.1),int(cover/0.1),Cor43Y,Cor33Z,Cor44Y,Cor44Z) #bottom
        # Outer rebars
        opsX.layer('straight',IDreinf,36,areaBar,Cor22Y,Cor22Z,Cor21Y,Cor21Z) #top
        opsX.layer('straight', IDreinf, 36,areaBar,Cor23Y,Cor23Z,Cor24Y,Cor24Z) #bottom
        opsX.layer('straight', IDreinf,17,areaBar,Cor22Y,Cor22Z-0.15,Cor23Y,Cor23Z+0.15) #left
        opsX.layer('straight', IDreinf,17,areaBar,Cor21Y,Cor21Z-0.15,Cor24Y,Cor24Z+0.15)  #right
        # Inner rebars
        opsX.layer('straight', IDreinf,27,areaBar,Cor32Y,Cor32Z,Cor31Y,Cor31Z)  #top
        opsX.layer('straight', IDreinf,27,areaBar,Cor33Y,Cor33Z,Cor34Y,Cor34Z)  #bottom
        opsX.layer('straight', IDreinf,11,areaBar,Cor32Y,Cor32Z-0.15,Cor33Y,Cor33Z+0.15)  #left
        opsX.layer('straight', IDreinf,11,areaBar,Cor31Y,Cor31Z-0.15,Cor34Y,Cor34Z+0.15)  #right
        #####---opensees内置纤维划分可视化
        # pltHandle = opsX.auxiliary_fiberSectionPlot(sectTag=1)
        # pltHandle.savefig("sect2.eps")
        # pltHandle.show()
        ################################################################################################################
        # =======================================
        #                 Pile
        # =======================================
        IDcover2=9923  # material ID tag -- unconfined cover concrete
        IDcore2=9924  # material ID tag -- confined core concrete
        # Cover Concrete 保护层混凝土材性参数
        coverfy_pile=-fc_pile * 1000
        coverdy_pile=-e_fc_pile
        coverdu_pile= -e_fcu_pile
        Ec_pile= 4700 * math.sqrt(fc_pile) * 1000 # Concrete Elastic Modulus [kPa]
        # Core Concrete 核心混凝土材性参数
        corefy_pile=-fcc_pile * 1000  # concrete core strength as uncertainty
        coredy_pile= -e_fcc_pile
        coredu_pile=-e_fccu_pile
        Ecc_pile=4700 * math.sqrt(fcc_pile)*1000
        # Establish concrete material 混凝土材料定义
        ops.uniaxialMaterial('Concrete04', 9933,coverfy_pile,coverdy_pile,coverdu_pile,Ec_pile)
        ops.uniaxialMaterial('Concrete04', 9934,corefy_pile,coredy_pile,coredu_pile,Ecc_pile)
        # Combined with minmax command
        ops.uniaxialMaterial('MinMax',IDcover2,9933, '-min',coverdu_pile)
        ops.uniaxialMaterial('MinMax',IDcore2,9934, '-min',coredu_pile)
        # Fiber section
        rho_l_pile=rho_pile
        numBarsPile=30  # Pile rebar number
        areaBar2=math.pi * 0.25 * (d - 2 * cover) ** 2 * rho_l_pile / numBarsPile
        intRad2=0.0
        extRad2=d / 2.0
        coreRad2=extRad2-cover # The distance from the the centre of the circle to the edge of the cover concrete
        nsubdivCirc2=10 # number of fibers for concrete in circle
        nsubdivRadcover2=2
        nsubdivRadcore2=8
        yCenter2=0
        zCenter2=0
        opsX.section('Fiber', 2, '-GJ', 1.0e10)
        opsX.patch('circ',IDcover2,nsubdivCirc2,nsubdivRadcover2,yCenter2,zCenter2,coreRad2,extRad2,0,360)
        opsX.patch('circ',IDcore2,nsubdivCirc2,nsubdivRadcore2,yCenter2,zCenter2,intRad2,coreRad2,0,360)
        opsX.layer('circ',IDreinf,numBarsPile,areaBar2,yCenter2,zCenter2,coreRad2)
        #####---opensees内置纤维划分可视化
        # pltHandle = opsX.auxiliary_fiberSectionPlot(sectTag=2)
        # pltHandle.savefig("sect2.eps")
        # pltHandle.show()
        ################################################################################################################
        ######################### Define nodes and mass for pile/cap/pier
        # ========== define mass ==========
        AgPier=4 * (Cor11Y * Cor11Z-Cor41Y * Cor41Z)
        AgPile=math.pi * (d / 2.0) ** 2
        arbiEleMassPier=AgPier*ySize*2.6 # define arbitrary element ($ySize=0.5m) mass for pier
        arbiEleMassPile=AgPile*ySize*2.6 # define arbitrary element ($ySize=0.5m) mass for pile
        superStrucMass=2*N*fc_pile*1000*AgPile*alpha/9.81-Mcap-AgPier*Hp*2.6
        # ========== define pile nodes ===========
        for row in range(N):
            pilenum_ctr=(2 *row + 1)*1000
            for i in range(1,(numPnode+1),1):
                xdim=row *p2pDis *d
                ydim=(i - 1) * ySize
                zdim=0.5 *p2pDis *d
                nodenum=pilenum_ctr +i
                opsX.node(nodenum,xdim,ydim,zdim)
                ops.mass(nodenum,arbiEleMassPile,arbiEleMassPile,arbiEleMassPile,0,0,0)
                ops.fix(nodenum,0,0,1,1,1,0)
            pilenum_ctr=(2 *row + 2)*1000
            for i in range(1, (numPnode + 1), 1):
                xdim=row *p2pDis *d
                ydim=(i - 1) * ySize
                zdim=- 0.5 *p2pDis *d
                nodenum=pilenum_ctr +i
                opsX.node(nodenum,xdim,ydim,zdim)
                ops.mass(nodenum,arbiEleMassPile,arbiEleMassPile,arbiEleMassPile,0,0,0)
                ops.fix(nodenum,0,0,1,1,1,0)
        # ==========================================================
        # =========== define cap nodes ==============
        # cap bottom
        CenterCorX=(N - 1) * p2pDis *d / 2.0
        xdim=CenterCorX
        ydim=numPele *ySize
        zdim=0
        nodenum=610 # Cap bottom
        opsX.node(nodenum,xdim,ydim,zdim)
        ops.mass(nodenum,0,0,0,0,0,0)
        ops.fix(nodenum,0,0,1,1,1,0)
        # cap center
        xdim=CenterCorX
        ydim=numPele *ySize + 0.5 *Hc *d
        zdim=0
        nodenum=611  # Cap center
        opsX.node(nodenum,xdim,ydim,zdim)
        ops.mass(nodenum,Mcap,Mcap,Mcap,0,0,0)
        ops.fix(nodenum,0,0,1,1,1,0)
        # cap top (pier bottom)
        xdim=CenterCorX
        ydim=numPele *ySize +Hc *d
        zdim=0
        nodenum=1000 +numPnode + 1
        opsX.node(nodenum,xdim,ydim,zdim)
        ops.mass(nodenum,0,0,0,0,0,0)
        ops.fix(nodenum,0,0,1,1,1,0)
        # ======== define pier nodes =============
        Pier_total_nodenum=round(numPnode + 1 +Hp /ySize)
        for i in range(numPnode+2,Pier_total_nodenum+1,1):
            xdim=CenterCorX
            ydim=(i - 2) * ySize +Hc *d
            zdim=0
            nodenum=1000 +i
            opsX.node(nodenum,xdim,ydim,zdim)
            ops.mass(nodenum,arbiEleMassPier,arbiEleMassPier,arbiEleMassPier,0,0,0)
            ops.fix(nodenum,0,0,1,1,1,0)
        # ========== define pierhead node ===========
        pierTopNodeNum=round(1000 +Pier_total_nodenum)  # define pierTopNodeNum, used for building bearing
        # ========== define superstructure node ===========
        superStrucNodeNum=round(pierTopNodeNum + 1) # define superStrucNodeNum, used for building bearing
        opsX.node(superStrucNodeNum,xdim,ydim,zdim)
        ops.mass(superStrucNodeNum,superStrucMass,superStrucMass,superStrucMass,0,0,0)
        ops.fix(superStrucNodeNum,0,0,1,1,1,0)
        # ===============================================
        # DEFINE BEARING
        # ===============================================
        # define bearing material
        ops.uniaxialMaterial('Elastic',9901,1E10)
        opsX.element('zeroLength',9999,pierTopNodeNum,superStrucNodeNum,'-mat',9901,9901,9901,9901,9901,9901,
                     '-dir',1,2,3,4,5,6)
        # ===============================================
        # DEFINE PILE ELEMENT
        # ===============================================
        # create coordinate-transformation object
        opsX.geomTransf('Corotational',1,0,0,-1) # Considering the P-△ effect
        opsX.geomTransf('Linear',2,0,1,0) # For pile-cap connections
        # define integration point number
        np=5
        # creat pile element
        for row in range(N):
            pilenum_ctr=(2 *row + 1)*1000
            for i in range(1,numPele+1,1):
                opsX.element('nonlinearBeamColumn',pilenum_ctr+i,pilenum_ctr+i,pilenum_ctr+i+1,np,2,1)
            pilenum_ctr=(2 *row + 2)*1000
            for i in range(1,numPele+1,1):
                opsX.element('nonlinearBeamColumn',pilenum_ctr+i,pilenum_ctr+i,pilenum_ctr+i+1,np,2,1)
        # creat pile-cap connection (rigid element)
        for row in range(N):
            pilenum_ctr=(2 *row + 1)*1000
            opsX.element('elasticBeamColumn',pilenum_ctr+numPnode,610,pilenum_ctr+numPnode,50,3.25e7,1.25e7,100,100,100,2)
            pilenum_ctr=(2 *row + 2)*1000
            opsX.element('elasticBeamColumn',pilenum_ctr+numPnode,610,pilenum_ctr+numPnode,50,3.25e7,1.25e7,100,100,100,2)
        opsX.element('elasticBeamColumn',10888,610,611,50,3.25e7,1.25e7,100,100,100,1)
        opsX.element('elasticBeamColumn',10999,611,1000+numPnode+1,50,3.25e7,1.25e7,100,100,100,1)
        # creat pier element
        for i in range(numPnode+1,Pier_total_nodenum,1):
            opsX.element('nonlinearBeamColumn',1000+i,1000+i,1000+i+1,np,1,1)
        ################################################################################################################
        ### define soil nodes
        for row in range(N):
            soilNum_ctr=(2 *row + 1)*1100000
            for i in range(1,numYnode+1,1):
                xdim=row *p2pDis *d
                ydim=(i - 1) * ySize
                zdim=0.5 *p2pDis *d  # $pp2pDis means the piele-to-pile distance (x times pile diameter)
                nodenum=soilNum_ctr +i
                opsX.node(nodenum,xdim,ydim,zdim)
                ops.fix(nodenum,1,1,1,1,1,1)
            soilNum_ctr=(2 *row + 2)*1100000
            for i in range(1,numYnode+1,1):
                xdim=row *p2pDis *d
                ydim=(i - 1) * ySize
                zdim=- 0.5 *p2pDis *d
                nodenum=soilNum_ctr +i
                opsX.node(nodenum,xdim,ydim,zdim)
                ops.fix(nodenum,1,1,1,1,1,1)
        ## soil properties (according to soil layer type)
        eleSize=ySize
        diameter=d
        # soil unit weight (kN/m^3)
        gamma=9.81 * (massDen-1.0)  # effective density
        # soil internal friction angle (degrees)
        phi=phi # friction angle DDS
        # soil shear modulus at pile tip (kPa)
        Gsoil=Gsoil;  # shear modulus DDS
        # select pult definition method for p-y curves
        # API (default) --> 1
        # Brinch Hansen --> 2
        puSwitch=1
        # variation in loading coefficent A for pult calculation
        # static loading, Reese, 2011   --> 1
        # cyclic loading, API           --> 2
        ASwitch=2
        # variation in coefficent of subgrade reaction with depth for p-y curves
        # API linear variation (default)   --> 1
        # modified API parabolic variation --> 2
        # ATC-1996                         --> 3
        kSwitch=1
        # effect of ground water on subgrade reaction modulus for p-y curves
        # above gwt --> 1
        # below gwt --> 2
        gwtSwitch=2
        # specify necessary procedures
        ###---number of pile rows
        if N==2:
            groupEffFactor_leading=P_trailing ### back
            groupEffFactor_trailing=P_trailing
            groupEffFactor_active=P_leading -P_trailing

            # p-y spring material
            for i in range(2,numYnode+1,1):
                # depth of current py node
                pyDepth=fullDepth -eleSize * (i-1)
                # procedure to define pult and y50
                pult0,y50=py_sand(pyDepth,gamma,phi,diameter,eleSize,puSwitch,ASwitch,kSwitch,gwtSwitch)
                pult_leadingPile=pult0*groupEffFactor_leading
                pult_trailingPile=pult0 *groupEffFactor_trailing
                ops.uniaxialMaterial('PySimple1',i+3000,2,pult_leadingPile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+1000,2,pult_trailingPile,y50,0.3)
                pult_active=pult0*groupEffFactor_active
                ops.uniaxialMaterial('QzSimple1',i+6000,2,pult_active,y50) ### parallel with p-y simulate pile group effect
                row=0
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                ####---element parallel to simulate pile group effect in y direction
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                ####---element parallel to simulate pile group effect in y direction
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)

                row=1
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0) ###---to consider pile group effect
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
        if N==3:
            groupEffFactor_leading=P_trailing
            groupEffFactor_middle=P_middle
            groupEffFactor_trailing=P_trailing
            groupEffFactor_active=P_leading-P_trailing
            # p-y spring material
            for i in range(2,numYnode+1,1):
                # depth of current py node
                pyDepth=fullDepth -eleSize * (i-1)
                # procedure to define pult and y50
                pult0, y50 = py_sand(pyDepth, gamma, phi, diameter, eleSize, puSwitch, ASwitch, kSwitch, gwtSwitch)
                pult_leadingPile=pult0 *groupEffFactor_leading
                pult_middlePile=pult0 *groupEffFactor_middle
                pult_trailingPile=pult0 *groupEffFactor_trailing
                ops.uniaxialMaterial('PySimple1',i+3000,2,pult_leadingPile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+2000,2,pult_middlePile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+1000,2,pult_trailingPile,y50,0.3)
                pult_active=pult0 *groupEffFactor_active
                ops.uniaxialMaterial('QzSimple1',i+6000,2,pult_active,y50)

                row=0
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)

                row=1
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+2000,'-dir',1)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+2000,'-dir',1)

                row=2
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
        if N==4:
            groupEffFactor_leading=P_trailing
            groupEffFactor_middle=P_trailing
            groupEffFactor_trailing=P_trailing
            groupEffFactor_active1=P_middle -P_trailing
            groupEffFactor_active2=P_leading -P_trailing
            # p-y spring material
            for i in range(2,numYnode+1,1):
                # depth of current py node
                pyDepth=fullDepth -eleSize * (i-1)
                # procedure to define pult and y50
                pult0, y50 = py_sand(pyDepth, gamma, phi, diameter, eleSize, puSwitch, ASwitch, kSwitch, gwtSwitch)
                pult_leadingPile=pult0 *groupEffFactor_leading
                pult_middlePile=pult0 *groupEffFactor_middle
                pult_trailingPile=pult0 *groupEffFactor_trailing

                ops.uniaxialMaterial('PySimple1',i+4000,2,pult_leadingPile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+3000,2,pult_middlePile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+2000,2,pult_middlePile,y50,0.3)
                ops.uniaxialMaterial('PySimple1',i+1000,2,pult_trailingPile,y50,0.3)

                pult_active1=pult0 *groupEffFactor_active1
                pult_active2=pult0 *groupEffFactor_active2

                ops.uniaxialMaterial('QzSimple1',i+6000,2,pult_active1,y50)
                ops.uniaxialMaterial('QzSimple1',i+7000,2,pult_active2,y50)

                row=0
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+7000,'-dir',1)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+1000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+7000,'-dir',1)

                row=1
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+2000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+2000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1)

                row=2
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+3000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+6000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)

                row=3
                soilNum_ctr=(2 *row + 1)*1100000
                pilenum_ctr=(2 *row + 1)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+4000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+7000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
                soilNum_ctr=(2 *row + 2)*1100000
                pilenum_ctr=(2 *row + 2)*1000
                opsX.element('zeroLength',soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+4000,'-dir',1)
                opsX.element('zeroLength',10*soilNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+7000,'-dir',1,
                             '-orient',-1,0,0,0,1,0)
        # t-z spring material
        for i in range(2,numYnode+1,1):
            # depth of current tz node
            pyDepth=fullDepth -eleSize * (i-1)
            tult,z50=tz_sand(phi,diameter,gamma,pyDepth,eleSize)
            ops.uniaxialMaterial('TzSimple1',i+10000,2,1.8*tult,z50)

        ops.uniaxialMaterial('ENT',10001,10000000)

        for row in range(N):
            eleNum_ctr=(2 *row + 1)*11000
            soilNum_ctr=(2 *row + 1)*1100000
            pilenum_ctr=(2 *row + 1)*1000
            for i in range(1,numYnode+1,1):
                opsX.element('zeroLength',eleNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+10000,'-dir',2)
            eleNum_ctr=(2 *row + 2)*11000
            soilNum_ctr=(2 *row + 2)*1100000
            pilenum_ctr=(2 *row + 2)*1000
            for i in range(1,numYnode+1,1):
                opsX.element('zeroLength',eleNum_ctr+i,soilNum_ctr+i,pilenum_ctr+i,'-mat',i+10000,'-dir',2)

        #### "=== Load Structure Gravity"
        GnumStep=1
        arbiEleWeightPile=arbiEleMassPile * 9.81
        arbiEleWeightPier=arbiEleMassPier * 9.81
        superStrucWeight=superStrucMass * 9.81
        ops.timeSeries('Constant', 1)
        ops.pattern('Plain', 1,1)
        # pile
        for row in range(N):
            pilenum_ctr=(2 *row + 1)*1000
            for i in range(1,numPnode+1,1):
                nodenum=pilenum_ctr +i
                ops.load(nodenum,0,-arbiEleWeightPile/GnumStep,0,0,0,0)
            pilenum_ctr=(2 *row + 2)*1000
            for i in range(1,numPnode+1,1):
                nodenum=pilenum_ctr +i
                ops.load(nodenum,0,-arbiEleWeightPile/GnumStep,0,0,0,0)
        # pile cap
        ops.load(611,0,-9.81*Mcap/GnumStep,0,0,0,0)
        # pier
        for i in range(numPnode+2,Pier_total_nodenum+1,1):
            nodenum=1000 +i
            ops.load(nodenum,0,-arbiEleWeightPier/GnumStep,0,0,0,0)
        # superstructure
        ops.load(superStrucNodeNum,0,-superStrucWeight/GnumStep,0,0,0,0)
        ## "--- Structural Gravity Analysis Begin..."
        recordELeList=[int(1000+numPnode+1)]
        opsX.integration_recorderElement(f"Data/{output}/Pile",'99000.out',recordELeList,'localForce')

        ops.test('NormDispIncr',1.0e-4,50)
        ops.algorithm('NewtonLineSearch',0.8)
        ops.integrator('Newmark',0.5,0.25)
        ops.system('BandGeneral')
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.analysis('Transient')
        ok = ops.analyze(1,1)
        if ok==0:
            print("--- Structure gravity analysis succeed!")
        else:
            print("--- Structure gravity analysis disconvergence!!!")
        ops.setTime(0.0)
        ops.wipeAnalysis()
        ops.remove("recorders")
        ########################################
        ### Beginning Modal Analysis..."
        ### determine Natural Period, Frequency & damping parameters for SDOF
        ### set xDamp 0.05; # damping ratio (0.02-0.05-typical)
        m_order=10 # model order
        lamiga=ops.eigen(m_order)
        periodList=[]
        for i in range(m_order):
            omegaa=(lamiga[i])**0.5
            Ti=2*3.1415926/omegaa
            periodList.append(Ti)
            print(f"{i+1}th period is: {Ti}")
        ################################################################################################################
        opsX.auxiliary_writeModelInformationToDB()  ###---将模型信息写入数据库，以便在SAPStructure中显示模型
        #################################################################################################################
        periodPath=os.path.join("Period", f"{caseNum}")
        os.makedirs(periodPath, exist_ok=True)
        opsX.auxiliary_writeDataIntoTxtFile(savePath=periodPath,filename='period',listData=[periodList],decimals=[6])

        xDamp=0.05
        MpropSwitch=1.0
        KcurrSwitch=1.0
        KcommSwitch=0.0
        KinitSwitch=0.0
        nEigenI=1
        nEigenJ=2
        lambdaI=lamiga[nEigenI-1]
        lambdaJ=lamiga[nEigenJ-1]
        omegaI=lambdaI**0.5
        omegaJ=lambdaJ**0.5
        alphaM=MpropSwitch*xDamp*(2*omegaI*omegaJ)/(omegaI+omegaJ)  # M-prop. damping; D = alphaM*M
        betaKcurr=KcurrSwitch * 2. *xDamp / (omegaI+omegaJ)  # current-K;      +beatKcurr*KCurrent
        betaKcomm=KcommSwitch * 2. *xDamp / (omegaI+omegaJ)  # last-committed K;   +betaKcomm*KlastCommitt
        betaKinit=KinitSwitch * 2. *xDamp / (omegaI+omegaJ)
        print("alphaM,betaK=",alphaM,betaKcurr)
        ### "*** Beginning Time History Analysis..."
        tsTag=10
        ops.timeSeries('Path',tsTag,'-dt',motionDT,'-filePath',f'Waves/{AccNumber}.acc','-factor',9.81*Scale)
        ops.pattern('UniformExcitation',400,1,'-accel',tsTag)
        # Load Recorder
        recDT=motionDT
        #### pile and cap response
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-10000.out','-timeSeries',tsTag,'-time',
                     '-dT',recDT,'-node',611,'-dof',1,'accel') ###cap acc. (-timeSeries: absolute = gm + relative)
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-11000.out','-time','-dT',recDT,
                     '-node',611,'-dof',1,2,6,'disp') #cap disp. and rotation
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-11010.out','-time','-dT',recDT,
                     '-node',610,'-dof',1,2,6,'disp') #cap disp. and rotation
        for row in range(N):
            pilenum_ctr=(2 *row + 1)*1000
            ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-{12000+row*20000}.out','-time',
                         '-dT',recDT,'-nodeRange',pilenum_ctr+1,pilenum_ctr+numPnode,'-dof',1,'disp')
            ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-{15000 + row * 20000}.out',
                         '-time','-dT',recDT,'-nodeRange',pilenum_ctr+1,pilenum_ctr+numPnode,'-dof',2,'disp')
            ops.recorder('Element','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-{13000 + row * 20000}.out',
                         '-time','-dT',recDT,'-eleRange',pilenum_ctr+1,pilenum_ctr+numPele,'localForce')
            ops.recorder('Element','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-{14000 +row * 20000}.out',
                         '-time','-dT',recDT,'-eleRange',pilenum_ctr+1,pilenum_ctr+numPele,'section',np,'deformation')
        ### pier and superstructure response
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-15000.out','-timeSeries',tsTag,
                     '-time','-dT',recDT,'-node',superStrucNodeNum,'-dof',1,'accel') ### superstruc. acc.
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-16000.out','-time','-dT',recDT,
                     '-node',superStrucNodeNum,'-dof',1,2,6,'disp')
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-17000.out','-time','-dT',recDT,
                     '-nodeRange',1000+numPnode+1,pierTopNodeNum,'-dof',1,'disp') ### pier node lateral disp.
        ops.recorder('Node','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-17500.out','-time','-dT',recDT,
                     '-nodeRange',1000+numPnode+1,pierTopNodeNum,'-dof',2,'disp')
        ops.recorder('Element','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-18000.out','-time','-dT',
                     recDT,'-eleRange',1000+numPnode+1,pierTopNodeNum-1,'localForce')
        ops.recorder('Element','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-19000.out','-time','-dT',recDT,
                     '-eleRange',1000+numPnode+1,pierTopNodeNum-1,'section',1,'deformation')
        ### reaction force and displacment in bearing
        ops.recorder('Element','-file',f'Data/{output}/Pile/{AccNumber}-{Scale}-19999.out','-time','-dT',recDT,
                     '-ele',10888,10999,'localForce')
        ###
        ops.constraints('Penalty',1.e20,1.e20)
        ops.test('NormDispIncr',1.0e-4,35,0)
        ops.algorithm('KrylovNewton')
        ops.numberer('RCM')
        ops.system('ProfileSPD')
        gamma=0.5
        beta=0.25
        ops.integrator('Newmark',gamma,beta)
        ops.rayleigh(alphaM,betaKcurr,betaKinit,betaKcomm)
        ops.analysis('Transient')

        tFinal=motionSteps *motionDT
        tFinal=10
        tCurrent=ops.getTime()
        ok=0
        step=1
        while (ok==0)and(tCurrent<tFinal):
            ok=ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print("Trying Broyden...")
                ops.algorithm('Broyden',8)
                ok=ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print("Trying Newton ..")
                ops.algorithm('Newton')
                ok=ops.analyze(1,motionDT)
            if ok!=0:
                ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print("Trying Newton With LineSearch ..")
                ops.algorithm('NewtonLineSearch')
                ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print("Trying Broyden-Fletcher-Goldfarb-Shanno (BFGS) ..")
                ops.algorithm('BFGS')
                ok=ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print('Trying SecantNewton ..')
                ops.algorithm('SecantNewton')
                ok=ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok!=0:
                print("Trying ModifiedNewton ..")
                ops.algorithm('ModifiedNewton','-initial')
                ok=ops.analyze(1,motionDT)
            if ok!=0:
                ok=ops.analyze(1,motionDT/2.0)
            if ok!=0:
                ok=ops.analyze(1,motionDT/4.0)
            if ok==0:
                print(f"Succeed! Step {step}/{motionSteps} of Acc {AccNumber} with Scalefactor {Scale} in Case {caseNum} ")
                ops.algorithm('KrylovNewton')
            if ok!=0:
                print(f"Disconvergence! Step {step}/{motionSteps} of Acc {AccNumber} with Scalefactor {Scale} in Case {caseNum}")
            step=step + 1
            tCurrent=ops.getTime()










































































        



