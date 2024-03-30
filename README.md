# SAPStructure (Version 0.1.0)
Seismic analysis platform for structures   
##########################################################################    
Author: Junjun Guo([HomePage](https://github.com/Junjun1guo))    
E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com      
##########################################################################
______
- [Tutorial-1:Auxiliary modules](#Tutorials-1)
- [Tutorial-2:Install SAPStructure and view OpenSeesPy model](#Tutorials-2)
- [Tutorial-3:Quickly generate model nodes and elements with AutoCAD dxf file](#Tutorial-3)
- [Tutorial-4:2D soil profile meshing  with AutoCAD dxf file](#Tutorial-4)
______
### Notes: compatible modules: records==0.5.3, sQLAlchemy==1.3.20    

## Tutorials-1
### Auxiliary modules
1. #### CalculateGroundMotionIMs module
   a class for calculate ground motion intensity measure. please use the command print(help(CalculateGroundMotionIMs)) to check the structure and the usage of the    class
2. #### GroundMotionProcess modulea 
   class for ground motion baseline correction ,fltering, and conversion among acceleration,velocity and displacement, please use the  command print(help(GroundMotionProcess)) to check the structure and the usage of the class
3. #### OpenSeesPyX module
   a class for the visualization and quick construction of OpenSeesPy model. please use the command print(help(OpenSeesPyX)) to check the structure    and the usage of    the class
4. #### SectionPropertyCalculate module
   a class for calculating the section properties. please use the command print(help(SectionPropertyCalculate)) to check the structure    and the usage of the class
5. #### SectMCAnalysis module
   a class for section moment curvature analysis. please use the command print(help(SectMCAnalysis)) to check the structure and the usage of the    class
6. #### ExciteAnyDirectionOpenSees module
   a class for horizontally rotate FE model,it is convinient to get the rotated node coordinates use this class. please use the        command  print(help(ExciteAnyDirectionOpenSees)) to check the structure and the usage of the class
7. #### PythonInteractSAP2000 module
   a class for python  interacting with the SAP2000 program.
8. #### ShakeTableTest module
   a class for processing shake table test data, such as calculating the periods and damping ratios of a structure
______
## Tutorials-2      
### Install SAPStructure and view OpenSeesPy model
1. Download the zip file
2. Run the example model (eg. exmple 1)
3. Download SAPStructure from https://fbs.sh/JunjunGuo/SAPStructure/SAPStructureSetup.exe, and install it
4. When encounter error after installation, just close the window, and reopen it.
5. Select SAPStructure and right click the mouse, then click the properties and choose running the program as an administrator.
6. Click loadResultDB button, and load the result database 
7. Then display the model and conduct post process.     
Prepare your own openseespy model by referring the examples       
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/model.JPG" width =100% height =100% div align="center">
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/modeShape.JPG" width =100% height =100% div align="center">
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/shellWall.jpg" width =100% height =100% div align="center">
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/timeHistory.JPG" width =100% height =100% div align="center">
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/hysteretic.JPG" width =100% height =100% div align="center">

______  
## Tutorial-3      
### Quickly generate model nodes and elements with AutoCAD dxf file                
1. Clicking the "ModelPath" button to select an folder to save the model database and the generated model nodes and elements.(There generates a template dxf file named "cadModelTemplate.dxf" with several specific layers,such as girder,pier,etc., plot each model part in the corresponding layer.), see Figure 2.1.
2. Setting the number of segments for arc and spline, the elements length of girder, pier and cap beam.
3. Clicking the "DXFModelLoad" button to load your model dxf file. (You can reference an example girder model dxf file constructed with the generated template file, called "girderBridgeExample.dxf", see Figure 2.2)
4. Clicking the "generateModel" button to automatically generate the model nodes and elements, and save them into the model database and txt files in the selected folder.
5. Go to the main windown of SAPStructure, clicking the button "loadModelDB" to load the generated model database, and display the model, see Figures 2.3 and 2.4. In addition, the model can be visualized with node and element tags, see Figures 2.5 and 2.6.
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/quicklyGenerateModel-1.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.1 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/dxfModel.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.2 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/quickModel.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.3 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/modelNodesElements.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.4 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/modelNodesElementsStr.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.5 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/modelNodesElementsInt.jpg" width =100% height =100% div align="center">
<p align="center">Figure 2.6 </p>      

______
## Tutorial-4      
### 2D soil profile meshing  with AutoCAD dxf file
1. Open the 2D soil profile divide window, see Figure 3.1.
2. Click "resultsSavePath" button to select an folder to save the generated nodes and elements files.
3. Click "dxfModelLoad" button to load the soil profile dxf file (see the example dxf file-soilProfile.dxf).
4. Enter the layer name of the rectangular border of the soil profile (Make sure the soil profile is within the rectangular border).
5. Specify the width and height of the divided elements.
6. Click the "meshSoilProfile" button to divide the rectangular region (see Figure 3.2).
7. Then, entering the layer name of the closed boundary of a local soil region, and enter the number of the soil region. 
8. Click "soilIdentify" button to identify the generated elements that within the specified soil region.
9. Repeat steps 7 and 8 until meshing all the soil regions (see Figure 3.3).
10. Finally, click "saveNodesEles" button to save the nodes and the elements of each soil region to the selected folder.

<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/soilProfileMeshPanel.jpg" width =100% height =100% div align="center">
<p align="center">Figure 3.1 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/soilProfileRectangleMesh.jpg" width =100% height =100% div align="center">
<p align="center">Figure 3.2 </p>
<img src="https://github.com/Junjun1guo/SAPStructure/blob/main/figures/soilProfileRegionMesh.jpg" width =100% height =100% div align="center">
<p align="center">Figure 3.3 </p>
