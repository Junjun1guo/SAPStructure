#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2023-04-03
########################################################################################################################
########################---import modules---#################################
import numpy as np
import os
import math
import pygmsh
import pyvista as pv
import matplotlib.pyplot as plt
import sys
sys.path.append("..")
########################################################################################################################
########################################################################################################################
from auxiliaryModules.mainMod import EleMeshPlotAndSelect
########################################################################################################################
########################################################################################################################
def example_1_mesh_polygon():
    """
    --------------------------------------------------------------------------------------------------------------------
    Flat shapes-polygon
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [0.0, 0.0],
                [1.0, -0.2],
                [1.1, 1.2],
                [0.1, 0.7],
            ],
            mesh_size=0.1,
        )
        mesh = geom.generate_mesh()
        meshProcessInstance = EleMeshPlotAndSelect()###---initialize class
        meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)###---add mesh to plot
        meshProcessInstance.meshProcess(showPointTag=False) ###---Process the mesh
        # meshProcessInstance.saveNodes(saveName="nodes") ###---Save nodes
        # meshProcessInstance.saveElements(saveName='elements') ###---save elements
        # meshProcessInstance.selectNodes_inLine(startNodeTag=1, endNodeTag=2, saveSetName='set1') ###---select nodes
        # meshProcessInstance.selectFaces_inPlane(planeNode1Tag=1, planeNode2Tag=2, planeNode3Tag=3, saveSetName='set4')
        meshProcessInstance.meshPlot() ###---visualize the meshes
########################################################################################################################
########################################################################################################################
def example_2_mesh_circle():
    """
    --------------------------------------------------------------------------------------------------------------------
    Flat shapes-circle
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_circle([0.0, 0.0], 1.0, mesh_size=0.2)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.selectNodes_betweenTwoConcentricCylinders(circleCenter=[0,0],radius1=0.99,radius2=1.01,Zmin=0,
                                                                  Zmax=0,saveSetName="circleNodes")
    meshProcessInstance.meshPlot()  ###---visualize the meshes
########################################################################################################################
########################################################################################################################
def example_3_mesh_BSpline():
    """
    --------------------------------------------------------------------------------------------------------------------
    Flat shapes-(B-)Splines
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        lcar = 0.1
        p1 = geom.add_point([0.0, 0.0], lcar)
        p2 = geom.add_point([1.0, 0.0], lcar)
        p3 = geom.add_point([1.0, 0.5], lcar)
        p4 = geom.add_point([1.0, 1.0], lcar)
        s1 = geom.add_bspline([p1, p2, p3, p4])

        p2 = geom.add_point([0.0, 1.0], lcar)
        p3 = geom.add_point([0.5, 1.0], lcar)
        s2 = geom.add_spline([p4, p3, p2, p1])

        ll = geom.add_curve_loop([s1, s2])
        pl = geom.add_plane_surface(ll)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes
########################################################################################################################
########################################################################################################################
def example_4_mesh_extrude():
    """
    --------------------------------------------------------------------------------------------------------------------
    Extrude
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [
                [0.0, 0.0],
                [1.0, -0.2],
                [1.1, 1.2],
                [0.1, 0.7],
            ],
            mesh_size=0.1,
        )
        geom.extrude(poly, [0.0, 0.3, 1.0], num_layers=5)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_5_mesh_revolve():
    """
    --------------------------------------------------------------------------------------------------------------------
    Revolve
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [
                [0.0, 0.2, 0.0],
                [0.0, 1.2, 0.0],
                [0.0, 1.2, 1.0],
            ],
            mesh_size=0.1,
        )
        geom.revolve(poly, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.8 * pi)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_6_mesh_twist():
    """
    --------------------------------------------------------------------------------------------------------------------
    Twist
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [
                [+0.0, +0.5],
                [-0.1, +0.1],
                [-0.5, +0.0],
                [-0.1, -0.1],
                [+0.0, -0.5],
                [+0.1, -0.1],
                [+0.5, +0.0],
                [+0.1, +0.1],
            ],
            mesh_size=0.05,
        )

        geom.twist(
            poly,
            translation_axis=[0, 0, 1],
            rotation_axis=[0, 0, 1],
            point_on_axis=[0, 0, 0],
            angle=pi / 3,
        )
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_7_mesh_BooleanIntersection():
    """
    --------------------------------------------------------------------------------------------------------------------
    Boolean_intersection
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi, cos

    with pygmsh.occ.Geometry() as geom: ####Gmsh also supports OpenCASCADE (occ), allowing for a CAD-style geometry specification.
        geom.characteristic_length_max = 0.1
        r = 0.5
        disks = [
            geom.add_disk([-0.5 * cos(7 / 6 * pi), -0.25], 1.0),
            geom.add_disk([+0.5 * cos(7 / 6 * pi), -0.25], 1.0),
            geom.add_disk([0.0, 0.5], 1.0),
        ]
        geom.boolean_intersection(disks)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_8_mesh_BooleanDifference():
    """
    --------------------------------------------------------------------------------------------------------------------
    Boolean_difference(subtract)
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = 0.1
        ellipsoid = geom.add_ellipsoid([0.0, 0.0, 0.0], [1.0, 0.7, 0.5])

        cylinders = [
            geom.add_cylinder([-1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.3),
            geom.add_cylinder([0.0, -1.0, 0.0], [0.0, 2.0, 0.0], 0.3),
            geom.add_cylinder([0.0, 0.0, -1.0], [0.0, 0.0, 2.0], 0.3),
        ]
        geom.boolean_difference(ellipsoid, geom.boolean_union(cylinders))

        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_9_mesh_BooleanUnion():
    """
    --------------------------------------------------------------------------------------------------------------------
    Boolean_difference and union
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.1
        geom.characteristic_length_max = 0.1

        rectangle = geom.add_rectangle([-1.0, -1.0, 0.0], 2.0, 2.0)
        disk1 = geom.add_disk([-1.2, 0.0, 0.0], 0.5)
        disk2 = geom.add_disk([+1.2, 0.0, 0.0], 0.5)

        disk3 = geom.add_disk([0.0, -0.9, 0.0], 0.5)
        disk4 = geom.add_disk([0.0, +0.9, 0.0], 0.5)
        flat = geom.boolean_difference(
            geom.boolean_union([rectangle, disk1, disk2]),
            geom.boolean_union([disk3, disk4]),
        )

        geom.extrude(flat, [0, 0, 0.3])

        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_10_mesh_refinement():
    """
    --------------------------------------------------------------------------------------------------------------------
    Mesh refinement
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [
                [0.0, 0.0],
                [2.0, 0.0],
                [3.0, 1.0],
                [1.0, 2.0],
                [0.0, 1.0],
            ],
            mesh_size=0.3,
        )

        field0 = geom.add_boundary_layer(
            edges_list=[poly.curves[0]],
            lcmin=0.05,
            lcmax=0.2,
            distmin=0.0,
            distmax=0.2,
        )
        field1 = geom.add_boundary_layer(
            nodes_list=[poly.points[2]],
            lcmin=0.05,
            lcmax=0.2,
            distmin=0.1,
            distmax=0.4,
        )
        geom.set_background_mesh([field0, field1], operator="Min")

        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_11_mesh_airfoil():
    """
    --------------------------------------------------------------------------------------------------------------------
    Airfoil
    --------------------------------------------------------------------------------------------------------------------
    """
    airfoil_coordinates = np.array(
        [
            [1.000000, 0.000000, 0.0],
            [0.999023, 0.000209, 0.0],
            [0.996095, 0.000832, 0.0],
            [0.991228, 0.001863, 0.0],
            [0.984438, 0.003289, 0.0],
            [0.975752, 0.005092, 0.0],
            [0.965201, 0.007252, 0.0],
            [0.952825, 0.009744, 0.0],
            [0.938669, 0.012538, 0.0],
            [0.922788, 0.015605, 0.0],
            [0.905240, 0.018910, 0.0],
            [0.886092, 0.022419, 0.0],
            [0.865417, 0.026096, 0.0],
            [0.843294, 0.029903, 0.0],
            [0.819807, 0.033804, 0.0],
            [0.795047, 0.037760, 0.0],
            [0.769109, 0.041734, 0.0],
            [0.742094, 0.045689, 0.0],
            [0.714107, 0.049588, 0.0],
            [0.685258, 0.053394, 0.0],
            [0.655659, 0.057071, 0.0],
            [0.625426, 0.060584, 0.0],
            [0.594680, 0.063897, 0.0],
            [0.563542, 0.066977, 0.0],
            [0.532136, 0.069789, 0.0],
            [0.500587, 0.072303, 0.0],
            [0.469022, 0.074486, 0.0],
            [0.437567, 0.076312, 0.0],
            [0.406350, 0.077752, 0.0],
            [0.375297, 0.078743, 0.0],
            [0.344680, 0.079180, 0.0],
            [0.314678, 0.079051, 0.0],
            [0.285418, 0.078355, 0.0],
            [0.257025, 0.077096, 0.0],
            [0.229618, 0.075287, 0.0],
            [0.203313, 0.072945, 0.0],
            [0.178222, 0.070096, 0.0],
            [0.154449, 0.066770, 0.0],
            [0.132094, 0.063005, 0.0],
            [0.111248, 0.058842, 0.0],
            [0.091996, 0.054325, 0.0],
            [0.074415, 0.049504, 0.0],
            [0.058573, 0.044427, 0.0],
            [0.044532, 0.039144, 0.0],
            [0.032343, 0.033704, 0.0],
            [0.022051, 0.028152, 0.0],
            [0.013692, 0.022531, 0.0],
            [0.007292, 0.016878, 0.0],
            [0.002870, 0.011224, 0.0],
            [0.000439, 0.005592, 0.0],
            [0.000000, 0.000000, 0.0],
            [0.001535, -0.005395, 0.0],
            [0.005015, -0.010439, 0.0],
            [0.010421, -0.015126, 0.0],
            [0.017725, -0.019451, 0.0],
            [0.026892, -0.023408, 0.0],
            [0.037880, -0.026990, 0.0],
            [0.050641, -0.030193, 0.0],
            [0.065120, -0.033014, 0.0],
            [0.081257, -0.035451, 0.0],
            [0.098987, -0.037507, 0.0],
            [0.118239, -0.039185, 0.0],
            [0.138937, -0.040493, 0.0],
            [0.161004, -0.041444, 0.0],
            [0.184354, -0.042054, 0.0],
            [0.208902, -0.042343, 0.0],
            [0.234555, -0.042335, 0.0],
            [0.261221, -0.042058, 0.0],
            [0.288802, -0.041541, 0.0],
            [0.317197, -0.040817, 0.0],
            [0.346303, -0.039923, 0.0],
            [0.376013, -0.038892, 0.0],
            [0.406269, -0.037757, 0.0],
            [0.437099, -0.036467, 0.0],
            [0.468187, -0.035009, 0.0],
            [0.499413, -0.033414, 0.0],
            [0.530654, -0.031708, 0.0],
            [0.561791, -0.029917, 0.0],
            [0.592701, -0.028066, 0.0],
            [0.623264, -0.026176, 0.0],
            [0.653358, -0.024269, 0.0],
            [0.682867, -0.022360, 0.0],
            [0.711672, -0.020466, 0.0],
            [0.739659, -0.018600, 0.0],
            [0.766718, -0.016774, 0.0],
            [0.792738, -0.014999, 0.0],
            [0.817617, -0.013284, 0.0],
            [0.841253, -0.011637, 0.0],
            [0.863551, -0.010068, 0.0],
            [0.884421, -0.008583, 0.0],
            [0.903777, -0.007191, 0.0],
            [0.921540, -0.005900, 0.0],
            [0.937637, -0.004717, 0.0],
            [0.952002, -0.003650, 0.0],
            [0.964576, -0.002708, 0.0],
            [0.975305, -0.001896, 0.0],
            [0.984145, -0.001222, 0.0],
            [0.991060, -0.000691, 0.0],
            [0.996020, -0.000308, 0.0],
            [0.999004, -0.000077, 0.0],
        ]
    )

    # Scale airfoil to input coord
    coord = 1.0
    airfoil_coordinates *= coord

    # Instantiate geometry object
    with pygmsh.geo.Geometry() as geom:
        # Create polygon for airfoil
        char_length = 1.0e-1
        airfoil = geom.add_polygon(airfoil_coordinates, char_length, make_surface=False)

        # Create surface for numerical domain with an airfoil-shaped hole
        left_dist = 1.0
        right_dist = 3.0
        top_dist = 1.0
        bottom_dist = 1.0
        xmin = airfoil_coordinates[:, 0].min() - left_dist * coord
        xmax = airfoil_coordinates[:, 0].max() + right_dist * coord
        ymin = airfoil_coordinates[:, 1].min() - bottom_dist * coord
        ymax = airfoil_coordinates[:, 1].max() + top_dist * coord
        domainCoordinates = np.array(
            [[xmin, ymin, 0.0], [xmax, ymin, 0.0], [xmax, ymax, 0.0], [xmin, ymax, 0.0]]
        )
        polygon = geom.add_polygon(domainCoordinates, char_length, holes=[airfoil])
        geom.set_recombined_surfaces([polygon.surface])

        ref = 10.525891646546
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_12_mesh_bsplines():
    """
    --------------------------------------------------------------------------------------------------------------------
    bsplines
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        lcar = 0.1
        p1 = geom.add_point([0.0, 0.0, 0.0], lcar)
        p2 = geom.add_point([1.0, 0.0, 0.0], lcar)
        p3 = geom.add_point([1.0, 0.5, 0.0], lcar)
        p4 = geom.add_point([1.0, 1.0, 0.0], lcar)
        s1 = geom.add_bspline([p1, p2, p3, p4])

        p2 = geom.add_point([0.0, 1.0, 0.0], lcar)
        p3 = geom.add_point([0.5, 1.0, 0.0], lcar)
        s2 = geom.add_bspline([p4, p3, p2, p1])

        ll = geom.add_curve_loop([s1, s2])
        pl = geom.add_plane_surface(ll)
        mesh = geom.generate_mesh(verbose=True) ###---verbose=True(print full mesh information)
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_13_mesh_circle():
    """
    --------------------------------------------------------------------------------------------------------------------
    circle
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_circle(
            [0.0, 0.0, 0.0],
            1.0,
            mesh_size=0.1,
            num_sections=4,
            # If compound==False, the section borders have to be points of the
            # discretization. If using a compound circle, they don't; gmsh can
            # choose by itself where to point the circle points.
            compound=True,
        )
        # geom.add_physical(c.plane_surface, "super disk")
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_14_mesh_circleTransform():
    """
    --------------------------------------------------------------------------------------------------------------------
    circle transform
    --------------------------------------------------------------------------------------------------------------------
    """
    radius=1.0
    with pygmsh.geo.Geometry() as geom:
        R = [
            pygmsh.rotation_matrix(np.eye(1, 3, d)[0], theta)
            for d, theta in enumerate(np.pi / np.array([2.0, 3.0, 5]))
        ]
        geom.add_circle([7.0, 11.0, 13.0], radius, 0.1, R[0] @ R[1] @ R[2])
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_15_mesh_cube():
    """
    --------------------------------------------------------------------------------------------------------------------
    cube
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_box(0, 1, 0, 1, 0, 1, 0.1)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_16_mesh_ellipsoid():
    """
    --------------------------------------------------------------------------------------------------------------------
    ellipsoid
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_ellipsoid([0.0, 0.0, 0.0], [1.0, 0.5, 0.75], 0.05)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_17_mesh_embed():
    """
    --------------------------------------------------------------------------------------------------------------------
    embed
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [
                [0, 0.3],
                [0, 1.1],
                [0.9, 1.1],
                [0.9, 0.3],
                [0.6, 0.7],
                [0.3, 0.7],
                [0.2, 0.4],
            ],
            mesh_size=[0.2, 0.2, 0.2, 0.2, 0.03, 0.03, 0.01],
        )
        geom.in_surface(poly.lines[4], poly)
        geom.in_surface(poly.points[6], poly)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_18_mesh_hex():
    """
    --------------------------------------------------------------------------------------------------------------------
    hex
    --------------------------------------------------------------------------------------------------------------------
    """
    lcar=1.0
    with pygmsh.geo.Geometry() as geom:
        lbw = [2, 3, 5]
        points = [geom.add_point([x, 0.0, 0.0], lcar) for x in [0.0, lbw[0]]]
        line = geom.add_line(*points)

        _, rectangle, _ = geom.extrude(
            line, translation_axis=[0.0, lbw[1], 0.0], num_layers=lbw[1], recombine=True
        )
        geom.extrude(
            rectangle,
            translation_axis=[0.0, 0.0, lbw[2]],
            num_layers=lbw[2],
            recombine=True,
        )
        # compute_volume only supports 3D for tetras, but does return surface area for
        # quads
        mesh = geom.generate_mesh()
        # mesh.remove_lower_dimensional_cells()
        # mesh.remove_orphaned_nodes()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=1, elementStartNumber=1)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_19_mesh_holeInSquareHighOrder():
    """
    --------------------------------------------------------------------------------------------------------------------
    Hole in square with higher order triangle
    --------------------------------------------------------------------------------------------------------------------
    """
    # Characteristic length
    lcar = 0.2

    # Coordinates of lower-left and upper-right vertices of a square domain
    xmin = 0.0
    xmax = 5.0
    ymin = 0.0
    ymax = 5.0

    # Vertices of a square hole
    squareHoleCoordinates = np.array([[1.0, 1.0], [4.0, 1.0], [4.0, 4.0], [1.0, 4.0]])

    with pygmsh.geo.Geometry() as geom:
        # Create square hole
        squareHole = geom.add_polygon(squareHoleCoordinates, lcar, make_surface=False)
        # Create square domain with square hole
        geom.add_rectangle(
            xmin, xmax, ymin, ymax, 0.0, lcar, holes=[squareHole.curve_loop]
        )
        mesh = geom.generate_mesh(order=2)
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_20_mesh_layers():
    """
    --------------------------------------------------------------------------------------------------------------------
    Layers
    --------------------------------------------------------------------------------------------------------------------
    """
    mesh_size=0.05
    with pygmsh.geo.Geometry() as geom:
        # Draw a cross with a circular hole
        circ = geom.add_circle(
            [0.0, 0.0, 0.0], 0.1, mesh_size=mesh_size, make_surface=False
        )
        poly = geom.add_polygon(
            [
                [+0.0, +0.5, 0.0],
                [-0.1, +0.1, 0.0],
                [-0.5, +0.0, 0.0],
                [-0.1, -0.1, 0.0],
                [+0.0, -0.5, 0.0],
                [+0.1, -0.1, 0.0],
                [+0.5, +0.0, 0.0],
                [+0.1, +0.1, 0.0],
            ],
            mesh_size=mesh_size,
            holes=[circ],
        )
        axis = [0, 0, 1.0]
        geom.extrude(poly, translation_axis=axis, num_layers=1)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_21_mesh_pacman():
    """
    --------------------------------------------------------------------------------------------------------------------
    Pacman
    --------------------------------------------------------------------------------------------------------------------
    """
    from numpy import cos, pi, sin
    lcar = 0.3
    with pygmsh.geo.Geometry() as geom:
        r = 1.25 * 3.4
        p1 = geom.add_point([0.0, 0.0, 0.0], lcar)
        # p2 = geom.add_point([+r, 0.0, 0.0], lcar)
        p3 = geom.add_point([-r, 0.0, 0.0], lcar)
        p4 = geom.add_point([0.0, +r, 0.0], lcar)
        p5 = geom.add_point([0.0, -r, 0.0], lcar)
        p6 = geom.add_point([r * cos(+pi / 12.0), r * sin(+pi / 12.0), 0.0], lcar)
        p7 = geom.add_point([r * cos(-pi / 12.0), r * sin(-pi / 12.0), 0.0], lcar)
        p8 = geom.add_point([0.5 * r, 0.0, 0.0], lcar)

        c0 = geom.add_circle_arc(p6, p1, p4)
        c1 = geom.add_circle_arc(p4, p1, p3)
        c2 = geom.add_circle_arc(p3, p1, p5)
        c3 = geom.add_circle_arc(p5, p1, p7)
        l1 = geom.add_line(p7, p8)
        l2 = geom.add_line(p8, p6)
        ll = geom.add_curve_loop([c0, c1, c2, c3, l1, l2])

        pacman = geom.add_plane_surface(ll)

        # test setting physical groups
        geom.add_physical(p1, label="c")
        geom.add_physical(c0, label="arc")
        geom.add_physical(pacman, "dummy")
        geom.add_physical(pacman, label="77")

        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_22_mesh_physical():
    """
    --------------------------------------------------------------------------------------------------------------------
    physical
    --------------------------------------------------------------------------------------------------------------------
    """
    lcar = 0.1
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], lcar)

        top, volume, lat = geom.extrude(poly, [0, 0, 2])

        geom.add_physical(poly, label="bottom")
        geom.add_physical(top, label="top")
        geom.add_physical(volume, label="volume")
        geom.add_physical(lat, label="lat")
        geom.add_physical(poly.lines[0], label="line")

        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_23_mesh_pipes():
    """
    --------------------------------------------------------------------------------------------------------------------
    pipes
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        sqrt2on2 = 0.5 * np.sqrt(2.0)
        R = pygmsh.rotation_matrix([sqrt2on2, sqrt2on2, 0], np.pi / 6.0)
        geom.add_pipe(
            inner_radius=0.3, outer_radius=0.4, length=1.0, R=R, mesh_size=0.04
        )

        R = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        geom.add_pipe(
            inner_radius=0.3,
            outer_radius=0.4,
            length=1.0,
            mesh_size=0.04,
            R=R,
            variant="circle_extrusion",
        )
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_24_mesh_quads():
    """
    --------------------------------------------------------------------------------------------------------------------
    quads
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        rectangle = geom.add_rectangle(0.0, 1.0, 0.0, 1.0, 0.0, 0.1)
        geom.set_recombined_surfaces([rectangle.surface])
        mesh = geom.generate_mesh(dim=2)

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_25_mesh_recombine():
    """
    --------------------------------------------------------------------------------------------------------------------
    recombine
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        pts = [
            geom.add_point((0.0, 0.0, 0.0), mesh_size=1.0),
            geom.add_point((2.0, 0.0, 0.0), mesh_size=1.0),
            geom.add_point((0.0, 1.0, 0.0), mesh_size=1.0),
            geom.add_point((2.0, 1.0, 0.0), mesh_size=1.0),
        ]
        lines = [
            geom.add_line(pts[0], pts[1]),
            geom.add_line(pts[1], pts[3]),
            geom.add_line(pts[3], pts[2]),
            geom.add_line(pts[2], pts[0]),
        ]
        ll0 = geom.add_curve_loop(lines)
        rs0 = geom.add_surface(ll0)

        geom.set_transfinite_curve(lines[3], 3, "Progression", 1.0)
        geom.set_transfinite_curve(lines[1], 3, "Progression", 1.0)
        geom.set_transfinite_curve(lines[2], 3, "Progression", 1.0)
        geom.set_transfinite_curve(lines[0], 3, "Progression", 1.0)
        geom.set_transfinite_surface(rs0, "Left", pts)
        geom.set_recombined_surfaces([rs0])

        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_26_mesh_rectangle():
    """
    --------------------------------------------------------------------------------------------------------------------
    rectangle
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        geom.add_rectangle(0.0, 1.0, 0.0, 1.0, 0.0, 0.1)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_27_mesh_regularExtrusion():
    """
    --------------------------------------------------------------------------------------------------------------------
    regularExtrusion
    --------------------------------------------------------------------------------------------------------------------
    """
    x = 5
    y = 4
    z = 3
    x_layers = 10
    y_layers = 5
    z_layers = 3

    with pygmsh.geo.Geometry() as geom:
        p = geom.add_point([0, 0, 0], 1)
        _, l, _ = geom.extrude(p, [x, 0, 0], num_layers=x_layers)

        _, s, _ = geom.extrude(l, [0, y, 0], num_layers=y_layers)
        geom.extrude(s, [0, 0, z], num_layers=z_layers)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_28_mesh_rotatedLayers():
    """
    --------------------------------------------------------------------------------------------------------------------
    rotated layers
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    mesh_size = 0.05
    with pygmsh.geo.Geometry() as geom:
        # Draw a square
        poly = geom.add_polygon(
            [
                [+0.5, +0.0, 0.0],
                [+0.0, +0.5, 0.0],
                [-0.5, +0.0, 0.0],
                [+0.0, -0.5, 0.0],
            ],
            mesh_size=mesh_size,
        )
        axis = [0, 0, 1.0]
        geom.twist(
            poly,
            translation_axis=axis,
            rotation_axis=axis,
            point_on_axis=[0.0, 0.0, 0.0],
            angle=0.5 * pi,
            num_layers=5,
            recombine=True,
        )
        mesh = geom.generate_mesh()
    print(mesh.cells_dict)

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_29_mesh_rotation():
    """
    --------------------------------------------------------------------------------------------------------------------
    rotation
    --------------------------------------------------------------------------------------------------------------------
    """
    angle = np.pi / 5

    # Generate reference geometry
    with pygmsh.geo.Geometry() as geom:
        rect = geom.add_rectangle(0.0, 2.0, 0.0, 1.0, 0.0, 0.1)
        mesh_unrot = geom.generate_mesh()
    vertex_index = mesh_unrot.cells_dict["vertex"]
    vertex_index = vertex_index.reshape((vertex_index.shape[0],))

    with pygmsh.geo.Geometry() as geom:
        # Generate rotated geometry
        geom = pygmsh.geo.Geometry()
        rect = geom.add_rectangle(0.0, 2.0, 0.0, 1.0, 0.0, 0.1)
        geom.rotate(rect.surface, (0, 0, 0), angle, (0, 0, 1))
        mesh = geom.generate_mesh()

    new_vertex_index = mesh.cells_dict["vertex"]
    new_vertex_index = new_vertex_index.reshape((new_vertex_index.shape[0],))

    # Generate rotation matrix and compare with rotated geometry
    Rm = pygmsh.helpers.rotation_matrix([0, 0, 1], angle)
    for v, v_new in zip(vertex_index, new_vertex_index):
        point = mesh_unrot.points[v, :]
        rot_point = np.dot(Rm, point)
        new_point = mesh.points[v, :]
        assert np.allclose(rot_point, new_point)

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_30_mesh_screw():
    """
    --------------------------------------------------------------------------------------------------------------------
    screw
    --------------------------------------------------------------------------------------------------------------------
    """
    mesh_size = 0.05
    with pygmsh.geo.Geometry() as geom:
        # Draw a cross with a circular hole
        circ = geom.add_circle([0.0, 0.0], 0.1, mesh_size=mesh_size)
        poly = geom.add_polygon(
            [
                [+0.0, +0.5],
                [-0.1, +0.1],
                [-0.5, +0.0],
                [-0.1, -0.1],
                [+0.0, -0.5],
                [+0.1, -0.1],
                [+0.5, +0.0],
                [+0.1, +0.1],
            ],
            mesh_size=mesh_size,
            holes=[circ],
        )

        geom.twist(
            poly,
            translation_axis=[0.0, 0.0, 1.0],
            rotation_axis=[0.0, 0.0, 1.0],
            point_on_axis=[0.0, 0.0, 0.0],
            angle=2.0 / 6.0 * np.pi,
        )
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_31_mesh_subdomains():
    """
    --------------------------------------------------------------------------------------------------------------------
    subdomains
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        lcar = 0.1
        circle = geom.add_circle([0.5, 0.5, 0.0], 1.0, lcar)
        triangle = geom.add_polygon(
            [[2.0, -0.5, 0.0], [4.0, -0.5, 0.0], [4.0, 1.5, 0.0]], lcar
        )
        rectangle = geom.add_rectangle(4.75, 6.25, -0.24, 1.25, 0.0, lcar)
        # hold all domain
        geom.add_polygon(
            [
                [-1.0, -1.0, 0.0],
                [+7.0, -1.0, 0.0],
                [+7.0, +2.0, 0.0],
                [-1.0, +2.0, 0.0],
            ],
            lcar,
            holes=[circle.curve_loop, triangle.curve_loop, rectangle.curve_loop],
        )
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_32_mesh_swissCheese():
    """
    --------------------------------------------------------------------------------------------------------------------
    swiss cheese
    --------------------------------------------------------------------------------------------------------------------
    """
    X0 = np.array(
        [[+0.0, +0.0, 0.0], [+0.5, +0.3, 0.1], [-0.5, +0.3, 0.1], [+0.5, -0.3, 0.1]]
    )
    R = np.array([0.1, 0.2, 0.1, 0.14])

    with pygmsh.geo.Geometry() as geom:
        holes = [
            geom.add_ball(x0, r, with_volume=False, mesh_size=0.2 * r).surface_loop
            for x0, r in zip(X0, R)
        ]
        # geom.add_box(
        #         -1, 1,
        #         -1, 1,
        #         -1, 1,
        #         mesh_size=0.2,
        #         holes=holes
        #         )
        geom.add_ball([0, 0, 0], 1.0, mesh_size=0.2, holes=holes)
        # geom.add_physical_volume(ball, label="cheese")
        mesh = geom.generate_mesh(algorithm=5)

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_33_mesh_symmetrize():
    """
    --------------------------------------------------------------------------------------------------------------------
    symmetrize
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [[0.0, 0.5], [1.0, 0.5], [1.0, 1.0], [0.0, 1.0]],
            mesh_size=0.05,
        )
        cp = geom.copy(poly)
        geom.symmetrize(cp, [0.0, 1.0, 0.0, -0.5])
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_34_mesh_tori():
    """
    --------------------------------------------------------------------------------------------------------------------
    tori
    --------------------------------------------------------------------------------------------------------------------
    """
    irad = 0.05
    orad = 0.6
    with pygmsh.geo.Geometry() as geom:
        R = pygmsh.rotation_matrix([1.0, 0.0, 0.0], np.pi / 2)
        geom.add_torus(irad=irad, orad=orad, mesh_size=0.03, x0=[0.0, 0.0, -1.0], R=R)

        R = pygmsh.rotation_matrix([0.0, 1.0, 0.0], np.pi / 2)
        geom.add_torus(
            irad=irad,
            orad=orad,
            mesh_size=0.03,
            x0=[0.0, 0.0, 1.0],
            variant="extrude_circle",
        )
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_35_mesh_torus():
    """
    --------------------------------------------------------------------------------------------------------------------
    torus
    --------------------------------------------------------------------------------------------------------------------
    """
    irad = 0.05
    orad = 0.6
    with pygmsh.geo.Geometry() as geom:
        R = pygmsh.rotation_matrix([1.0, 0.0, 0.0], np.pi / 2)
        geom.add_torus(irad=irad, orad=orad, mesh_size=0.03, x0=[0.0, 0.0, -1.0], R=R)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_36_mesh_torusCrowd():
    """
    --------------------------------------------------------------------------------------------------------------------
    torus corwd
    --------------------------------------------------------------------------------------------------------------------
    """
    # internal radius of torus
    irad = 0.15
    # external radius of torus
    orad = 0.27

    Z_pos = (irad + orad) * np.concatenate(
        [+np.ones(8), -np.ones(8), +np.ones(8), -np.ones(8)]
    )

    Alpha = np.concatenate(
        [
            np.arange(8) * np.pi / 4.0,
            np.arange(8) * np.pi / 4.0 + np.pi / 16.0,
            np.arange(8) * np.pi / 4.0,
            np.arange(8) * np.pi / 4.0 + np.pi / 16.0,
        ]
    )

    A1 = (
            (irad + orad)
            / np.tan(np.pi / 8.0)
            * np.concatenate(
        [1.6 * np.ones(8), 1.6 * np.ones(8), 1.9 * np.ones(8), 1.9 * np.ones(8)]
    )
    )

    with pygmsh.geo.Geometry() as geom:
        for alpha, a1, z in zip(Alpha, A1, Z_pos):
            # Rotate torus to the y-z-plane.
            R1 = pygmsh.rotation_matrix([0.0, 1.0, 0.0], 0.5 * np.pi)
            R2 = pygmsh.rotation_matrix([0.0, 0.0, 1.0], alpha)
            x0 = np.array([a1, 0.0, 0.0])
            x1 = np.array([0.0, 0.0, z])
            # First rotate to y-z-plane, then move out to a1, rotate by angle
            # alpha, move up by z.
            #
            #    xnew = R2*(R1*x+x0) + x1
            #
            geom.add_torus(
                irad=irad,
                orad=orad,
                mesh_size=0.1,
                R=np.dot(R2, R1),
                x0=np.dot(R2, x0) + x1,
            )
        geom.add_box(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, mesh_size=0.3)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_37_mesh_transfinite():
    """
    --------------------------------------------------------------------------------------------------------------------
    transfinite
    --------------------------------------------------------------------------------------------------------------------
    """
    lcar = 0.1
    with pygmsh.geo.Geometry() as geom:
        poly = geom.add_polygon(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], lcar
        )
        geom.set_transfinite_surface(poly, "Left", corner_pts=[])
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_38_mesh_unorderedUnoriented():
    """
    --------------------------------------------------------------------------------------------------------------------
    unordered unoriented
    --------------------------------------------------------------------------------------------------------------------
    """
    import random
    with pygmsh.geo.Geometry() as geom:
        # Generate an approximation of a circle
        t = np.arange(0, 2.0 * np.pi, 0.05)
        x = np.column_stack([np.cos(t), np.sin(t), np.zeros_like(t)])
        points = [geom.add_point(p) for p in x]

        # Shuffle the orientation of lines by point order
        o = [0 if k % 3 == 0 else 1 for k in range(len(points))]

        lines = [
            geom.add_line(points[k + o[k]], points[k + (o[k] + 1) % 2])
            for k in range(len(points) - 1)
        ]
        lines.append(geom.add_line(points[-1], points[0]))

        # Shuffle the order of lines
        random.seed(1)
        random.shuffle(lines)

        oriented_lines = pygmsh.orient_lines(lines)
        ll = geom.add_curve_loop(oriented_lines)
        geom.add_plane_surface(ll)

        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_39_mesh_BallWithStick():
    """
    --------------------------------------------------------------------------------------------------------------------
    ball with stick
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.1
        geom.characteristic_length_max = 0.1

        ball = geom.add_ball([0.0, 0.0, 0.0], 1.0)
        box1 = geom.add_box([0, 0, 0], [1, 1, 1])
        box2 = geom.add_box([-2, -0.5, -0.5], [1.5, 0.8, 0.8])

        cut = geom.boolean_difference(ball, box1)
        frag = geom.boolean_fragments(cut, box2)

        # The three fragments are:
        # frag[0]: The ball with two cuts
        # frag[1]: The intersection of the stick and the ball
        # frag[2]: The stick without the ball
        geom.add_physical([frag[0], frag[1]], label="Sphere cut by box 1")
        geom.add_physical(frag[2], label="Box 2 cut by sphere")

        mesh = geom.generate_mesh(algorithm=6)

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_40_mesh_logo():
    """
    --------------------------------------------------------------------------------------------------------------------
    logo
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        # test setters, getters
        print(geom.characteristic_length_min)
        print(geom.characteristic_length_max)
        geom.characteristic_length_min = 2.0
        geom.characteristic_length_max = 2.0

        rect1 = geom.add_rectangle([10.0, 0.0, 0.0], 20.0, 40.0, corner_radius=5.0)
        rect2 = geom.add_rectangle([0.0, 10.0, 0.0], 40.0, 20.0, corner_radius=5.0)
        disk1 = geom.add_disk([14.5, 35.0, 0.0], 1.85)
        disk2 = geom.add_disk([25.5, 5.0, 0.0], 1.85)

        rect3 = geom.add_rectangle([10.0, 30.0, 0.0], 10.0, 1.0)
        rect4 = geom.add_rectangle([20.0, 9.0, 0.0], 10.0, 1.0)

        r1 = geom.add_rectangle([9.0, 0.0, 0.0], 21.0, 20.5, corner_radius=8.0)
        r2 = geom.add_rectangle([10.0, 00.0, 0.0], 20.0, 19.5, corner_radius=7.0)
        diff1 = geom.boolean_difference(r1, r2)
        r22 = geom.add_rectangle([9.0, 10.0, 0.0], 11.0, 11.0)
        inter1 = geom.boolean_intersection([diff1, r22])

        r3 = geom.add_rectangle([10.0, 19.5, 0.0], 21.0, 21.0, corner_radius=8.0)
        r4 = geom.add_rectangle([10.0, 20.5, 0.0], 20.0, 20.0, corner_radius=7.0)
        diff2 = geom.boolean_difference(r3, r4)
        r33 = geom.add_rectangle([20.0, 19.0, 0.0], 11.0, 11.0)
        inter2 = geom.boolean_intersection([diff2, r33])

        geom.boolean_difference(
            geom.boolean_union([rect1, rect2]),
            geom.boolean_union([disk1, disk2, rect3, rect4, inter1, inter2]),
        )

        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_41_mesh_ball():
    """
    --------------------------------------------------------------------------------------------------------------------
    ball
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.add_ball([0.0, 0.0, 0.0], 1.0, mesh_size=0.1)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_42_mesh_booleanUnion():
    """
    --------------------------------------------------------------------------------------------------------------------
    boolean union
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.1
        geom.characteristic_length_max = 0.1
        rectangle = geom.add_rectangle([-1.0, -1.0, 0.0], 2.0, 2.0)
        disk_w = geom.add_disk([-1.0, 0.0, 0.0], 0.5)
        disk_e = geom.add_disk([+1.0, 0.0, 0.0], 0.5)
        geom.boolean_union([rectangle, disk_w, disk_e])
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_43_mesh_booleanIntersection():
    """
    --------------------------------------------------------------------------------------------------------------------
    boolean intersection
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        angles = [math.pi * 3 / 6, math.pi * 7 / 6, math.pi * 11 / 6]
        disks = [
            geom.add_disk([math.cos(angles[0]), math.sin(angles[0]), 0.0], 1.5),
            geom.add_disk([math.cos(angles[1]), math.sin(angles[1]), 0.0], 1.5),
            geom.add_disk([math.cos(angles[2]), math.sin(angles[2]), 0.0], 1.5),
        ]
        geom.boolean_intersection(disks)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_44_mesh_booleanDifference():
    """
    --------------------------------------------------------------------------------------------------------------------
    boolean difference
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.1
        geom.characteristic_length_max = 0.1
        rectangle = geom.add_rectangle([-1.0, -1.0, 0.0], 2.0, 2.0)
        disk_w = geom.add_disk([-1.0, 0.0, 0.0], 0.5)
        disk_e = geom.add_disk([+1.0, 0.0, 0.0], 0.5)
        geom.boolean_union([disk_w, disk_e])
        geom.boolean_difference(rectangle, geom.boolean_union([disk_w, disk_e]))
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_45_mesh_booleanAll():
    """
    --------------------------------------------------------------------------------------------------------------------
    boolean all
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.1
        geom.characteristic_length_max = 0.1

        rectangle = geom.add_rectangle([-1.0, -1.0, 0.0], 2.0, 2.0)
        disk1 = geom.add_disk([-1.0, 0.0, 0.0], 0.5)
        disk2 = geom.add_disk([+1.0, 0.0, 0.0], 0.5)
        union = geom.boolean_union([rectangle, disk1, disk2])

        disk3 = geom.add_disk([0.0, -1.0, 0.0], 0.5)
        disk4 = geom.add_disk([0.0, +1.0, 0.0], 0.5)
        geom.boolean_difference(union, geom.boolean_union([disk3, disk4]))
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_46_mesh_booleansSquareCircleSlice():
    """
    --------------------------------------------------------------------------------------------------------------------
    boolean all
    --------------------------------------------------------------------------------------------------------------------
    """

    def square_loop(geom):
        """Construct square using built in geometry."""
        points = [
            geom.add_point([-0.5, -0.5], 0.05),
            geom.add_point([-0.5, 0.5], 0.05),
            geom.add_point([0.5, 0.5], 0.05),
            geom.add_point([0.5, -0.5], 0.05),
        ]
        lines = [
            geom.add_line(points[0], points[1]),
            geom.add_line(points[1], points[2]),
            geom.add_line(points[2], points[3]),
            geom.add_line(points[3], points[0]),
        ]
        return geom.add_curve_loop(lines)

    def circle_loop(geom):
        """construct circle using geo geometry module."""
        points = [
            geom.add_point([+0.0, +0.0], 0.05),
            geom.add_point([+0.0, +0.1], 0.05),
            geom.add_point([-0.1, +0.0], 0.05),
            geom.add_point([+0.0, -0.1], 0.05),
            geom.add_point([+0.1, +0.0], 0.05),
        ]
        quarter_circles = [
            geom.add_circle_arc(points[1], points[0], points[2]),
            geom.add_circle_arc(points[2], points[0], points[3]),
            geom.add_circle_arc(points[3], points[0], points[4]),
            geom.add_circle_arc(points[4], points[0], points[1]),
        ]
        return geom.add_curve_loop(quarter_circles)


    with pygmsh.occ.Geometry() as geom:
        square = square_loop(geom)
        curve_loop = circle_loop(geom)
        surf1 = geom.add_plane_surface(square)
        surf2 = geom.add_plane_surface(curve_loop)
        geom.boolean_fragments(surf1, surf2)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_47_mesh_box():
    """
    --------------------------------------------------------------------------------------------------------------------
    box
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.add_box([0.0, 0.0, 0.0], [1, 2, 3], mesh_size=0.1)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_48_mesh_builtinMix():
    """
    --------------------------------------------------------------------------------------------------------------------
    builtin mix
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = 0.1
        p0 = geom.add_point([-0.5, -0.5, 0], 0.01)
        p1 = geom.add_point([+0.5, -0.5, 0], 0.01)
        p2 = geom.add_point([+0.5, +0.5, 0], 0.01)
        p3 = geom.add_point([-0.5, +0.5, 0], 0.01)
        l0 = geom.add_line(p0, p1)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p3)
        l3 = geom.add_line(p3, p0)
        ll0 = geom.add_curve_loop([l0, l1, l2, l3])
        square_builtin = geom.add_plane_surface(ll0)
        square_occ = geom.add_rectangle([0, 0, 0], 1.0, 1.0)
        geom.boolean_difference(square_occ, square_builtin)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_49_mesh_cone():
    """
    --------------------------------------------------------------------------------------------------------------------
    cone
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.occ.Geometry() as geom:
        geom.add_cone(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            1.0,
            0.3,
            mesh_size=0.1,
            angle=1.25 * pi,
        )
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_50_mesh_cylinder():
    """
    --------------------------------------------------------------------------------------------------------------------
    cylinder
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.occ.Geometry() as geom:
        geom.add_cylinder(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.5, 0.25 * pi, mesh_size=0.1
        )
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_51_mesh_ellipsoid():
    """
    --------------------------------------------------------------------------------------------------------------------
    ellipsoid
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.occ.Geometry() as geom:
        geom.add_ellipsoid([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], mesh_size=0.1)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_52_mesh_extrude():
    """
    --------------------------------------------------------------------------------------------------------------------
    extrude
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = 0.05

        rectangle = geom.add_rectangle([-1.0, -1.0, 0.0], 2.0, 2.0, corner_radius=0.2)
        disk1 = geom.add_disk([-1.2, 0.0, 0.0], 0.5)
        disk2 = geom.add_disk([+1.2, 0.0, 0.0], 0.5, 0.3)

        disk3 = geom.add_disk([0.0, -0.9, 0.0], 0.5)
        disk4 = geom.add_disk([0.0, +0.9, 0.0], 0.5)
        flat = geom.boolean_difference(
            geom.boolean_union([rectangle, disk1, disk2]),
            geom.boolean_union([disk3, disk4]),
        )
        geom.extrude(flat, [0, 0, 0.3])
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_53_mesh_extrude():
    """
    --------------------------------------------------------------------------------------------------------------------
    extrude
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = 1.0

        mesh_size = 1
        h = 25
        w = 10
        length = 100
        # x_fin = -0.5 * length
        cr = 1

        f = 0.5 * w
        y = [-f, -f + cr, +f - cr, +f]
        z = [0.0, h - cr, h]
        f = 0.5 * cr
        x = [-f, f]
        points = [
            geom.add_point((x[0], y[0], z[0]), mesh_size=mesh_size),
            geom.add_point((x[0], y[0], z[1]), mesh_size=mesh_size),
            geom.add_point((x[0], y[1], z[1]), mesh_size=mesh_size),
            geom.add_point((x[0], y[1], z[2]), mesh_size=mesh_size),
            geom.add_point((x[0], y[2], z[2]), mesh_size=mesh_size),
            geom.add_point((x[0], y[2], z[1]), mesh_size=mesh_size),
            geom.add_point((x[0], y[3], z[1]), mesh_size=mesh_size),
            geom.add_point((x[0], y[3], z[0]), mesh_size=mesh_size),
        ]

        lines = [
            geom.add_line(points[0], points[1]),
            geom.add_circle_arc(points[1], points[2], points[3]),
            geom.add_line(points[3], points[4]),
            geom.add_circle_arc(points[4], points[5], points[6]),
            geom.add_line(points[6], points[7]),
            geom.add_line(points[7], points[0]),
        ]

        curve_loop = geom.add_curve_loop(lines)
        surface = geom.add_plane_surface(curve_loop)
        geom.extrude(surface, translation_axis=[length, 0, 0])

        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_54_mesh_regularExtrusion():
    """
    --------------------------------------------------------------------------------------------------------------------
    regular extrusion
    --------------------------------------------------------------------------------------------------------------------
    """
    x = 5
    y = 4
    z = 3
    x_layers = 10
    y_layers = 5
    z_layers = 3
    with pygmsh.occ.Geometry() as geom:
        p = geom.add_point([0, 0, 0], 1)
        _, l, _ = geom.extrude(p, [x, 0, 0], num_layers=x_layers)
        _, s, _ = geom.extrude(l, [0, y, 0], num_layers=y_layers)
        geom.extrude(s, [0, 0, z], num_layers=z_layers)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_55_mesh_torus():
    """
    --------------------------------------------------------------------------------------------------------------------
    torus
    --------------------------------------------------------------------------------------------------------------------
    """
    from math import pi
    with pygmsh.occ.Geometry() as geom:
        geom.add_torus([0.0, 0.0, 0.0], 1.0, 0.3, 1.25 * pi, mesh_size=0.1)
        mesh = geom.generate_mesh()
    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_56_mesh_wedge():
    """
    --------------------------------------------------------------------------------------------------------------------
    wedge
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.add_wedge([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], top_extent=0.4, mesh_size=0.1)
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_57_mesh_translations():
    """
    --------------------------------------------------------------------------------------------------------------------
    translations
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.05
        geom.characteristic_length_max = 0.05
        disk = geom.add_disk([0, 0, 0], 1)
        disk2 = geom.add_disk([1.5, 0, 0], 1)
        geom.translate(disk, [1.5, 0, 0])
        geom.boolean_union([disk2, disk])
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes

def example_58_mesh_translations():
    """
    --------------------------------------------------------------------------------------------------------------------
    translations
    --------------------------------------------------------------------------------------------------------------------
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = 0.2
        geom.characteristic_length_max = 0.2
        ball = geom.add_ball([0, 0, 0], 1)
        ball2 = geom.add_ball([1.5, 0, 0], 1)
        geom.translate(ball, [1.5, 0, 0])
        geom.boolean_union([ball2, ball])
        mesh = geom.generate_mesh()

    meshProcessInstance = EleMeshPlotAndSelect()  ###---initialize class
    meshProcessInstance.addMesh(mesh, nodeStartNumber=0, elementStartNumber=0)  ###---add mesh to plot
    meshProcessInstance.meshProcess(showPointTag=False)  ###---Process the mesh
    meshProcessInstance.meshPlot()  ###---visualize the meshes
########################################################################################################################
########################################################################################################################
if __name__ == "__main__":
    #############################################################################
    ####---example_1_mesh
    example_1_mesh_polygon()
    #############################################################################
    ####---example_2_mesh
    # example_2_mesh_circle()
    #############################################################################
    ####---example_3_mesh
    # example_3_mesh_BSpline()
    #############################################################################
    ####---example_4_mesh
    # example_4_mesh_extrude()
    #############################################################################
    ####---example_5_mesh
    # example_5_mesh_revolve()
    #############################################################################
    ####---example_6_mesh
    # example_6_mesh_twist()
    #############################################################################
    ####---example_7_mesh
    # example_7_mesh_BooleanIntersection()
    #############################################################################
    ####---example_8_mesh
    # example_8_mesh_BooleanDifference()
    #############################################################################
    ####---example_9_mesh
    # example_9_mesh_BooleanUnion()
    #############################################################################
    ####---example_10_mesh
    # example_10_mesh_refinement()
    #############################################################################
    ####---example_11_mesh
    # example_11_mesh_airfoil()
    #############################################################################
    ####---example_12_mesh
    # example_12_mesh_bsplines()
    #############################################################################
    ####---example_13_mesh
    # example_13_mesh_circle()
    #############################################################################
    ####---example_14_mesh
    # example_14_mesh_circleTransform()
    #############################################################################
    ####---example_15_mesh
    # example_15_mesh_cube()
    #############################################################################
    ####---example_16_mesh
    # example_16_mesh_ellipsoid()
    #############################################################################
    ####---example_17_mesh
    # example_17_mesh_embed()
    #############################################################################
    ####---example_18_mesh
    # example_18_mesh_hex()
    #############################################################################
    ####---example_19_mesh
    # example_19_mesh_holeInSquareHighOrder()
    #############################################################################
    ####---example_20_mesh
    # example_20_mesh_layers()
    #############################################################################
    ####---example_21_mesh
    # example_21_mesh_pacman()
    #############################################################################
    ####---example_22_mesh
    # example_22_mesh_physical()
    #############################################################################
    ####---example_23_mesh
    # example_23_mesh_pipes()
    #############################################################################
    ####---example_24_mesh
    # example_24_mesh_quads()
    #############################################################################
    ####---example_25_mesh
    # example_25_mesh_recombine()
    #############################################################################
    ####---example_26_mesh
    # example_26_mesh_rectangle()
    #############################################################################
    ####---example_27_mesh
    # example_27_mesh_regularExtrusion()
    #############################################################################
    ####---example_28_mesh
    # example_28_mesh_rotatedLayers()
    #############################################################################
    ####---example_29_mesh
    # example_29_mesh_rotation()
    #############################################################################
    ####---example_30_mesh
    # example_30_mesh_screw()
    #############################################################################
    ####---example_31_mesh
    # example_31_mesh_subdomains()
    #############################################################################
    ####---example_32_mesh
    # example_32_mesh_swissCheese()
    #############################################################################
    ####---example_33_mesh
    # example_33_mesh_symmetrize()
    #############################################################################
    ####---example_34_mesh
    # example_34_mesh_tori()
    #############################################################################
    ####---example_35_mesh
    # example_35_mesh_torus()
    #############################################################################
    ####---example_36_mesh
    # example_36_mesh_torusCrowd()
    #############################################################################
    ####---example_37_mesh
    # example_37_mesh_transfinite()
    #############################################################################
    ####---example_38_mesh
    # example_38_mesh_unorderedUnoriented()
    #############################################################################
    ####---example_39_mesh
    # example_39_mesh_BallWithStick()
    #############################################################################
    ####---example_40_mesh
    # example_40_mesh_logo()
    #############################################################################
    ####---example_41_mesh
    # example_41_mesh_ball()
    #############################################################################
    ####---example_42_mesh
    # example_42_mesh_booleanUnion()
    #############################################################################
    ####---example_43_mesh
    # example_43_mesh_booleanIntersection()
    #############################################################################
    ####---example_44_mesh
    # example_44_mesh_booleanDifference()
    #############################################################################
    ####---example_45_mesh
    # example_45_mesh_booleanAll()
    #############################################################################
    ####---example_46_mesh
    # example_46_mesh_booleansSquareCircleSlice()
    #############################################################################
    ####---example_47_mesh
    # example_47_mesh_box()
    #############################################################################
    ####---example_48_mesh
    # example_48_mesh_builtinMix()
    #############################################################################
    ####---example_49_mesh
    # example_49_mesh_cone()
    #############################################################################
    ####---example_50_mesh
    # example_50_mesh_cylinder()
    #############################################################################
    ####---example_51_mesh
    # example_51_mesh_ellipsoid()
    #############################################################################
    ####---example_52_mesh
    # example_52_mesh_extrude()
    #############################################################################
    ####---example_53_mesh
    # example_53_mesh_extrude()
    #############################################################################
    ####---example_54_mesh
    # example_54_mesh_regularExtrusion()
    #############################################################################
    ####---example_55_mesh
    # example_55_mesh_torus()
    #############################################################################
    ####---example_56_mesh
    # example_56_mesh_wedge()
    #############################################################################
    ####---example_57_mesh
    # example_57_mesh_translations()
    #############################################################################
    ####---example_58_mesh
    # example_58_mesh_translations()
