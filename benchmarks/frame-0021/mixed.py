# Test of the warping DOF
#
# Linear 7-DOF analysis of a cantilever subjected to a torque
#
import os

import veux
from veux.motion import Motion
from xsection.analysis import SaintVenantSectionAnalysis
from xsection._benchmarks import load_shape
from xsection.library import WideFlange, HollowRectangle, Channel, Rectangle, Circle

import xara
from xara import Section, Material

# External libraries
import numpy as np
import matplotlib.pyplot as plt
try:
    plt.style.use("veux-web")
except:
    pass

from plots import (
    PlotTwist2, 
    PlotConvergence,
    plot_resultants, 
    plot_kinematics
)
from model import create_cantilever, analyze
from shaft import Shaft

if __name__ == "__main__":
    th = 0.1 #0.05
    depth = 20
    width = 8 # depth*0.4
    Failed = []

    material = Material(
        G = 11.2e3,
        E = 29e3
    )

    print(f"G = {material['G']}, E = {material['E']}")
    # Mmax   = 1.2e3

    element = os.environ.get("Element", "ExactFrame")
    section = os.environ.get("Section", "MixedFiber")

    WarpTypes = os.environ.get("WarpType", "NT,NR").split(",")
    Boundary  = os.environ.get("Boundary", "b,c").split(",")


    Shapes = os.environ.get("Shape", "w,h,r").split(",")
    mesh_scale = 1/int(os.environ.get("Mesh", "200"))

    for shape_name in Shapes:
        print(f"Shape {shape_name.upper()}")
        if shape_name == "c":
            shape = Channel(
                        tf=th*depth,
                        tw=th*depth,
                        d=depth,
                        b=depth*0.4,
                        material=material,
                        mesh_type="T6",
                        mesher="gmsh",
                        mesh_scale=1/15 #800
            )
            # shape = shape.translate(-shape._analysis.shear_center())

        elif shape_name == "r":
            shape = Rectangle(d=depth, 
                              b=0.4*depth, 
                              material=material,
                              mesh_scale=1/20, 
                              mesh_type="T6", 
                              mesher="gmsh")

        elif shape_name == "o":
            shape = Circle(depth/2, 
                           divisions=8, 
                           mesh_scale=1/200, 
                           material=material)

        elif shape_name == "w":
            # W21x93
            shape = WideFlange(
                        tf = 0.93,
                        tw = 0.58,
                        d  = 21.62,
                        b  = 8.42,
                        material=material,
                        mesher="gmsh",
                        mesh_type="T6",
                        mesh_scale=1/4
                    )
            # shape = WideFlange(
            #             tf = th*depth,
            #             tw = th*depth,
            #             d  = depth,
            #             b  = width,
            #             mesh_scale=1/10,
            #             mesher="gmsh"
            #         )
        elif shape_name == "h":
            shape = HollowRectangle(
                        tf = th*depth*2,
                        tw = th*depth,
                        d  = depth,
                        b  = width,
                        material=material,
                        mesher="gmsh",
                        mesh_type="T6",
                        mesh_scale=1/8
                    )
        elif shape_name == "h02":
            shape = HollowRectangle(
                        tf = th*depth,
                        tw = th*depth,
                        d  = depth,
                        b  = width,
                        material=material,
                        mesher="gmsh",
                        mesh_scale=1/8#0
                    )
        else:
            shape = load_shape(shape_name, mesh_scale=1/10, mesher="gmsh", material=material)
            shape = shape.translate(-shape._analysis.shear_center())
            # shape = shape.translate(-shape.centroid)


        # print(shape.summary())

        # veux.serve(veux.render(shape.model))

        sv = SaintVenantSectionAnalysis(shape)

        GJ = sv.twist_rigidity()
        Mmax = GJ/(depth/2)*np.pi*2*1e-5


        for boun in Boundary:

            key = f"{element[:5]}-{section[:5]}-{shape_name}-{boun}"

            for slenderness in [2,1,0.5]:#, 1, 0.5]:

                p1 = PlotTwist2(Mmax, GJ, boun,
                                title=f"Shape {shape_name.upper()}, Case {boun.upper()}, slenderness {slenderness}",
                                skip="Batch" in os.environ)

                for warp_type in WarpTypes:

                    print(f"Shape {shape_name.upper()} Case {boun}, slenderness {slenderness}, warp {warp_type}")
                    p1.reset(label=f"warp = {warp_type}")


                    model = create_cantilever(
                        slenderness,
                        shape,
                        material,
                        boun,
                        ne=8, #16,
                        warp_type=warp_type,
                        nen=3 if element == "ExactFrame" else 2,
                        section = Section(type=section, 
                                          shape=shape,
                                          material=material,
                                          mixed_type=warp_type),
                        element = element)

                    analyze(model, Mmax, element,tol=1e-12)

                    model.reactions()
                    p1.update(model)

                    end = model.getNodeTags()[-1]
                    L = model.nodeCoord(end, 1)
                    shaft = Shaft(boun, L, shape, warp_type=warp_type, sv=sv)
                    # p1.draw(L, shaft, Mmax)
                    try:
                        model.eval(f"verify error [nodeDisp {end} 4] {shaft.twist(L, Mmax):.12f} 1e-3 \"{boun}\"")
                    except:
                        Failed.append(f"{key} slenderness {slenderness} warp {warp_type}")
                        pass
                        # raise

            
                p1.finalize()


    if "Batch" not in os.environ:
        plt.show()

    print("Failed cases:")
    for i,f in enumerate(Failed):
        print(" ", i+1, " ", f)
