# Tip-loaded cantilever beam with shear deformation, plane section.
#
# May 30 2025
#
import os
import xara
from xara.benchmarks import Prism
import xara.units.iks as units
from xsection.library import from_aisc
import numpy as np


def analyze(model, P):
    nn = len(model.getNodeTags())

    model.pattern("Plain", 1, "Linear", load={nn: [0,0,1, 0,0,0]})
    model.integrator("LoadControl", P)
    model.analysis("Static")
    model.test("Residual", 1e-10, 1)

    assert model.analyze(1) == 0



if __name__ == "__main__":
#   shape = Rectangle(d=14, b=10, mesh_scale=1/100) #


    E  = 29e3*units.ksi
    G  = 11.2e3*units.ksi
    material = xara.Material(E=E, G=11.2e3*units.ksi)
    shape = from_aisc("W14x48", units=units, mesh_scale=1/100, material=material)

    section_type = os.environ.get("Section", "MultiaxialFiber")
    section  = xara.FrameSection(section_type, material, shape, mixed=False)

    A  = shape.cnn()[0,0]
    GA = A*G
    EI = E*shape.cmm()[1,1]

    L = shape.d*4
    for shear in [0,1]:
        prism = Prism(shape=shape,
                      length=L,
                      boundary=((1,1,1,  1,1,1),
                                (0,0,0,  0,0,0)),
                      material=material,
                      element=os.environ.get("Element", "ForceFrame"),
                      section=section,
                      shear=shear,
                      divisions=1,
                      order=1
                )
        model = prism.create_model()

        P = 10
        analyze(model, P)

        end = len(model.getNodeTags())
        uz = model.nodeDisp(end, 3)
        u_euler = P*L**3/(3*EI)
        u_shear = P*L/GA * shear

        print(f"Uz = {uz:.6f}, Uz theory = {u_euler+u_shear:.6f} ({u_euler:.6f} + {u_shear:.6f})")

        model.eval(f"verify value [nodeDisp {end} 3] {u_euler+u_shear:.12f} 1e-6")
