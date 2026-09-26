# Test of the warping DOF
#
# Linear 7-DOF analysis of a cantilever subjected to a torque
#
import os
import sys

import veux
from veux.motion import Motion
from shps.shapes import WideFlange, HollowRectangle
import opensees.openseespy as ops

# External libraries
import numpy as np

def create_cantilever(shape, case="a", element="ExactFrame", section="Elastic"):

    E = 29e3 # ksi
    v = 0.30 #0.5*E/G - 1
    G = 0.5*E/(1+v) # 11,200 ksi

    L  = 240
    ne =  40
    nen = 2
    nn = ne*(nen-1)+1

    model = ops.Model(ndm=3, ndf=7)

    mat = 1
    sec = 1
    model.material('ElasticIsotropic', mat, E, v)


    if section == "Elastic":
        cmm = shape.cmm()
        cnn = shape.cnn()
        cnv = shape.cnv()
        cnm = shape.cnm()
        cmw = shape.cmw()
        A = cnn[0,0]
        model.section("ElasticFrame", sec,
                        E=E,
                        G=G,
                        A=A,
                        Ay=1*A,
                        Az=1*A,
                        Qy=cnm[0,1],
                        Qz=cnm[2,0],
                        Iy=cmm[1,1],
                        Iz=cmm[2,2],
                        J =shape.torsion.torsion_constant(),
                        Ry= cnv[1,0],
                        Rz=-cnv[2,0],
                        Sy= cmw[1,0],
                        Sz=-cmw[2,0]
        )
    else:
        model.section("ShearFiber", 1, GJ=0)
        for fiber in shape.create_fibers():
            model.fiber(**fiber, material=mat, section=1)


    model.geomTransf("Linear", 1, (0,0,1))

    for i,x in enumerate(np.linspace(0, L, nn)):
        model.node(i, (x,0,0))

    for i in range(ne):
        start = i * (nen - 1)
        nodes = list(range(start, start + nen))
        model.element(element, i+1, nodes, section=sec, transform=1, shear=1)

    wi = int(case in "cb")
    wj = int(case in "c")

    model.fix(0,     (1,1,1,  1,1,1, wi))
    model.fix(nn-1,  (0,0,0,  0,0,0, wj))
    return model, shape


if __name__ == "__main__":
    case = "b" # "c" "a"

    # W21x93
    shape = WideFlange(
                tf = 0.93,
                tw = 0.58,
                d  = 21.62,
                b  = 8.42,
                mesh_scale=1/100
    )

    print(shape.summary())

#   shape = HollowRectangle(
#               tf = 0.4,
#               tw = 0.4,
#               d  = 21.62,
#               b  = 8.42,
#               mesh_scale=1/200
#           )


    speed = 1/10

    model, shape = create_cantilever(shape, case,
                                     section = os.environ.get("Section", "ShearFiber"),
                                     element = os.environ.get("Element", "ExactFrame"))
    end = len(model.getNodeTags()) - 1

    # Render
    artist = veux.create_artist(model,
                                model_config=dict(
                                    extrude_outline=shape,
                                    section_warping=shape.model.create_handle(shape._analysis.solution())
                                )
    )
    motion = Motion(artist)

    # Apply torsional moment
    nsteps =  15
    Mmax   = 1.2e3
    model.pattern("Plain", 1, "Linear")
#       model.load(0  , (0,0,0,  0,0,0, -1e-3), pattern=1)
    model.load(end, (0,0,0,  1,0,0,  0), pattern=1)

    model.system('Umfpack')
    model.integrator("LoadControl", Mmax/nsteps)
    model.test("NormDispIncr",1e-12,5,1)
    model.algorithm("Newton")
    model.analysis("Static")

    u = []
    P = []
    while model.getTime() < Mmax:
        u.append(model.nodeDisp(end, 4))
        P.append(model.getTime())
        if model.analyze(1) != 0:
            print(f"Failed at time = {model.getTime()}")
            break


        if motion is not None:

            motion.advance(time=model.getTime()*speed)

            motion.draw_sections(rotation=model.nodeRotation,
                                 position=model.nodeDisp,
            )
            motion.set_field(lambda tag: model.nodeDisp(tag, 7))



    if motion is not None:
        motion.add_to(artist.canvas)
        if len(sys.argv) > 1:
            artist.save(sys.argv[1])
        else:
            veux.serve(artist)


