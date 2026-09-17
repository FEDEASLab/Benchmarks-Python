#
# Cantilever
#
import os
import time
import xara

# External libraries
import numpy as np
import matplotlib.pyplot as plt


def create_cantilever(shape,
                      material,
                      L,
                      element,
                      transform="Corotational02",
                      shear=1,
                      ne=10):

    nen = 2
    nn = ne*(nen-1)+1

    model = xara.Model(ndm=3, ndf=6)

    sec = 1


    model.section("ElasticFrame", sec, **material, **shape)

    model.geomTransf(transform, 1, (0,0,1))

    for i,x in enumerate(np.linspace(0, L, nn)):
        model.node(i, (x,0,0))


    for i in range(ne):
        start = i * (nen - 1)
        nodes = tuple(range(start, start + nen))
        model.element(element, i+1, nodes, section=sec, transform=1, shear=shear)


    model.fix(0,     (1,1,1,  1,1,1))
    model.fix(nn-1,  (0,0,0,  0,0,0))
    return model


def analyze(model, Mmax, Fmax, nsteps=1):

    end = len(model.getNodeTags()) - 1
    #
    # Apply loading
    #
    model.pattern("Plain", 1, "Linear")
    model.load(end, (0,0,Fmax,  0,0,Mmax), pattern=1)

    model.system('BandGen')
    model.integrator("LoadControl", 1/nsteps)#, 6, 1/nsteps/10, 1/nsteps)
    model.test("Energy", 1e-16, 10, 0) #50 #20
    model.algorithm("Newton")
    model.analysis("Static")

    ux = []
    uy = []
    P  = []
    ni = []
    L  = model.nodeCoord(end, 1)

    ux.append(1+model.nodeDisp(end, 1)/L)
    uy.append(  model.nodeDisp(end, 2)/L)
    P.append(model.getTime())

    while model.getTime() < 1.0:
        if model.analyze(1) != 0:
            print(f"Failed at time = {model.getTime()}")
            break
        ni.append(model.numIter()-1)
        ux.append(1+model.nodeDisp(end, 1)/L)
        uy.append(  model.nodeDisp(end, 2)/L)
        P.append(model.getTime())


    print(sum(ni)/len(ni), np.array(model.nodeDisp(end))[:3]) # + model.nodeCoord(end))
    return ux, uy, P, ni



if __name__ == "__main__":
    length = 10.0

    material = dict(
        E = 1e4,
        G = 1e4
    )

    shape = dict(
        A   = 1,
        Ay  = 1,
        Az  = 1,
        Iy  = 1e-2,
        Iz  = 1e-2,
        J   = 1e-2
    )

    M0 = 20*np.pi # 2 pi EI /L

    analyze(create_cantilever(
                    shape,
                    material,
                    length,
                    ne=10,
                    shear=1,
                    element = "ExactFrame",
                    transform="Linear"),
              Mmax=M0/8, Fmax=1/16, nsteps=1)


    Transforms = [
        "Corotational",
        "Corotational02",
        "Corotational03",
        "Corotational04",
        "Corotational05"
    ]


    for transform in Transforms: #
        print(f"Transform: {transform}")
        model = create_cantilever(shape,
                                material,
                                length,
                                ne=10,
                                element = "ForceFrame",
                                transform=transform)

        analyze(model, Mmax=M0/8, Fmax=1/16, nsteps=1)


