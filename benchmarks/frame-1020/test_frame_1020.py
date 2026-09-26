#
# Cantilever beam subjected to follower end load.
#
# Reference
# ---------
#
# Simo, J. C., and L. Vu-Quoc. 
#   "A Three-Dimensional Finite-Strain Rod Model. Part II: Computational Aspects."
#   Computer Methods in Applied Mechanics and Engineering 58, no. 1 (1986): 79–116. 
#   https://doi.org/10.1016/0045-7825(86)90079-4.
#
import numpy as np
import xara


def create_cantilever(ne, element, transform, rotations=None):

    model = xara.Model(ndm=3, ndf=6)

    L  = 100

    nen = 2
    nmn = ne*(nen-1)+1
    sec = 1
    A = 1.61538e8
    I = 3.5e7
    model.section("ElasticFrame", sec,
                    E=1,
                    G=1,
                    A=A,
                    Ay=A,
                    Az=A,
                    Iy=I*100,
                    Iz=I,
                    J =I
    )

    model.geomTransf(transform, 1, (0,0,1))

    for i,x in enumerate(np.linspace(0, L, nmn)):
        model.node(i, (x,0,0))

    for i in range(ne):
        start = i * (nen - 1)
        nodes = list(range(start, start + nen))
        model.element(element, i+1, nodes,
                      section=sec, transform=1, shear=1)

    model.fix(0,  (1,1,1,  1,1,1))
    for i in range(nmn):
        model.nodeRotation(i)

    return model


def analyze(element, transform):
    ne = 10

    model = create_cantilever(ne, element=element, transform=transform)

    #
    # Apply vertical load
    #
    Pmax   = 150e3 # N
    model.pattern("Plain", 1, "Linear")


    model.eleLoad("Frame", # TODO: use xara.FrameLoad
                  "Point",
                  force = [0, 1, 0],
                  basis = "director",
                  offset=[1.0,0,0],
                  pattern=1,
                  elements=[ne]
    )

    model.system('Umfpack')
    model.integrator("LoadControl", Pmax/500)
    model.test("Energy", 1e-16, 15, 0)
#   model.test('NormUnbalance',1e-6,100,1)
    model.algorithm("Newton")
    model.analysis("Static")

    u = []
    v = []
    w = []
    P = []
    while model.state.time < Pmax:
        assert model.analyze(1) == 0
        u.append(-model.nodeDisp(ne, 1))
        v.append( model.nodeDisp(ne, 2))
        w.append( model.nodeDisp(ne, 3))
        P.append( model.getTime())

    return u, v, w, P


def test_exactframe():
    analyze(element = "ExactFrame", transform = "Linear")


def test_corotational():
    analyze(element = "ForceFrame", transform = "Corotational02")


def test_cosserat():

    analyze(element = "CosseratFrame", transform = "Identity")

    analyze(element = "CosseratFrame", transform = "Spherical")

    analyze(element = "CosseratFrame", transform = "Corotational02")


if __name__ == "__main__":
    test_exactframe()
    test_corotational()
    test_cosserat()
