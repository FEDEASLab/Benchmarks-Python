# Verification of joint offsets in nonlinear frame transformations
#
# Cantilever beam with offset axial load
#
#   e
#  | |
#  o-.
#    |
#    |
#    |
#    |
#   _#_
#
import xara
from xara.helpers import find_node
import pytest 
import numpy as np


L = 120
e = L/40 #L/24
P = 2000
E = 29e3
I = 6*12**3/12
A = 1e3
Element = "ElasticFrame"
ShearFlag = 0


Pcr = E*I*np.pi**2/(L*2)**2


def create_model_offset(transform, element, offset=True, nen=2, ne=6, ndm=3):
    """
    offset:
       True:  offset at joint
    """

    nn = ne*(nen-1)+1

    section = xara.FrameSection("Elastic",
      E  = E,
      G  = 11e3,
      J  = 2*12**4/3,
      Iy = I,
      Iz = 6*12**3/12,
      A  = A
    )

    model = xara.Model(ndm=ndm, ndf=6 if ndm == 3 else 3)
    model.section(section)

    for i, y in enumerate(np.linspace(0, L, nn)):
        if offset and i == ne:
            x = -e
        else:
            x = 0
        model.node(i, (x, y, 0))

    model.geomTransf(transform, 1, (0, 0, 1))
    model.geomTransf(transform, 2, (0, 0, 1), 
                    jntOffset=(0, 0, 0,    
                               e if offset else 0, 0,0))


    for i in range(ne):
        transform = 2 if i == ne-1 else 1

        start = i * (nen - 1)
        nodes = list(range(start, start + nen))

        model.element(element, i, nodes, 
                      section=section, 
                      transform=transform,
                      shear=ShearFlag if "Cosserat" not in element else 1)


    tip = nn-1
    base = 0
    model.fix(0,  (1, 1, 1, 1, 1, 1))
    model.fix(tip, (0, 0, 1, 1, 1, 0))
    for i in range(1, nn-1):
        model.fix(i, (0, 0, 1, 1, 1, 0))


    #
    # Load
    #
    model.pattern(
        xara.StaticPattern([
            xara.NodalLoad(model, {tip:  (0,-P, 0,  0, 0, P*e if not offset else 0),
                                   base: (0, 0, 0,  0, 0, 0)})
        ])
    )

    return model



def analyze(transform, element, offset=True, model=None, n=40, dlam=1/10, ne=6, tol=1e-12, verbose=False):

    # Create model
    if model is not None:
        pass
        case = ""
    elif isinstance(offset, str):
        case = "Load"
        model = create_model_follower(transform=transform, element=element, ne=ne)
    else:
        case = "Joint" if offset else "Moment"
        model = create_model_offset(transform=transform, element=element, offset=offset, ne=ne)

    analysis = xara.StaticAnalysis(model, 
                                   test=("NormDispIncr", tol, 100, 1 if verbose else 9),
                                #    test=("NormDispIncr", 1e-9, 50,  9),
                                   system="BandGeneral",
                                   integrator=("LoadControl", dlam))

    tip = find_node(model, y=L)
    u, N, ic = [], [], []
    if not verbose:
        print(f"{transform:<12} \t {element:<16} {case:<8} \t ", end="")

    for i in range(n):
        assert analysis.analyze(1) == 0
        ic.append(model.getIterationCount())
        u.append(model.nodeDisp(tip, 6))
        N.append(P*model.getLoadFactor(1)/Pcr)
    else:
        if not verbose:
            print(model.nodeDisp(tip, 1), model.nodeDisp(tip, 2), u[-1])

    return model



def test_linear():
    model = analyze("Linear", "PrismFrame", 
            offset=True, 
            tol=1e-14, 
            n=1, # 1 step
            ne=1, 
            verbose=False, 
            dlam=1)

    tip = find_node(model, y=L)
    assert model.nodeDisp(tip, 1) == pytest.approx(-1.724137931034483, abs=1e-10)
    assert model.nodeDisp(tip, 2) == pytest.approx(-0.094482758620690, abs=1e-10)
    assert model.nodeDisp(tip, 6) == pytest.approx( 0.028735632183908, abs=1e-10)



def test_pdelta():
    model = analyze("PDelta02", "PrismFrame", 
            offset=True, 
            tol=1e-14, 
            n=1, # 1 step
            ne=1, 
            verbose=False, 
            dlam=1)

    tip = find_node(model, y=L)
    assert model.nodeDisp(tip, 1) == pytest.approx(-2.795031055900625, abs=1e-10)
    assert model.nodeDisp(tip, 2) == pytest.approx(-0.134641250803170, abs=1e-10)
    assert model.nodeDisp(tip, 6) == pytest.approx( 0.042121796244735, abs=1e-10)


# def test_C3():
#     model = analyze("PDelta02", "PrismFrame", 
#             offset=True, 
#             tol=1e-14, 
#             n=1, # 1 step
#             ne=1, 
#             verbose=False, 
#             dlam=1)

#     tip = find_node(model, y=L)
#     assert model.nodeDisp(tip, 1) == pytest.approx(-2.795031055900625, abs=1e-10)
#     assert model.nodeDisp(tip, 2) == pytest.approx(-0.134641250803170, abs=1e-10)
#     assert model.nodeDisp(tip, 6) == pytest.approx( 0.042121796244735, abs=1e-10)

