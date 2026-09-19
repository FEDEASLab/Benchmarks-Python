"""
Test distributed element loads.
"""

import xara
import pytest
from xara.load import FrameLoad
import math
isclose = lambda a,b,rel_tol=1e-9: math.isclose(a,b,rel_tol=rel_tol)

L = 48
E = 29000
G = 12000
h = 6
b = h
A = h*b
I = h**3*b/12
J = 2*I

w = 0.1
m = (h/2)*w

F = 1
offset = b/2
T = F*offset
Elements = [
    "ExactFrame", "ForceFrame"
]

def _run_axial(element, load_type="center"):
    """
    A cantilever is subjected to a uniformly distributed axial
    force at an offset from the neutral axis.

    Two equivalent representations are tested:
        1) A uniformly distributed couple and axial force
        2) A uniformly distributed axial force with an offset
    """
    model = xara.Model(ndm=3, ndf=6)

    model.node(1, (0,0,0))
    model.node(2, (L,0,0))

    model.fix(1, (1,1,1, 1,1,1))

    model.geomTransf('Linear', 1, (0,0,1))

    model.section('ElasticFrame',1,
                  E=E,A=A,Iy=I,Iz=I,
                  G=G,
                  J=J,
                  Ay=A*500,
                  Az=A*500)

    model.element(element, 1, [1,2], section=1, transform=1, 
                  shear=1 if element == "ExactFrame" else 0)

    if load_type == "offset":
        load = FrameLoad(model,
                        shape='Uniform', 
                        basis='local',
                        elements=[1],
                        couple=[0,0,-m],
                        offset=[0,h/2,0])
    elif load_type == "center":
        load = FrameLoad(model,
                        shape='Uniform', 
                        basis='local',
                        elements=[1],
                        force=[w, 0, 0],
                        couple=[0,0,-m],
                        offset=[0,0,0])

    model.pattern(xara.StaticPattern(load))

    assert model.getEleLoadClassTags()  == 141414
    assert model.getEleLoadClassTags(1) == 141414
    assert model.getEleLoadClassTags(pattern=1) == 141414

    analysis = xara.StaticAnalysis(model)
    analysis.analyze()

    assert model.getLoadFactor(1) == 1.0

    model.reactions()
    return model


def test_axial():
    for load_type in ["offset", "center"]:
        for element in ["ForceFrame", "ExactFrame"]:
            if load_type == "center" and element == "ExactFrame":
                continue
            tol = 3e-1 if element == "ExactFrame" else 1e-9

            print(element, load_type)
            model = _run_axial(element, load_type=load_type)

            # Root moment
            assert model.nodeReaction(1,6) == pytest.approx(m*L, rel=1e-6)
            # Tip rotation
            assert model.nodeDisp(2,6) == pytest.approx(-m*L**2/(2*E*I), rel=1e-6)
            # Tip deflection
            assert model.nodeDisp(2,2) == pytest.approx(-m*L**3/(3*E*I), rel=tol)
            # assert isclose(-m*L**3/(3*E*I), u22),   (-m*L**3/(3*E*I), u22)


def _run_twist(element, load_type="center"):
    model = xara.Model(ndm=3, ndf=6)

    model.node(1, (0,0,0))
    model.node(2, (L,0,0))
    model.fix(1, (1,1,1, 1,1,1))

    # beam axis = local/global x
    model.geomTransf("Linear", 1, (0,0,1))

    model.section(
        "ElasticFrame", 1,
        E=E, A=A, Iy=I, Iz=I,
        G=G, J=J,
        Ay=A*500, 
        Az=A*500,
    )

    model.element(element, 1, [1,2], section=1, transform=1, shear=1)

    if load_type == "offset":
        # +z force with +y eccentricity -> +x torque
        load = FrameLoad(
            model,
            shape="Point",
            basis="local",
            elements=[1],
            force=[0,0,F],
            offset=[1,b/2, 0],
        )
    elif load_type == "center":
        load = FrameLoad(
            model,
            shape="Point",
            basis="local",
            elements=[1],
            force=[0,0,F],
            couple=[T,0,0],
            offset=[1,0,0],
        )


    model.pattern(xara.StaticPattern(load))
    analysis = xara.StaticAnalysis(model)
    analysis.analyze()
    model.reactions()
    return model


def _run_twist_column(element, load_type="center", basis="local"):
    # Member axis is global z.
    model = xara.Model(ndm=3, ndf=6)

    model.node(1, (0, 0, 0))
    model.node(2, (0, 0, L))

    model.fix(1, (1, 1, 1, 1, 1, 1))

    model.geomTransf("Linear", 1, (1, 0, 0))

    model.section(
        "ElasticFrame", 1,
        E=E, A=A, Iy=I, Iz=I,
        G=G, J=J,
        Ay=A*500,
        Az=A*500,
    )

    model.element(element, 1, [1, 2], section=1, transform=1, shear=1)

    if load_type == "offset":
        if basis == "local":
            # In local coordinates, use r = (0, b/2, 0), F = (0, 0, F)
            # so that r x F = (+T, 0, 0), i.e. positive torsion.
            load = FrameLoad(
                model,
                shape="Point",
                basis=basis,
                elements=[1],
                force=[0, 0, F],
                offset=[1, b/2, 0],
            )

    elif load_type == "center":
        if basis == "local":
            load = FrameLoad(
                model,
                shape="Point",
                basis=basis,
                elements=[1],
                force=[0,0,F],
                couple=[T,0,0],
                offset=[1,0,0],
            )
        elif basis == "global":
                load = FrameLoad(
                model,
                shape="Point",
                basis="global",
                elements=[1],
                force=[0, 0, 0],
                couple=[0, 0, T],   # global z = member axis
                offset=[1, 0, 0],
            )
    else:
        raise ValueError(f"Unknown load_type {load_type!r}")

    model.pattern(xara.StaticPattern(load))
    analysis = xara.StaticAnalysis(model)
    analysis.analyze()
    model.reactions()
    return model


def test_twist_couple():
    model = _run_twist("ForceFrame", load_type="center")

    twist = model.nodeDisp(2, 4)
    root_torque = model.nodeReaction(1, 4)

    assert isclose(-T, root_torque), (root_torque, -T)
    assert isclose(T*L/(G*J), twist), (twist, T*L/(G*J))


def test_twist_offset():
    model = _run_twist("ForceFrame", load_type="offset")

    assert model.nodeDisp(2, 4) == pytest.approx(T*L/(G*J), rel=1e-9)
    assert model.nodeReaction(1, 4) == pytest.approx(-T, rel=1e-9)


def test_twist_column_local_couple():
    for element in Elements:
        model = _run_twist_column(element, load_type="center", basis="local")

        print("Displacements ", model.nodeDisp(2))
        print("Reactions ", model.nodeReaction(1))

        assert model.nodeReaction(1,6) == pytest.approx(-T, rel=1e-4)
        assert model.nodeDisp(2,6) == pytest.approx(T*L/(G*J), rel=1e-4)


# def test_twist_column_global_couple():
#     for element in Elements:
#         model = _run_twist_column(element, load_type="center", basis="global")

#         print("Displacements ", model.nodeDisp(2))
#         print("Reactions ", model.nodeReaction(1))
#         twist = model.nodeDisp(2, 6)
#         root_torque = model.nodeReaction(1, 6)

#         assert isclose(-T, root_torque), (root_torque, -T)
#         assert isclose(T*L/(G*J), twist), (twist, T*L/(G*J))


if __name__ == "__main__":
    test_axial()
    test_twist_couple()
    test_twist_offset()
    test_twist_column_local_couple()
