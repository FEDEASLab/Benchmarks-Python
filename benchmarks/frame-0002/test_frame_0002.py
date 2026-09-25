#
# Vertical cantilever with terminal loads implemented with FrameLoad
#
import xara
from xara.load import FrameLoad, NodalLoad
import pytest
import numpy as np

L = 140
ne = 2


E = 29e3
G = 11e3
d = 12
b = 6
Iy = b*d**3/12
Iz = d*b**3/12
J  = Iy + Iz
A  = 1e3
Ay  = 5/6*b*d
Az  = Ay

#     V.    V.      N
F = [1000, 2000,   60]

M = [10.,30.,15.0]
e = (b/3, d/3)

# Global basis vectors
e_major, e_minor, e_axial = np.eye(3)
# Local basis vectors
i_major = np.array([0., 1., 0.])
i_minor = np.array([0., 0., 1.])
i_axial = np.array([1., 0., 0.])


def create_prism(shape, shear=0, load_type="global"):

    model = xara.Model(ndm=3, ndf=6)

    model.node(1, (0, 0, 0))
    model.node(2, (0, 0, L))

    model.fix(1, (1, 1, 1, 1, 1, 1))

    model.geomTransf("Linear", 1, (0, 1, 0))

    section = xara.FrameSection("Elastic", **shape)
    model.section(section)

    model.element("ForceFrame", 1, (1, 2), 
                  section=section, 
                  shear=shear,
                  transform=1)

    if load_type == "node":
        f = F
        r = e[0]*e_major + e[1]*e_minor
        m = M + np.cross(r, F)
        model.pattern(xara.StaticPattern(
            NodalLoad(model, {
                2: (*f, *m)
            })                   
        ))

    else:
        if load_type == "global":
            f = F
            m = M
        elif load_type == "local":
            f = ( np.dot(F, e_axial)*i_axial
                + np.dot(F, e_minor)*i_minor
                + np.dot(F, e_major)*i_major).tolist()
            m = ( np.dot(M, e_axial)*i_axial
                + np.dot(M, e_minor)*i_minor
                + np.dot(M, e_major)*i_major).tolist()

        model.pattern(xara.StaticPattern(
            FrameLoad(model,
                    shape="Point",
                    elements=[1],
                    offset=(1, *e),
                    basis = load_type,
                    force=f,
                    couple=m,
            )                   
        ))

    analysis = xara.StaticAnalysis(model)
    analysis.analyze()
    return model


def solution(shape, shear=False):

    FX, FY, N = F
    r = e[0]*e_major + e[1]*e_minor
    MY, MZ, T = M + np.cross(r, F)

    EIy = shape["E"]*shape["Iy"]
    EIz = shape["E"]*shape["Iz"]
    EA  = shape["E"]*shape["A"]
    GJ  = shape["G"]*shape["J"]

    UX = FX*L**3/(3*EIz) + MZ*L**2/(2*EIz)
    UY = FY*L**3/(3*EIy) - MY*L**2/(2*EIy)
    uz = N*L/EA

    rx = MY*L/(EIy) - FY*L**2/(2*EIy)
    ry = MZ*L/(EIz) + FX*L**2/(2*EIz)
    rz = T*L/GJ

    if shear:
        GAy = shape["G"]*shape["Ay"]
        GAz = shape["G"]*shape["Az"]

        UX += FX*L/GAy
        UY += FY*L/GAz

    return list(map(float, [UX, UY, uz, rx, ry, rz]))


@pytest.mark.parametrize("load_type", ["global", "local", "node"])
def test_terminal_loads_euler(load_type):
    shape = {
        "E":  E,
        "G":  G,
        "A":  A,
        "Iy": Iy,
        "Iz": Iz,
        "J":  J
    }
    model = create_prism(shape, load_type=load_type, shear=0)
    Umodel = model.nodeDisp(2)
    Usoln = solution(shape)
    assert Umodel[0] == pytest.approx(Usoln[0], abs=1e-10)
    assert Umodel[1] == pytest.approx(Usoln[1], abs=1e-10)
    assert Umodel[2] == pytest.approx(Usoln[2], abs=1e-10)

    assert Umodel[3] == pytest.approx(Usoln[3], abs=1e-10)
    assert Umodel[4] == pytest.approx(Usoln[4], abs=1e-10)
    assert Umodel[5] == pytest.approx(Usoln[5], abs=1e-10)


    print("Solution: ", Usoln)
    print("Result:   ", Umodel)


@pytest.mark.parametrize("load_type", ["global", "local", "node"])
def test_terminal_loads_shear(load_type):
    shape = {
        "E":  E,
        "G":  G,
        "A":  A,
        "Ay": Ay,
        "Az": Az,
        "Iy": Iy,
        "Iz": Iz,
        "J":  J
    }
    model = create_prism(shape, load_type=load_type, shear=1)
    Umodel = model.nodeDisp(2)
    Usoln = solution(shape, shear=True)
    assert Umodel[0] == pytest.approx(Usoln[0], abs=1e-10)
    assert Umodel[1] == pytest.approx(Usoln[1], abs=1e-10)
    assert Umodel[2] == pytest.approx(Usoln[2], abs=1e-10)

    assert Umodel[3] == pytest.approx(Usoln[3], abs=1e-10)
    assert Umodel[4] == pytest.approx(Usoln[4], abs=1e-10)
    assert Umodel[5] == pytest.approx(Usoln[5], abs=1e-10)

    print("Solution: ", Usoln)
    print("Result:   ", Umodel)
    for i in range(6):
        print(f"    U[{i}]: ", (Umodel[i] - Usoln[i])/Usoln[i]*100, "%")



if __name__ == "__main__":
    for load_type in ["node", "global", "local"]:
        print(f"\nRunning load_type={load_type}")
        test_terminal_loads_euler(load_type)
        print()
        test_terminal_loads_shear(load_type)
