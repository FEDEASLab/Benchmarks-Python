#
# Vertical cantilever with terminal loads
#
# Adapted from frame-0002 to check sections in non-principal coordinates
#
import xara
from xara.load import NodalLoad
from xsection.library import Rectangle
import pytest
import numpy as np
from shps.rotor import exp 

L  = 140.0
ne = 2


E = 29e3
G = 11e3
nu = E/(2*G) - 1.0
ElasticMaterial = xara.MultiaxialMaterial("ElasticIsotropic", E=E, G=G)

d = 12.0
b =  6.0

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


def _create_prism(element: str, section, shear, eccentricity="moment", rotation=0.0):

    model = xara.Model(ndm=3, ndf=6)

    model.node(1, (0, 0, 0))
    model.node(2, (0, 0, L))

    model.fix(1, (1, 1, 1, 1, 1, 1))

    R = exp([0, 0, rotation])

    model.geomTransf("Linear", 1, (0, 1, 0))

    model.section(section)

    model.element(element, 1, (1, 2), 
                  section=section, 
                  shear=shear,
                  transform=1)

    f = F
    m = M 
    if eccentricity == "moment":
        r = e[0]*e_major + e[1]*e_minor
        m +=  np.cross(r, F)
    else:
        # TODO: Implement version with eccentricity implemented by 
        # translating shape
        raise NotImplementedError(f"Unsupported eccentricity type: {eccentricity}")

    model.pattern(xara.StaticPattern(
        NodalLoad(model, {
            2: (*f, *m)
        })                   
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



def run_terminal_loads(element, section, shape: dict, shear: int):
    model = _create_prism(element, section, shear=shear)
    Umodel = model.nodeDisp(2)

    Usoln = solution(shape, shear=bool(shear))

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



@pytest.mark.parametrize("element", ["ForceFrame", "ElasticFrame"])
@pytest.mark.parametrize("shear", [0, 1])
def test_shape_rectangle(shear: int, element: str):

    shape = Rectangle(d=d, b=b, 
                      mesh_scale=1/20,
                      mesh_type="T6",
                      material=ElasticMaterial)
    A  = b*d
    GA = G*A
    assert shape.area == pytest.approx(A, abs=1e-12)

    assert shape.poisson_average_area() == pytest.approx(nu, rel=1e-8)

    section = xara.FrameSection("Elastic", shape)

    Cnn = section.cnn()
    Cmm = section.cmm()
    C   = section.constants
    assert C._validate_resultants(["N", "Vy", "Vz", "T", "My", "Mz"])

    assert Cnn[0][0]/E == pytest.approx(A, abs=1e-12)
    assert Cnn[0][0]/E == pytest.approx(C["EA"]/E, abs=1e-12)

    assert Cmm[1][1] == C["EIy"]
    assert Cmm[2][2] == C["EIz"]

    assert Cmm[0][1] == Cmm[0][2] == 0.0
    assert Cmm[1][2]/E == pytest.approx(0.0, abs=1e-10)
    assert Cmm[2][1]/E == pytest.approx(0.0, abs=1e-10)
    assert Cmm[1][1]/E == pytest.approx(b*d**3/12, rel=1e-8)
    assert Cmm[2][2]/E == pytest.approx(d*b**3/12, rel=1e-8)

    assert 0.5 < C["GAy"]/(GA) < 1.0
    assert 0.5 < C["GAz"]/(GA) < 1.0


    shape_dict = {
        "E":  E,
        "G":  G,
        "A":  section.constants["EA"]/E,
        "Ay": section.constants["GAy"]/G,
        "Az": section.constants["GAz"]/G,
        "J":  section.constants["GJ"]/G,
        "Iy": section.constants["EIy"]/E,
        "Iz": section.constants["EIz"]/E,
    }
    run_terminal_loads(element, section, shape_dict, shear)


# @pytest.mark.parametrize("element", ["ForceFrame"]) #, "ElasticFrame"
# @pytest.mark.parametrize("shear", [0, 1])
def _test_shape_rotated(shear: int, element: str):
    angle = np.pi/4

    shape_principal = Rectangle(d=d, b=b, 
                      mesh_scale=1/20,
                      mesh_type="T6",
                      material=ElasticMaterial
    )
    shape_rotated = shape_principal.rotate(angle)

    A  = b*d
    GA = G*A
    assert shape_rotated.area == pytest.approx(A, abs=1e-12)

    assert shape_rotated.poisson_average_area() == pytest.approx(nu, rel=1e-8)

    section = xara.FrameSection("Elastic", shape_rotated)
    section_principal = xara.FrameSection("Elastic", shape_principal)

    Cnn = section.cnn()
    Cmm = section.cmm()
    C   = section.constants
    assert C._validate_resultants(["N", "Vy", "Vz", "T", "My", "Mz"])

    assert Cnn[0][0]/E == pytest.approx(A, abs=1e-12)
    assert Cnn[0][0]/E == pytest.approx(C["EA"]/E, abs=1e-12)

    assert Cmm[1][1] == C["EIy"]
    assert Cmm[2][2] == C["EIz"]

    assert Cmm[0][1] == Cmm[0][2] == 0.0
    assert Cmm[1][2]/E == pytest.approx(0.0, abs=1e-10)
    assert Cmm[2][1]/E == pytest.approx(0.0, abs=1e-10)
    assert Cmm[1][1]/E == pytest.approx(b*d**3/12, rel=1e-8)
    assert Cmm[2][2]/E == pytest.approx(d*b**3/12, rel=1e-8)

    assert 0.5 < C["GAy"]/(GA) < 1.0
    assert 0.5 < C["GAz"]/(GA) < 1.0

    # For computing the reference solution
    shape_dict = {
        "E":  E,
        "G":  G,
        "A":  section_principal.constants["EA"]/E,
        "Ay": section_principal.constants["GAy"]/G,
        "Az": section_principal.constants["GAz"]/G,
        "J":  section_principal.constants["GJ"]/G,
        "Iy": section_principal.constants["EIy"]/E,
        "Iz": section_principal.constants["EIz"]/E,
    }
    run_terminal_loads(element, section, shape_dict, shear)





if __name__ == "__main__":
    for element in ["ForceFrame", "ElasticFrame"]:
        for shear in [0, 1]:
            print(f"\nRunning element={element}, shear={shear}")

            test_shape_rectangle(shear, element)
