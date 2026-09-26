# Cantilever subjected to uniform torsion (6-DOF)
#
# This test ensures that torsional stiffness is correctly captured in frame elements
# when sections are defined by shape objects.
#
import xara
import xsection.library as xs
from xsection.analysis import WarpingAnalysis
import xara.units.iks as units
import pytest

def create_model(element, section, shape, material):
    print(f"Element = {element}, Section = {section}")
    model = xara.Model(ndm=3, ndf=6)

    L = 10*units.ft

    model.node(1, 0,0,0)
    model.node(2, L,0,0)
    model.fix(1, 1,1,1, 1,1,1)


    model.material(material)

    model.section(xara.FrameSection(section, shape))

    model.geomTransf("Linear", 1, (0,0,1))

    model.element(element, 1, (1,2), section=1, transform=1, shear=0)

    return model


def test_frame_0011():
    T = 1000
    L = 10*units.ft
    shape = "W14x48"

    material = xara.MultiaxialMaterial(
        "ElasticIsotropic",
        E = 29e3*units.ksi,
        nu = 0.3
    )

    shape = xs.from_aisc(shape, units=units, mesh_scale=1/10, material=material)

    G = material["E"]/(2.0*(1+material["nu"]))
    GJ = WarpingAnalysis(shape).twist_rigidity()
    print(GJ/G)


    for element in "CubicFrame", "ForceFrame", "PrismFrame":

        for section in "Elastic", "Fiber", "MultiaxialFiber":
            model = create_model(element, section, shape, material)

            # Apply a torque loading
            loads = xara.NodalLoad(model, {2: [0,0,0, T,0,0]})
            model.pattern(xara.StaticPattern(loads))

            # Analyze
            xara.StaticAnalysis(model, test=("Residual", 1e-8, 2)).analyze()

            print(T*L/(model.state.u(2,4)))
            assert T*L/model.state.u(2,4) == pytest.approx(GJ, rel=1e-5)


if __name__ == "__main__":
    test_frame_0011()

