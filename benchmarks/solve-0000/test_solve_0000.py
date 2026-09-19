import xara
import pytest
import numpy as np


def create_model(system):
    E = 0.5e8
    A = 1.0
    L = 144.0

    model = xara.Model(ndm=2, ndf=2)
    model.node(1, (0.0, 0.0))
    model.node(2, ( L,  L))
    model.node(3, ( L, 0.0))

    model.fix(1, (1, 1))
    model.fix(3, (1, 1))

    material = xara.UniaxialMaterial("Elastic", 1, E)
    model.material(material)

    section = xara.TrussSection("Truss", material=material, area=A)
    model.section(section)

    model.element("Truss", 1, (1, 2), section=1)
    model.element("Truss", 2, (2, 3), section=1)

    analysis = xara.StaticAnalysis(model, system=system)
    model.pattern("Plain", 1, "Linear")
    model.load(2, (0.0, -1.0), pattern=1)
    assert analysis.analyze(1) == 0

    return model


@pytest.mark.parametrize("system", ["FullGeneral", "ProfileSPD", "BandSPD", "FullGeneral"])
def test_solve(system):
    model = create_model(system)
    A = model.getTangent()

    rng = np.random.default_rng(seed=42)
    B = rng.random(A.shape[0])

    Xref = np.linalg.solve(A, B)

    X = model.solveA(B)

    print(Xref)
    print(X)
    for i in range(len(X)):
        assert X[i] == pytest.approx(Xref[i], rel=1e-6, abs=1e-6)


if __name__ == "__main__":
    test_solve("FullGeneral")
    test_solve("ProfileSPD")
    test_solve("BandSPD")
    test_solve("FullGeneral")
