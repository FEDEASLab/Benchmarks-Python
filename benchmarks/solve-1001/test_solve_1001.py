# Problem from Section 9.9.2 of [1]
#
# [1] M. A. Crisfield, Non-linear finite element analysis of solids and structures. 
#     Volume 1: Essentials, Repr. Chichester: Wiley, 2001.
#

import xara
import pytest

reference = -1.0097631086926619

Reference = [
   {"Iter":     1, "EnergyIncr":   0.0508013, "Residual":    0.106934, "Correction":    0.950143},
   {"Iter":     2, "EnergyIncr": 0.000328586, "Residual":   0.0122893, "Correction":    0.053475},
   {"Iter":     3, "EnergyIncr": 4.46957e-06, "Residual":  0.00145456, "Correction":  0.00614558}
]

def setup(test_name, tol,
          element = "ExactTruss", # "CorotTruss" #
          strain = "RE",
          search_type="InitialInterpolated", 
          algorithm="NewtonLineSearch",
          mnr=True):
    print(algorithm)
    print(test_name, search_type, "(Modified Newton)" if mnr else "(Newton-Raphson)")


    E = 0.5e8

    model = xara.Model("basic", ndm=2, ndf=2)

    model.node(1,    0.0, 0.0)
    model.node(2, 2500.0, 25.0)

    model.fix(1, 1, 1)
    model.fix(2, 1, 0)

    model.uniaxialMaterial("Elastic", 1, E)
    model.section("Truss", 1, material=1, area=1.0)

    Strains = {
        "Linear": 0,
        "GL":  1,
        "RE":  2
    }
    if "Exact" in element:
        model.element(element, 1, 1, 2, section=1, 
                      strain=Strains[strain])
    else:
        model.element(element, 1, 1, 2, section=1)

    pattern = xara.StaticPattern([
        xara.NodalLoad(model, {
            2: (0.0, -1.0)
        })
    ])
    model.pattern(pattern)

    # model.pattern("Plain", 1, "Linear")
    # model.load(2, 0.0, -1.0)

    model.system("FullGeneral")
    model.integrator("LoadControl", 1.9)

    if algorithm == "NewtonLineSearch":
        model.algorithm(
            "NewtonLineSearch",
            0.8,
            "-minEta", 0.01,
            "-maxEta", 25,
            "-prediction-tangent", "current", #"initial", #
            "-correction-tangent", "predictor" if mnr else  "current",
            pFlag=0,
            type=search_type
        )
    else:
        model.algorithm(algorithm)

    model.test(test_name, tol, 20, 1)
    model.analysis("Static")

    return model


def verify(value, reference, tol):
    error = abs(value - reference)
    assert error <= tol, (
        f"value = {value:.16g}, "
        f"reference = {reference:.16g}, "
        f"error = {error:.3e}"
    )



def test_solve_1001_criteria():
    Tests = [
        ("EnergyIncr", 1.0e-5),
        ("Correction", 1.0e-2),
        ("RelativeNormUnbalance", 1.0e-3),
        ("Residual", 1.0e-2),
    ]

    for test_name, tol in Tests:
        model = setup(test_name, tol)
        #
        model.analyze(1)

        verify(model.nodeDisp(2, 2), reference, 1.0e-8)

        norms = model.testNorms()
        for i, norm in enumerate(norms):
            if test_name == "RelativeNormUnbalance":
                verify(norm, Reference[i]["Residual"]/1.9, 1e-6)
            else:
                verify(norm, Reference[i][test_name], 1e-6)


@pytest.mark.parametrize("search_type", [
    "InitialInterpolated",
    # "RegulaFalsi"
])
def test_solve_1001_line_search_modified_newton(search_type):
    # model = setup("Residual", 1.0e-2)
    model = setup("RelativeNormUnbalance", 1.0e-3, search_type=search_type, mnr=True)
    assert model.analyze(6) == 0


@pytest.mark.parametrize("search_type", [
    "InitialInterpolated",
    # "RegulaFalsi"
])
def test_solve_1001_line_search_newton_raphson(search_type):
    model = setup("RelativeNormUnbalance", 1.0e-3, search_type=search_type, mnr=False)
    assert model.analyze(6) == 0



def test_solve_1001_quasi(algorithm="BFGS"):
    model = setup("RelativeNormUnbalance", 1.0e-3, algorithm=algorithm, mnr=False)
    assert model.analyze(5) == 0


if __name__ == "__main__":
    test_solve_1001_criteria()

    print("\n\n")
    test_solve_1001_line_search_modified_newton(search_type="InitialInterpolated")

    test_solve_1001_quasi("BFGS")

    test_solve_1001_quasi("KrylovNewton")

    test_solve_1001_line_search_newton_raphson(search_type="InitialInterpolated")

    test_solve_1001_line_search_newton_raphson(search_type="Bisection")

    test_solve_1001_line_search_newton_raphson(search_type="Armijo")

    # test_solve_1001_line_search_newton_raphson(search_type="RegulaFalsi")

    # test_solve_1001_line_search_modified_newton(search_type="Bisection")

    # test_solve_1001_line_search_newton_raphson(search_type="Secant")

    # test_solve_1001_line_search_modified_newton(search_type="RegulaFalsi")
