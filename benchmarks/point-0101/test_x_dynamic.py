"""Verification of Newmark integrators against Chopra's hand calculations.

References
----------
Chopra, A.K. "Dynamics of Structures: Theory and Applications",
Prentice Hall, 4th Edition, 2012, Chapter 5.
  - Linear:    Examples 5.3 / 5.4
  - Nonlinear: Examples 5.5 / 5.6
"""

from math import sqrt, pi
import pytest
import xara

# tolerances

TOL = 1.0e-3
DT  = 0.1

# model properties (shared by all cases)

M          = 0.2533
K          = 10.0
DAMP_RATIO = 0.05

# reference data

LINEAR_CASES = {
    "average_acceleration": {
        "integrator_cmd": "integrator Newmark 0.5 0.25 -alpha 1",
        "disp": [0.0437, 0.2326, 0.6121, 1.0825, 1.4309, 1.4231, 0.9622, 0.1908, -0.6044, -1.1442],
        "vel":  [0.8733, 2.9057, 4.6833, 4.7260, 2.2421, -2.3996, -6.8192, -8.6092, -7.2932, -3.5026],
        "accel": [17.4666, 23.1801, 12.3719, -11.5175, -38.1611, -54.6722, -33.6997, -2.1211, 28.4423, 47.3701],
    },
    "linear_acceleration": {
        "integrator_cmd": "integrator Newmark 0.5 [expr 1.0/6.0] -alpha 1",
        "disp": [0.0300, 0.2193, 0.6166, 1.1130, 1.4782, 1.4625, 0.9514, 0.1273, -0.6954, -1.2208],
        "vel":  [0.8995, 2.9819, 4.7716, 4.7419, 2.1802, -2.6911, -7.1468, -8.7758, -7.1539, -3.0508],
        "accel": [17.9904, 23.6566, 12.1372, -12.7305, -39.9425, -56.0447, -33.0689, 0.4892, 31.9491, 50.1114],
    },
}

NONLINEAR_CASES = {
    "Newton": {
        "algorithm_cmd": "algorithm Newton",
        "disp": [0.0437, 0.2326, 0.6121, 1.1143, 1.6214, 1.9891, 2.0951, 1.9240, 1.5602, 1.415],
        "vel":  [0.8733, 2.9057, 4.6833, 5.3624, 4.7792, 2.5742, -0.4534, -2.960, -4.3075, -4.0668],
        "accel": [17.4666, 23.1801, 12.3719, 1.2103, -12.8735, -31.2270, -29.3242, -20.9876, -5.7830, 10.5962],
    },
    "ModifiedNewton": {
        "algorithm_cmd": "algorithm ModifiedNewton",
        "disp": [0.0437, 0.2326, 0.6121, 1.1143, 1.6214, 1.9891, 2.0951, 1.9240, 1.5602, 1.414],
        "vel":  [0.8733, 2.9057, 4.6833, 5.3623, 4.7791, 2.5741, -0.4534, -2.960, -4.3076, -4.0668],
        "accel": [17.4666, 23.1801, 12.3719, 1.2095, -12.8734, -31.2270, -29.3242, -20.9879, -5.7824, 10.5969],
    },
}


# helpers

def _build_model(yield_disp: float = 0.0) -> xara.Model:
    """Build a 1-DOF SDOF model (linear when *yield_disp* is 0)."""
    wn = sqrt(K / M)

    model = xara.Model(ndm=1, ndf=1)
    model.node(1, 0.0)
    model.node(2, 0.0, mass=M)

    if yield_disp == 0.0:
        model.uniaxialMaterial("Elastic", 1, K)
    else:
        model.uniaxialMaterial("ElasticPP", 1, K, yield_disp)

    model.element("zeroLength", 1, 1, 2, mat=1, dir=1)
    model.fix(1, 1)

    a0 = 2.0 * wn * DAMP_RATIO
    model.rayleigh(a0, 0.0, 0.0, 0.0)
    return model


def _apply_trig_load(model: xara.Model):
    """Apply the trigonometric force time-series used in all Chopra examples."""
    model.timeSeries("Trig", 1, 0.0, 0.6, 1.2, factor=10.0)
    model.pattern("Plain", 1, 1, load={2: [1.0]})


def _build_analysis(model, integrator_cmd, algo_cmd, *, linear=False):
    """Configure a transient analysis on *model*."""
    model.constraints("Plain")
    model.numberer("Plain")
    model.eval(integrator_cmd)
    max_iter = 2 if linear else 6
    model.test("NormDispIncr", 1.0e-4, max_iter, 0)
    model.eval(algo_cmd)
    model.system("ProfileSPD")
    model.analysis("Transient")


def _run_and_collect(model, n_steps):
    """Run *n_steps* and return lists of (disp, vel, accel) at node 2."""
    disps, vels, accels = [], [], []
    for _ in range(n_steps):
        model.analyze(1, DT)
        disps.append(model.nodeDisp(2, 0))
        vels.append(model.nodeVel(2, 0))
        accels.append(model.nodeAccel(2, 0))
    return disps, vels, accels


def _assert_close(computed, expected, label, tol=TOL):
    """Check every entry in *computed* vs *expected* within TOL."""
    for i, (c, e) in enumerate(zip(computed, expected)):
        assert abs(c - e) <= tol, (
            f"{label} mismatch at step {i + 1}: "
            f"computed {c:.6f}, expected {e:.6f}, "
            f"diff {abs(c - e):.2e} > tol {tol:.2e}"
        )


# tests

class TestNewmarkLinear:
    """Chopra §5.4 – Newmark on a linear SDOF system."""

    @pytest.mark.parametrize("case_name", LINEAR_CASES)
    def test_displacement(self, case_name):
        case = LINEAR_CASES[case_name]
        model = _build_model(yield_disp=0.0)
        _apply_trig_load(model)
        _build_analysis(model, case["integrator_cmd"], "algorithm Newton", linear=True)
        disps, _, _ = _run_and_collect(model, len(case["disp"]))
        _assert_close(disps, case["disp"], "displacement")

    @pytest.mark.parametrize("case_name", LINEAR_CASES)
    def test_velocity(self, case_name):
        case = LINEAR_CASES[case_name]
        model = _build_model(yield_disp=0.0)
        _apply_trig_load(model)
        _build_analysis(model, case["integrator_cmd"], "algorithm Newton", linear=True)
        _, vels, _ = _run_and_collect(model, len(case["vel"]))
        _assert_close(vels, case["vel"], "velocity", tol=1e-1)

    @pytest.mark.parametrize("case_name", LINEAR_CASES)
    def test_acceleration(self, case_name):
        case = LINEAR_CASES[case_name]
        model = _build_model(yield_disp=0.0)
        _apply_trig_load(model)
        _build_analysis(model, case["integrator_cmd"], "algorithm Newton", linear=True)
        _, _, accels = _run_and_collect(model, len(case["accel"]))
        _assert_close(accels, case["accel"], "acceleration", tol=5e-2)


class TestNewmarkNonlinear:
    """Chopra §5.7 – Newmark on a nonlinear (elastic-perfectly-plastic) SDOF."""

    INTEGRATOR_CMD = "integrator Newmark 0.5 0.25"
    YIELD_DISP = 0.75
    N_STEPS = 9  # original test ran 9 steps

    @pytest.mark.parametrize("case_name", NONLINEAR_CASES)
    def test_displacement(self, case_name):
        case = NONLINEAR_CASES[case_name]
        model = _build_model(yield_disp=self.YIELD_DISP)
        _apply_trig_load(model)
        _build_analysis(model, self.INTEGRATOR_CMD, case["algorithm_cmd"], linear=False)
        disps, _, _ = _run_and_collect(model, self.N_STEPS)
        _assert_close(disps, case["disp"][:self.N_STEPS], "displacement")

    @pytest.mark.parametrize("case_name", NONLINEAR_CASES)
    def test_velocity(self, case_name):
        case = NONLINEAR_CASES[case_name]
        model = _build_model(yield_disp=self.YIELD_DISP)
        _apply_trig_load(model)
        _build_analysis(model, self.INTEGRATOR_CMD, case["algorithm_cmd"], linear=False)
        _, vels, _ = _run_and_collect(model, self.N_STEPS)
        _assert_close(vels, case["vel"][:self.N_STEPS], "velocity", tol=5e-2)

    @pytest.mark.parametrize("case_name", NONLINEAR_CASES)
    def test_acceleration(self, case_name):
        case = NONLINEAR_CASES[case_name]
        model = _build_model(yield_disp=self.YIELD_DISP)
        _apply_trig_load(model)
        _build_analysis(model, self.INTEGRATOR_CMD, case["algorithm_cmd"], linear=False)
        _, _, accels = _run_and_collect(model, self.N_STEPS)
        _assert_close(accels, case["accel"][:self.N_STEPS], "acceleration", tol=5e-2)

