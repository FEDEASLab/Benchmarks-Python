#
# Beam with distributed loads
#

import pytest
from xara_prism import AISC_360_C1, check_span


rtol = {"u_tran": 1e-2, "u_long": 1e-12}

analysis_options={
    "test": ("Energy", 1e-18, 2, 0)
}


@pytest.mark.parametrize("shear", [0, 1])
def test_aisc(shear):

    column = AISC_360_C1(shear=shear)


    for basis in ["local", "global"]:
        loading = dict(basis=basis) if basis is not None else None

        check_span(
            column.analyze(3, "ForceFrame", "Linear", ne=4, steps=2, loading_options=loading),
            column.solution("linear"),
            time=1,
            rtol=rtol
        )
        if not shear:
            continue

        # Lagrange-interpolated element; needs at least nen=5 (quartic displacement)
        # for uniform load solution
        check_span(
            column.analyze(3, "ShearFrame", "Linear", ne=4, nen=5, steps=2, loading_options=loading),
            column.solution("linear"),
            time=1,
            rtol=rtol
        )
