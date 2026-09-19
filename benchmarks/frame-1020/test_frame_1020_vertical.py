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
from xara_prism import Span1, BasicShape, check_span
from xara_prism.loading import ShearForce
from xara_prism.solutions.reissner_terminal import ReissnerTerminalCantilever
import numpy as np


L = 100.0
A = 1.61538e8
I = 3.5e7
Pmax = 150e3
shape = BasicShape(
    E=1,
    G=1,
    A=A,
    Ay=A,
    Az=A,
    Iy=I,#*100,
    Iz=I,
    J =I
)

span = Span1(length=L, shape=shape, shear=1, smax=1.0, loads=[
    ShearForce(Pmax, x=L, scale="linear", follower=True)
])

solution = ReissnerTerminalCantilever(span)

steps = 500
slice = 50
time = np.linspace(1/steps, 1.0, steps)


def test_cosserat():
    rtol = None
    atol = {"u_tran": 5, "u_long": 6}
    check_span(
        span.analyze(3, "ExactFrame", "Linear", ne=10, nen=2, steps=steps, 
                    analysis_options={"test": ("NormDispIncr", 1e-12, 10)}),
        solution,
        atol=atol,
        space=[span.length],
        time=time[::slice]
    )

def test_identity():
    atol = {"u_tran": 5, "u_long": 6}
    check_span(
        span.analyze(3, "CosseratFrame01", "Identity", ne=10, nen=2, steps=steps, 
                     analysis_options={"test": ("NormDispIncr", 1e-12, 10)}),
        solution,
        atol=atol,
        space=[span.length],
        time=time[::slice]
    )


def test_spherical():
    atol = {"u_tran": 5, "u_long": 6}
    check_span(
        span.analyze(3, "CosseratFrame01", "Spherical", ne=10, nen=2, steps=steps, 
                    analysis_options={"test": ("NormDispIncr", 1e-10, 15)}),
        solution,
        atol=atol,
        space=[span.length],
        time=time[::slice]
    )


def test_corotational():
    rtol = None
    atol = {"u_tran": 5.5, "u_long": 6}
    check_span(
        span.analyze(3, "ForceFrame", "Corotational02", ne=10, nen=2, steps=steps, 
                    analysis_options={"test": ("NormDispIncr", 1e-12, 10)}),
        solution,
        atol=atol,
        space=[span.length],
        time=time[::slice]
    )
    check_span(
        span.analyze(3, "ForceFrame", "Corotational03", ne=10, nen=2, steps=steps, 
                    analysis_options={"test": ("NormDispIncr", 1e-12, 10, 0)}),
        solution,
        atol=atol,
        space=[span.length],
        time=time[::slice]
    )

if __name__ == "__main__":
    test_cosserat()
    test_corotational()