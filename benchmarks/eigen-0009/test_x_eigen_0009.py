# ArpackSOE::addM with Transformation-condensed equalDOF.
#
# 2D plane-strain quad column (H=10, N=10), base fixed, left/right nodes at
# every level joined by equalDOF so the mesh is a pure 1D shear beam. Closed
# form: T1 = 4H/Vs with Vs=sqrt(G/rho). G=1e5, rho=2 => Texact = 0.178885 s.
#
# https://github.com/peer-open-source/xara/pull/113
#
import xara
import math
import pytest

H   = 10.0
N   = 10
dy  = H/N
w   = 1.0
nu  = 0.0
G   = 100000.0
E   = 2.0*G*(1.0+nu)
rho = 2.0
Vs     = math.sqrt(G/rho)
Texact = 4.0*H/Vs
pi     = math.pi

# ARPACK

def _create_model():
    model = xara.Model(ndm=2, ndf=2)
    model.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for j in range(N+1):
        y = -H + j*dy
        model.node(2*j+1,  0.0, y)
        model.node(2*j+2,   w , y)

    model.fix( 1, (1, 1))
    model.fix( 2, (1, 1))

    for j in range(1, N+1):
        model.element("quad", j,
                      2*(j-1)+1, 2*(j-1)+2, 2*j+2, 2*j+1, 
                      1.0, "PlaneStrain", 1
        )
        model.equalDOF(2*j+1, 2*j+2, 1, 2)

    model.constraints("Transformation")
    model.numberer("Plain")
    model.system("UmfPack")
    return model



def _run_arpack():
    model = _create_model()
    return 2.0*pi/math.sqrt(model.eigen(3)[0])


# fullGenLapack control
def _run_lapack():
    model = _create_model()
    return 2.0*pi/math.sqrt(model.eigen("-fullGenLapack", 3)[0])


def test():
    T_lapack = _run_lapack()
    T_arpack = _run_arpack()

    rExact  = T_arpack/Texact
    rLapack = T_arpack/T_lapack
    sqrt2   = math.sqrt(2.0)

    print("T_arpack=%.6f  T_lapack=%.6f  Texact=%.6f  T_arpack/Texact=%.6f  T_arpack/T_lapack=%.6f" \
        .format(T_arpack, T_lapack, Texact, rExact, rLapack))

    # verify value $rLapack  1.0 1.0e-6

    if abs(rLapack - 1.0) > 1.0e-6:
        assert False

    if abs(rExact - 1.0) > 1.03e-3:
        assert False, rExact-1.0



