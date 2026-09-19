#
# Adapted from:
#  https://openseesdigital.com/2025/02/23/two-node-links-awakening/
#
import xara
from math import isclose
import pytest
import numpy as np


kt = 10 # Translational stiffness
kr = 20 # Rotational stiffness
P =   5 # Load
L = 10  # Link length
c = 0.5 # Shear distance

def test_orient_2D():
    model = xara.Model(ndm=2, ndf=3)

    model.uniaxialMaterial("Elastic",1,kt)
    model.uniaxialMaterial("Elastic",2,kr)

    model.node(0,0,0)
    model.node(1,0,0)
    model.node(2,0,L)
    model.fix(0,1,1,1)
    model.fix(1,0,1,0)
    model.fix(2,0,1,0)

    # Zero-length link between nodes 0 and 1
    model.element("zeroLength", 1, (0,1), mat=(1,2), dir=(2,3), orient=(0,1,0))
    # Two-node link between nodes 0 and 2
    model.element("twoNodeLink",2, (0,2), mat=(1,2), dir=(2,3), shearDist=c)

    model.pattern("Plain",1,"Constant")
    model.load(1, (P,0,0), pattern=1)
    model.load(2, (P,0,0), pattern=1)

    model.analysis("Static")
    model.analyze(1)

    # X-displacement
    u1 = model.nodeDisp(1,1)
    u2 = model.nodeDisp(2,1)

    assert isclose(u1,P/kt)
    # assert isclose(u2,P/kt + P*L*(1-c)/kr*L*(1-c))
    assert u2 == pytest.approx(P/kt + P*L*(1-c)/kr*L*(1-c), abs=1e-10)

    # Rotation
    u1 = model.nodeDisp(1,3)
    u2 = model.nodeDisp(2,3)

    assert isclose(u1,0,abs_tol=1e-10)
    # assert isclose(u2,-P*L*(1-c)/kr)
    assert u2 == pytest.approx(-P*L*(1-c)/kr, abs=1e-10)



EX = np.array([1.0, 0.0, 0.0])
EY = np.array([0.0, 1.0, 0.0])


def _as_tuple(v):
    a = np.asarray(v, dtype=float)
    return tuple(float(x) for x in a)


def _HatSO3(v):
    x, y, z = v
    return np.array([
        [0.0, -z,  y],
        [z,   0.0, -x],
        [-y,  x,   0.0],
    ])


def _ExpSO3(v):
    """
    Return the rotation matrix corresponding to the given rotation vector v.
    """
    v = np.asarray(v, dtype=float)
    theta = np.linalg.norm(v)
    if theta < 1.0e-14:
        return np.eye(3)

    k = v / theta
    K = _HatSO3(k)
    return np.eye(3) + np.sin(theta)*K + (1.0 - np.cos(theta))*(K @ K)


def _LogSO3(R):
    """Return the rotation vector corresponding to the given rotation matrix R. 
    This is just a special case of the general matrix logarithm like MATLAB's logm. 
    Equivalent to scipy.spatial.transform.Rotation.from_matrix(R).as_rotvec().
    """
    R = np.asarray(R, dtype=float)

    cos_theta = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    if theta < 1.0e-14:
        return np.zeros(3)

    if np.pi - theta < 1.0e-8:
        A = 0.5 * (R + np.eye(3))
        axis = np.array([
            np.sqrt(max(A[0, 0], 0.0)),
            np.sqrt(max(A[1, 1], 0.0)),
            np.sqrt(max(A[2, 2], 0.0)),
        ])

        if axis[0] > 1.0e-8:
            axis[1] = np.copysign(axis[1], R[0, 1])
            axis[2] = np.copysign(axis[2], R[0, 2])
        elif axis[1] > 1.0e-8:
            axis[2] = np.copysign(axis[2], R[1, 2])

        axis /= np.linalg.norm(axis)
        return theta * axis

    axis = np.array([
        R[2, 1] - R[1, 2],
        R[0, 2] - R[2, 0],
        R[1, 0] - R[0, 1],
    ]) / (2.0 * np.sin(theta))

    return theta * axis


def _localize_state(model, node, R):
    d = np.asarray(model.nodeDisp(node), dtype=float)
    u = (R.T @ d[:3]).tolist()
    r = (R.T @ d[3:]).tolist()
    return u, r


def _create_model(R):
    R = np.asarray(R, dtype=float)
    if R.shape != (3, 3):
        raise ValueError("R must be a 3x3 matrix.")
    if not np.allclose(R.T @ R, np.eye(3), atol=1.0e-12):
        raise ValueError("R must be orthogonal.")
    if not np.isclose(np.linalg.det(R), 1.0, atol=1.0e-12):
        raise ValueError("R must be a proper rotation matrix.")

    rotvec = _LogSO3(R).tolist()

    # Reference embedding of the original 2D problem:
    #   link axis  = +Y
    #   load dir   = +X
    #
    # choose local x = +Y and local y = +X, so local z = -Z.
    # This preserves the sign convention of the original 2D rotation check.
    link_axis = R @ EY
    load_dir  = R @ EX

    p0 = np.zeros(3)
    p2 = L * EY

    model = xara.Model(ndm=3, ndf=6)

    model.uniaxialMaterial("Elastic", 1, kt)
    model.uniaxialMaterial("Elastic", 2, kr)

    # Physical nodes
    model.node(0, _as_tuple(R @ p0))
    model.node(1, _as_tuple(R @ p0))
    model.node(2, _as_tuple(R @ p2))

    # Support nodes used to realize the rotated support directions
    model.node(101, _as_tuple(R @ p0))
    model.node(102, _as_tuple(R @ p2))

    model.fix(0,   (1, 1, 1, 1, 1, 1))
    model.fix(101, (0, 1, 1, 1, 1, 0))
    model.fix(102, (0, 1, 1, 1, 1, 0))

    model.constrain(101, 1, rotate=rotvec)
    model.constrain(102, 2, rotate=rotvec)

    # Zero-length link between nodes 0 and 1
    model.element(
        "zeroLength", 1, (0, 1),
        mat=(1, 2), dir=(2, 6),
        x=_as_tuple(link_axis),
        y=_as_tuple(load_dir),
    )

    # Two-node link between nodes 0 and 2
    model.element(
        "twoNodeLink", 2, (0, 2),
        mat=(1, 2), dir=(2, 6),
        y=_as_tuple(load_dir),
        shearDist=(c, 0.5),
    )

    Fx, Fy, Fz = map(float, P * load_dir)
    model.pattern("Plain", 1, "Constant")
    model.load(1, (Fx, Fy, Fz, 0.0, 0.0, 0.0), pattern=1)
    model.load(2, (Fx, Fy, Fz, 0.0, 0.0, 0.0), pattern=1)

    model.constraints("Transformation")
    model.numberer("Plain")
    model.system("BandGeneral")
    model.integrator("LoadControl", 1.0)
    model.algorithm("Linear")
    model.analysis("Static")

    return model


def test_xy_3D():
    R = _ExpSO3(np.array([0.31, -0.27, 0.19]))

    model = _create_model(R)
    model.analyze(1)

    arm = L * (1.0 - c)

    u1, r1 = _localize_state(model, 1, R)
    u2, r2 = _localize_state(model, 2, R)

    assert u1 == pytest.approx([P / kt, 0.0, 0.0], abs=1e-10)
    assert u2 == pytest.approx([P / kt + P * arm**2 / kr, 0.0, 0.0], abs=1e-10)

    assert r1 == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)
    assert r2 == pytest.approx([0.0, 0.0, -P * arm / kr], abs=1e-10)


    u11 = model.eleResponse(1, "material", 1, "strain")

    assert u11 == pytest.approx(P/kt, abs=1e-10)


if __name__ == "__main__":
    test_xy_3D()
