# 
import xara
import pytest
from xara.units.english import inch, kip, ksi
from math import atan2
import numpy as np
from shps.rotor import exp


# Geometry
H = 144.0 * inch
W = 144.0 * inch

# Section properties
width = 12.0 * inch
height = 12.0 * inch
Iz = (width * height**3) / 12.0
Iy = (height * width**3) / 12.0
A = width * height
J = (width * height**3) / 3.0

# Material
E = 29000 * ksi
G = 11500 * ksi
nu = 0.30

k_penalty = 1.0e10 * kip / inch

# Load
P = -10.0 * kip
angle = atan2(0.6, 0.8) # 0.6435011088 # radians



Reference = {
    (1, 2): 5.811,
    (4, 2): 5.236
}

expected_Fz_jt1 = 5.811
expected_F3_jt4 = 5.236


def error(value, reference):
    return abs(value - reference) / reference * 100.0



def test_penalty_2d():
    # Model
    model = xara.Model('basic', ndm=2, ndf=3)

    # Nodes
    model.node(1, (0.0, 0.0))      # Fixed support
    model.node(2, (0.0, H))        # Top-left
    model.node(3, (W/2, H))        # Load point
    model.node(4, (W,   H))        # Top-right
    model.node(5, (W, 0.0))        # Skewed support

    # Fixed support at node 1
    model.fix(1, (1, 1, 1))

    # Skewed support at node 5 using penalty method
    cos_theta = 0.8
    sin_theta = 0.6

    R = exp([0.0, 0.0, angle])
    ground_node = 100
    model.node(ground_node, (W, 0.0))
    model.fix(ground_node, (1, 1, 1))

    k_penalty = 1.0e10 * kip / inch
    model.uniaxialMaterial('Elastic', 1001, k_penalty)

    # Normal direction for constraint
    nx, ny, = R[:2,0]
    tx, ty = cos_theta, sin_theta

    model.element('zeroLength', 1001, (ground_node, 5), 
                  mat=1001, 
                  dir=2,
                  x=(nx, ny, 0))#, 0.0, tx, ty, 0.0))

    # Frame section and elements
    model.material('ElasticIsotropic', 1, E, nu)
    model.section("ElasticFrame", 1, E=E, G=G, A=A, Iy=Iy, Iz=Iz, J=J)
    model.geomTransf("Linear", 1)

    model.element('PrismFrame', 1, (1, 2), section=1, transform=1)
    model.element('PrismFrame', 2, (2, 3), section=1, transform=1)
    model.element('PrismFrame', 3, (3, 4), section=1, transform=1)
    model.element('PrismFrame', 4, (4, 5), section=1, transform=1)

    # Load
    model.timeSeries('Constant', 1)
    model.pattern('Plain', 1, 1)
    model.load(3, (0.0, P, 0.0))

    # Analysis
    model.system('BandGeneral')
    model.numberer('RCM')
    model.constraints('Transformation')
    model.integrator('LoadControl', 1.0)
    model.algorithm('Newton')
    model.test('Energy', 1e-10, 10)
    model.analysis('Static')
    model.analyze(1)


    # Results

    print("\n" + "="*60)
    print("Model E: Skewed Support - Example 1-005e")
    print("="*60)

    # Reactions
    model.reactions()
    Fz_jt1 = model.nodeReaction(1, 2)/kip
    R4 = model.nodeReaction(ground_node)[:2]
    r4 = R[:2,:2].T @ R4 / kip
    F3_jt4 = r4[1]


    # Verification
    print("\n" + "-"*60)
    print("Verification against Example 1-005e:")
    print("-"*60)

    error_Fz = error(Fz_jt1, expected_Fz_jt1)*100 
    error_F3 = error(F3_jt4, expected_F3_jt4)*100 

    print(f"  Fz (jt. 1):  {Fz_jt1:8.3f} kip  (expected {expected_Fz_jt1:.3f})  Error: {error_Fz:.2f}%")
    print(f"  F3 (jt. 4):  {F3_jt4:8.3f} kip  (expected {expected_F3_jt4:.3f})  Error: {error_F3:.2f}%")

    if error_Fz < 1.0 and error_F3 < 1.0:
        print("\nPASSED")
        assert True
    else:
        print("\nFAILED")
        assert False

    print("="*60)
    

def test_penalty_3d():
    # Model

    model = xara.Model('basic', ndm=3, ndf=6)

    # Nodes
    model.node(1, (0.0, 0.0, 0.0))      # Fixed support
    model.node(2, (0.0, 0.0, H))        # Top-left
    model.node(3, (W/2, 0.0, H))        # Load point
    model.node(4, (W,   0.0, H))        # Top-right
    model.node(5, (W,   0.0, 0.0))      # Skewed support

    # Fixed support at node 1
    model.fix(1, (1, 1, 1, 1, 1, 1))

    # Skewed support at node 5 using penalty method

    R = exp([0.0, -angle, 0.0])
    ground_node = 100
    model.node(ground_node, (W, 0.0, 0.0))
    model.fix(ground_node,  (1, 1, 1, 1, 1, 1))

    model.uniaxialMaterial('Elastic', 1001, k_penalty)

    # Normal direction for constraint
    nx, ny, nz = R[:,0] # == R @ [1, 0, 0] == R*Ex
    tx, ty, tz = R[:,1]

    model.element('zeroLength', 1001, (ground_node, 5), 
                  mat=1001, 
                  dir=3,
                  orient=(nx, ny, nz, tx, ty, tz))

    # Frame section and elements
    model.material('ElasticIsotropic', 1, E, nu)
    model.section("ElasticFrame", 1, E=E, G=G, A=A, Iy=Iy, Iz=Iz, J=J)
    model.geomTransf("Linear", 1, (0, 1, 0)) # Columns
    model.geomTransf("Linear", 2, (0, 1, 0)) # Beams

    model.element('PrismFrame', 1, (1, 2), section=1, transform=1)
    model.element('PrismFrame', 2, (2, 3), section=1, transform=2)
    model.element('PrismFrame', 3, (3, 4), section=1, transform=2)
    model.element('PrismFrame', 4, (4, 5), section=1, transform=1)

    # Load
    model.timeSeries('Constant', 1)
    model.pattern('Plain', 1, 1)
    model.load(3, (0.0, 0.0, P, 0.0, 0.0, 0.0))

    # Analysis
    model.system('BandGeneral')
    model.numberer('RCM')
    model.integrator('LoadControl', 1.0)
    model.algorithm('Newton')
    model.test('Energy', 1e-8, 10)
    model.analysis('Static')
    model.analyze(1)


    # Results

    print("\n" + "="*60)
    print("Model E: Skewed Support - Example 1-005e")
    print("="*60)

    # Reactions
    model.reactions()
    Fz_jt1 = model.nodeReaction(1, 3)/kip
    R4 = model.nodeReaction(ground_node)[:3]
    r4 = R.T @ R4 / kip
    F3_jt4 = r4[2]


    # Verification
    print("\n" + "-"*60)
    print("Verification against Example 1-005e:")
    print("-"*60)

    error_Fz = error(Fz_jt1, expected_Fz_jt1)*100 
    error_F3 = error(F3_jt4, expected_F3_jt4)*100 

    print(f"  Fz (jt. 1):  {Fz_jt1:8.3f} kip  (expected {expected_Fz_jt1:.3f})  Error: {error_Fz:.2f}%")
    print(f"  F3 (jt. 4):  {F3_jt4:8.3f} kip  (expected {expected_F3_jt4:.3f})  Error: {error_F3:.2f}%")

    if error_Fz < 1.0 and error_F3 < 1.0:
        print("\nPASSED")
        assert True
    else:
        print("\nFAILED")
        assert False



@pytest.mark.parametrize("constraint_solver", ["Transformation", "Auto"])
def test_constrain_3d(constraint_solver):
    # Model

    model = xara.Model('basic', ndm=3, ndf=6)

    # Nodes
    model.node(1, (0.0, 0.0, 0.0))      # Fixed support
    model.node(2, (0.0, 0.0, H))        # Top-left
    model.node(3, (W/2, 0.0, H))        # Load point
    model.node(4, (W,   0.0, H))        # Top-right
    model.node(5, (W,   0.0, 0.0))      # Skewed support

    # Fixed support at node 1
    model.fix(1, (1, 1, 1, 1, 1, 1))

    # Skewed support at node 5 using penalty method

    R = exp([0.0, -angle, 0.0])
    ground_node = 100
    model.node(ground_node, (W, 0.0, 0.0))
    model.fix(ground_node,  (0, 1, 1, 0, 0, 0))


    model.constrain((ground_node, 5), rotate=[0, -angle, 0])

    # Frame section and elements
    model.material('ElasticIsotropic', 1, E, nu)
    model.section("ElasticFrame", 1, E=E, G=G, A=A, Iy=Iy, Iz=Iz, J=J)
    model.geomTransf("Linear", 1, (0, 1, 0)) # Columns
    model.geomTransf("Linear", 2, (0, 1, 0)) # Beams

    model.element('PrismFrame', 1, (1, 2), section=1, transform=1)
    model.element('PrismFrame', 2, (2, 3), section=1, transform=2)
    model.element('PrismFrame', 3, (3, 4), section=1, transform=2)
    model.element('PrismFrame', 4, (4, 5), section=1, transform=1)

    # Load
    model.timeSeries('Constant', 1)
    model.pattern('Plain', 1, 1)
    model.load(3, (0.0, 0.0, P, 0.0, 0.0, 0.0))

    # Analysis
    model.system('BandGeneral')
    model.numberer('RCM')
    model.constraints(constraint_solver)
    model.integrator('LoadControl', 1.0)
    model.algorithm('Newton')
    model.test('Energy', 1e-8, 10)
    model.analysis('Static')
    model.analyze(1)


    # Results

    print("\n" + "="*60)
    print("Model E: Skewed Support - Example 1-005e")
    print("="*60)

    # Reactions
    model.reactions()
    Fz_jt1 = model.nodeReaction(1, 3)/kip
    R4 = model.nodeReaction(5)[:3]
    r4 = R.T @ R4 / kip
    F3_jt4 = r4[2]
    print(f"R4: {R4}, r4: {r4}")


    # Verification
    print("\n" + "-"*60)
    print("Verification against Example 1-005e:")
    print("-"*60)

    error_Fz = error(Fz_jt1, expected_Fz_jt1)*100 
    error_F3 = error(F3_jt4, expected_F3_jt4)*100 

    print(f"  Fz (jt. 1):  {Fz_jt1:8.3f} kip  (expected {expected_Fz_jt1:.3f})  Error: {error_Fz:.2f}%")
    print(f"  F3 (jt. 4):  {F3_jt4:8.3f} kip  (expected {expected_F3_jt4:.3f})  Error: {error_F3:.2f}%")

    if error_Fz < 1.0 and error_F3 < 1.0:
        print("\nPASSED")
        assert True
    else:
        print("\nFAILED")
        assert False

    print("="*60)



if __name__ == "__main__":
    test_penalty_2d()
    test_penalty_3d()
    test_constrain_3d('Transformation')
    test_constrain_3d('Auto')