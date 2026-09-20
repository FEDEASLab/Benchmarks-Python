
# - https://openseesdigital.com/2025/08/10/a-model-of-inconsistency/
# - https://web.archive.org/web/20240818123710/https://portwooddigital.com/2024/08/18/secant-accelerated-newton-algorithm/
# - https://openseesdigital.com/2024/08/18/secant-accelerated-newton-algorithm/

import xara
import pytest

VERBOSITY = 9

def create_model():

    ops = xara.Model(ndm=1, ndf=1)
    
    ops.node(0,0); 
    ops.fix(0,1)
    ops.node(1,0)
    ops.node(2,0)

    ops.uniaxialMaterial('Steel01',1,10,10,0.1)
    ops.uniaxialMaterial('Steel01',2,4,2,0.5)
    ops.uniaxialMaterial('Steel01',3,7,7,0)

    # Spring elements
    ops.element('zeroLength',1,0,1,'-mat',1,'-dir',1)
    ops.element('zeroLength',2,1,2,'-mat',2,'-dir',1)
    ops.element('zeroLength',3,0,2,'-mat',3,'-dir',1)

    # Dummy elastic material
    ops.uniaxialMaterial('Elastic',0,0)

    # Diagonal '-1' stiffness
    ops.uniaxialMaterial('Penalty',4,0,-1,'-noStress')
    ops.element('zeroLength',4,0,1,'-mat',4,'-dir',1)

    # Off-diagonal '+0.5' stiffness
    ops.uniaxialMaterial('Penalty',5,0,-0.5,'-noStress')
    ops.element('zeroLength',5,1,2,'-mat',5,'-dir',1)
    ops.uniaxialMaterial('Penalty',6,0,0.5,'-noStress')
    ops.element('zeroLength',6,0,1,'-mat',6,'-dir',1)
    ops.element('zeroLength',7,0,2,'-mat',6,'-dir',1)

    ops.pattern('Plain',1,"Constant")
    ops.load(1,6)
    ops.load(2,12)

    return ops


def test_krylov():
    ops = create_model()

    ops.test('RelativeNormUnbalance',1e-4,8,VERBOSITY)
    ops.algorithm('KrylovNewton')
    ops.analysis('Static')
    
    ops.analyze(1)

    assert ops.nodeDisp(1) == pytest.approx(2.0, rel=1e-8)
    assert ops.nodeDisp(2) == pytest.approx(5.0, rel=1e-8)


def test_secant():
    ops = create_model()

    ops.test('RelativeNormUnbalance',1e-6,20,VERBOSITY)
    ops.algorithm('SecantNewton')
    ops.analysis('Static')
    
    ops.analyze(1)

    assert ops.nodeDisp(1) == pytest.approx(2.0, rel=1e-4)
    assert ops.nodeDisp(2) == pytest.approx(5.0, rel=1e-4)


if __name__ == "__main__":
    VERBOSITY = 1
    test_krylov()
    test_secant()
