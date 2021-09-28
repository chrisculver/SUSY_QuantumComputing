from src.pauli_matrices import *

import sympy as sp

def test_paulis():
    qubit = 3
    rule = pauliQubitToMatrix(qubit)
    keys=[]
    for s in ['X','Y','Z','I']:
        keys.append(sp.Symbol(s+'^'+str(qubit)))
    
    assert (list(rule.keys())==keys)
    
    assert (rule[sp.Symbol('X^3')]==np.array([[0,1],[1,0]])).all()
    assert (rule[sp.Symbol('Y^3')]==np.array([[0,-1j],[1j,0]])).all()
    assert (rule[sp.Symbol('Z^3')]==np.array([[1,0],[0,-1]])).all()
    assert (rule[sp.Symbol('I^3')]==np.array([[1,0],[0,1]])).all()