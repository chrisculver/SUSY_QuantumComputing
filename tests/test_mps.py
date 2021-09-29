from src.MatrixToPauliString import *
from src.BinaryEncodings import *
from src.sympy_utilities import *

import numpy as np
import sympy as sp

def test_matrix_to_pauli_string():
    #test for number matrix
    matrix = np.array([[0,0,0,0],[0,1,0,0],[0,0,2,0],[0,0,0,3]])
    pauli_strings = matrix_to_pauli_strings(matrix, standard_encode)
    expected = 1.5*sp.Symbol('I^0')*sp.Symbol('I^1')
    expected -= 1.0*sp.Symbol('I^0')*sp.Symbol('Z^1')
    expected -= 0.5*sp.Symbol('I^1')*sp.Symbol('Z^0')
    assert sp.expand(pauli_strings)

    assert (getMatrix(sp.expand(pauli_strings))==matrix).all()
    #test for matrix with off-diagonal elements
    
    matrix = np.array([[1,0,2],[0,3,0],[4,0,5]])
    ps = matrix_to_pauli_strings(matrix, standard_encode)
    ps = sp.expand(ps)
    sym = lambda p,q: sp.Symbol(p+'^'+str(q))
    expected = 2.25*sym('I',0)*sym('I',1)
    expected += 1.5*sym('I',0)*sym('X',1)
    expected -= 0.5*sp.I*sym('I',0)*sym('Y',1)
    expected -= 0.25*sym('I',0)*sym('Z',1)
    expected += 0.75*sym('Z',0)*sym('I',1)
    expected += 1.5*sym('Z',0)*sym('X',1)
    expected -= 0.5*sp.I*sym('Z',0)*sym('Y',1)
    expected -= 1.75*sym('Z',0)*sym('Z',1)
    
    assert sp.expand(ps)==expected
    
    assert (getMatrix(sp.expand(ps)) == np.array([[1,0,2,0],[0,3,0,0],[4,0,5,0],[0,0,0,0]])).all()