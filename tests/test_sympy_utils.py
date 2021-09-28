from src.sympy_utilities import *

def test_ps_list():
    expr = 1.2*sp.Symbol('X^0')*sp.Symbol('Y^1')
    expr -= 2.1j*sp.Symbol('I^0')*sp.Symbol('Y^1')
    
    expected = defaultdict(complex)
    expected['XY']=1.2
    expected['IY']=-2.1j
    
    assert pauli_strings_to_list(expr) == expected