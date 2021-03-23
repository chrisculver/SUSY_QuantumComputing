from sympy import symbols,expand
from src.BinaryEncodings import *
from src.MatrixToPauliString import *

def test_standard_encoding():
    assert standard_encode(7,3) == '111'
    assert standard_encode(7,5) == '00111'
    
    
def test_gray_code():
    assert gray_code(7,3) == '100'
    assert gray_code(7,5) == '00100'


    
def test_standard_encoding_diagonal_matrix():
    # equation 25 and 26 of https://arxiv.org/pdf/1909.12847.pdf
    # are being tested here
    matrix = [[0,0,0],[0,1,0],[0,0,2]]
    mtops = MatrixToPauliString(matrix)
    mtops.convert(standard_encode)
    result = 0.75*symbols('I^0')*symbols('I^1')+0.25*symbols('Z^0')*symbols('I^1')-0.25*symbols('I^0')*symbols('Z^1')-0.75*symbols('Z^0')*symbols('Z^1')
            
    assert expand(mtops.sum_of_paulis) == result
    
    matrix = [[0,0,0,0],[0,1,0,0],[0,0,2,0],[0,0,0,3]]
    mtops = MatrixToPauliString(matrix)
    mtops.convert(standard_encode)
    result = 1.5*symbols('I^0')*symbols('I^1')+0.5*symbols('Z^0')*symbols('I^1')-1.0*symbols('I^0')*symbols('Z^1')
            
    assert expand(mtops.sum_of_paulis) == result
