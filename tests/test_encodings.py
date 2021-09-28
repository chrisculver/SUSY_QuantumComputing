from src.BinaryEncodings import *

def test_standard():
    assert standard_encode(7,3) == '111'
    assert standard_encode(7,5) == '00111'
    assert standard_encode(6,3) == '110'
    
def test_gray():
    assert gray_code(7,3) == '100'
    assert gray_code(7,5) == '00100'
    assert gray_code(6,3) == '101'
