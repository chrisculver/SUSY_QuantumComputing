from src.HamiltonianTerms import *

import math

def test_matrices():
    adag4 = adag(4)
    
    assert (adag4 == np.array([[0,0,0,0],[1,0,0,0],[0,math.sqrt(2),0,0],[0,0,math.sqrt(3),0]])).all()
    
    a4 = a(4)
    
    assert (a4 == np.array([[0,1,0,0],[0,0,math.sqrt(2),0],[0,0,0,math.sqrt(3)],[0,0,0,0]])).all()
    
    n4 = number_op(4)
    
    assert (n4 == np.array([[0,0,0,0],[0,1,0,0],[0,0,2,0],[0,0,0,3]])).all()