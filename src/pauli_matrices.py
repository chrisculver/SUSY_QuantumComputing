import numpy as np
import sympy as sp

pauliX = np.array([[0,1],[1,0]])
pauliY = np.array([[0,-1j],[1j,0]])
pauliZ = np.array([[1,0],[0,-1]])


pauliSymbolToMatrix = {'X': pauliX, 'Y': pauliY, 'Z': pauliZ, 'I': np.eye(2)}

def pauliQubitToMatrix(n):
    nf=str(n)
    return {sp.Symbol('X^'+nf): pauliX, sp.Symbol('Y^'+nf): pauliY, sp.Symbol('Z^'+nf): pauliZ, sp.Symbol('I^'+nf): np.eye(2)}