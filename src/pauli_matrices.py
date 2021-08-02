import numpy as np

pauliX = np.array([[0,1],[1,0]])
pauliY = np.array([[0,-1j],[1j,0]])
pauliZ = np.array([[1,0],[0,-1]])

pauliSymbolToMatrix = {'X': pauliX, 'Y': pauliY, 'Z': pauliZ, 'I': np.eye(2)}