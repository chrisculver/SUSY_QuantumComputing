from MatrixToPauliString import MatrixToPauliString
from PauliStringToQISKITCircuit import PauliStringsToQISKITCircuit

class MatrixToQISKITCircuit:
    def __init__(self, matrix):
        self.matrix = matrix
        
    def convert(self, encoding):
        mtops = MatrixToPauliString(self.matrix)
        mtops.convert(encoding)
        pstoc = PauliStringsToQISKITCircuit(mtops.sum_of_paulis)
        pstoc.string_to_terms()
        pstoc.make_qiskit_circuit( len(self.matrix) )
        self.circuit = pstoc.circuit