from src.MatrixToPauliString import MatrixToPauliString
from src.PauliStringToQISKITCircuit import PauliStringsToQISKITCircuit

class MatrixToQISKITCircuit:
    def __init__(self, matrix, evo_time):
        self.matrix = matrix
        self.evo_time = evo_time
        
    def convert(self, encoding):
        mtops = MatrixToPauliString(self.matrix)
        mtops.convert(encoding)
        pstoc = PauliStringsToQISKITCircuit(mtops.sum_of_paulis)
        pstoc.string_to_terms()
        pstoc.make_qiskit_circuit( self.evo_time )
        self.quantum_register = pstoc.quantum_register
        self.circuit = pstoc.circuit
