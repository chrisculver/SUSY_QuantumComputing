from sympy import expand, preorder_traversal
from collections import defaultdict
from qiskit.aqua.operators import WeightedPauliOperator
from qiskit import transpile, QuantumRegister, QuantumCircuit
from src.utilities import max_sympy_exponent
import math

def sympy_expr_to_list(expr):
    arg_list=[]

    for a in preorder_traversal(expr):
        if(a!=expr):
            arg_list.append(a)

    return arg_list


class PauliStringsToQISKITCircuit:
    def __init__(self, pauli_strings):
        self.pauli_strings = expand(pauli_strings)

        #pauli_terms is a dictionary with keys corresponding to unique pauli strings,
        #with values corresponding to the coefficient
        self.pauli_terms = defaultdict(complex)


    def string_to_terms(self):
        #goes through all terms in pauli_strings (aka all a_i in sum_i a_i)
        for arg in self.pauli_strings.args:
            #each term is a product of a number * paulis, convert each elem to list
            arg_list=sympy_expr_to_list(arg)


            paulis = [ '' for i in range(0,len(arg_list)-1) ]
            start = 1
            coef = complex(arg_list[0])

            if(len(str(arg_list[1]))==1):
                coef*=1j
                start=2
                paulis.remove('')

            for ai in range(start, len(arg_list)):
                symbol = str(arg_list[ai])
                parts = symbol.split('^')
                paulis[ int(parts[1]) ] = parts[0]

            key = ''
            for i in paulis:
                key += i

            self.pauli_terms[key] += coef


    def make_qiskit_circuit(self, time):
        pauli_dict={'paulis': []}
        for k,v in self.pauli_terms.items():
            pauli_dict['paulis'].append({"coeff": {"real": v.real, "imag": v.imag}, "label":k})

        op = WeightedPauliOperator.from_dict(pauli_dict)
        
        nqubits = max_sympy_exponent(self.pauli_strings)+1
        print(nqubits)
        self.quantum_register = QuantumRegister(nqubits, 'q')
        self.circuit = QuantumCircuit(self.quantum_register)
        self.circuit = op.evolve(quantum_registers=self.quantum_register, evo_time = time)
        self.circuit = transpile(self.circuit, basis_gates = ['cx', 'u1', 'u2', 'u3', 'H', 'X', 'Y', 'Z'])
