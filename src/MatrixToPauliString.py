from sympy import symbols, expand
from src.utilities import *
from collections import defaultdict


def basis_to_pauli_string(i, j, q):
    if(i=='0' and j=='0'):
        return 0.5*symbols('I^'+q) + 0.5*symbols('Z^'+q)

    elif(i=='0' and j=='1'):
        return 0.5*symbols('X^'+q) + 0.5*1j*symbols('Z^'+q)

    elif(i=='1' and j=='0'):
        return 0.5*symbols('X^'+q) - 0.5*1j*symbols('Z^'+q)

    else:
        return 0.5*symbols('I^'+q) - 0.5*symbols('Z^'+q)


class MatrixToPauliString:
    def __init__(self, matrix):
        self.matrix=matrix #row-col ordering, m[i][j] i-th row, j-th col
        self.N = len(self.matrix)
        self.pauli_strings=0
        self.pauli_terms = defaultdict(complex)
    
    def convert(self, encoding):
        for i in range(0, self.N):
            for j in range(0, self.N):
                self.pauli_strings += self.convert_element(i, j, encoding)
    
    def convert_element(self, i, j, encoding):
        nBinDigits = len(format(self.N-1,'b'))
        i_bin = encoding(i, nBinDigits)
        j_bin = encoding(j, nBinDigits)
        qubit_product = 1
        
        for q in range(0,len(i_bin)):
            qubit_product *= basis_to_pauli_string(i_bin[q], j_bin[q], str(len(i_bin)-1-q))
            
        return self.matrix[i][j]*qubit_product
    
    
    def pauli_strings_as_list(self):
        self.pauli_strings = expand(self.pauli_strings)
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
        
        return self.pauli_terms
    

    
