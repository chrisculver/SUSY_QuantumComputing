from sympy import symbols


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
        self.sum_of_paulis=0
    
    def convert(self, encoding):
        for i in range(0, self.N):
            for j in range(0, self.N):
                self.sum_of_paulis += self.convert_element(i, j, encoding)
    
    def convert_element(self, i, j, encoding):
        nBinDigits = len(format(self.N-1,'b'))
        i_bin = encoding(i, nBinDigits)
        j_bin = encoding(j, nBinDigits)
        qubit_product = 1
        
        for q in range(0,len(i_bin)):
            qubit_product *= basis_to_pauli_string(i_bin[q], j_bin[q], str(len(i_bin)-1-q))
            
        return self.matrix[i][j]*qubit_product
    

    
