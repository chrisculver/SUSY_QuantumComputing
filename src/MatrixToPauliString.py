from sympy import symbols, expand
from collections import defaultdict


def basis_to_pauli_string(i, j, q):
    if(i=='0' and j=='0'):
        return 0.5*symbols('I^'+q) + 0.5*symbols('Z^'+q)

    elif(i=='0' and j=='1'):
        return 0.5*symbols('X^'+q) + 0.5*1j*symbols('Y^'+q)

    elif(i=='1' and j=='0'):
        return 0.5*symbols('X^'+q) - 0.5*1j*symbols('Y^'+q)

    else:
        return 0.5*symbols('I^'+q) - 0.5*symbols('Z^'+q)


def matrix_to_pauli_strings(matrix, encoding):
    #row-col ordering, m[i][j] i-th row, j-th col
    N = len(matrix)
    pauli_strings=0
    for i in range(0, N):
        for j in range(0, N):
            #print('convert_element check:')
            #print('  args: m[i][j]={}, i={}, j={}, N={}, encoding={}'.format(matrix[i][j],i,j,N,encoding))
            #print('  res: {}'.format(convert_element(matrix[i][j],i,j,N,encoding)))
            pauli_strings += convert_element(matrix[i][j], i, j, N, encoding)
    
    return pauli_strings


def convert_element(elem, i, j, N, encoding):
    nBinDigits = len(format(N-1,'b'))
    i_bin = encoding(i, nBinDigits)
    j_bin = encoding(j, nBinDigits)
    qubit_product = 1
       
    for q in range(0,len(i_bin)):
        qubit_product *= basis_to_pauli_string(i_bin[q], j_bin[q], str(len(i_bin)-1-q))
            
    return elem*qubit_product