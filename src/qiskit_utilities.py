from qiskit.opflow import I,X,Y,Z,Zero,PauliTrotterEvolution
from qiskit.circuit import Parameter

from src.sympy_utilities import sympy_expr_to_list
from sympy import expand
from collections import defaultdict

char_to_op={'I': I, 'X': X, 'Y': Y, 'Z': Z}

def pauli_string_to_trotter_step(ps, time):
    return op_to_trotter(pauli_string_to_op(ps),time)

def pauli_string_to_op(ps):
    return pauli_dict_to_op(pauli_string_to_dict(ps))


def pauli_string_to_dict(ps):
    pauli_dict = defaultdict(complex)
    ps = expand(ps)
    #goes through all terms in pauli_strings (aka all a_i in sum_i a_i)
    for arg in ps.args:
        #print(arg)
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

        pauli_dict[key] += coef
        
    return pauli_dict


def pauli_dict_to_op(pd):
    expr = 0
    for pstring, coef in pd.items():
        term = char_to_op[pstring[0]]
        for char in pstring[1:]:
            term = term^char_to_op[char]
        expr += coef.real*term
    return expr

def op_to_trotter(op, time):
    evo_time = Parameter('t')
    evolution_op = (evo_time*op).exp_i()
    fixed_time = evolution_op.bind_parameters({evo_time: time})
    return PauliTrotterEvolution(trotter_mode='suzuki').convert(fixed_time)
    
    