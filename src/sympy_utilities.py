import src.HamiltonianTerms as hmats

import sympy as sp


p, q = sp.symbols('p, q', commutative=False)
a, adag = sp.symbols('a, ad', commutative=False)
b, bdag = sp.symbols('b, bd', commutative=False)

m = sp.Symbol('m')
g = sp.Symbol('g')

qp_to_ada = {q: 0.5*sp.sqrt(2/m)*(a + adag), 
            p: complex(0,1)*sp.sqrt(2*m)*(adag - a)/2}



ada_matrices = {a*adag: sp.Matrix(hmats.a(5)*hmats.adag(5)),
                adag*a: sp.Matrix(hmats.adag(5)*hmats.a(5))}



adaga_sub = {a: sp.Matrix(hmats.a(9)), adag: sp.Matrix(hmats.adag(9))}


def convert(expr):
    new_expr = sp.Matrix(sp.zeros(9))
    for elem in expr.args:
        new_expr += convert_term_to_matrix(elem)

    return new_expr


def convert_term_to_matrix(term):
    new_elem = 1
    for elem in term.args:
        is_operator = False
        for i in sp.preorder_traversal(elem):
            if( (i is a) or (i is adag) ):
                is_operator = True
                
        if is_operator:
            lst = elem.args

            if(len(lst)>1):
                for i in range(lst[1]):
                    new_elem*=lst[0].subs(adaga_sub)
            else:
                new_elem*=elem.subs(adaga_sub)
        else:
            new_elem*=elem
        
    return new_elem


    
    
    
    
    
    
    
    
# TODO: Some(all?) of these need to be moved to the MatrixToPauliString function
# as they are not generic and for sympy expressions of a certain type.
# Maybe I should make a PauliStringSympy class to hold all of this stuff safely.
def unique_sympy_symbols(expr):
    return set([ str(elem).split('^')[0] for elem in expr.free_symbols ])
    
def max_sympy_exponent(expr):
    return max(set([ int(str(elem).split('^')[1]) for elem in expr.free_symbols ]))

def sympy_expr_to_list(expr):
    arg_list=[]

    for a in sp.preorder_traversal(expr):
        if(a!=expr):
            arg_list.append(a)

    return arg_list

def identity_qubit_padded_H(ham):
    new_h = 0
    
    qubits = [i for i in range(max_sympy_exponent(ham)+1)]
    
    for elem in ham.args:
        qubits_with_pauli = []
        for e in elem.free_symbols:
            qubits_with_pauli.append(int(str(e).split('^')[1]))
            
        #np.setdiff1d doesn't seem to work
        ids_to_pad = []
        for qi in qubits:
            if qi not in qubits_with_pauli:
                ids_to_pad.append(qi)
    
        for iden in ids_to_pad:
            elem*=sp.Symbol('I^'+str(iden))
        
        new_h+=elem
        
    return new_h