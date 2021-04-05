from sympy import preorder_traversal, Symbol


# TODO: Some(all?) of these need to be moved to the MatrixToPauliString function
# as they are not generic and for sympy expressions of a certain type.
# Maybe I should make a PauliStringSympy class to hold all of this stuff safely.
def unique_sympy_symbols(expr):
    return set([ str(elem).split('^')[0] for elem in expr.free_symbols ])
    
def max_sympy_exponent(expr):
    return max(set([ int(str(elem).split('^')[1]) for elem in expr.free_symbols ]))

def sympy_expr_to_list(expr):
    arg_list=[]

    for a in preorder_traversal(expr):
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
            elem*=Symbol('I^'+str(iden))
        
        new_h+=elem
        
    return new_h