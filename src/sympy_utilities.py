from sympy import preorder_traversal

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