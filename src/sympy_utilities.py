import src.HamiltonianTerms as hmats
import src.MatrixToPauliString as mps
from src.pauli_matrices import *
from collections import defaultdict

import sympy as sp
import numpy as np

import copy


p, q = sp.symbols('p, q', commutative=False)
a, adag = sp.symbols('a, ad', commutative=False)
b, bdag = sp.symbols('b, bd', commutative=False)

m = sp.Symbol('m')
g = sp.Symbol('g')
mu = sp.Symbol('mu')

qp_to_ada = {q: 0.5*sp.sqrt(2/m)*(a + adag), 
            p: complex(0,1)*sp.sqrt(2*m)*(adag - a)/2}



def adaga_sub(cutoff): 
    return {a: sp.Matrix(hmats.a(cutoff)), adag: sp.Matrix(hmats.adag(cutoff))}



def convert_to_matrix(expr, cutoff, buffer):
    
    new_expr = np.zeros([cutoff+buffer,cutoff+buffer])

    if type(expr)==sp.core.add.Add:
        for elem in expr.args:
            new_expr = new_expr + convert_term_to_matrix(elem,cutoff+buffer)
            
    elif type(expr)==sp.core.mul.Mul:
        tmp=convert_term_to_matrix(expr,cutoff+buffer)
        new_expr = new_expr+tmp
        
    elif type(expr)==sp.core.numbers.Float:
        new_expr = new_expr + convert_term_to_matrix(expr,cutoff+buffer)
    else:
        raise ValueError('Cannot convert type {} to matrix'.format(type(expr)))
        
        
    return np.array(new_expr.tolist())[:cutoff,:cutoff]


def convert_term_to_matrix(term, cutoff):

    new_elem = np.eye(cutoff)
    has_aadag = False
    for elem in term.args:
        for i in sp.preorder_traversal(elem):
            if( (i is a ) or (i is adag) ):
                has_aadag=True
    
    if has_aadag and type(term)==sp.core.mul.Mul:
        for elem in term.args:
            is_operator = False
            for i in sp.preorder_traversal(elem):
                if( (i is a) or (i is adag) ):
                    is_operator = True
                
            if is_operator:
                lst = elem.args

                if(len(lst)>1): # this should handle pow?
                    for i in range(lst[1]):
                        new_elem=np.matmul(new_elem,lst[0].subs(adaga_sub(cutoff)))
                else:

                    new_elem=np.matmul(new_elem,elem.subs(adaga_sub(cutoff))) 
            else:

                new_elem=elem*new_elem
        
        return new_elem
    
    elif has_aadag and type(term)!=sp.core.mul.Mul:
        raise ValueError('Havent implemented doing convert_to_matrix for type {}'.format(type(term)))
    
    else:
        if(term.args==()):
            new_elem=new_elem*term
        else:
            for elem in term.args:
                new_elem=new_elem*elem
        return new_elem*np.eye(cutoff)
    


def commutation_rhs(fq):
    return sp.Symbol('Z^'+str(fq))
    
class Hamiltonian():
    def __init__(self, bosonic, fermionic, params, cutoff, encoding):
        self.bosonic = bosonic
        self.cutoff = cutoff
        #TODO: calculate the buffer...ensor charg
        self.buffer = 16

        #let's just store the bosonic hamiltonian in all interesting forms
        self.harmonic = sp.expand(self.bosonic.subs(qp_to_ada)).subs(params)
        self.bmatrix = convert_to_matrix(self.harmonic, cutoff, self.buffer)
        self.bosonPauliStrings = sp.expand(sp.N(mps.matrix_to_pauli_strings(self.bmatrix, encoding)))
        
        #now lets do the fermionic part
        self.fermionic = sp.expand(fermionic.subs(qp_to_ada).subs(params))

        self.fq = max_sympy_exponent(self.bosonPauliStrings)+1 
        #sub in the bdag and b, not sure this always happens
        
        self.get_fermionic_matrix()
        self.fermionPauliStrings = sp.expand(sp.N(mps.matrix_to_pauli_strings(self.fmatrix, encoding)))
    
        
        self.hamMatrix=np.kron(np.eye(2),self.bmatrix)+self.fmatrix
        self.pauliStrings=sp.N(sp.expand(mps.matrix_to_pauli_strings(self.hamMatrix,encoding)))
        self.pauliStrings = self.pauliStrings.xreplace(dict([(n,0) for n in self.pauliStrings.atoms(sp.Float) if abs(n) < 1e-12]))
        self.qubitMatrix=getMatrix(self.pauliStrings)


    
    
    def get_fermionic_matrix(self):
        mat=convert_to_matrix(self.fermionic,self.cutoff,self.buffer)
        #print(mat)
        self.fmatrix=-np.kron(pauliZ,mat)
        
    def solve_exactly(self):
        sol = np.linalg.eig(np.array(self.hamMatrix,dtype=float))
        self.eigenvalues = sol[0]
        self.eigenvectors = sol[1]
        
        idx = self.eigenvalues.argsort()
        self.eigenvalues=self.eigenvalues[idx]
        self.eigenvectors=self.eigenvectors[idx]
        

        
    
    
    
    
    
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

def identity_qubit_padded(ham,nqubits):    

    qubits = [i for i in range(nqubits+1)]
    
    if type(ham)==sp.core.add.Add:
        new_h = 0

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

    
    
    elif type(ham)==sp.core.mul.Mul:
        new_h = ham
        qubits_with_pauli=[]
        for e in ham.free_symbols:
            qubits_with_pauli.append(int(str(e).split('^')[1]))
            
        ids_to_pad=[]
        for qi in qubits:
            if qi not in qubits_with_pauli:
                ids_to_pad.append(qi)
        
        for iden in ids_to_pad:
            new_h*=sp.Symbol('I^'+str(iden))

        return new_h

    
    
    

    
def pauli_strings_to_list(pauli_strings):
    if(type(pauli_strings)!=sp.core.add.Add):
        raise ValueError('not implemented for non-sp-add types')
    pauli_strings = sp.expand(pauli_strings)
    ps_list = defaultdict(complex)
    #goes through all terms in pauli_strings (aka all a_i in sum_i a_i)
    for arg in pauli_strings.args:
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

        ps_list[key] += coef
        
    return ps_list



def getMatrix(pauliString):
    Ndim=2**(max_sympy_exponent(pauliString.args[0])+1)
    hamMat=np.zeros([Ndim,Ndim],dtype=complex)
    for elem in pauliString.args:
        arg_list = sympy_expr_to_list(elem)
        
        coef = complex(arg_list[0])
        paulis = [ '' for i in range(0,len(arg_list)-1) ]    
        start = 1
        if(len(str(arg_list[1]))==1):
            coef*=1j
            start=2
            paulis.remove('')
        
        for ai in range(start, len(arg_list)):
            symbol = str(arg_list[ai])
            parts = symbol.split('^')
            paulis[ int(parts[1]) ] = parts[0]
        
        res = pauliSymbolToMatrix[paulis[len(paulis)-1]]
        for p in reversed(paulis[:-1]):
            res = np.kron(res,pauliSymbolToMatrix[p])
        
        hamMat+=coef*res
        
    return hamMat