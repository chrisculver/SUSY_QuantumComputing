import math
import numpy as np

def a_element(i,j):
    if(i-1 == j):
        return math.sqrt(i)
    else:
        return 0

def adag_element(i,j):
    if(i+1 == j):
        return math.sqrt(i+1)
    else:
        return 0
    
def n_element(i,j):
    if(i == j):
        return i
    else:
        return 0

def a(n):
  return np.array( [[ a_element(i,j) for i in range(n)] for j in range(n)] )

def adag(n):
  return np.array( [[ adag_element(i,j) for i in range(n)] for j in range(n)] )

def number_op(n):
  return np.array( [[ n_element(i,j) for i in range(n)] for j in range(n)] )
