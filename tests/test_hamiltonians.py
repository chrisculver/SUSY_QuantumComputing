from src.HamiltonianTerms import *

def test_terms():
  assert (adag(4)*a(4)).all() == number_op(4).all()
