import numpy as np

def norm(x,mpi=None):
  return np.sqrt(np.dot(x,x)) # linalg.norm(x)

def dot(x,y,mpi=None):
  return np.dot(x,y)

def multtranspose(A,x,mpi=None):
# y = A^T x = (x^T A)^T
  y = np.matmul(x.transpose(),A).transpose()
  return y
