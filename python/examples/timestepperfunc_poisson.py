from aa4py import DotDict
from pyioutils.hdfio import load, save
import numpy as np
from collections.abc import Iterable

def timestepperfunc_poisson(xl, yl, fetchOutput, iter_, A, b, calcNorms=False, mpi=None):

#   See andacc.py for a description of what the different values of fetchOutput mean

  if mpi is None:
    x=xl
    y=yl
  else:
    from pympiutils import xloc2glob, xglob2loc
    x=xloc2glob(xl,layout=mpi.layout, MPI=mpi.MPI)

# Run the model (if using MPI, only the root rank does this)
  if mpi is None or (mpi is not None and mpi.rank==0):
    sh=x.shape
    x.shape=(x.shape[0],1)
    
    gy=None

    if x is not None:
      iniFile = "vini"

    if (fetchOutput==-1) or (fetchOutput==0):
      if x is not None:
# 		we save this to compute norms and convergence later
        if calcNorms:
          v = x          
          if fetchOutput==0:
            save({"v": v}, iniFile)

    M = np.diag(1 / (A.diagonal()))

    t1 = np.matmul(M, b)
    t2 = np.matmul(np.matmul(M, A.toarray()), x)
    t3 = 6 * np.matmul(M, np.exp(x))
    gv = t1 - t2 + t3 + x    
    
    if isinstance(fetchOutput, Iterable) and len(fetchOutput) == 0:
        print("Running Model ...")

    else:
        print("Submitting Model ...")

    if (isinstance(fetchOutput, Iterable) and len(fetchOutput) == 0) or fetchOutput:
        print("Finished running model")

        if calcNorms:
          if fetchOutput==1:
            v = load(iniFile).v
          vnorms = np.linalg.norm(gv - v)
        else:
          vnorms = None

        gx = gv

        externalConvergence = DotDict()
        externalConvergence.convdata = DotDict()
        externalConvergence.convdata.testScalar=3.5
    else:
        return None

    x.shape=sh
    gx.shape=sh

  else:
    gx = None
    gy = None
    vnorms = None
    externalConvergence = DotDict()

# Data to send back
  if mpi is None:  
    gxl = gx
    gyl = gy
    vnorms_out = vnorms
    externalConvergence_out = externalConvergence
  else:
    gxl=xglob2loc(gx,layout=mpi.layout, MPI=mpi.MPI)
    gyl = None
#   Broadcast other data to all ranks
    vnorms_out = mpi.comm.bcast(vnorms, root = 0)
    externalConvergence_out = mpi.comm.bcast(externalConvergence, root = 0)
    mpi.comm.Barrier()
    
  return gxl, gyl, vnorms_out, externalConvergence_out  # SPK why is this here and not above in the if clause?
