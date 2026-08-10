from aa4py import DotDict
from pyioutils.hdfio import load, save
import numpy as np
from collections.abc import Iterable

def timestepperfunc_poisson_xy(xl, yl, fetchOutput, iter_, A, b, calcNorms=False, mpi=None):

#   See andacc.py for a description of what the different values of fetchOutput mean

  # First gather global vectors and assemble combined vector for evaluation function below
  if mpi is None:
    x=np.concatenate((xl,yl))
  else:
    # Should this import be here?
    from pympiutils import xloc2glob, xglob2loc
    xg=xloc2glob(xl,layout=mpi.layout, MPI=mpi.MPI)
    yg=xloc2glob(yl,layout=mpi.layouty, MPI=mpi.MPI)
    if mpi.rank==0:
      x=np.concatenate((xg,yg))

# Evaluate the function (if using MPI, only the root rank does this)
  if mpi is None or (mpi is not None and mpi.rank==0):
    sh=x.shape
    x.shape=(x.shape[0],1)
    
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
    vnorms = None
    externalConvergence = DotDict()

  if mpi is not None: mpi.comm.Barrier()
    
# Data to send back
  if mpi is None:  
    gxl = gx[0:800]
    gyl = gx[800:]
    vnorms_out = vnorms
    externalConvergence_out = externalConvergence
  else:
    # split global solution gx to (x,y)
    if mpi.rank == 0:
      gxg=gx[0:800]
      gyg=gx[800:]
    else:
      gxg=None
      gyg=None

    mpi.comm.Barrier()
    gxl=xglob2loc(gxg,layout=mpi.layout, MPI=mpi.MPI)
    gyl=xglob2loc(gyg,layout=mpi.layouty, MPI=mpi.MPI)

#   Broadcast other data to all ranks
    vnorms_out = mpi.comm.bcast(vnorms, root = 0)
    externalConvergence_out = mpi.comm.bcast(externalConvergence, root = 0)
    mpi.comm.Barrier()
    
  return gxl, gyl, vnorms_out, externalConvergence_out  # SPK why is this here and not above in the if clause?
