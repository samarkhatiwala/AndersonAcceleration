# export OMP_NUM_THREADS=1
# export MKL_NUM_THREADS=1
# export OPENBLAS_NUM_THREADS=1

from math import ceil
from aa4py import (
    DotDict,
    andacc,
    AndersonAccelerationObj,
    CONVERGED_REASON,
    saveaarestart, loadaarestart
)
from pyioutils.hdfio import load, save
import numpy as np
import sys
import os
from time import time

from poisson_nonlinear_fun import poisson_nonlinear_fun
from timestepperfunc_poisson_xy import timestepperfunc_poisson_xy

useMPI=False

if useMPI:
  from mpi4py import MPI
  from pympiutils import generate_mpi_layout, xglob2loc

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

f = lambda x, y: x * y * 0

bx0 = lambda y: y * 0
bxf = lambda y: y * 0

by0 = lambda x: x * 0
byf = lambda x: x * 0

x0 = 0
xf = 1
y0 = 0
yf = 1
D = np.array([x0, xf, y0, yf])

nx = 32
ny = 32

hx = (xf - 0) / (nx + 1)
x = np.linspace(x0, xf, nx + 2)
y = np.linspace(y0, yf, ny + 2)

X, Y = np.meshgrid(x, y)
X = np.transpose(X)
Y = np.transpose(Y)

Iint = np.arange(1, nx + 1)
Jint = np.arange(1, ny + 1)

A, b = poisson_nonlinear_fun(f, bx0, bxf, by0, byf, D, nx, ny)

v0 = np.zeros(nx * ny)
x0 = v0[0:800]
y0 = v0[800:]

AAparams = DotDict()
AAparams.mMax = 30
AAparams.itmax = 100
AAparams.beta = 1.0
AAparams.restartAANormStagnation = 0
AAparams.restartAANormDiff = 0.05
AAparams.restartAASuccessiveIters = 3

histParams = DotDict()
histParams.nhistfreq = 20
histParams.nhistmax = ceil(AAparams.itmax / histParams.nhistfreq) + 5

startTime = time()

if useMPI:
  mpiaa=DotDict()
  mpiaa.MPI=MPI
  mpiaa.comm=comm
  mpiaa.rank=rank
  mpiaa.size=size
  n=len(x0)
  mpiaa.layout=generate_mpi_layout(nglob=n, nprocs=size, MPI=mpiaa.MPI)
  n=len(y0)
  mpiaa.layouty=generate_mpi_layout(nglob=n, nprocs=size, MPI=mpiaa.MPI)
else:
  mpiaa=None

doPrint = (not useMPI) or (useMPI and mpiaa.rank==0)

fn="aa_restart"
if os.path.isfile(fn + ".h5"):
  if doPrint: print("Reading restart file")
  s=loadaarestart(fn,mpi=mpiaa)
  xl0=s.x
  yl0=s.y
else:
  s = AndersonAccelerationObj({})
  if useMPI:
    xl0=xglob2loc(x0,layout=mpiaa.layout, MPI=mpiaa.MPI)
    yl0=xglob2loc(y0,layout=mpiaa.layouty, MPI=mpiaa.MPI)
  else:
    xl0=x0
    yl0=y0

def checkpointingFunction(i,s,mpi=None):
  if (i>0) and (i % 10 == 0):
    fn="spinup_solution_{:04d}".format(i)
    res_hist = s.res_hist[i, :]
    to_save = {'iter_':i, 'x': s.x, 'res_hist': res_hist}
    mpix=None
    mpiy=None
    if mpi is not None:
      mpix=DotDict({'comm':mpi.comm,'rank':mpi.rank,'layout':mpi.layout,'globalVars':["res_hist"]})
      if s.y is not None:
        mpiy=DotDict({'comm':mpi.comm,'rank':mpi.rank,'layout':mpi.layouty})
    save(to_save,fn,groupname="x",mpi=mpix)
    if s.y is not None:
      to_save = {'y': s.y}
      save(to_save,fn,groupname="y",append=True,mpi=mpiy)
  return 0

histParams.checkpointingFunction=checkpointingFunction  
histParams.trimHistory=False

g = lambda x, y, fetchOutput, iter_: timestepperfunc_poisson_xy(x, y, fetchOutput, iter_, A, b, calcNorms=True, mpi=mpiaa)
suff='_test'
ret = andacc(g, xl0, s, AAparams=AAparams, histParams=histParams, doRestart=False, fileSuff=suff, y=yl0, debug=True, mpi=mpiaa)

xsol, iter_, ysol, converged = ret

endTime = time()

if doPrint: print("Total time: ", endTime-startTime)

fn="aa_restart.h5"
saveaarestart(s,fn,mpi=mpiaa)
