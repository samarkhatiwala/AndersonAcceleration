from aa4py import (
    DotDict,
    andacc,
    AndersonAccelerationObj,
    CONVERGED_REASON,
    saveaarestart, loadaarestart
)
from pyioutils.hdfio import load, save
import numpy as np
from math import ceil

# Define fixed point function
def g(x, y, fetchOutput, iter_):
  gx=np.array([(x[0]**2 + x[1]**2 + 8.)/10., (x[0]*x[1]**2 + x[0] + 8.)/10.])
  return gx, None, None, {}

AAparams = DotDict()
AAparams.mMax = 2
AAparams.itmax = 10
AAparams.beta = 1.0

s = AndersonAccelerationObj({})
x0=np.array([0.,0.])

def checkpointingFunction(i,s,mpi=None):
  if (i % 2 == 0):
    print(f"Solution at iter {i}: {s.x}")
  if (i>0) and (i % 5 == 0):
    fn="aa_checkpoint_{:04d}".format(i)
    save(s,fn,mpi=mpi,debug=False)
  return 0

histParams = DotDict()
# Uncomment the following to store history and switch on checkpointing
histParams.nhistfreq = 2
histParams.nhistmax = ceil(AAparams.itmax / histParams.nhistfreq) + 5
histParams.checkpointingFunction=checkpointingFunction  

suff='_test'
ret = andacc(g, x0, s, AAparams=AAparams, histParams=histParams, doRestart=False, fileSuff=suff, y=None, debug=True)
  
xsol, iter_, ysol, converged = ret
  
fn="aa_restart.h5"
saveaarestart(s,fn)
  
