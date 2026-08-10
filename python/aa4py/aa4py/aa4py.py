'''
This file implements Anderson Acceleration (AA). It is specifically designed for the 
computation of equilibrium solutions of periodically-forced ocean and land biogeochemical 
models (see Khatiwala (2023; https://doi.org/10.1029/2022MS003447) and Khatiwala (2024; 
https://doi.org/10.1126/sciadv.adn2839)) but can be used for any fixed point iteration 
problem. It is based on AndAcc.m originally written by Homer F. Walker (listed in H. Walker, 
Anderson acceleration: Algorithms and implementations, Worcester Polytechnic Institute 
Mathematical Sciences Department Research Report MS-6-15-50, June, 2011). The original 
code has been reimplemented in an object oriented manner with extensive modifications 
as per Khatiwala (2023, 2024). These include algorithm checkpointing; enforcement of positivity 
constraints on the solution vector; the ability to apply AA to a user-specified subset of 
the model state vector (called "x/y split" in the code); and MPI-based parallelization. 
Notes:
1) Saving/loading the AA state for restarts, as well as checkpointing, requires installing 
the pyioutils modules ('pip install pyioutils'). This module provides hdfio.py, a set of 
wrapper functions around the h5py module (built against your system HDF5 library).
2) Parallelization is performed via the mpi4py module, which you must build against your 
system MPI library. It also requires installing the pympiutils module ('pip install pyioutils'). 
You must ensure that the h5py module is installed against a MPI-enabled HDF5 library, and 
that the mpi4py and h5py modules are built using the same MPI library and associated compilers.
3) The qrdelete and related functions have been translated from Reichel and Gragg (1990; Algorithm 
686: FORTRAN subroutines for updating the QR decomposition, ACM Transactions on Mathematical 
Software (TOMS), Vol. 16 (4), https://dl.acm.org/doi/10.1145/98267.98291) and parallelized. 
4) The positivity constraint uses the cvxopt module (https://cvxopt.org/) and the (included) 
lsqlin.py wrapper written by Valeriy Vishnevskiy and released under a WTFPL license 
(http://maggotroot.blogspot.ch/2013/11/constrained-linear-least-squares-in.html). The version 
included here has been modified slightly from the original. This feature is only practical for 
relatively small problems (e.g., land models) and currently only works in serial.
'''

__authors__ = 'Samar Khatiwala and Jamie Carr'
__email__ = 'samkat6@gmail.com'

import numpy as np
from math import ceil
import scipy.linalg as linalg
import os.path
from time import time
from pyioutils.hdfio import load, save

np.set_printoptions(precision=4, suppress=True)

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        return self.__dict__.update(d)


class AndersonAccelerationObj(DotDict):
    """Stores the current state of the program. Currently just a dictionary."""

# Negative values indicate some kind of failure
CONVERGED_REASON=(
  {'tol': 1, 
   'itmax': 2,
   'maxNumItersPerRun': 3,
   'wallTime': 4,
   'killFile': 5,
   'modelWallTime': 6,
   'IncorrectParameter': -1,
   'modelFailed': -2,
   'modelPreprocessFailed': -3,
   'modelPostprocessFailed': -4,
   'modelRestartWriteFailed': -5,
   'modelRestartReadFailed': -6,
   'checkpointingFunctionFailed': -7,
   'OtherFailure': -99,
  }
)
CONVERGED_REASON=DotDict(CONVERGED_REASON)


def add_param_to_s(s, paramStr, params):
    if paramStr not in params:
        return
    paramVal = params.get(paramStr)
    # if the default is an int, but we're setting it to a float, then we probably want an int.
    # the behavior can be adjusted by changing the default value types. e.g 0.0 vs 0.
    s[paramStr] = (
        int(paramVal)
        if isinstance(s.get(paramStr), int) and isinstance(paramVal, float)
        else paramVal
    )


def set_default_state(s, x, y, AAparams, histParams, doPrint=True):
    """
    Program state is stored in state dictionary. If restarting, load it, else, initialize.
    """

    if not isinstance(x,np.ndarray):
        raise Exception("Error: x must be a NumPy array")
    
    s.fetchOutput = 0

    s.mMax = min(50, len(x))
    s.itmax = 100
    s.atol = 1.0e-10
    s.rtol = 1.0e-10
    s.droptol = 1.0e10
    s.beta = 1.0
    s.AAstart = 0
    s.restartAANormStagnation = 0
    s.restartAANormDiff = 0.01
    s.restartAASuccessiveIters = 2
    s.restartAAPeriodic = 0
    s.applyAAtoY = True
    s.enforcePositive = False
    s.minNegativeValue = -1.e-10

    # return variables
    # variables
    s.x = x
    s.y = y
    s.restartCount = 0
    s.flushAA = 0
    s.mAA = 0
    s.ithist = 0
    s.nfeval = 0
    s.iter0 = 0
    s.res_norm = -1.0
    
    if AAparams is not None:
        AAKeys = [
            "mMax",
            "itmax",
            "atol",
            "rtol",
            "droptol",
            "beta",
            "AAstart",
            "restartAANormStagnation",
            "restartAANormDiff",
            "restartAASuccessiveIters",
            "restartAAPeriodic",
            "applyAAtoY",
            "enforcePositive",
            "minNegativeValue",
        ]

        for AAKey in AAKeys:
            add_param_to_s(s, AAKey, AAparams)

    # %   Set history defaults and overwrite with values from histParams
    # %   nhistfreq = frequency with which to save history
    # %   ncheckpointfreq = frequency with which to write checkpoint (if < 0
    # %                     then no checkpoints are written; this is the default)

    s.nhistfreq = -1
    s.ncheckpointfreq = -1
    s.savexhist = True
    s.saveyhist = True

    for histKey in ["nhistfreq", "ncheckpointfreq","savexhist","saveyhist"]:
        add_param_to_s(s, histKey, histParams)

    s.res_hist = None

    if s.y is None:
        s.applyAAtoY = False

    # Ideally, these would be created with order='F' but HDF5 stores data 
    # in C row-major order and the order is lost when saving/reading restarts
    s.DG = np.zeros((len(x), s.mMax))
    s.R = np.zeros((s.mMax, s.mMax))
    s.Q = np.zeros((len(x), s.mMax))

    s.f_old = 0.0
    s.g_old = 0.0

    if s.applyAAtoY:
        s.DGy = np.zeros((len(y), s.mMax))
        s.gy_old = 0.0

    if s.nhistfreq > 0:
        if histParams.get("nhistmax"):
            if histParams.nhistmax < (ceil(s.itmax / s.nhistfreq) + 5):
                if doPrint: print("ERROR: histParams.nhistmax is too small!")
                return 1
            else:
                s.nhistmax = histParams.nhistmax
        else:
            s.nhistmax = ceil(s.itmax / s.nhistfreq) + 5  # space to store history

        s.aahist = DotDict({})
        s.aahist.iterhist = np.zeros((s.nhistmax, 1))
        s.aahist.nfevalhist = np.zeros((s.nhistmax, 1))
        if s.savexhist:
            s.aahist.xhist = [[] for _ in range(s.nhistmax)]
        if s.y is not None and s.saveyhist:
            s.aahist.yhist = [[] for _ in range(s.nhistmax)]
        s.aahist.vnormshist = [[] for _ in range(s.nhistmax)]
        s.aahist.convdatahist = [[] for _ in range(s.nhistmax)]
    else:
      # To be safe we override values set in histParams
      s.savexhist = False
      s.saveyhist = False

    return 0


def run_checks(s, AAparams, doPrint=True):
    """Run checks before launching anderson acceleration iteration"""
    aaHasChanged = 0
    if s.iter0 > 0:
        if AAparams.get("beta"):
            if s.beta != AAparams.beta:
                if doPrint: print(f"WARNING: beta has changed! Setting new beta value to {AAparams.beta}")
                s.beta = AAparams.beta
        if AAparams.get("itmax"):
            if s.itmax != AAparams.itmax:
                if doPrint: print("WARNING: itmax has changed!")
                if s.itmax > AAparams.itmax:
                    if doPrint: print("ERROR: new itmax < old itmax")
                    return 1
                else:
                    aaHasChanged = 1
                    if doPrint: print("Resizing res_hist")
                    tmp = s.res_hist
                    s.res_hist = np.zeros((int(AAparams.itmax) + 1, tmp.shape[1]))
                    s.res_hist[: s.itmax + 1, :] = tmp
                    s.itmax = AAparams.itmax
                    del tmp
        if s.res_hist.shape[0] != s.itmax + 1:
            if doPrint: print(f"shape={s.res_hist.shape}")
            if doPrint: print(f"itmax: {s.itmax}")
            aaHasChanged = 1
            if doPrint: print("Resizing res_hist")
            tmp = s.res_hist
            s.res_hist = np.zeros((int(s.itmax) + 1, tmp.shape[1]))
            s.res_hist[: tmp.shape[0], :] = tmp
            s.itmax = AAparams.itmax
            del tmp
        if aaHasChanged:
            if s.nhistfreq > 0:
                nnew=ceil(s.itmax / s.nhistfreq) + 5
                if s.nhistmax < (ceil(s.itmax / s.nhistfreq) + 5) or (len(s.aahist.iterhist) != s.nhistmax):
                    if doPrint: print("Resizing aahist")
                    ncurr=len(s.aahist.vnormshist)
                    nappend=nnew-ncurr
                    s.aahist.iterhist=np.append(s.aahist.iterhist,np.zeros(nappend))
                    s.aahist.nfevalhist=np.append(s.aahist.nfevalhist,np.zeros(nappend))
                    listext=[[] for _ in range(nappend)]
                    if s.savexhist:
                        s.aahist.xhist.extend(listext)
                    if s.y is not None and s.saveyhist:
                        s.aahist.yhist.extend(listext)
                    s.aahist.vnormshist.extend(listext)
                    s.aahist.convdatahist.extend(listext)
#                     errmsg = """
#                     ERROR: Option to resize aahist on restart has not been implemented yet. 
#                     If you intend to restart a run with a different itmax then it would be 
#                     best to set histParams.nhistfreq=-1. 
#                     If you really want to restart this run then set s.nhistfreq=-1 (but keep 
#                     in mind that subsequent history will not be stored).
#                     """
#                     raise Exception(errmsg)
    return 0

def andacc(g, x, s, AAparams=DotDict({}), histParams=DotDict({}), doRestart=False, fileSuff="", y=None, debug=False, mpi=None):
    """
    Begin Anderson acceleration on function g, input x with a AndersonAcceleration state s.
    and other parameters.
    """

    # Generally None is used to denote undefined variables rather than an empty list
    # as is convention in MATLAB.

    # andacc always returns a tuple (x, iter_, y, converged), even when there is a failure 
    # of some kind. A positive converged value indicates andacc terminated without an error. 
    # A negative value signals some kind of failure.
    
    # andacc calls the function g with arguments (x, y, fetchOutput, iter_)
    # fetchOutput: This is an integer which takes on different values depending on 
    #              the 'mode' in which andacc is being used (doRestart 1/0 or True/False)
    #  if doRestart=False or 0:
    #   -1: run the model and read and return the output. Valid return arguments from g are 
    #       a tuple (if everything was ok) or an integer error code.
    #  if doRestart=True or 1:
    #    1: fetch model output (this assumes the model has been submitted/run either via a 
    #          previous call to g or separately by the user). Valid return arguments from g 
    #          are a tuple (if everything was ok) or an integer error code.
    #    0: submit/run the model and return without reading the output. Valid return arguments 
    #          from g are None (if everything was ok) or an integer error code.

    doPrint = (mpi is None) or (mpi is not None and mpi.rank==0)

    if not callable(g):
        raise Exception("Error: g must be a callable function")

    # Hack to check whether s can accessed with dot notation. Checking whether s is of type 
    # DotDict fails with restarts since the type of s is hdfio.DotDict.
    try:
        s.abcdefgh
    except:
        raise Exception("Error: s must be a dictionary of type DotDict")

    # Note: DotDict seems to make a copy, which is not ideal. Should this be changed?
    if not isinstance(AAparams,DotDict):
        AAparams = DotDict(AAparams)

    if not isinstance(histParams,DotDict):
        histParams = DotDict(histParams)
        
    # Use None if a value is not defined, not an empty list here.
    if s.x is None:
        rr = set_default_state(s, x, y, AAparams, histParams, doPrint=doPrint)
        if rr != 0:
            if doPrint: print(f"Terminating AndAcc because of a problem in set_default_state")
            iter_ = s.iter0
            converged = CONVERGED_REASON['IncorrectParameter']
            s.converged = converged
            return s.x, iter_, s.y, converged
    else:
        rr = run_checks(s, AAparams, doPrint=doPrint)
        if rr != 0:
            if doPrint: print(f"Terminating AndAcc because of a problem in run_checks")
            iter_ = s.iter0
            converged = CONVERGED_REASON['IncorrectParameter']
            s.converged = converged
            return s.x, iter_, s.y, converged

    end = int(s.itmax + 1)
    iter_ = s.iter0
    
    if mpi is None:
      try:
          from aa4py.linalg import norm, dot, multtranspose
      except ImportError as e:
          if doPrint: print(f"There was a problem importing linear algebra stub routines from aa4py.linalg: {e}")
          if doPrint: print("Terminating AndAcc because of an import error: The linear algebra stub stub routines could not be importted from from aa4py.linalg")
          converged = CONVERGED_REASON['OtherFailure']
          s.converged = converged
          return s.x, iter_, s.y, converged 
    else:
      try:
          from pympiutils import norm, dot, multtranspose
      except ImportError as e:
          if doPrint: print("Terminating AndAcc because of an import error: The pympiutils module is needed for parallel linear algebra. Make sure it is in your path")
          converged = CONVERGED_REASON['OtherFailure']
          s.converged = converged
          return s.x, iter_, s.y, converged 
    
    if mpi is not None: mpi.comm.Barrier()
    
    if (s.iter0 >= end):
      if doPrint: print(f"Terminating AndAcc because maximum number of iterations {s.itmax} has already been reached")
      converged = CONVERGED_REASON['itmax']
      s.converged = converged
      return s.x, iter_, s.y, converged

    if (not isinstance(doRestart,bool)) and (doRestart not in [0, 1]):
        if doPrint: print(f"Terminating AndAcc because doRestart (={doRestart}) must be True/False or 1/0")
        converged = CONVERGED_REASON['IncorrectParameter']
        s.converged = converged
        return s.x, iter_, s.y, converged 

    killFile = AAparams.get("killFile","killaa")

    trimHistory = histParams.get("trimHistory",True)

    if s.y is not None and s.applyAAtoY is False:
      if doPrint: print(f"Warning: Problem has x/y split but AA is switched off for y")

    if s.enforcePositive:
        if mpi is not None:
            if doPrint: print("Terminating AndAcc because MPI is not supported with enforcePositive=True!")
            converged = CONVERGED_REASON['IncorrectParameter']
            s.converged = converged
            return s.x, iter_, s.y, converged 
        try:
            from aa4py.lsqlin import lsqlin, cvxopt_to_numpy_matrix
        except ImportError as e:
            if doPrint: print("Terminating AndAcc because of an import error: The lsqlin module is needed to solve the non-negative least-squares problem. Make sure it is in your path")
            converged = CONVERGED_REASON['OtherFailure']
            s.converged = converged
            return s.x, iter_, s.y, converged 
      
    doDamping = callable(s.beta) or (s.beta > 0 and s.beta != 1)
    if doDamping and s.applyAAtoY:
        converged = CONVERGED_REASON['IncorrectParameter']
        s.converged = converged
        if doPrint: print("Terminating AndAcc because damping with x/y split is not supported!")
        return s.x, iter_, s.y, converged 

    if ("wallTime" in AAparams) and ("startTime" not in AAparams):
        converged = CONVERGED_REASON['IncorrectParameter']
        s.converged = converged
        if doPrint: print("Terminating AndAcc because wallTime is specified but startTime is not given!")
        return s.x, iter_, s.y, converged 
        
    s.endOfRun = 0
    converged = 0

    if mpi is not None: mpi.comm.Barrier()

    # Beginning of iteration.
    for iter_ in range(s.iter0, end):
        if "maxNumItersPerRun" in AAparams:
            if (iter_-s.iter0) == AAparams.maxNumItersPerRun:
                if doPrint: print(f"Terminating AndAcc because maximum number of iterations per run reached at iter {iter_}")
                s.iter0 = iter_
                converged = CONVERGED_REASON['maxNumItersPerRun']
                s.converged = converged
                return s.x, iter_, s.y, converged

        if "wallTime" in AAparams:
            if time() - AAparams.startTime > AAparams.wallTime:
                if doPrint: print(f"Terminating AndAcc because of wall time limit at iter {iter_}")
                s.iter0 = iter_
                converged = CONVERGED_REASON['wallTime']
                s.converged = converged
                return s.x, iter_, s.y, converged

        if os.path.isfile(killFile):
            if doPrint: print(f"Terminating AndAcc because of external kill signal at iter {iter_}")
            s.iter0 = iter_
            converged = CONVERGED_REASON['killFile']
            s.converged = converged
            return s.x, iter_, s.y, converged

        if not doRestart:
            # Apply g and compute the current residual norm.
            args = (s.x, s.y, -1, iter_)
            if mpi is not None: mpi.comm.Barrier()
            rr = g(*args)
            if isinstance(rr,int):
                key = next((k for k, v in CONVERGED_REASON.items() if v == rr), None)
                if key is not None:
                    if doPrint: print(f"Terminating AndAcc because of converged reason {key} at iter {iter_}")
                    converged = CONVERGED_REASON[key]
                else:
                    if doPrint: print(f"Terminating AndAcc because of model exit code {rr} at iter {iter_}")
                    converged = CONVERGED_REASON['OtherFailure']
                s.iter0 = iter_
                s.converged = converged
                return s.x, iter_, s.y, converged
            gval, gyval, vnorms, externalconvergence = rr
            if mpi is not None: mpi.comm.Barrier()
        else:
            if s.fetchOutput:
                args = (s.x, s.y, 1, iter_)
                if mpi is not None: mpi.comm.Barrier()
                rr = g(*args)
                if isinstance(rr,int):
                    key = next((k for k, v in CONVERGED_REASON.items() if v == rr), None)
                    if key is not None:
                        if doPrint: print(f"Terminating AndAcc because of converged reason {key} at iter {iter_}")
                        converged = CONVERGED_REASON[key]
                    else:
                        if doPrint: print(f"Terminating AndAcc because of model exit code {rr} at iter {iter_}")
                        converged = CONVERGED_REASON['OtherFailure']
                    s.iter0 = iter_
                    s.converged = converged
                    return s.x, iter_, s.y, converged
                gval, gyval, vnorms, externalconvergence = rr
                if mpi is not None: mpi.comm.Barrier()
                s.fetchOutput = 0
            else:
                # We save the restart first before submitting the next job.
                # fetchOutput is set to 1 for the subsequent (restarted) iteration but note that
                # the function g is called below with the fetchOutput argument set to 0
                s.fetchOutput = 1
                s.iter0 = iter_
                if mpi is not None: mpi.comm.Barrier()
                rr = g(s.x, s.y, 0, iter_)  # Submit Job
                # We always expect None or an integer back (if there was an error)
                if rr==None: # all okay
                  break
                else:
                  if isinstance(rr,int):
                      key = next((k for k, v in CONVERGED_REASON.items() if v == rr), None)
                      if key is not None:
                          if doPrint: print(f"Terminating AndAcc because of converged reason {key} at iter {iter_}")
                          converged = CONVERGED_REASON[key]
                      else:
                          if doPrint: print(f"Terminating AndAcc because of model exit code {rr} at iter {iter_}")
                          converged = CONVERGED_REASON['OtherFailure']
                      s.iter0 = iter_
                      s.converged = converged
                      return s.x, iter_, s.y, converged
                  else:
                      if doPrint: print("Terminating AndAcc because of unknown return argument in call to g; epected None or an integer!")
                      converged = CONVERGED_REASON['OtherFailure']
                      s.iter0 = iter_
                      s.converged = converged
                      return s.x, iter_, s.y, converged

        s.nfeval = s.nfeval + 1
        fval = gval - s.x

        if mpi is not None: mpi.comm.Barrier()
          
        res_norm = norm(fval,mpi=mpi)
        s.res_norm = res_norm

        # We calculate the norm for monitoring even if applyAAtoY is False
        if s.y is not None:
            fyval = gyval - s.y
            res_normy = norm(fyval,mpi=mpi)

        if doPrint:
            if s.y is not None:
                print(f"Iter {iter_}: {res_norm}, {res_normy}")
            else:
                print(f"Iter {iter_}: {res_norm}")

        # Ensure vnorms is always an array.
        if vnorms is not None:
            isArr = isinstance(vnorms, np.ndarray)
            vnorms = vnorms if isArr else np.array([vnorms])

        if iter_ == 0:
            if vnorms is not None:
              s.res_hist = np.zeros((int(s.itmax) + 1, 3 + vnorms.size))
            else:
              s.res_hist = np.zeros((int(s.itmax) + 1, 3))

            s.tol = max(s.atol, s.rtol * res_norm)

            if s.nhistfreq > 0:
                # NOTE: we increment ithist after populating aahist because of 0-based indexing
                s.aahist.iterhist[s.ithist] = iter_
                s.aahist.nfevalhist[s.ithist] = s.nfeval
                if s.savexhist:
                    s.aahist.xhist[s.ithist] = s.x
                if s.y is not None and s.saveyhist:
                    s.aahist.yhist[s.ithist] = s.y
                if vnorms is not None:
                    s.aahist.vnormshist[s.ithist] = np.append(vnorms, res_norm)
                else:
                    s.aahist.vnormshist[s.ithist] = np.array([res_norm])
                if convdata := externalconvergence.get("convdata"):
                    s.aahist.convdatahist[s.ithist] = convdata
                s.ithist = s.ithist + 1

        row = np.array([iter_, s.nfeval])
        if vnorms is not None:
            row = np.append(row, vnorms)
        s.res_hist[iter_, :] = np.append(row, np.array([res_norm]))

        # Check external convergence criteria
        if extconverged := externalconvergence.get("converged"):
            if extconverged == 0:
                if res_norm <= s.tol:
                    msg = "Overwriting default convergence reason as external convergence criteria have not been met "
                    if doPrint: print(msg, res_norm, extconverged)

            elif extconverged > 0:
                msg = "External convergence criteria have been met "
                if doPrint: print(msg, res_norm, extconverged)
                s.endOfRun = 1
                converged = extconverged
                break

            else:
                msg = "Divergence indicated by external convergence check "
                if doPrint: print(msg, res_norm, extconverged)
                s.endOfRun = 1
                converged = extconverged
                break
        # Use default convergence criterion
        else:
            if res_norm <= s.tol:
                if doPrint: print("Terminate with residual norm = ", res_norm)
                s.endOfRun = 1
                converged = CONVERGED_REASON['tol']
                break

        if s.mMax == 0 or iter_ < s.AAstart:
            # Without acceleration, update x <- g(x) to obtain the next approximate solution.
            s.x = gval
            if s.y is not None:
                s.y = gyval

        else:
            # Apply Anderson Acceleration

            if iter_ > 0 and s.mAA > 0:
                if os.path.isfile("flushaa"):
                    msg = f"Restarting AA at iter {iter_} mAA={s.mAA} due to external signal"
                    if doPrint: print(msg)
                    s.flushAA = 1
                    os.remove("flushaa")
                elif flushAA := externalconvergence.get("flushAA"):
                    if flushAA in [0, 1]:
                        s.flushAA = flushAA
                        if s.flushAA:
                            msg = f"Restart of AA is being triggered at iter {iter_}, mAA={s.mAA} due to external convergence signal"
                            if doPrint: print(msg)
                elif s.restartAAPeriodic > 0 and iter_ % s.restartAAPeriodic == 0:
                    s.flushAA = 1
                    msg = f"Restart of AA is being triggered at iter {iter_} due to periodic restart condition"
                    if doPrint: print(msg)
                elif s.restartAANormStagnation:
                    s.itc = iter_ + 1
                    s.itp = max(1, s.itc - 1)
                    # Not sure
                    rres = s.res_hist[s.itp : s.itc + 1, -1]
                    rd = rres[2] - rres[1]
                    if rd > 0.0:
                        s.restartCount = s.restartCount + 1
                        if s.restartCount >= s.restartAASuccessiveIters:
                            s.flushAA = 1
                            if doPrint:
                              print(
                                  f"Restart of AA is being triggered at iter {iter_}, mAA = {s.mAA} due to stagnation:",
                                  rres(1),
                                  rres(2),
                                  rd,
                              )
                    else:
                        s.restartCount = 0

            if s.flushAA:
                s.AAstart = iter_
                s.mAA = 0
                s.DG = np.zeros((len(s.x), s.mMax))
                s.R = np.zeros((s.mMax, s.mMax))
                s.Q = np.zeros((len(s.x), s.mMax))
                s.restartCount = 0
                s.flushAA = 0
                if s.applyAAtoY:
                    s.DGy = np.zeros((len(s.y), s.mMax))

            # Update the df vector and DG array

            if iter_ > s.AAstart:
                df = fval - s.f_old
                if s.mAA < s.mMax:
                    s.DG[:, s.mAA] = gval - s.g_old
                    if s.applyAAtoY:
                        s.DGy[:, s.mAA] = gyval - s.gy_old
                else:
                    # Roll all columns back one, then add newest dg on the last.
#                     s.DG = np.roll(s.DG, -1, axis=1)
                    for j in range(s.mAA-1):
                      s.DG[:,j] = s.DG[:,j+1]
                    s.DG[:, -1] = gval - s.g_old
                    if s.applyAAtoY:
#                         s.DGy = np.roll(s.DGy, -1, axis=1)
                        for j in range(s.mAA-1):
                          s.DGy[:,j] = s.DGy[:,j+1]
                        s.DGy[:, -1] = gyval - s.gy_old

                s.mAA = s.mAA + 1
              
            s.f_old = fval
            s.g_old = gval
            if s.applyAAtoY:
                s.gy_old = gyval

            if s.mAA == 0:
                # TODO sure we want this?
                if not s.flushAA:
                    s.x = gval
                    if s.y is not None:
                        s.y = gyval

            else:
                if mpi is not None: mpi.comm.Barrier()
                if s.mAA == 1:
                    # If mAA == 1, form the initial QR decomposition.
                    dfn=norm(df,mpi=mpi)
                    s.R[0, 0] = dfn
                    s.Q[:, 0] = df/dfn
                else:
                    if s.mAA > s.mMax:
                        _qrdelete(s.Q, s.R, 0)
                        s.mAA = s.mAA - 1
                    if mpi is not None: mpi.comm.Barrier()
                    
#                   Update the QR decomposition to incorporate the new column
                    for j in range(0, s.mAA - 1):
#                         a = s.Q[:, j] # a = np.matrix.transpose(s.Q[:, j : j + 1])
                        s.R[j, s.mAA - 1] = dot(df, s.Q[:, j],mpi=mpi)
                        df = df - s.R[j, s.mAA - 1] * s.Q[:, j]
#                         if mpi is not None: mpi.comm.Barrier()

                    dfn=norm(df,mpi=mpi)
                    s.R[s.mAA - 1, s.mAA - 1] = dfn
                    s.Q[:, s.mAA - 1] = df/dfn

                if s.droptol > 0:
                    # view R are the right size
                    condDF = np.linalg.cond(s.R[: s.mAA, : s.mAA], p=2)
                    while condDF > s.droptol and s.mAA > 2:
                        if doPrint: print(f" cond(D) = {condDF}, reducing mAA to {s.mAA-1} \n")
                        _qrdelete(s.Q, s.R, 0)
#                       Note: we don't change the size of DG as that 
#                       will throw an error in the next iteration
#                         s.DG = np.roll(s.DG,-1,axis=1)
                        for j in range(s.mAA-1):
                          s.DG[:,j] = s.DG[:,j+1]
                        if s.applyAAtoY:
#                             s.DGy = np.roll(s.DGy,-1,axis=1)
                            for j in range(s.mAA-1):
                              s.DGy[:,j] = s.DGy[:,j+1]

                        s.mAA = s.mAA - 1
                        condDF = np.linalg.cond(s.R[: s.mAA, : s.mAA], p=2)

                if mpi is not None: mpi.comm.Barrier()                   

                # b is temporary variable for cleaner code
                # We only use the first mAA rows of the product below so really only need to 
                # pass the first mAA columns of Q
                b = multtranspose(s.Q[:,:s.mAA],fval,mpi=mpi)
                # view R at right size
                if mpi is None:
                  gamma = linalg.solve(s.R[: s.mAA, : s.mAA], b)
                else:  
                  if mpi.rank==0:
                    gamma_in = linalg.solve(s.R[: s.mAA, : s.mAA], b)
                  else:
                    gamma_in=None
                  gamma = mpi.comm.bcast(gamma_in, root = 0)
                  mpi.comm.Barrier()

#               Note: all of these are (I think) local calculations
                s.x = gval - np.matmul(s.DG[:, : s.mAA], gamma)
                if s.y is not None:
                    if s.applyAAtoY:
                        s.y = gyval - np.matmul(s.DGy[:, : s.mAA], gamma)
                    else:
                        s.y = gyval

                if s.enforcePositive:
                    s.x[:]=-1.e-9
                    if np.any(s.x < s.minNegativeValue):
                        if debug: print(f"Minimum x is {np.min(s.x)}; enforcing positive constraint")
                        # Solve least squares problem enforcing nonnegativity constraint on solution x (the constraint is NOT on gamma)
                        gamma=lsqlin(s.R[: s.mAA, : s.mAA], b, A=s.DG[:, : s.mAA],b=gval, opts={'show_progress': False})['x']
                        # gamma=cvxopt_to_numpy_matrix(gamma)
                        gamma=np.array(gamma) # this is actually what is done in cvxopt_to_numpy_matrix but without the squeeze()
                        gamma.shape = (gamma.shape[0],)
                        # Note: all of these are (I think) local calculations
                        s.x = gval - np.matmul(s.DG[:, : s.mAA], gamma)
                        if s.y is not None:
                            if s.applyAAtoY:
                                s.y = gyval - np.matmul(s.DGy[:, : s.mAA], gamma)

                if doDamping:
#                   I think this should be okay with MPI but I haven't tested it
                    args = (s.Q[:, : s.mAA], np.matmul(s.R[: s.mAA, : s.mAA], gamma))
                    dx = fval - np.matmul(*args)
                    if callable(s.beta):
                        s.x = s.x - (1 - s.beta(iter_)) * dx
#                         if s.applyAAtoY:
                          # Damping with x/y split is not currently supported
                    else:
                        if s.beta > 0 and s.beta != 1:
                            s.x = s.x - (1 - s.beta) * dx
#                             if s.applyAAtoY:
                              # Damping with x/y split is not currently supported

        if s.nhistfreq > 0 and (iter_ > 0 and iter_ % s.nhistfreq == 0):
            # NOTE: we increment ithist after populating aahist because of 0-based indexing
            s.aahist.iterhist[s.ithist] = iter_
            s.aahist.nfevalhist[s.ithist] = s.nfeval
            if s.savexhist:
                s.aahist.xhist[s.ithist] = s.x
            if s.y is not None and s.saveyhist:
                s.aahist.yhist[s.ithist] = s.y
            if vnorms is not None:
                s.aahist.vnormshist[s.ithist] = np.append(vnorms, res_norm)
            else:
                s.aahist.vnormshist[s.ithist] = np.array([res_norm])
            if convdata := externalconvergence.get("convdata"):
                s.aahist.convdatahist[s.ithist] = convdata

            s.ithist = s.ithist + 1

        if mpi is not None: mpi.comm.Barrier()

        if s.ncheckpointfreq > 0 and (iter_ > 0 and iter_ % s.ncheckpointfreq == 0):
            fn = f"aa_checkpoint_{iter_}{fileSuff}.h5"
            # Hack to make this restart usable for actual restarts
            # Save iter0 before setting it to current iter
            iter0_ = s.iter0
            # Careful here: we're at the end of the loop so we want to restart
            # with the next iteration number
            s.iter0 = iter_ + 1
            saveaarestart(s,fn,mpi=mpi)
            # Restore iter0
            s.iter0 = iter0_

        if histParams.get("checkpointingFunction"):
            rr = histParams.checkpointingFunction(iter_, s, mpi=mpi)
            if rr != 0:
              if doPrint: print(f"Terminating AndAcc because of converged reason checkpointingFunctionFailed at iter {iter_}")
              converged = CONVERGED_REASON['checkpointingFunctionFailed']
              s.iter0 = iter_
              s.converged = converged
              return s.x, iter_, s.y, converged
            
        if iter_ == s.itmax:
            s.endOfRun = 1
            converged = CONVERGED_REASON['itmax']

#       End of loop
        if mpi is not None: mpi.comm.Barrier()

    if mpi is not None: mpi.comm.Barrier() # this deals with the break statements above
    
    if s.endOfRun:
        if res_norm > s.tol and iter_ == s.itmax:
            if doPrint:
              print(f"Terminate after itmax = {s.itmax} iterations.")
              print(f"Residual norm = {res_norm}")

        # store final iteration
        if s.nhistfreq > 0:
            # NOTE: we increment ithist after populating aahist because of 0-based indexing
            s.aahist.iterhist[s.ithist] = iter_
            s.aahist.nfevalhist[s.ithist] = s.nfeval
            if s.savexhist:
                s.aahist.xhist[s.ithist] = s.x
            if s.y is not None and s.saveyhist:
                s.aahist.yhist[s.ithist] = s.y
            if vnorms is not None:
                s.aahist.vnormshist[s.ithist] = np.append(vnorms, res_norm)
            else:
                s.aahist.vnormshist[s.ithist] = np.array([res_norm])
            if convdata := externalconvergence.get("convdata"):
                s.aahist.convdatahist[s.ithist] = convdata

            s.ithist = s.ithist + 1

        # This is to allow continuation of a previous solve
        s.iter0 = iter_ + 1

        if s.ncheckpointfreq > 0:
            if mpi is not None: mpi.comm.Barrier()
            fn = "aa_checkpoint_final.h5"
            saveaarestart(s,fn,mpi=mpi)

        if trimHistory:
            s.res_hist = s.res_hist[: iter_ + 1, :]
            if s.nhistfreq > 0:
                s.aahist.iterhist = s.aahist.iterhist[: s.ithist]
                s.aahist.nfevalhist = s.aahist.nfevalhist[: s.ithist]
                if s.savexhist:
                    s.aahist.xhist = s.aahist.xhist[: s.ithist]
                if s.y is not None and s.saveyhist:
                  s.aahist.yhist = s.aahist.yhist[: s.ithist]
                s.aahist.vnormshist = s.aahist.vnormshist[: s.ithist]
                s.aahist.convdatahist = s.aahist.convdatahist[: s.ithist]

        xsol = s.x
        ysol = s.y
    else:
        xsol = np.array([])
        ysol = np.array([])

    s.converged = converged

    if mpi is not None: mpi.comm.Barrier()
    
    return xsol, iter_, ysol, converged

# The following are translations of the Fortran routines of L. Reichel and W. B. Gragg (1990), 
# Algorithm 686: FORTRAN subroutines for updating the QR decomposition, ACM Transactions 
# on Mathematical Software (TOMS), Vol. 16 (4), https://dl.acm.org/doi/10.1145/98267.98291. 
# All these routines modify Q and R in-place. Only _qrdelete should be called by the end-user. 
# Currently only deletion of a single column K (0-based indexing) is supported. It can be 
# used in serial or parallel (this assumes that Q and R are NumPy arrays with Q distributed 
# row-wise, i.e., each process holds a contiguous subset of the rows of Q, and that the full 
# (much smaller) R matrix is identical and available on all processes). 
# The sign of the returned Q and R may be different from that computed by scipy.linalg.qr_delete

def _qrdelete(Q,R,K):

  M=Q.shape[0]
  N=R.shape[1]
  LDR=N

# C     UPDATE R AND Q.
  N1=N-1
  for L in range(K,N1): # K->N-2
    R[L,L+1],R[L+1,L+1],CS,SN = _DADFGR(R[L,L+1],R[L+1,L+1])
    _DAPLRR(CS,SN,R,LDR,N,L,L+2)
    _DAPLRC(CS,SN,M,Q[:,L],Q[:,L+1])

#   C     MOVE COLUMNS OF R.
  for L in range(K,N1):
    R[:,L]=R[:,L+1]

#   C     SET NTH COLUMNS TO ZERO.
  R[:,-1]=0.0
  Q[:,-1]=0.0

def _DADFGR(X,Y):
# C     DADFGR DEFINES AND APPLIES THE SYMMETRIC GIVENS REFLECTOR
  from math import copysign
  if (Y==0.0):
    C=1.0
    S=0.0
#     return C, S
  else:
# C
    MU=max(np.abs(X),np.abs(Y))
    X1=X/MU
    Y1=Y/MU
    T=copysign(MU*np.sqrt(X1*X1+Y1*Y1),X)
    C=X/T
    S=Y/T
    X=T
    Y=0.0
  return X,Y,C,S

# C DAPLRR(CS,SN,R,LDR,N,L,L+2) L: K->N-1, L+2: K+2->N+1
def _DAPLRR(C,S,A,LDA,N,IRW,L): #IRW: K->N-1, L: K+2->N+1
# C     DAPLRR APPLIES A GIVENS REFLECTOR G DEFINED BY C AND S, WHERE
# C     G HAS BEEN DETERMINED BY SUBROUTINE DADFGR. WHEN CALLED DAPLRR
# C     REPLACES THE TWO ROW MATRIX ( A(IRW,J), A(IRW+1,J) )', L<=J<=N,
# C     BY G * ( A(IRW,J), A(IRW+1,J) )', L<=J<=N.

  for I in range(L,N):
    TMP=C*A[IRW,I]+S*A[IRW+1,I]
    A[IRW+1,I]=S*A[IRW,I]-C*A[IRW+1,I]
    A[IRW,I]=TMP

#   DAPLRC(CS,SN,M,Q(1,L),Q(1,L+1)) L: K->N-1
def _DAPLRC(C,S,N,X,Y):
# C     DAPLRC APPLIES A GIVENS REFLECTOR G DEFINED BY C AND S, WHERE
# C     G HAS BEEN DETERMINED BY SUBROUTINE DADFGR. WHEN CALLED DAPLRC
# C     REPLACES THE TWO COLUMN MATRIX ( X, Y )  BY ( X, Y ) * G.

    TMP=C*X[:]+S*Y[:]
    Y[:]=S*X[:]-C*Y[:]
    X[:]=TMP

def saveaarestart(state, filename, mpi=None, debug=False):
  mpix=None
  mpiy=None
  if mpi is not None:
    mpix=DotDict()
    mpix.comm=mpi.comm
    mpix.rank=mpi.rank
    mpix.layout=mpi.layout
    mpix.globalVars=['R','res_hist','iterhist','nfevalhist','vnormshist']
    if state.y is not None:
      mpiy=DotDict()
      mpiy.comm=mpi.comm
      mpiy.rank=mpi.rank
      mpiy.layout=mpi.layouty
  save(state,filename,groupname="x",skipvars=['y','DGy','gy_old','yhist'],mpi=mpix,debug=debug)
  if state.y is not None:
    save(state,filename,groupname="y",savevars=['y','DGy','gy_old','yhist'],skipvars=['convdatahist', 'vnormshist', 'xhist'],append=True,mpi=mpiy,debug=debug)
  
def loadaarestart(filename, mpi=None, debug=False):
  mpix=None
  mpiy=None
  if mpi is not None:
    mpix=DotDict()
    mpix.comm=mpi.comm
    mpix.rank=mpi.rank
    mpix.layout=mpi.layout
    mpix.globalVars=['R','res_hist','iterhist','nfevalhist','vnormshist']
    if mpi.layouty is not None:
      mpiy=DotDict()
      mpiy.comm=mpi.comm
      mpiy.rank=mpi.rank
      mpiy.layout=mpi.layouty
  s=load(filename,groupname="x",mpi=mpix,debug=debug)
  if s is None:
    print(f"Group x not found in {filename}. Make sure it is a valid restart file")
  sy=load(filename,groupname="y",mpi=mpiy,debug=debug)
  if sy is not None:
    s.y=sy.y
    s.DGy=sy.DGy
    s.gy_old=sy.gy_old
    if s.saveyhist:
      s.aahist.yhist=sy.aahist.yhist
      
  return s