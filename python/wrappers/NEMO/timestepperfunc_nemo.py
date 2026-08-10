from aa4py import DotDict, CONVERGED_REASON
from aa4py.linalg import norm
from pyioutils.hdfio import load, save
import numpy as np
import os
from nemo_andacc_utils import writexyloc2NEMOGlobalRestart, writexy2NEMORestart, readFromNEMOGlobalRestart2xyloc, readFromNEMORestart2xy

def timestepperfunc_nemo(x, y, fetchOutput, iter_, data, jobId, doSubmit, mpi=None):

#   See andacc.py for a description of what the different values of fetchOutput mean

    # Any failure in the first fetchOutput block causes this function to terminate 
    # and return an integer error code. If there is no error and the second fetchOutput 
    # block below is NOT executed, this function returns a None. If it is executed, 
    # it returns a tuple (gx, gy, vnorms, externalConvergence)

    # The following block is always executed and thus requires a valid data argument
    gy=None

    useMPI = (mpi is not None)
    doPrint = (not useMPI) or (useMPI and mpi.rank==0)

    if isinstance(data,list):
      useGlobalFile=data[0].useGlobalFile
      calcNorms=False
      restartFileStart=[d.restartFileStart for d in data]
      restartFileEnd=[d.restartFileEnd for d in data]
      writeTRB=data[0].writeTRB
      run_model_cmd=data[0].run_model_cmd
      if "preprocess_cmd" in data[0]:
        preprocess_cmd=data[0].preprocess_cmd
      else:
        preprocess_cmd=None
      if "postprocess_cmd" in data[0]:
        postprocess_cmd=data[0].postprocess_cmd
      else:
        postprocess_cmd=None
    else:
      useGlobalFile=data.useGlobalFile
      calcNorms = ("calcNorms" in data and data.calcNorms)
      # Serial (+global or + tile)
      # MPI+global & rank 0
      # MPI+tile
      calcNorms = calcNorms and ((not useMPI) or ((useMPI and useGlobalFile) and mpi.rank==0) or (useMPI and (not useGlobalFile)))
      restartFileStart=data.restartFileStart
      restartFileEnd=data.restartFileEnd
      writeTRB=data.writeTRB
      run_model_cmd=data.run_model_cmd
      if "preprocess_cmd" in data:
        preprocess_cmd=data.preprocess_cmd
      else:
        preprocess_cmd=None
      if "postprocess_cmd" in data:
        postprocess_cmd=data.postprocess_cmd
      else:
        postprocess_cmd=None

    if x is not None:
      iniFile=f"trini_{jobId:04d}"
      
    # First fetchOutput block
    if (fetchOutput==-1) or (fetchOutput==0):
      if preprocess_cmd is not None:
        if isinstance(preprocess_cmd,str):
          os.system(preprocess_cmd)
        elif callable(preprocess_cmd):
          rr=preprocess_cmd(iter_,mpi)
          if rr != 0:
            if doPrint: print(f"Non-zero exit code {rr} returned by preprocess function. Returning to andacc")
            return rr
        else:
          if doPrint: print("ERROR: preprocess_cmd is neither callable nor a string. Returning to andacc")
          rr = CONVERGED_REASON['IncorrectParameter']
          return rr

      if x is not None:
        if useGlobalFile and useMPI: # MPI+global
          # This is a collective call but only root writes. If there is a failure in writing 
          # an exception will only be triggered on root. I think the (convoluted) logic below 
          # handles this correctly but I'm not 100% sure!
          rr = None
          try:
            writexyloc2NEMOGlobalRestart(x, y, data, restartFileStart, mpi, vbuf=data.vbuf, writeTRB=writeTRB)
            rr = 0
          except Exception as e:
            if doPrint: print(f"Writing to restart file failed with exception {e}. Returning to andacc")
            rr = CONVERGED_REASON['modelRestartWriteFailed']
          returncode = mpi.comm.bcast(rr, root = 0)
        else: # All other cases (MPI+tiles, serial+tiles, serial+global)
          try:
            writexy2NEMORestart(x, y, data, restartFileStart, singleFile=(useGlobalFile or useMPI), vbuf=data.vbuf, writeTRB=writeTRB)
            returncode = 0
          except Exception as e:
            if doPrint: print(f"Writing to restart file failed with exception {e}. Returning to andacc")
            returncode = CONVERGED_REASON['modelRestartWriteFailed']
          if useMPI: mpi.comm.barrier()
        if returncode != 0:
          return returncode

        # We save this to compute norms and convergence later
        # Note: For MPI+global only rank 0 (which holds the full tracer field) writes; for MPI+tile each processor 
        #       writes its local piece to its own file. See above for situations in which calcNorms=T.
        if calcNorms:
          if fetchOutput==0:
            save({"v": data.vbuf}, iniFile)
          else:
            # If we're here then fetchOutput=-1 and we keep a copy of v for use below
            v = data.vbuf.copy()

      if doSubmit:
        # Run the model (if using MPI, only the root rank does this)
        if doPrint: print("Running model ...")
        if isinstance(run_model_cmd,str):
          os.system(run_model_cmd)
        elif callable(run_model_cmd):
          rr=run_model_cmd(mpi)
          if rr != 0:
            if doPrint: print(f"Non-zero exit code {rr} returned by model. Returning to andacc")
            return rr
        else:
          if doPrint: print("ERROR: run_model_cmd is neither callable nor a string. Returning to andacc")
          rr = CONVERGED_REASON['IncorrectParameter']
          return rr

        if mpi is not None: mpi.comm.Barrier()
        
        # SPK Would it better to move this to the fetchOutput block below?
        if postprocess_cmd is not None:
          if isinstance(postprocess_cmd,str):
            os.system(postprocess_cmd) # SPK why is this here and not in the block below?
          elif callable(postprocess_cmd):
            rr=postprocess_cmd(iter_,mpi)
            if rr != 0:
              if doPrint: print(f"Non-zero exit code {rr} returned by postprocess function. Returning to andacc")
              return rr
          else:
            if doPrint: print("ERROR: postprocess_cmd is neither callable nor a string. Returning to andacc")
            rr = CONVERGED_REASON['IncorrectParameter']
            return rr

          if mpi is not None: mpi.comm.Barrier()

    # Second fetchOutput block
    if (fetchOutput==-1) or (fetchOutput==1):
      if doPrint: print("Finished running model")

      # read gv from restart
      if useGlobalFile and useMPI:
        # This is a collective call but only root reads. If there is a failure in reading 
        # an exception will only be triggered on root. I think the (convoluted) logic below 
        # handles this correctly but I'm not 100% sure!
        rr = None
        try:
          gx,gy=readFromNEMOGlobalRestart2xyloc(data, restartFileEnd, mpi, vbuf=data.vbuf)
          rr = 0
        except Exception as e:
          if doPrint: print(f"Reading from restart file failed with exception {e}. Returning to andacc")
          rr = CONVERGED_REASON['modelRestartReadFailed']
        returncode = mpi.comm.bcast(rr, root = 0)
      else:
        try:
          gx,gy=readFromNEMORestart2xy(data, restartFileEnd, singleFile=(useGlobalFile or useMPI), vbuf=data.vbuf)
          returncode = 0
        except Exception as e:
          if doPrint: print(f"Reading from restart file failed with exception {e}. Returning to andacc")
          returncode = CONVERGED_REASON['modelRestartReadFailed']
        if useMPI: mpi.comm.barrier()
      if returncode != 0:
        return returncode

      if calcNorms:
        if fetchOutput==1:
          v = load(iniFile).v
        gv = data.vbuf
        # calculate norms
        # This is stupid. gv=data.vbuf is scaled in place to x units when reading 
        # file above and mapping to x. We need to scale it back to model units 
        # when computing the norm below
        vnorms = np.zeros(data.numTracers)
        for itr in range(data.numTracers):
          vnorms[itr]=norm(gv[:,itr]*data.XTovScaleFac[itr,itr]-v[:,itr])

        if doPrint:
          with np.printoptions(suppress=False):
            print(vnorms)
      else:
        vnorms = None

      if useMPI:
        vnorms_out = mpi.comm.bcast(vnorms, root = 0)
        mpi.comm.Barrier()
      else:
        vnorms_out = vnorms

      externalConvergence = DotDict()
#       if "externalConvergence" in data:
#         externalConvergence = data.externalConvergence(iter_,tmparr)
#         if not isinstance(externalConvergence,dict):
#           raise Exception(f"ERROR: function {data.externalConvergence} must return a dictionary!")

      return gx, gy, vnorms_out, externalConvergence
    else:
      return None

    
