from math import ceil, floor
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
from nemo_andacc_utils import readFromNEMOGlobalRestart2xyloc, readFromNEMORestart2xy, xy2nemo, xyloc2nemoglob
from time import time, sleep
from timestepperfunc_nemo import timestepperfunc_nemo
import netCDF4 as nc
import shlex, subprocess

startTime = time()

modelRunTime=470 # [seconds] Set this to just a bit (~10-20 s) longer than the time it takes to run the model for a year
wallTime=1*60*60

# Uncomment one of these blocks
# Important: This script doesn't currently work on grouped tiles

# MPI+tiles:
# useMPI=True
# from mpi4py import MPI
# from pympiutils import generate_mpi_layout, xloc2glob
# useGlobalFile=False # do NOT change
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
# datatiles=load("data_for_nemo_pisces_tiled.h5").data
# datainactive=load("data_for_nemo_pisces_tiled.h5").datainactive
# dataglob=load("data_for_nemo_pisces.h5")
# comm.Barrier()
# if size != len(datatiles):
#   raise Exception(f"Number of ranks ({size}) must be equal to the number of tiles or groups of tiles ({len(datatiles)})")
# data=datatiles[rank]
# vbuf=np.zeros(shape=(data.nb,data.numTracers))
# mpiaa=DotDict()
# mpiaa.MPI=MPI
# mpiaa.comm=comm
# mpiaa.rank=rank
# mpiaa.size=size
# nloc=data.x2vIb.size
# mpiaa.layout=generate_mpi_layout(nloc=nloc,MPI=mpiaa.MPI)
# nloc=data.y2vIb.size
# mpiaa.layouty=generate_mpi_layout(nloc=nloc,MPI=mpiaa.MPI)
# restartFileStart = f"restart_trc_{data.iFile:04d}.nc"
# restartFileEnd = f"ORCA2_00005840_restart_trc_{data.iFile:04d}.nc"
# iJob=data.iFile # Tag for temporary files
# doPrint=(rank==0)
# diagnosticsFile = "ORCA2_1y_00010101_00011231_bioscalar.nc"
# comm.Barrier()

# MPI+global:
# useMPI=True
# from mpi4py import MPI
# from pympiutils import generate_mpi_layout, xloc2glob
# useGlobalFile=True
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
# data=load("data_for_nemo_pisces.h5")
# if rank==0:
#   vbuf=np.zeros(shape=(data.nb,data.numTracers))
# else:
#   vbuf=None
# comm.Barrier()
# mpiaa=DotDict()
# mpiaa.MPI=MPI
# mpiaa.comm=comm
# mpiaa.rank=rank
# mpiaa.size=size
# n=data.x2vIb.size
# mpiaa.layout=generate_mpi_layout(nglob=n, nprocs=size, MPI=mpiaa.MPI)
# n=data.y2vIb.size
# mpiaa.layouty=generate_mpi_layout(nglob=n, nprocs=size, MPI=mpiaa.MPI)
# restartFileStart = f"restart_trc.nc"
# restartFileEnd = f"ORCA2_00005840_restart_trc.nc"
# iJob=rank # Tag for temporary files
# doPrint=(rank==0)
# diagnosticsFile = "ORCA2_1y_00010101_00011231_bioscalar.nc"
# comm.Barrier()

# Serial+tiles:
useMPI=False
useGlobalFile=False
mpiaa=None
data=load("data_for_nemo_pisces.h5")
vbuf=np.zeros(shape=(data.nb,data.numTracers))
# For distributed restart files give the file prefix only below (everything before ".nc")
numTiles=144 # If useGlobalFile=False, substitute with the actual number of tiles
restartFileStart = "restart_trc"
restartFileEnd = "ORCA2_00005840_restart_trc"
iJob=1 # Tag for temporary files
doPrint=True
diagnosticsFile = "ORCA2_1y_00010101_00011231_bioscalar.nc"

# Serial+global:
# useMPI=False
# useGlobalFile=True
# mpiaa=None
# data=load("data_for_nemo_pisces.h5")
# vbuf=np.zeros(shape=(data.nb,data.numTracers))
# restartFileStart = f"restart_trc.nc"
# restartFileEnd = f"ORCA2_00005840_restart_trc.nc"
# iJob=1 # Tag for temporary files
# doPrint=True
# diagnosticsFile = "ORCA2_1y_00010101_00011231_bioscalar.nc"

data.vbuf=vbuf

AAparams = DotDict()
AAparams.mMax = 50
AAparams.itmax = 500
AAparams.beta = 1.0 
AAparams.restartAANormStagnation = 0
AAparams.restartAANormDiff = 0.05
AAparams.restartAASuccessiveIters = 3
AAparams.restartAAPeriodic = 0
AAparams.startTime = startTime
AAparams.wallTime=wallTime
# Preconditioned AA
# AAparams.mPAA=50
# AAparams.PAAfreq=60
# AAparams.PAAstart=10

histParams = DotDict()
histParams.nhistfreq = -1
# histParams.nhistmax = ceil(AAparams.itmax / histParams.nhistfreq) + 5

doRestart = 1
doSubmit=0

# Note: We swap the restart files as a preprocessing step. Thus, the ancillary data that 
# are in the End restart file are automatically available for the next model run. 

def runModel(mpi=None):
  elapTime = time() - startTime
  remTime = wallTime - elapTime
  runCmd=shlex.split("mpiexec -n 144 -ppn 48 ./nemo > log")
  if remTime > 1.5*modelRunTime:
    maxTime = floor(remTime)
    if (mpi is None) or (mpi is not None and mpi.rank==0):
      try:
        runout = subprocess.run(runCmd,capture_output=True,text=True,timeout=maxTime)
        rr=0 #runout.returncode
        sleep(2) # allow time for model/file system to finish writing
      except subprocess.TimeoutExpired as e:
        print(f"Model call timedout")
        rr = CONVERGED_REASON['modelWallTime']
      except Exception as e:
        print(f"Model call failed with exception {e}")
        rr = CONVERGED_REASON['modelFailed']
    else:  
      sleep(modelRunTime)
      sleep(2) # allow time for model/file system to finish writing
      rr=None
    if mpi is not None:
      mpi.comm.Barrier()
      returncode = mpi.comm.bcast(rr, root = 0)
    else:
      returncode = rr
  else:
    # Insufficient time left to run model
    returncode = CONVERGED_REASON['modelWallTime']
  return returncode

def preProcess(i,mpi=None):
  # Note: we copy rather than move the files here to make restarts more robust
  if useGlobalFile:
    if (mpi is None) or (mpi is not None and mpi.rank==0):
      try:
        os.system(f"cp -p {restartFileEnd} {restartFileStart}")
        rr=0
      except Exception as e:
        print(f"Preprocessing call failed with exception {e}")
        rr = CONVERGED_REASON['modelPreprocessFailed']
    else:
      rr=None  
    if mpi is not None:
      mpi.comm.Barrier()
      returncode = mpi.comm.bcast(rr, root = 0)
    else:
      returncode = rr
  else:
    try:
      if mpi is not None:
        os.system(f"cp -p {restartFileEnd} {restartFileStart}")
        returncode=0
        mpi.comm.Barrier()
      else:
        os.system(f"printf \"%04d\n\" `seq 0 {numTiles-1}` | xargs -IDD cp -p {restartFileEnd}_DD.nc {restartFileStart}_DD.nc")
        returncode=0
    except Exception as e:
      print(f"Preprocessing call failed with exception {e}")
      returncode = CONVERGED_REASON['modelPreprocessFailed']
  return returncode

# def postProcess(i,mpi=None):
#   os.system("sleep 20")
#   if mpi is not None: mpi.comm.Barrier()

data.restartFileStart=restartFileStart
data.restartFileEnd=restartFileEnd
data.useMPI=useMPI
data.useGlobalFile=useGlobalFile
data.writeTRB=True # If using Euler time step set this to True
data.run_model_cmd=runModel
data.preprocess_cmd=preProcess
# data.postprocess_cmd=postProcess
data.calcNorms=True

fn="aa_restart.h5"
if os.path.isfile(fn):
  if doPrint: print("Reading restart file")
  s=loadaarestart(fn,mpi=mpiaa)
  x0=s.x
  y0=s.y
  # On restarting we want to use the ancillary BGC variables from the previous run. However, if 
  # for some reason the restart file is not available we use the initial condition.
  # Actually, I think this script will crash at the preprocessing stage if any of the files are missing.
  # We only do this when not in doRestart mode. In doRestart mode, restartFileEnd should already 
  # be present from the completed run and will be read (if s.fetchOutput==1). If restartFileEnd is 
  # not present, an error should be generated when timestepper attempts to read it.
  if not doRestart:
#   if (doRestart==0) or (doRestart==1 and s.fetchOutput != 1):
    if useGlobalFile:
      if (not useMPI) or (useMPI and mpiaa.rank==0):
        if not os.path.isfile(restartFileEnd):
          os.system(f"cp -p InitialCondition/{restartFileEnd} {restartFileEnd}")
    else:
      if useMPI: # each process copies its corresponding tile
        if not os.path.isfile(restartFileEnd):
           os.system(f"cp -p InitialCondition/{restartFileEnd} {restartFileEnd}")
      else:
        for iFile in range(numTiles):
          fn=f"{restartFileEnd}_{iFile:04d}.nc"
          if not os.path.isfile(fn):
            os.system(f"cp -p InitialCondition/{fn} {fn}")
      if useMPI and len(datainactive)>0:
        if mpiaa.rank==0:
          for d in datainactive:
            iFile=d.iFile
            restartFileStartInactive = f"restart_trc_{iFile:04d}.nc"
            restartFileEndInactive = f"ORCA2_00005840_restart_trc_{iFile:04d}.nc"
            # Note: The model will look for the start file
            if not os.path.isfile(restartFileStartInactive):
              os.system(f"cp -p InitialCondition/{restartFileEndInactive} {restartFileStartInactive}")
    if useMPI: mpiaa.comm.Barrier()
else:
  s = AndersonAccelerationObj({})
  # On first job submission need to copy over restart file
  if useGlobalFile:
    if (not useMPI) or (useMPI and mpiaa.rank==0):
      os.system(f"cp -p InitialCondition/{restartFileEnd} {restartFileEnd}")
  else:
    if useMPI: # each process copies its corresponding tile
      os.system(f"cp -p InitialCondition/{restartFileEnd} {restartFileEnd}")
    else:
      os.system(f"printf \"%04d\n\" `seq 0 {numTiles-1}` | xargs -IDD cp -p InitialCondition/{restartFileEnd}_DD.nc {restartFileEnd}_DD.nc")
    if useMPI and len(datainactive)>0:
      if mpiaa.rank==0 :
        for d in datainactive:
          iFile=d.iFile
          restartFileStartInactive = f"restart_trc_{iFile:04d}.nc"
          restartFileEndInactive = f"ORCA2_00005840_restart_trc_{iFile:04d}.nc"
          # Note: The model will look for the start file
          os.system(f"cp -p InitialCondition/{restartFileEndInactive} {restartFileStartInactive}")
  if useMPI: mpiaa.comm.Barrier()
  # Read initial condition from same file
  if useGlobalFile and useMPI: # MPI+global
    x0,y0=readFromNEMOGlobalRestart2xyloc(data, restartFileEnd, mpiaa, vbuf=data.vbuf)
  else: # All other cases (MPI+grouped tiles, MPI+tiles, serial+tiles, serial+global)
    x0,y0=readFromNEMORestart2xy(data, restartFileEnd, singleFile=(useGlobalFile or useMPI), vbuf=data.vbuf)
  if doPrint: print("Finished reading IC")

def checkpointingFunction(i,s,mpi=None):
  # TBD: error handling code
  if mpi is None or (mpi is not None and mpi.rank==0):
    fn=f"aa_co2flux_{i:04d}"
    ds=nc.Dataset(diagnosticsFile)
    tcflx = ds['tcflx'][:].data
    ds.close()
    to_save = {'iter_':i, 'tcflx': tcflx}
    save(to_save,fn)
    with np.printoptions(suppress=False):
      print("CO2flux: ",tcflx)
  if mpi is not None: mpi.comm.Barrier()
  if (i % 25 == 0):
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
#     os.system(f"cp -p {restartFileEnd} {restartFileEnd}_{i:04d}")
#   if (i>0) and (i % 25 == 0):
#     if mpi is not None: mpi.comm.Barrier()
#     fn="aa_restart_checkpoint"
# #   Hack to make this restart usable for actual restarts
# #   Save iter0 before setting it to current iter
#     iter0_=s.iter0
# #   Careful here: this function is called at the end of the loop so we want to restart
# #   with the next iteration number
#     s.iter0=i+1
#     saveaarestart(s,fn,mpi=mpiaa)
# #   Restore iter0
#     s.iter0 = iter0_
# #     os.system(f"cp -p {restartFileEnd} {restartFileEnd}_checkpoint")
#     if mpi is not None: mpi.comm.Barrier()
  return 0

histParams.checkpointingFunction=checkpointingFunction  

if useMPI: mpiaa.comm.Barrier()

g = lambda x, y, fetchOutput, iter_: timestepperfunc_nemo(x, y, fetchOutput, iter_, data, iJob, doSubmit, mpi=mpiaa)
suff='_test'
xsol, iter_, ysol, converged = andacc(g, x0, s, AAparams, histParams, doRestart, fileSuff=suff, y=y0, debug=False, mpi=mpiaa)

if converged in (CONVERGED_REASON['tol'],CONVERGED_REASON['itmax']):
# AA finished
  if doPrint: print("Saving solution")
  if useMPI:
    xg0=xloc2glob(x0,layout=mpiaa.layout, MPI=mpiaa.MPI)
    yg0=xloc2glob(y0,layout=mpiaa.layouty, MPI=mpiaa.MPI)
    if useGlobalFile: # MPI+global
      xsolg=xloc2glob(xsol,layout=mpiaa.layout, MPI=mpiaa.MPI)
      ysolg=xloc2glob(ysol,layout=mpiaa.layouty, MPI=mpiaa.MPI)
      if mpiaa.rank==0:
        tmparr=xy2nemo(xsolg, ysolg, data, vbuf=data.vbuf)
      mpiaa.comm.Barrier()
    else: # MPI+grouped tiles or MPI+tiles
      tmparr,xsolg,ysolg=xyloc2nemoglob(xsol, ysol, dataglob, datatiles, mpiaa)
    if mpiaa.rank==0:
      to_save = {'x0': xg0, 'y0':yg0, 'xsol':xsolg, 'iter_':iter_, 'ysol':ysolg, 'arr':tmparr, 'converged':converged}
      fn="andacc_solution.h5"
      save(to_save,fn)
    mpiaa.comm.Barrier()
  else: # serial+global, serial+tiles
    tmparr = xy2nemo(xsol, ysol, data, vbuf=data.vbuf)
    to_save = {'x0': x0, 'y0':y0, 'xsol':xsol, 'iter_':iter_, 'ysol':ysol, 'arr':tmparr, 'converged':converged}
    fn="andacc_solution.h5"
    save(to_save,fn)
else:
  if converged in (CONVERGED_REASON['wallTime'],CONVERGED_REASON['modelWallTime']):
    if doPrint: print(f"Wall time exceeded")

fn="aa_restart.h5"
saveaarestart(s,fn,mpi=mpiaa)

if converged >= 0:
  sys.exit(0)
else:
  sys.exit(1)

# os.system("./signal_model 0")
# sleep(10)
