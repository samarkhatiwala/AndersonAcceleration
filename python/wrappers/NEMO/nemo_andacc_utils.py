import numpy as np
import netCDF4 as nc
import os.path
# I've moved this import to functions below to remove dependency on MPI
# from pympiutils import xloc2glob, xglob2loc

def writeToNEMORestart(trNames, arr, restartFile, singleFile=True, writeTRB=False):
  """
  Write global tracer arrays in list arr to restart file
  """

### NOTE: Set writeTRB=True ONLY IF using Euler time stepping; this will write a copy of 
### the TRN* arrays (which should be in trNames) to the corresponding TRB* arrays in the file. 
### It is unclear if or why NEMO needs these but just to be safe!

  numTracers=len(trNames)
  if len(arr) != numTracers:
    raise Exception("Number of elements in arr does not agree with number of tracer names")

  if singleFile:
    ncfile = nc.Dataset(restartFile,mode='r+',format='NETCDF4')
    ic=0
    for itr in range(numTracers):
      if len(arr[itr].shape)==4:
        ncfile[trNames[itr]][ic,:,:,:]=arr[itr][0,:,:,:]
        if writeTRB:
          tr='TRB'+trNames[itr][3:]
          ncfile[tr][ic,:,:,:]=arr[itr][0,:,:,:]
      elif len(arr[itr].shape)==3:
        # This is specific to MEDUSA sediment tracers
        ncfile[trNames[itr]][ic,:,:]=arr[itr][0,:,:]
        if writeTRB:
          tr='B_'+trNames[itr][2:]
          ncfile[tr][ic,:,:]=arr[itr][0,:,:]
      else:
        raise Exception("Unknown dimensions!")
    ncfile.close()

  else:
    dFile=restartFile + "_{:04d}.nc".format(0)
    ncfile = nc.Dataset(dFile)
    numDomains=ncfile.DOMAIN_number_total
    ncfile.close()
  
    # write data to subdomain files from global array
    for iFile in range(numDomains):
      dFile=restartFile + "_{:04d}.nc".format(iFile)
      ncfile = nc.Dataset(dFile,mode='r+',format='NETCDF4')
      nx,ny=ncfile.DOMAIN_size_local
      i1,j1=ncfile.DOMAIN_position_first-1
      i2,j2=ncfile.DOMAIN_position_last-1
      for itr in range(numTracers):
        if len(arr[itr].shape)==4:
          ncfile[trNames[itr]][0,:,:ny,:nx]=arr[itr][0,:,j1:j1+ny,i1:i1+nx]
          if writeTRB:
            tr='TRB'+trNames[itr][3:]
            ncfile[tr][0,:,:ny,:nx]=arr[itr][0,:,j1:j1+ny,i1:i1+nx]
        elif len(arr[itr].shape)==3:
          ncfile[trNames[itr]][0,:ny,:nx]=arr[itr][0,j1:j1+ny,i1:i1+nx]
#         This is specific to MEDUSA sediment tracers
          if writeTRB:
            tr='B_'+trNames[itr][2:]
            ncfile[tr][0,:ny,:nx]=arr[itr][0,j1:j1+ny,i1:i1+nx]
        else:
          raise Exception("Unknown dimensions!")
      ncfile.close()

def readFromNEMORestart(trNames, restartFile, tmask, tmaskId, singleFile=True):

  """
  Read tracer arrays from restart file and return as global arrays in list arr
  """

  numTracers=len(trNames)

  if len(tmaskId) != numTracers:
    raise Exception("Number of elements in tmaskId does not agree with number of tracer names!")

  arr=[tmask[tmaskId[itr]].copy() for itr in range(numTracers)]
  
  if singleFile:
    ncfile = nc.Dataset(restartFile)
    ic=0
    for itr in range(numTracers):
      if len(arr[itr].shape)==4:
        arr[itr][0,:,:,:]=ncfile[trNames[itr]][ic,:,:,:].data
      elif len(arr[itr].shape)==3:
        arr[itr][0,:,:]=ncfile[trNames[itr]][ic,:,:].data
      else:
        raise Exception("Unknown dimensions!")
    ncfile.close()

  else:
    dFile=restartFile + "_{:04d}.nc".format(0)
    ncfile = nc.Dataset(dFile)
    numDomains=ncfile.DOMAIN_number_total
    ncfile.close()

    # read data from subdomain files into global array
    for iFile in range(numDomains):
      dFile=restartFile + "_{:04d}.nc".format(iFile)
      ncfile = nc.Dataset(dFile)
      nx,ny=ncfile.DOMAIN_size_local
      i1,j1=ncfile.DOMAIN_position_first-1
      i2,j2=ncfile.DOMAIN_position_last-1
      for itr in range(numTracers):
        if len(arr[itr].shape)==4:
          arr[itr][0,:,j1:j1+ny,i1:i1+nx]=ncfile[trNames[itr]][0,:,:ny,:nx]
        elif len(arr[itr].shape)==3:
          arr[itr][0,j1:j1+ny,i1:i1+nx]=ncfile[trNames[itr]][0,:ny,:nx]
        else:
          raise Exception("Unknown dimensions!")
      ncfile.close()
    
  return arr

def xy2nemo(x, y, data, vbuf=None):

  if len(data.tmaskId) != data.numTracers:
    raise Exception("Number of elements in tmaskId does not agree with number of tracer names")

  if vbuf is None:
    v = np.zeros(shape=(data.nb,data.numTracers))
  else:
    v = vbuf
    if v.shape != (data.nb,data.numTracers):
      raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")

  v=_xy2v(x,y,data,vbuf=v)

  arr=[data.tmask[data.tmaskId[itr]].copy() for itr in range(data.numTracers)]

  # insert v into global array
  for itr in range(data.numTracers):
    tId=data.tmaskId[itr]
    nbt=len(data.m2vIb[tId])
    arr[itr].flat[data.m2vIb[tId]]=v[0:nbt,itr]

  return arr

def nemo2xy(arr, data, vbuf=None):

  if len(arr) != data.numTracers:
    raise Exception("Number of elements in arr does not agree with expected number of tracers")

  if vbuf is None:
    v = np.zeros(shape=(data.nb,data.numTracers))
  else:
    v = vbuf
    if v.shape != (data.nb,data.numTracers):
      raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")

  # insert global array data into v
  for itr in range(data.numTracers):
    tId=data.tmaskId[itr]
    nbt=len(data.m2vIb[tId])
    v[0:nbt,itr]=arr[itr].flat[data.m2vIb[tId]]

  x,y=_v2xy(v,data)

  return x, y

def _xy2v(x, y, data, vbuf=None):

  if vbuf is None:
    v = np.zeros(shape=(data.nb,data.numTracers))
  else:
    v = vbuf
    if v.shape != (data.nb,data.numTracers):
      raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")
    
  v.flat[data.x2vIb]=x

  if y is not None:
    v.flat[data.y2vIb]=y
  
  if "XTovScaleFac" in data:
    for itr in range(data.numTracers):
      v[:,itr] = v[:,itr]*data.XTovScaleFac[itr,itr]

  return v

def _v2xy(v, data):

  if v.shape[1] != data.numTracers:
    raise Exception("Number of columns in v does not agree with expected number of tracers")

  if "vToXScaleFac" in data:
    for itr in range(data.numTracers):
      v[:,itr] = v[:,itr]*data.vToXScaleFac[itr,itr]
  gtr=v
          
  x=gtr.flat[data.x2vIb]

  if np.array_equal(data.y2vIb,np.array([])):
    y=None
  else:
    y=gtr.flat[data.y2vIb]

  return x, y

def read2dNEMOField(fldName, fileName, singleFile=True):
  if singleFile:
    ncfile = nc.Dataset(fileName)
#     nxg,nyg=ncfile.DOMAIN_size_global
    nxg=ncfile.dimensions["x"].size
    nyg=ncfile.dimensions["y"].size
    arr=np.zeros(shape=(1,nyg,nxg))
    arr[0,:,:]=ncfile[fldName][:,:,:].data
    ncfile.close()
  else:
    dFile=fileName + "_{:04d}.nc".format(0)
    ncfile = nc.Dataset(dFile)
    nxg,nyg=ncfile.DOMAIN_size_global
    arr=np.zeros(shape=(1,nyg,nxg))
    numDomains=ncfile.DOMAIN_number_total
    ncfile.close()
#   read data from subdomain files into global array
    for iFile in range(numDomains):
      dFile=fileName + "_{:04d}.nc".format(iFile)
      ncfile = nc.Dataset(dFile)
      nx,ny=ncfile.DOMAIN_size_local
      i1,j1=ncfile.DOMAIN_position_first-1
      i2,j2=ncfile.DOMAIN_position_last-1
      arr[0,j1:j1+ny,i1:i1+nx]=ncfile[fldName][0,:ny,:nx]
      ncfile.close()
  return arr

def read3dNEMOField(fldName, fileName, singleFile=True):
  if singleFile:
    ncfile = nc.Dataset(fileName)
    nxg,nyg=ncfile.DOMAIN_size_global
    nz=ncfile["nav_lev"].size
    arr=np.zeros(shape=(1,nz,nyg,nxg))
    arr[0,:,:,:]=ncfile[fldName][:,:,:,:].data
    ncfile.close()
  else:
    dFile=fileName + "_{:04d}.nc".format(0)
    ncfile = nc.Dataset(dFile)
    nxg,nyg=ncfile.DOMAIN_size_global
    nz=ncfile["nav_lev"].size
    arr=np.zeros(shape=(1,nz,nyg,nxg))
    numDomains=ncfile.DOMAIN_number_total
    ncfile.close()
#   read data from subdomain files into global array
    for iFile in range(numDomains):
      dFile=fileName + "_{:04d}.nc".format(iFile)
      ncfile = nc.Dataset(dFile)
      nx,ny=ncfile.DOMAIN_size_local
      i1,j1=ncfile.DOMAIN_position_first-1
      i2,j2=ncfile.DOMAIN_position_last-1
      arr[0,:,j1:j1+ny,i1:i1+nx]=ncfile[fldName][0,:,:ny,:nx]
      ncfile.close()
  return arr

def writexyloc2NEMOGlobalRestart(x, y, data, restartFile, mpi, vbuf=None, writeTRB=False):
  # MPI+global  
  # Collective

  from pympiutils import xloc2glob
  
  # First gather all local pieces to rank 0
  xg=xloc2glob(x,layout=mpi.layout, MPI=mpi.MPI)
  yg=xloc2glob(y,layout=mpi.layouty, MPI=mpi.MPI)
  
  if mpi.rank==0:
    if vbuf is None:
      v = np.zeros(shape=(data.nb,data.numTracers))
    else:
      v = vbuf
      if v.shape != (data.nb,data.numTracers):
        raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")
  
    v=_xy2v(xg,yg,data,vbuf=v)

    if "setNegativesToZero" in data:
      if data.setNegativesToZero:
        if "tracersToClip" in data:
          for tr in data.tracersToClip:
            itr=data.trNames.index(tr)
            v[:,itr]=np.clip(v[:,itr],0.,None,out=v[:,itr])
        else:
          # clip all tracers
          v=np.clip(v,0.,None,out=v)

    # Note we only create an array for each type (3-d, 2-d etc). arr will thus be a list of 
    # length of data.tmask
    arr=[data.tmask[tId] for tId in range(len(data.tmask))]
    # We write the fields one at a time
    for itr in range(data.numTracers):
      tId=data.tmaskId[itr]
      # insert v into array
      nbt=len(data.m2vIb[tId])
      arr[tId].flat[data.m2vIb[tId]]=v[0:nbt,itr]
      # Be careful to pass everything as a list
      writeToNEMORestart([data.trNames[itr]], [arr[tId]], restartFile, singleFile=True, writeTRB=writeTRB)
    
  mpi.comm.Barrier()
  
def writexy2NEMORestart(x, y, data, restartFile, singleFile=True, vbuf=None, writeTRB=False):
  # MPI+grouped tiles
  # MPI+tiles
  # serial+tiles
  # serial+global
  # For serial+tiles, pass singleFile=False
  # For MPI+tiles and serial+global, pass singleFile=True
  # For MPI+grouped tiles, singleFile is ignored and hardwired to True below

  if isinstance(data,list): # MPI+grouped tiles
    nx=[d.x2vIb.size for d in data]
    iStartx=0
    if y is not None:
      ny=[d.y2vIb.size for d in data]
      iStarty=0
    for i, d in enumerate(data):
      iEndx=iStartx+nx[i]
      xt=x[iStartx:iEndx]
      iStartx=iEndx
      if y is not None and ny[i]>0:
        iEndy=iStarty+ny[i]
        yt=y[iStarty:iEndy]
        iStarty=iEndy
      else:
        yt=None
      v=_xy2v(xt,yt,d,vbuf=None)

      if "setNegativesToZero" in d:
        if d.setNegativesToZero:
          if "tracersToClip" in d:
            for tr in d.tracersToClip:
              itr=d.trNames.index(tr)
              v[:,itr]=np.clip(v[:,itr],0.,None,out=v[:,itr])
          else:
            # clip all tracers
            v=np.clip(v,0.,None,out=v)

      # Note we only create an array for each type (3-d, 2-d etc). arr will thus be a list of 
      # length of data.tmask
      arr=[d.tmask[tId] for tId in range(len(d.tmask))]
      # We write the fields one at a time
      for itr in range(d.numTracers):
        tId=d.tmaskId[itr]
        # insert v into array
        nbt=len(d.m2vIb[tId])
        arr[tId].flat[d.m2vIb[tId]]=v[0:nbt,itr]
        # Be careful to pass everything as a list
        writeToNEMORestart([d.trNames[itr]], [arr[tId]], restartFile[i], singleFile=True, writeTRB=writeTRB)
  
  else: # MPI+tiles, serial+tiles or serial+global
    if vbuf is None:
      v = np.zeros(shape=(data.nb,data.numTracers))
    else:
      v = vbuf
      if v.shape != (data.nb,data.numTracers):
        raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")

    v=_xy2v(x,y,data,vbuf=v)

    if "setNegativesToZero" in data:
      if data.setNegativesToZero:
        if "tracersToClip" in data:
          for tr in data.tracersToClip:
            itr=data.trNames.index(tr)
            v[:,itr]=np.clip(v[:,itr],0.,None,out=v[:,itr])
        else:
          # clip all tracers
          v=np.clip(v,0.,None,out=v)

    # Note we only create an array for each type (3-d, 2-d etc). arr will thus be a list of 
    # length of data.tmask
    arr=[data.tmask[tId] for tId in range(len(data.tmask))]
    # We write the fields one at a time
    for itr in range(data.numTracers):
      tId=data.tmaskId[itr]
      # insert v into array
      nbt=len(data.m2vIb[tId])
      arr[tId].flat[data.m2vIb[tId]]=v[0:nbt,itr]
      # Be careful to pass everything as a list
      writeToNEMORestart([data.trNames[itr]], [arr[tId]], restartFile, singleFile=singleFile, writeTRB=writeTRB)

def readFromNEMOGlobalRestart2xyloc(data, restartFile, mpi, vbuf=None):
  # MPI+global  
  # Collective

  from pympiutils import xglob2loc

  if mpi.rank==0:
    arr = readFromNEMORestart(data.trNames, restartFile, data.tmask, data.tmaskId, singleFile=True)
    xg,yg = nemo2xy(arr, data, vbuf=vbuf)
  else:
    xg=None
    yg=None

  mpi.comm.Barrier()
  
  # Now scatter global vector to other ranks
  x=xglob2loc(xg,layout=mpi.layout, MPI=mpi.MPI)
  y=xglob2loc(yg,layout=mpi.layouty, MPI=mpi.MPI)
  
  return x, y

def readFromNEMORestart2xy(data, restartFile, singleFile=True, vbuf=None):
  # MPI+grouped tiles
  # MPI+tiles
  # serial+tiles
  # serial+global
  # For serial+tiles, pass singleFile=False
  # For MPI+tiles and serial+global, pass singleFile=True
  # For MPI+grouped tiles, singleFile is ignored and hardwired to True below

  if isinstance(data,list): # MPI+grouped tiles
    nx=[d.x2vIb.size for d in data]
    x=np.zeros(sum(nx))
    iStartx=0
    y=None
    iStarty=None
    # Assemble local x from multiple tiles
    for i, d in enumerate(data):
      v = np.zeros(shape=(d.nb,d.numTracers))
      # We read the fields one at a time
      for itr in range(d.numTracers):
        tId=d.tmaskId[itr]
        # Be careful to pass everything as a list; note that arr is returned as a list with a single element
        arr = readFromNEMORestart([d.trNames[itr]], restartFile[i], [d.tmask[tId]], tmaskId=[0], singleFile=True)
        # insert array data into v
        nbt=len(d.m2vIb[tId])
        v[0:nbt,itr]=arr[0].flat[d.m2vIb[tId]]    
      xt,yt=_v2xy(v,d)
      iEndx=iStartx+nx[i]
      x[iStartx:iEndx]=xt
      iStartx=iEndx
      if yt is not None:
        if iStarty is None:
          ny=[d.y2vIb.size for d in data]
          if sum(ny)>0:
            y=np.zeros(sum(ny))
            iStarty=0
        if ny[i]>0:
          iEndy=iStarty+ny[i]
          y[iStarty:iEndy]=yt
          iStarty=iEndy
  else: # MPI+tiles, serial+tiles or serial+global
    if vbuf is None:
      v = np.zeros(shape=(data.nb,data.numTracers))
    else:
      v = vbuf
      if v.shape != (data.nb,data.numTracers):
        raise Exception(f"Buffer array has the wrong shape {v.shape}; expected shape: {(data.nb,data.numTracers)}")
   
    # We read the fields one at a time
    for itr in range(data.numTracers):
      tId=data.tmaskId[itr]
      # Be careful to pass everything as a list; note that arr is returned as a list with a single element
      arr = readFromNEMORestart([data.trNames[itr]], restartFile, [data.tmask[tId]], tmaskId=[0], singleFile=singleFile)
      # insert array data into v
      nbt=len(data.m2vIb[tId])
      v[0:nbt,itr]=arr[0].flat[data.m2vIb[tId]]
    x,y=_v2xy(v,data)

  return x, y
  
def xyloc2nemoglob(x, y, dataglob, datatiles, mpi):
  # MPI+grouped tiles, MPI+tiles
  # Collective

  from pympiutils import xloc2glob
  
  # nemoglob2xyloc and xyloc2nemoglob are provided to map to and from global fields and AA vectors 
  # distributed according to tiles. You would only use them when running AA with MPI and tiled 
  # restart files because in that case processes handle individual tiles independently of others 
  # and there is no information available (or needed) as to how each tile fits into the global field. 
  # However, a user may want to view the complete tracer field from the AA solutions, for example 
  # at the end of a spin-up. Or read initial conditions from a global restart file at the beginning 
  # of a spin-up. These routines provide that functionality (somewhat clunkily but I couldn't come 
  # up with a better way).
  
  # check we have the right data
  if mpi.size != len(datatiles):
    raise Exception("Number of elements in datatiles does not agree with number of MPI ranks")
    
  # First gather all local pieces to rank 0
  xg=xloc2glob(x,layout=mpi.layout, MPI=mpi.MPI)
  yg=xloc2glob(y,layout=mpi.layouty, MPI=mpi.MPI)

  # Now assemble the full model fields from local pieces
  if mpi.rank==0:
    # global tracer field
    arrg=[dataglob.tmask[dataglob.tmaskId[itr]].copy() for itr in range(dataglob.numTracers)]

    for i in range(mpi.size):
      d=datatiles[i]
      # Get local piece of global vector for rank i
      iStart=mpi.layout.displ[i]
      iEnd=iStart+mpi.layout.nlocal_sizes[i]
      xl=xg[iStart:iEnd]
      if yg is not None:
        iStart=mpi.layouty.displ[i]
        iEnd=iStart+mpi.layouty.nlocal_sizes[i]
        yl=yg[iStart:iEnd]
      else:
        yl=None
        
      # Insert xl,yl into tile field(s)
      if isinstance(d,list): # MPI+grouped tiles
        nx=[dt.x2vIb.size for dt in d]
        iStartx=0
        if y is not None:
          ny=[dt.y2vIb.size for dt in d]
          iStarty=0
        for i, dt in enumerate(d):
          iEndx=iStartx+nx[i]
          xt=xl[iStartx:iEndx]
          iStartx=iEndx
          if y is not None and ny[i]>0:
            iEndy=iStarty+ny[i]
            yt=yl[iStarty:iEndy]
            iStarty=iEndy
          else:
            yt=None
          arrt=xy2nemo(xt, yt, dt, vbuf=None)
          # Insert local tile into global field
          [nxt,nyt,i1,j1,i2,j2] = dt.locsize
          for itr in range(dataglob.numTracers):
            if len(arrg[itr].shape)==4:
              arrg[itr][0,:,j1:j1+nyt,i1:i1+nxt]=arrt[itr][0,:,:nyt,:nxt]
            elif len(arrg[itr].shape)==3:
              arrg[itr][0,j1:j1+nyt,i1:i1+nxt]=arrt[itr][0,:nyt,:nxt]
            else:
              raise Exception("Unknown dimensions!")
      else: # MPI+tiles
        arrt=xy2nemo(xl, yl, d, vbuf=None)
        # Insert local tile into global field
        [nxt,nyt,i1,j1,i2,j2] = d.locsize
        for itr in range(dataglob.numTracers):
          if len(arrg[itr].shape)==4:
            arrg[itr][0,:,j1:j1+nyt,i1:i1+nxt]=arrt[itr][0,:,:nyt,:nxt]
          elif len(arrg[itr].shape)==3:
            arrg[itr][0,j1:j1+nyt,i1:i1+nxt]=arrt[itr][0,:nyt,:nxt]
          else:
            raise Exception("Unknown dimensions!")
  else:
    arrg = None

  mpi.comm.Barrier()

  return arrg, xg, yg

def nemoglob2xyloc(arrg, datatiles, mpi):
  # MPI+grouped tiles, MPI+tiles
  # Collective (arrg should be None on all but rank 0)

  from pympiutils import xglob2loc

  # nemoglob2xyloc and xyloc2nemoglob are provided to map to and from global fields and AA vectors 
  # distributed according to tiles. You would only use them when running AA with MPI and tiled 
  # restart files because in that case processes handle individual tiles independently of others 
  # and there is no information available (or needed) as to how each tile fits into the global field. 
  # However, a user may want to view the complete tracer field from the AA solutions, for example 
  # at the end of a spin-up. Or read initial conditions from a global restart file at the beginning 
  # of a spin-up. These routines provide that functionality (somewhat clunkily but I couldn't come 
  # up with a better way).

  # check we have the right data
  if mpi.size != len(datatiles):
    raise Exception("Number of elements in datatiles does not agree with number of MPI ranks")

  if mpi.rank==0:
    xg = np.zeros(mpi.layout.n)
    if mpi.layouty is not None:
      yg = np.zeros(mpi.layouty.n)
    else:
      yg = None
      
    for i in range(mpi.size):
      d=datatiles[i]
      if isinstance(d,list):
        # assemble local vector from multiple tiles
        nx=[d.x2vIb.size for d in data]
        xl=np.zeros(sum(nx))
        iStartx=0
        yl=None
        iStarty=None
        for i, dt in enumerate(d):
          # extract tile field from global field
          [nxt,nyt,i1,j1,i2,j2] = dt.locsize
          arrt=[dt.tmask[dt.tmaskId[itr]].copy() for itr in range(dt.numTracers)]
          for itr in range(numTracers):
            if len(arrg[itr].shape)==4:
              arrt[itr][0,:,:nyt,:nxt]=arrg[itr][0,:,j1:j1+nyt,i1:i1+nxt]
            elif len(arrg[itr].shape)==3:
              arrt[itr][0,:nyt,:nxt]=arrg[itr][0,j1:j1+nyt,i1:i1+nxt]
            else:
              raise Exception("Unknown dimensions!")
          xt,yt=nemo2xy(arrt, dt, vbuf=None)
          iEndx=iStartx+nx[i]
          xl[iStartx:iEndx]=xt
          iStartx=iEndx
          if yt is not None:
            if iStarty is None:
              ny=[d.y2vIb.size for d in data]
              if sum(ny)>0:
                yl=np.zeros(np.sum(ny))
                iStarty=0
            if ny[i]>0:
              iEndy=iStarty+ny[i]
              yl[iStarty:iEndy]=yt
              iStarty=iEndy
      else:
        # extract tile field from global field
        [nxt,nyt,i1,j1,i2,j2] = d.locsize
        arrt=[d.tmask[d.tmaskId[itr]].copy() for itr in range(d.numTracers)]
        for itr in range(numTracers):
          if len(arrg[itr].shape)==4:
            arrt[itr][0,:,:nyt,:nxt]=arrg[itr][0,:,j1:j1+nyt,i1:i1+nxt]
          elif len(arrg[itr].shape)==3:
            arrt[itr][0,:nyt,:nxt]=arrg[itr][0,j1:j1+nyt,i1:i1+nxt]
          else:
            raise Exception("Unknown dimensions!")
      
        xl,yl=nemo2xy(arrt, d, vbuf=None)
      
      # Insert local vector for rank i into global vector
      iStart=mpi.layout.displ[i]
      iEnd=iStart+mpi.layout.nlocal_sizes[i]
      xg[iStart:iEnd]=xl
      if yg is not None:
        iStart=mpi.layouty.displ[i]
        iEnd=iStart+mpi.layouty.nlocal_sizes[i]
        yg[iStart:iEnd]=yl
  else:
    xg = None
    yg = None

  mpi.comm.Barrier()
  
  # finally scatter global vector to local vectors
  xl=xglob2loc(xg,layout=mpi.layout, MPI=mpi.MPI)
  yl=xglob2loc(yg,layout=mpi.layouty, MPI=mpi.MPI)

  return xl, yl
