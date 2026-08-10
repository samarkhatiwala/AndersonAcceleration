import numpy as np
from aa4py import DotDict
from pyioutils.hdfio import load, save
import netCDF4 as nc

useSED=True

fileName="mesh_mask"

dFile=fileName + "_{:04d}.nc".format(0)
ncfile = nc.Dataset(dFile)
numDomains=ncfile.DOMAIN_number_total
ncfile.close()

tmaskl=[[] for _ in range(numDomains)]
locsize=[[] for _ in range(numDomains)]
for iFile in range(numDomains):
  dFile=fileName + "_{:04d}.nc".format(iFile)
  ncfile = nc.Dataset(dFile)
  tmaskl[iFile]=ncfile['tmask'][:]*1.
  nx,ny=ncfile.DOMAIN_size_local
  i1,j1=ncfile.DOMAIN_position_first-1
  i2,j2=ncfile.DOMAIN_position_last-1
  locsize[iFile]=[nx,ny,i1,j1,i2,j2]
  ncfile.close()

if useSED: 
  sedmaskl=[tmaskl[iFile][:,0,:,:] for iFile in range(numDomains)]

# In the following tmask generally refers to the land/sea mask 
# and Idx (m2vIb) maps the wet points of the model to "v"

# Idxtl=[np.flatnonzero(tmaskl[iFile]) for iFile in range(numDomains)]
# Idxsl=[np.flatnonzero(sedmaskl[iFile]) for iFile in range(numDomains)]

# nbl=[len(Idxl[iFile]) for iFile in range(numDomains)] # number of wet points

tr = ['CHN','CHD','PHN','PHD','ZMI','ZME','DIN','SIL','FER','DET','PDS','DTC','DiC','ALK','OXY']
wtr = ['None']*len(tr)*2
wtr[0:len(wtr):2] = ['TRB'+tr[i] for i in range(len(tr))]
wtr[1:len(wtr):2] = ['TRN'+tr[i] for i in range(len(tr))]

if useSED:
  tr = ['SED_N','SED_FE','SED_SI','SED_C','SED_CA']
  sedtr = ['None']*len(tr)*2
  sedtr[0:len(sedtr):2] = ['B_'+tr[i] for i in range(len(tr))]
  sedtr[1:len(sedtr):2] = ['N_'+tr[i] for i in range(len(tr))]

data=[DotDict() for _ in range(numDomains)]
for iFile in range(numDomains):
  data[iFile].iFile=iFile
  if useSED:
    data[iFile].tmask=[tmaskl[iFile], sedmaskl[iFile]]
    data[iFile].tmaskId=[0 for itr in range(len(wtr))] + [1 for itr in range(len(sedtr))]
    data[iFile].trNames=wtr+sedtr
  else:
    data[iFile].tmask=[tmaskl[iFile]]
    data[iFile].tmaskId=[0 for itr in range(len(wtr))]
    data[iFile].trNames=wtr
  data[iFile].numTracers=len(data[iFile].trNames)
  Idxt=np.flatnonzero(tmaskl[iFile])
  if useSED:
    Idxs=np.flatnonzero(sedmaskl[iFile])
    Idx=[Idxt, Idxs]
  else:
    Idx=[Idxt]
  data[iFile].trnb=[len(Idx[data[iFile].tmaskId[itr]]) for itr in range(data[iFile].numTracers)]
  data[iFile].trnb=np.array(data[iFile].trnb)
  data[iFile].nb=max(data[iFile].trnb)
  data[iFile].m2vIb=Idx
  data[iFile].locsize=locsize[iFile] # (nx,ny,i1,j1,i2,j2)
  vscaleFac=np.ones(data[iFile].numTracers)
  data[iFile].XTovScaleFac=np.diag(vscaleFac)
  data[iFile].vToXScaleFac=np.diag(1./vscaleFac)
  vib=np.zeros(shape=(data[iFile].nb,data[iFile].numTracers))
  for itr in range(data[iFile].numTracers):
    nbt=data[iFile].trnb[itr]
    vib[0:nbt,itr]=1
  data[iFile].x2vIb=np.flatnonzero(vib==1)
  data[iFile].y2vIb=np.flatnonzero(vib==2)

# Delete tiles without wet points
data=[d for d in data if d.nb>0]

save({'data': data}, "data_for_nemo_medusa_with_sed_tiled")
