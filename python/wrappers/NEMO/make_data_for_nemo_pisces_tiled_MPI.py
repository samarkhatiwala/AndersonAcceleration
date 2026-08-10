import numpy as np
from aa4py import DotDict
from pyioutils.hdfio import load, save
import netCDF4 as nc

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

# In the following tmask generally refers to the land/sea mask 
# and Idx (m2vIb) maps the wet points of the model to "v"

wtr=["TRNDIC","TRNAlkalini","TRNO2","TRNCaCO3","TRNPO4","TRNPOC","TRNSi","TRNPHY","TRNZOO","TRNDOC","TRNPHY2","TRNZOO2","TRNDSi","TRNFer","TRNBFe","TRNGOC","TRNSFe","TRNDFe","TRNGSi","TRNNFe","TRNNCHL","TRNDCHL","TRNNO3","TRNNH4"]

data=[DotDict() for _ in range(numDomains)]
for iFile in range(numDomains):
  data[iFile].iFile=iFile
  data[iFile].tmask=[tmaskl[iFile]]
  data[iFile].tmaskId=[0 for itr in range(len(wtr))]
  data[iFile].trNames=wtr
  data[iFile].numTracers=len(data[iFile].trNames)
  Idxt=np.flatnonzero(tmaskl[iFile])
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
datainactive=[d for d in data if d.nb==0]
dataactive=[d for d in data if d.nb>0]

save({'data': dataactive, 'datainactive': datainactive}, "data_for_nemo_pisces_tiled")
