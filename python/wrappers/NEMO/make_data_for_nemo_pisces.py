import numpy as np
from aa4py import DotDict
from pyioutils.hdfio import load, save
from nemo_andacc_utils import read3dNEMOField, read2dNEMOField
import netCDF4 as nc

tmask=read3dNEMOField("tmask", "mesh_mask", singleFile=False)

# In the following tmask generally refers to the land/sea mask 
# and Idx (m2vIb) maps the wet points of the model to "v"

Idxt=np.flatnonzero(tmask)  #np.nonzero(np.ravel(tmask != 0))[0]

wtr=["TRNDIC","TRNAlkalini","TRNO2","TRNCaCO3","TRNPO4","TRNPOC","TRNSi","TRNPHY","TRNZOO","TRNDOC","TRNPHY2","TRNZOO2","TRNDSi","TRNFer","TRNBFe","TRNGOC","TRNSFe","TRNDFe","TRNGSi","TRNNFe","TRNNCHL","TRNDCHL","TRNNO3","TRNNH4"]

Idx=[Idxt]

data=DotDict()
data.tmask=[tmask]
data.tmaskId=[0 for itr in range(len(wtr))]
data.trNames=wtr
data.numTracers=len(data.trNames)
data.trnb=[len(Idx[data.tmaskId[itr]]) for itr in range(data.numTracers)]
data.trnb=np.array(data.trnb)
data.nb=max(data.trnb)
data.m2vIb=Idx

vscaleFac=np.ones(data.numTracers)
data.XTovScaleFac=np.diag(vscaleFac)
data.vToXScaleFac=np.diag(1./vscaleFac)

# Add surface area to data file for getting global co2 flux
e1t=read2dNEMOField("e1t", "mesh_mask", singleFile=False)
e2t=read2dNEMOField("e2t", "mesh_mask", singleFile=False)
dA=e1t*e2t
data.dA=dA

vib=np.zeros(shape=(data.nb,data.numTracers))
for itr in range(data.numTracers):
  nbt=data.trnb[itr]
  vib[0:nbt,itr]=1

data.x2vIb=np.flatnonzero(vib==1)
data.y2vIb=np.flatnonzero(vib==2)

save(data, "data_for_nemo_pisces.h5")
