# Software for solving fixed point problems with Anderson Acceleration

This repository contains MATLAB and Python software for solving fixed point problems using Anderson Acceleration. While it can be used for any fixed point problem, the primary motivation in developing it was for applications to computing equilibrium solutions of seasonally-forced ocean and land biogeochemical models. As such, the code is specifically designed to run on the batch HPC systems on which such models run. Furthermore, the Python version is fully parallelized making it feasible to apply it to large problems.

The software is based on the MATLAB routine AndAcc.m originally written by Homer F. Walker (Anderson acceleration: Algorithms and implementations, Worcester Polytechnic Institute Mathematical Sciences Department Research Report MS-6-15-50, 10/14/2011) with extensive modifications as per Khatiwala (2023), "Fast spin-up of geochemical tracers in ocean circulation and climate models, J. Adv. Model. Earth Sys., https://doi.org/10.1029/2022MS003447. A huge thanks to [Jamie Carr](https://github.com/colourfulLanguage) for helping port the MATLAB version to Python.

Note: I've currently paused developing the MATLAB version further while I focus on the  Python version.

**IMPORTANT**: Please do NOT post this code on your own github or other website. See LICENSE.txt for licensing information.

Feel free to email if you have any questions: <samkat6@gmail.com>

**Citing this software:**

Please cite this software as: Khatiwala, S. (2023). Anderson Acceleration software for solving fixed point problems and spin-up of seasonally forced ocean and land biogeochemical models, https://doi.org/10.5281/zenodo.17070247.

Additionally, please cite:

Khatiwala, S. (2023). Fast spin-up of geochemical tracers in ocean circulation and climate models, J. Adv. Model. Earth Sys., https://doi.org/10.1029/2022MS003447.

Khatiwala, S. (2024). Efficient spin-up of Earth System Models using Sequence Acceleration, Science Advances, https://doi.org/10.1126/sciadv.adn2839.

**Acknowledgements**:

This work has been funded in part by UK NERC grant NE/W004976/1 as part of the Agile Initiative at the Oxford Martin School, University of Oxford, and Waseda University’s Global Research Center (GRC) Support Program (Project No. GRCkojin-2503).

**Installation (Python version)**

```bash
pip install aa4py
```

(This will install various dependencies such my [pyioutils](https://github.com/samarkhatiwala/pyioutils) module, `numpy` and `scipy`.)

To use the software in parallel, you will additionally need to install my [pympiutils](https://github.com/samarkhatiwala/pympiutils) package:

```bash
pip install pympiutils
```

Note that parallelization depends on [mpi4py](https://mpi4py.readthedocs.io/en/stable/) and [h5py](https://www.h5py.org/) built against a MPI-enabled HDF5 library (instructions coming soon ...).

**Usage (Python version)**

The following is a minimal working example. (See `run_simple_example.py` in the examples directory for a more elaborate version, and`run_poisson_example.py` and `run_poisson_xy_example.py` for more complex examples, including running in parallel.) 

```python 
from aa4py import (
    DotDict,
    andacc,
    AndersonAccelerationObj,
    saveaarestart
)
import numpy as np

# Define the fixed point function
def g(x, y, fetchOutput, iter_):
  gx=np.array([(x[0]**2 + x[1]**2 + 8.)/10., (x[0]*x[1]**2 + x[0] + 8.)/10.])
  return gx, None, None, {}

# AA parameters dictionary
AAparams = DotDict()
AAparams.mMax = 2
AAparams.itmax = 10
AAparams.beta = 1.0

# Instantiate the AA solver object
s = AndersonAccelerationObj({})

# Initial guess
x0=np.array([0.,0.])

# Call the solver
xsol, iter_, ysol, converged = andacc(g, x0, s, AAparams=AAparams, doRestart=False, fileSuff='_test', y=None, debug=True)

print(f"Solution after {iter_} iterations is {xsol}")

# Save the solver state to file (can be loaded with saveaarestart)
fn="aa_restart.h5"
saveaarestart(s,fn)


```

**Spin-up of ocean and land biogeochemical models**

As described in [Khatiwala (2023)](https://doi.org/10.1029/2022MS003447) and [Khatiwala (2024)](https://doi.org/10.1126/sciadv.adn2839), Anderson Acceleration (AA) can be applied to the spin-up of seasonally-forced biogeochemical models by posing it as a fixed point iteration, x<sub>n+1</sub>  = g(x<sub>n</sub>), where x<sub>n</sub> is the model state vector at iteration n and g the model time-stepper that integrates the the model for one forcing period (typically a year) and returns the solution. AA thus treats the model as a "black box". To apply this software to the spin-up problem, you must provide a model-specific interface ("wrapper") that writes the vector x (mapping appropriately between the vector and 2-d/3-d fields) to the model's restart file, runs the model, and reads from the model's restart file and returns the solution vector g(x).

Wrappers have been developed for several ocean (NEMO-PISCES, NEMO-MEDUSA, CESM-MARBL-MOM6, MITgcm-BLING, MITgcm-dic, UVic-MOBI and MOPS, NORESM-HAMOCC, MIROC, MRI, etc) and land (JULES, CLM) models. I intend to make these available here (see the`python/wrappers` directory). For now I've provided the [NEMO wrapper](https://github.com/samarkhatiwala/AndersonAcceleration/python/wrappers/NEMO) as an example (see the User Guide cited below) and I'm happy to share the others (just email me).

For more information on applying this software to the spin-up problem and the structure of the wrapper code see this (slightly dated) [User Guide](https://doi.org/10.5281/zenodo.14741088) (many thanks to Andrew Wilson for putting it together!).
