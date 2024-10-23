# BSD-3-Clause 
#
# Copyright 2024 S.A Gilchrist
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors 
# may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A !PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT !LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#

import numpy as np
from numpy import ctypeslib as nct
import ctypes

# Float/int 
CRUMP_FP  = np.float64
CRUMP_IT  = np.int64
CRUMP_ITC = ctypes.c_int64

def solve(uin,tol=1e-13,nmax=None,libpath='crumpf.so'):
  """
    Compute potential magnetic field and scalar potential. 
    Calculation is performed in place.

    Parameters:
    -----------
      usol:(4,nz,ny,nx)
         Solution array that combines the magnetic field and space for
         the scalar potential

         b = uin[1:3,:,:,:] 
         u = uin[3  ,:,:,:]

         The contents of u is not used, but the solution for the potential
         is placed there. Boundary conditions for b are extracted from b.

      tol: float
        Convergence tolerance. Should be close to machine eps
      nmax: int
        Max number of iterations

      Returns:
      --------
        uin: (4,nz,ny,nx)

         b = uin[1:3,:,:,:] 
         u = uin[3  ,:,:,:]

        res: dict
            success: bool that indicates if tol was reached
          tol: Calculation ends when max diff. between iterations is less than this
            delta: Change in scalar potential at last iteration. Should be less than tol if converged
          nmax: Max. number of iterations set
            liter: Iteration number of last iteration. Iterations go from 1:nmax. If converged, liter < nmax
  """

  # Setting nmax<0 will cause CRUMP to use HUGE
  if(nmax is None):
    nmax = -1

  # Get external lib
  libc = ctypes.cdll.LoadLibrary(libpath)

  # Get length of option vector
  OPT_LEN = libc.get_opt_length()

  # Set argument types
  arg1 = nct.ndpointer(CRUMP_IT,ndim=1,flags=('C','A','W')) 
  arg3 = nct.ndpointer(CRUMP_FP,ndim=1,flags=('C','A','W')) 
  arg4 = nct.ndpointer(CRUMP_FP,ndim=4,flags=('C','A','W')) 

  # [nshape,iopt,ropt,b]
  libc.crump_solve.argtypes = [arg1,arg1,arg3,arg4]

  # Setup
  nshape = np.array(uin.shape[::-1]).astype(CRUMP_IT)
  iopt   = np.zeros(OPT_LEN,dtype=CRUMP_IT)
  ropt   = np.zeros(OPT_LEN,dtype=CRUMP_FP)

  # Set options 
  ropt[libc.get_ropt_tol()]  = tol
  iopt[libc.get_iopt_nmax()] = nmax

  # Call external. Convert CRUMP error code into bool
  ierr = libc.crump_solve(nshape,iopt,ropt,uin)
  ierr = (ierr==libc.get_noerror())

  # Pack parameters into a dict
  keys=["get_ropt_tol","get_ropt_delta","get_iopt_nmax","get_iopt_liter"]
  res = {"success":ierr}
  for key in keys:
    ic        = getattr(libc,key)()
    key2,vec  = key.split("_")[-1:-3:-1]
    if("ropt" in vec):
      res[key2] = ropt[ic]
    else:
      res[key2] = iopt[ic]

  return uin,res

