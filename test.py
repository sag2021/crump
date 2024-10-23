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

#
# Simple test case for CRUMP code. 
# Goal of the test is to confirm that the method achieves a truncation error scaling of E ~ h^2  
# The code is applied to a known test case and the truncation error is estimated
# as the max. abs. difference between the numerical and analytic solutions. The 
# run time is also recorded. For SOR, the run time should have scaling t ~ N^4,
# where N is the number of points in a given dimension.
#

# Standard
import numpy as np
import time

# Local
import crump

def test_case(X,Y,Z,wn=2*np.pi):
  """
    Simple test case. Picking wn = 2*PI*n gives zero normal
    component on the boundary

    Parameters:
    -----------
      X,Y,Z: (nz,ny,nx) 
        Cartesian coordinate mesh. Must span [0,1] in X direction

    Returns:
    -------- 
      U: (4,nz,ny,nx)
        Magnetic field and scalar potential in one array. Non-dimensional units.
        scalar pot: U[3,:,:,:]
        b-field: U[1:3,:,:,:]
     
  """

  # Set parameters
  l  = np.sqrt(2)*wn
  nf = np.sinh(l)*l

  # Output
  U = np.zeros((4,)+X.shape,dtype=crump.CRUMP_FP)

  # Magnetic field
  U[0,:,:,:] = -wn*np.sin(wn*X)*np.cos(wn*Y)*np.cosh(l*(Z-1))/nf
  U[1,:,:,:] = -wn*np.cos(wn*X)*np.sin(wn*Y)*np.cosh(l*(Z-1))/nf
  U[2,:,:,:] =  +l*np.cos(wn*X)*np.cos(wn*Y)*np.sinh(l*(Z-1))/nf 
             
  # Scalar potential 
  U[3,:,:,:] = +np.cos(wn*X)*np.cos(wn*Y)*np.cosh(l*(Z-1))/nf 

  # Ensure mean(U) = 0. For Neumann boundary conditions potential
  # is unique up to a constant. Computing u such that mean(u) = 0 is 
  # how CRUMP resolves this ambiguity, so the test case needs to 
  # match this assumption 
  meanu      = U[3,:,:,:].mean()
  U[3,:,:,:] = U[3,:,:,:] - meanu
  
  return U

# ---------------------------------------------------------------------

def power_law_fit(x,y):
  """
    Fit a power law y = A*x^gamma
      
    Parameters:
    -----------
      x: array
        Array of x values
      y: array
        Array of y values

    Returns:
    --------
      gamma: float
        power-law index 
      A: float
        Coefficent of x^gamma
      ev: lambda
        Lambda that returns y = A*x^gamma for the given parameters

  """
  Lx    = np.log10(x)
  Ly    = np.log10(y)

  p  = np.polyfit(Lx,Ly,1)
  A  = 10.**p[1]
  ev = lambda x : A*x**p[0]
 
  return p[0],A,ev 

# ---------------------------------------------------------------------

# Set coarsest mesh
nshape_base = np.array([21,18,23],dtype=crump.CRUMP_IT)

# Set factor by which number of mesh points increases
sfactors = np.array([1,2,2.5,3,4,5,8],dtype=crump.CRUMP_IT)

#
# Compute solution at different resolutions
# 
dxs = np.zeros(len(sfactors))
Eb  = np.zeros(len(sfactors))
Eu  = np.zeros(len(sfactors))
dt  = np.zeros(len(sfactors))
print("# dx,Emax(B),Emax(u),time[s]")
for i,sfac in enumerate(sfactors):

  # Scale up mesh size 
  nshape = np.array(nshape_base*sfac,dtype=crump.CRUMP_IT)

  # Uniform mesh. Only works when dx=dy=dz
  nz,ny,nx = nshape[::]
  x        = np.linspace(0.,1.,nx) # Must be non-dimensional [0,1] in x 
  dx       = x[1] - x[0]
  y        = np.arange(ny)*dx 
  z        = np.arange(nz)*dx
  Z,Y,X    = np.meshgrid(z,y,x,indexing='ij')
  dxs[i]   = dx

  # Get test case
  U = test_case(X,Y,Z,wn=2.3*np.pi)

  # Solve using CRUMP
  U1      = U.copy()
  t1      = time.time()
  U2,res  = crump.solve(U)
  t2      = time.time()
  dt[i]   = t2-t1

  # Check convergence
  if(not res["success"]):
    print("Warning: not converged",nshape)

  # Get magnetic field
  b1 = U1[:3,::]
  b2 = U2[:3,::]

  # Get scalar potential
  u1 = U1[3,::]
  u2 = U2[3,::]

  # Get max error
  Eb[i] = np.linalg.norm(b1-b2,axis=0).max()
  Eu[i] = np.abs(u1-u2).max()
  
  outl = "{:g}\t{:g}\t{:g}\t{:g}".format(dx,Eb[i],Eu[i],dt[i])
  print(outl)

# Power law fit
gamma_b,Ab,evb = power_law_fit(dxs,Eb)
gamma_u,Au,evu = power_law_fit(dxs,Eu)
gamma_t,At,evt = power_law_fit(1./dxs,dt)

# Indices
print("Index B: ",gamma_b)
print("Index u: ",gamma_u)
print("Index t: ",gamma_t)

