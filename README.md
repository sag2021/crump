# Overview 

Checkerboard Relaxation Method Used for Potential field reconstructions (CRUMP). The code computes the current-free (potential) magnetic field in a domain that matches the normal component of a specified magnetic field
on the boundary. The domain is a Cartesian box and the boundary is the six planar faces of the box. For a fixed normal component, the potential field represents a minimum energy state 
for the magnetic field. The potential magnetic field is found from the gradient of a scalar potential using a checkerboard Successive Over-Relaxation (SOR) method for the solution of Laplace's equation. The code has a Python frontend and Fortran backend. This code was written as a quick way to perform a specific calculation, so it's not particularly versatile. 

A previous version of this code was used in Mastrano, Wheatland, and Gilchrist (2018)[^1].
The difference between this version and that one is superficial. The original was designed to be called from IDL, whereas this version
has a Python frontend. 

More information on the method of SOR in general can be found in most text books on numerical methods e.g. Press et al. [^2].

[^1]: Mastrano, A., Wheatland, M.S., and Gilchrist, S.A.:2018, Solar Physics,293,130 (doi:10.1007/s11207-018-1351-0).
[^2]: Press, W.H. Teukolsky, S.A., Vetterling, W.T., Flannery, B.P.: 2007,Numerical Recipes 3rd Edition: The Art of Scientific Computing,Cambridge University Press,ISBN:9780521880688

# Neumann boundary conditions 

When the normal component is specified on all boundaries, the scalar potential, u, is only unique up to a constant. Hence, CRUMP enforces the additional
condition $\mbox{mean}(u) = 0$ in the volume, which resolves this ambiguity. This choice doesn't affect the magnetic field. 

In addition, when B normal is fixed on all the boundaries, an additional compatibility condition comes into play. Specifically, the net magnetic 
flux must balance when integrated over all six faces of the box. If this isn't the case, then $\nabla\cdot\mathbf B = 0$ can't be achieved in the volume. 

# Code details

The code takes a given magnetic field, B,in a Cartesian box and computes a potential field based on the value of the normal component of B on each of the six boundaries.  
It returns both the magnetic field and the magnetic scalar potential. Everything is computed in non-dimensional units. 

The calculation is performed using a second-order-finite-difference scheme. Successive Over-Relaxation (SOR) is used to solve the finite-difference system. 
Hence, the truncation error should scaled as $E \sim h^2$, where $h$ is the mesh spacing. The run time will scale as $t \sim N^4$, where $N$ is the number of mesh points in a given dimension (and therefore N^3 in total). 

The Fortran backend is written in Fortran 2003 and parallelized for shared memory parallel computers using OpenMP. 

# Usage

First compile the shared library. This has been tested under gfortran. 

See the test.py script for how to use the Python frontend and the docstring for crump.solve. 

# Mesh

The mesh spacing is assumed to be uniform, i.e. $dx=dy=dz$. The mesh spacing is not passed as an argument, instead dx is inferred from the array shape
as $dx = 1/(nx-1)$, where nx = shape[-1]. In other words, the domain is always assumed to have a non-dimensional length of 1 in the x direction. This
assumption is pretty restrictive, but the code was designed to solve one particular problem on a uniform mesh, not be a general method.

# REAL and INTEGER types

There are three basic types used by CRUMP: REAL(FP),INTEGER(IT), and INTEGER(C_INT). 
The former two are used in calculations. They may be set in the CRUMP_MOD module. By default they 
are C_DOUBLE and C_INT64_T. Setting these to different types (e.g. single) won't break the Fortran module, but
CRUMP can no longer be called by the Python module: crump.py expects double and int64.

Some integers are also C_INT. These correspond to positions in the options vector and function return types indended to be
called from Python. 




