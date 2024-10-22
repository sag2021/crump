# Overview 

Checkerboard Relaxation Method Used for Potential field reconstructions (CRUMP). Module computes a magnetic potential field in a Cartesian domain from the normal component of B on the six planar boundaries.
The magnetic field is found from the gradient of a scalar potential using a checkerboard successive over-relaxation method for the solution of Laplace's equation. The code has a Python frontend
and Fortran backend. 

A previous version of this code was used in Mastrano, A., Wheatland, M.S., and Gilchrist, S.A.:2018, Solar Physics,293,130 (doi:10.1007/s11207-018-1351-0).
The differences between this version and that one are superficial. The original was designed to be called from IDL, whereas this version
has a Python frontend. 

This code was written as a quick way to perform a specific calculation, so it's not particularly versatile. 

# Code details

The code takes a given magnetic field, B,in a Cartesian box and computes a potential field based on the value of the normal component of B on each of the six boundaries.  
It returns both the magnetic field and the magnetic scalar potential. Everything is computed in non-dimensional units. 

The calculation is performed using a second-order-finite-difference scheme. Successive OverRelaxation (SOR) is used to solve the finite-difference system. 
Hence, the truncation error should scaled as E ~ h^2, where h is the mesh spacing. The run time will scale as t ~ N^4, where N is the number of mesh
points in a given dimension (and therefore N^3 in total). 

The Fortran backend is written in Fortran 2003 and parallelized for shared memory parallel computers using OpenMP. 

# Usage

First compile the shared library. This has been tested under gfortran. 

See the test.py script for how to use the Python frontend. 

# Mesh

The mesh must span the unit box [0,1]x[0,1]x[0,1] and be uniform, i.e. dx=dy=dz. 


