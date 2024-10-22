# Overview 

Checkerboard Relaxation Method Used for Potential field reconstructions (CRUMP).
Module computes a magnetic potential field in a Cartesian domain from the normal component of B on the six planar boundaries.
The magnetic field is found from the gradient of a scalar potential using a checkerboard successive over-relaxation method for the 
solution of Laplace's equation. 

# Code details

The code takes a given magnetic field, B,in a Cartesian box and computes a potential field based on the value of the
normal component of B on each of the six boundaries. The output is written to file (except when called using IDL).  

The part of the code responsible for the computation is written in Fortran2003, the main program, which acts as an interface, is written in C to facilitate command 
line input and compatibility with IDL. 

The code is parallelised for shared memory parallel computers using OpenMP. 

# Compilation 

