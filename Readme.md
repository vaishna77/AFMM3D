"testFMM3Dsolve.cpp" solves the Fredholm integral equation of second kind with 1/r kernel function, using GMRES. Wherein each iteration of GMRES, invoves a matrix-vector product which is performed using AFMM3D that is built using NCA.

"testFMM3DRAMeff2.cpp" performs matrix vector product with the kernel function defined as per the input choice

Both files takes these inputs at run time:

cubeRootN: It determines the system size, that is N=pow(cubeRootN,3)

nParticlesInLeafAlong1D: it determines the maximum number of particles in a leaf node, which is pow(nParticlesInLeafAlong1D,3).

L: half side length of the cube which is the computational domain

TOL_POW: tolerance set for NCA

Qchoice: for various choices look at the kernel.hpp file; Choose Qchoice 17 for "testFMM3Dsolve.cpp".
