# MATH233_FinalProject
Diffusion on a Y-network using CG and BiCGSTAB in Matlab and C++

To run code in Matlab, start with SciCOmpFinalProjectNorms.m. If all the matlab files are downloaded, the code should be setup to run. To change from CG to BiCGSTAB, go into the SciCompFunc.m file and change first line in the while loop.

To run in C++, be sure to change the location of the classes to be correct before running. The code does not use parallelization, but is built to include OpenMP in case you want to parallelize something down the road. Start with the main file. 
