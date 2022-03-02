# Balas-Algorithm

Currently, the program accepts inputs as an output of the R program only. The format of this file is as follows

nvar

ncons

nnz (number of nonzeros in A)

c1 c2 ... c_nvar (this part also needs to have the indices of the nonzero objective coefficients, will change it)

i1 j1 A[i1][j1] i2 j2 A[i2][j2] ... i_nnz j_nnz A[i_nnz][j_nnz]

In the algorithm, we must have c[i] >= 0 for all i and all constraints must be of type "<=". This is not implemented yet.


Usage:

./balas outfile.txt
