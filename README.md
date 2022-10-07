# Balas-Algorithm

Currently, the program accepts inputs as an output of the R program only. The format of this file is as follows

N (number of variables)

M (number of constraints)

nnz (number of nonzeros in A)

obj_nnz (number of nonzero objective coefficients)

k1 c[k1] k2 c[k2] ... obj_nnz c[obj_nnz] (indices and values for nonzeros coefficients of c)

b[1] b[2] ... b[M]

i1 j1 A[i1][j1] i2 j2 A[i2][j2] ... i_nnz j_nnz A[i_nnz][j_nnz] (indices and values for nonzeros coefficients of A)

In the algorithm, we must have c[i] >= 0 for all i and all constraints must be of type "<=". This is not implemented yet.


Usage:

./balas outfile.txt
