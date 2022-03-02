# Install package to read in problems in CPLEX's LP format
#install.packages("Rglpk")

library("Rglpk")

# Set working directory
setwd("/Users/selinbayramoglu/Desktop/ComputationalMethods/Program")

problem <- Rglpk_read_file("small-example.lp", type="CPLEX_LP")

nvar = problem$objective$nrow   # number of variables

ncons = length(problem$constraints[[3]])    # number of rows

nnz = length(problem$constraints[[1]]$v)  # number of nonzero entries in A

obj_nnz = length(problem$objective$i)   # number of nonzero entries in c

rhs = problem$constraints[[3]]    # rhs vector b

mymatA = rbind(problem$constraints[[1]]$i, problem$constraints[[1]]$j, problem$constraints[[1]]$v)

newmatA = as.vector(t(matrix(mymatA, ncol=1)))

myc = rbind(problem$objective$i, problem$objective$v)

newc = as.vector(myc)

# Create file
cat(nvar, file="outfile.txt",sep="\n")
cat(ncons, file="outfile.txt", sep="\n", append = "TRUE")
cat(nnz, file="outfile.txt", sep="\n", append = "TRUE")
cat(obj_nnz, file="outfile.txt", sep="\n", append = "TRUE")
cat(newc, file="outfile.txt", append = "TRUE") ##nonzero objective coefficients
cat(file="outfile.txt", sep="\n", append = "TRUE")
cat(rhs,file="outfile.txt", append = "TRUE") ##rhs
cat(file="outfile.txt", sep="\n", append = "TRUE")
cat(newmatA,file="outfile.txt", append = "TRUE")


