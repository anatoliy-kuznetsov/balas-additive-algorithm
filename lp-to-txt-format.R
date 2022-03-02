install.packages("Rglpk")
library("Rglpk")

setwd("/Users/selinbayramoglu/Desktop/ComputationalMethods/Program")

problem <- Rglpk_read_file("small-example.lp", type="CPLEX_LP")

nvar = problem$objective$nrow

ncons = length(problem$constraints[[3]])

nnz = length(problem$constraints[[1]]$v)

mymat = rbind(problem$constraints[[1]]$i, problem$constraints[[1]]$j, problem$constraints[[1]]$v)

newmat = as.vector(t(matrix(mymat, ncol=1)))

cat(nvar, file="outfile.txt",sep="\n")
cat(ncons, file="outfile.txt", sep="\n", append = "TRUE")
cat(nnz, file="outfile.txt", sep="\n", append = "TRUE")

cat(problem$objective$v,file="outfile.txt", append = "TRUE")
cat(file="outfile.txt", sep="\n", append = "TRUE")
cat(newmat,file="outfile.txt", append = "TRUE")


