#install.packages("Rglpk")
library("Rglpk")

setwd("/Users/selinbayramoglu/Desktop/ComputationalMethods/Program")

#args <- commandArgs()
#probname = args[7]

probname = "MPS-format/academictimetablebig.mps"
problem <- Rglpk_read_file(probname, type="MPS_free")


# ========= Transform into Balas form if needed ======

offset = 0

nvar = problem$constraints[[1]]$ncol

ncons = problem$constraints[[1]]$nrow

nnz = length(problem$constraints[[1]]$v)

obj_nnz = length(problem$objective$i)



# Variable transformation
if (obj_nnz > 0) {
  for( n in 1:obj_nnz)
  {
    if( problem$objective$v[n] < 0 )
    {
      myvar = problem$objective$i
      # varlist = c(varlist, myvar)
      offset = offset + problem$objective$v[n]
      problem$objective$v[n] = -1 * problem$objective$v[n] # make cj >= 0
      
      locations = which(problem$constraints[[1]]$j == n) 
      rowindices = problem$constraints[[1]]$i[locations]
      
      problem$constraints[[3]][rowindices] = problem$constraints[[3]][rowindices] - problem$constraints[[1]]$v[locations] # update rhs
      problem$constraints[[1]]$v[locations] = -1 * problem$constraints[[1]]$v[locations] # update A.j (A.j *= -1)
    }
  }
}


# Transformation of ">=" constraints
for( m in 1:ncons )
{
  if(problem$constraints[[2]][m] == ">=") {
    # conslist1 = c(conslist1, m)
    
    problem$constraints[[2]][m] = "<="
    
    locations = which(problem$constraints[[1]]$i == m)
    
    problem$constraints[[3]][m] = -1 * problem$constraints[[3]][m]
    problem$constraints[[1]]$v[locations] = -1 * problem$constraints[[1]]$v[locations]
  }
}

# Transformation of "=" constraints

ctr = 0

for( m in 1:ncons )
{
  if(problem$constraints[[2]][m] == "==") {
    # conslist2 = c(conslist2, m)
    ctr = ctr + 1
    problem$constraints[[2]][m] = "<="
    
    locations = which(problem$constraints[[1]]$i == m)
    start_loc = which(problem$constraints[[1]]$i == m)[1]
    
    problem$constraints[[2]] = c(problem$constraints[[2]], "<=") # append new constraint at the bottom of A
    
    nlocations = sum( problem$constraints[[1]]$i == m )
    locations = which( problem$constraints[[1]]$i == m )
    
    if(nlocations > 0) {
      problem$constraints[[1]]$i = c( problem$constraints[[1]]$i, rep(0, nlocations) ) #extend vector for new row inds
      problem$constraints[[1]]$j = c( problem$constraints[[1]]$j, rep(0, nlocations) ) #extend vector for new col inds
      problem$constraints[[1]]$v = c( problem$constraints[[1]]$v, rep(0, nlocations) ) #extend vector for new nonzero values
      
      problem$constraints[[3]] = c(problem$constraints[[3]],-1 * problem$constraints[[3]][m]) # append new rhs to the end
      problem$constraints[[1]]$i[nnz + 1:nlocations] = rep(ncons + ctr, nlocations)
      problem$constraints[[1]]$j[nnz + 1:nlocations] = problem$constraints[[1]]$j[locations]
      problem$constraints[[1]]$v[nnz + 1:nlocations] = -1 * problem$constraints[[1]]$v[locations]
      nnz = nnz + nlocations # update number of nonzeros, i.e. the lengths of i,j and v vectors
    }
  }
}



# ======== End transform =======



nvar = problem$constraints[[1]]$ncol

ncons = problem$constraints[[1]]$nrow + ctr

nnz = length(problem$constraints[[1]]$v)

obj_nnz = length(problem$objective$i)

rhs = problem$constraints[[3]]

mymatA = rbind(problem$constraints[[1]]$i, problem$constraints[[1]]$j, problem$constraints[[1]]$v)

newmatA = as.vector(t(matrix(mymatA, ncol=1)))

myc = rbind(problem$objective$i, problem$objective$v)

newc = as.vector(myc)

outfile = paste(probname,".txt", sep="")

cat(nvar, file = outfile, sep="\n")
cat(ncons, file = outfile, sep="\n", append = "TRUE")
cat(nnz, file = outfile, sep="\n", append = "TRUE")
cat(obj_nnz, file = outfile, sep="\n", append = "TRUE")
cat(newc, file = outfile, append = "TRUE") ##nonzero objective coefficients
cat(file = outfile, sep="\n", append = "TRUE")
cat(rhs, file = outfile, append = "TRUE") ##rhs
cat(file = outfile, sep="\n", append = "TRUE")
cat(newmatA, file = outfile, append = "TRUE")

