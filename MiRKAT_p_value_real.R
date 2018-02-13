#!/usr/bin/R


library(MiRKAT, quietly=TRUE)
#> This is vegan 2.4-3
library(GUniFrac, quietly=TRUE)


#myArgs[1] = "/Users/faezeh/Projects/python/notebook/extended_dist_mat_single.txt"
#myArgs[2] = "/Users/faezeh/Projects/python/notebook/health_status.txt"

#p_val_file = "/Users/faezeh/Projects/python/notebook/p-value.txt"
#matrix_file = "/Users/faezeh/Projects/python/notebook/extended_dist_mat_single.txt"
#health_status_file = "/Users/faezeh/Projects/python/notebook/health_status.txt"
#info_file = "/Users/faezeh/Projects/python/notebook/info.txt"
myArgs <- commandArgs(trailingOnly = TRUE)

###################################################
### code chunk number 4: data4
###################################################
set.seed(123)

#Smoker =(HealthStatus == "Smoker") **2
#K.weighted = D2K(D.weighted)
#print(length(myArgs))
matrix_file = myArgs[1]
health_status_file = myArgs[2]
#info_file = myArgs[3]
p_val_file = myArgs[3]
single_p_val_file = myArgs[4]

#health_status = c(0,0,0,1,1,1)
#print(matrix_file)
mcf_distance_matrix = read.table(file = matrix_file, sep = '\t', header = TRUE, row.names = 1)
health_status = scan(file = health_status_file, sep = '\t')

#print(mcf_distance_matrix)
#print(health_status)
D.mcf = as.matrix(mcf_distance_matrix * 0.0001)

D2K <- function(D){
  n <- nrow(D)
  centerM <- diag(n) - 1/n
  K <- -0.5*centerM %*% (D*D) %*% centerM
  eK <- eigen(K, symmetric=TRUE)
  # to ensure that K is positive semi definite we replace negative eigen values with zero
  K <- eK$vector %*% diag(pmax(0,eK$values)) %*% t(eK$vector)
  return(K)
}

K.mcf <- D2K(D.mcf)
#print("D2K end")
nrow(K.mcf)
nrow(t(health_status))
###################################################
### code chunk number 5: data5
###################################################
p_value_permutaion_MiRKAT = MiRKAT(y = health_status, Ks = K.mcf, X = NULL, out_type = "D", method = "permutation")
#p_value_permutation_MMiRKAT = MMiRKAT(Y = t(t(health_status)), K = K.mcf, X = NULL)

#p_value_moment = MiRKAT(y = health_status, Ks = K.mcf, X = NULL, out_type = "D", method = "moment")
#p_value_davies = MiRKAT(y = health_status, Ks = K.mcf, X = NULL, out_type = "D", method = "davies")
#print("p-value = " )
#print(p_value_permutaion)
#write(c(p_value_permutaion_MiRKAT, p_value_permutation_MMiRKAT, DMR), file = p_val_file, ncolumns = 3, append = TRUE, sep = "\t")



write(c(p_value_permutaion_MiRKAT), file = p_val_file, ncolumns = 1, append = TRUE, sep = "\t")
write(p_value_permutaion_MiRKAT, file = single_p_val_file, ncolumns = 1, append = FALSE, sep = "\t")
