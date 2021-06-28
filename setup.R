n.counties = 50
N = 100
library(mvtnorm)
library(ecoreg)
coeffs = c(0.05,1,2)
#rho = c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1)
rho = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)

#creates list of matrices, each matrix is the biases for a given covariance rho.
#Row of matrix is a vector of c(muhat, alphahat, betahat) - c(mu, alpha, beta) for a single imulation 
createlist = function(list.length,nsims,cols){
  list.matrix = list()
  for (i in 1:list.length){
    list.matrix[[i]] = matrix(NA,nrow = nsims, ncol = cols)
  }
  return(list.matrix)
}

#takes a list of matrices of the form above and takes means for plotting.
bindmeans = function(list, index){
  vec = c()
  for (i in 1:length(list)){
    vec = c(vec, mean(list[[i]][,index]))
  }
  return(vec)
}

# Simulation function
#Takes in the means of the cts covariate,rho, and creates individual dataframe, as well as county dataframe
simfunc = function(mean1, mean2, N, mu, alpha, beta, rho){
  n.counties = length(mean1)
  aggy = rep(NA, n.counties)
  aggx1 = rep(NA, n.counties)
  aggx2 = rep(NA, n.counties)
  mat = matrix(data = NA, nrow = n.counties*N, ncol = 4)
  
  for (i in 1:n.counties){
    mean.x1 = mean1[i]
    mean.x2 = mean2[i]
    cases = 0
    sum.x1 = 0
    sum.x2 = 0
    
    for (j in 1:N){
      xs = rmvnorm(n = 1, mean = c(mean.x1,mean.x2), sigma = rbind(c(1, rho),c(rho, 1)))
      linear.pred = mu[i] + alpha*xs[1] + beta*xs[2] 
      pij = exp(linear.pred)/(1 + exp(linear.pred))
      case = rbinom(n = 1, size = 1, prob = pij)
      cases = cases + case
      sum.x1 = sum.x1 + xs[1] 
      sum.x2 = sum.x2 + xs[2]
      mat[N*(i-1)+j,1] = i
      mat[N*(i-1)+j,2] = xs[1]
      mat[N*(i-1)+j,3] = xs[2]
      mat[N*(i-1)+j,4] = case
    }
    
    aggy[i] = cases
    aggx1[i] = sum.x1/N
    aggx2[i] = sum.x2/N
  }
  df.idata = as.data.frame(cbind(group = mat[,1], case = mat[,4], indiv.airpollution = mat[,2], indiv.income = mat[,3]))
  aggdata = as.data.frame(cbind(y = aggy, group.airpollution = aggx1, group.income = aggx2))
  return(c(df.idata,aggdata))
}