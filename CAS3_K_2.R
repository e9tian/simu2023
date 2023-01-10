# load libraries

library(MASS)

library("Matrix")

lib_loc <- '/data/user95/app/R/lib'

library("expm", lib.loc = lib_loc)

#setwd("/Users/Victor/Desktop/Simu1")
# functions for ReSEM

source("ReSEM.R")

# for reproducibility
set.seed(424)
K=2
N=1000+2000
n=200+400
n1=100+200
n0=n-n1
Stra.N=c(1000,2000)
Stra.n=c(200,400)
Stra.n1=c(100,200)


##compute related things
Stra.Pi=Stra.N/N
Stra.pi=Stra.n/n
Stra.weight=Stra.Pi^2/Stra.pi
Stra.e1=Stra.n1/Stra.n
Stra.e0=1-Stra.e1
Stra.f=Stra.n/Stra.N
##

Stra.ind = list()
for(i in 1:K){
  Stra.ind[[i]]=c(1:Stra.N[i])+sum(Stra.N[1:i])-Stra.N[i]
}

# generate covariate
W = matrix(NA, nrow = N, ncol = 1)
X = matrix(NA, nrow = N, ncol = 2)
E = matrix(NA, nrow = N, ncol = 3)
C = matrix(NA, nrow = N, ncol = 4)

C = covariate(rh,N)

W = as.matrix(C[,1])
X = C[,1:2]
E = C[,1:3]
Stra.D = rt(K,3)

# generate potential outcomes


Y1 = rep(NA,N)
Y0 = rep(NA,N)


for (k in 1:K) {
  b11 <- rt(4,3) 
  b12 <- 0.1*rt(4,3)
  b01 <- b11 + rt(4,3)
  b02 <- b12 + 0.1*rt(4,3)
  for (i in 1:Stra.N[k]){
    j = sum(Stra.N[1:k])-Stra.N[k]+i
    Y1[j] <- t(C[j,])%*%b11+exp(t(C[j,])%*%b12)+Stra.D[k]
    Y0[j] <- t(C[j,])%*%b01+exp(t(C[j,])%*%b02)+Stra.D[k]
  }
}


SNR=5

Y1=Y1+rnorm(N,var(Y1)/SNR)

Y0=Y0+rnorm(N,var(Y0)/SNR) 


tau = Y1-Y0 

# asymptotic acceptance probabilities
p.S = 0.01
p.T = 0.01

# generate constrained Gaussian variables to compute quantile range
num.simu = 10^6

L.J.aS <-   generate.constrained.Gaussian(num.simu, ncol(W), qchisq(p.S, df = ncol(W)))

L.K.aT <-   generate.constrained.Gaussian(num.simu, ncol(X), qchisq(p.T, df = ncol(X)))


epsilon = rnorm(num.simu, mean = 0, sd = 1)

# get treatment assignment from CRSE
iter.max = 10000
assign.CRSE.all <- assign.Stra.ReSEM.sampling(num.assign = iter.max, N, n, W, prob.S=1,K,Stra.N,Stra.n)$ind.assign
for (k in 1:nrow(assign.CRSE.all)) {
  sampled = (assign.CRSE.all[k,] == 0)
  treatment.ind = assign.Stra.ReSEM.treatment(n, n1, X = X[sampled,], prob.T=1,K,Stra.n,Stra.n1)$ind.assign
  assign.CRSE.all[k, sampled] = treatment.ind
}

# get treatment assignment from ReSEM


ind.ReSEM.all.list <- assign.Stra.ReSEM.sampling(num.assign = iter.max, N, n, W, prob.S=p.T,K,Stra.N,Stra.n)$ind.assign
for (k in 1:nrow(ind.ReSEM.all.list)) {
  sampled = (ind.ReSEM.all.list[k,] == 0)
  treatment.ind = assign.Stra.ReSEM.treatment(n, n1, X = X[sampled,], prob.T=p.S,K,Stra.n,Stra.n1)$ind.assign
  ind.ReSEM.all.list[k, sampled] = treatment.ind
}
assign.ReSEM.ST.all <-  ind.ReSEM.all.list



ind.ReSEM.all.list <- assign.Stra.ReSEM.sampling(num.assign = iter.max, N, n, W, prob.S=p.S,K,Stra.N,Stra.n)$ind.assign
for (k in 1:nrow(ind.ReSEM.all.list)) {
  sampled = (ind.ReSEM.all.list[k,] == 0)
  treatment.ind = assign.Stra.ReSEM.treatment(n, n1, X = X[sampled,], prob.T=1,K,Stra.n,Stra.n1)$ind.assign
  ind.ReSEM.all.list[k, sampled] = treatment.ind
}
assign.ReSEM.S.all <- ind.ReSEM.all.list


ind.ReSEM.all.list <- assign.Stra.ReSEM.sampling(num.assign = iter.max, N, n, W, prob.S=1,K,Stra.N,Stra.n)$ind.assign
for (k in 1:nrow(ind.ReSEM.all.list)) {
  sampled = (ind.ReSEM.all.list[k,] == 0)
  treatment.ind = assign.Stra.ReSEM.treatment(n, n1, X = X[sampled,], prob.T=p.T,K,Stra.n,Stra.n1)$ind.assign
  ind.ReSEM.all.list[k, sampled] = treatment.ind
}
assign.ReSEM.T.all <-  ind.ReSEM.all.list




# different estimators

tau.Stra.diff.CRSE=rep(NA,iter.max)
tau.Stra.diff.ReSEM.S=rep(NA,iter.max)
tau.Stra.diff.ReSEM.T=rep(NA,iter.max)
tau.Stra.diff.ReSEM.ST=rep(NA,iter.max)

tau.Stra.diff.CRSE.reg=rep(NA,iter.max)
tau.Stra.diff.ReSEM.S.reg=rep(NA,iter.max)
tau.Stra.diff.ReSEM.T.reg=rep(NA,iter.max)
tau.Stra.diff.ReSEM.ST.reg=rep(NA,iter.max)


for(iter in 1:iter.max){
  tau.Stra.diff.CRSE[iter] = Weighted.diff.in.means(assign.CRSE.all[iter,], Y = obs.outcome(assign.CRSE.all[iter,], Y1, Y0), Stra.N, N, K)
  tau.Stra.diff.ReSEM.S[iter] = Weighted.diff.in.means(assign.ReSEM.S.all[iter,], Y = obs.outcome(assign.ReSEM.S.all[iter,], Y1, Y0), Stra.N, N, K)
  tau.Stra.diff.ReSEM.T[iter] = Weighted.diff.in.means(assign.ReSEM.T.all[iter,], Y = obs.outcome(assign.ReSEM.T.all[iter,], Y1, Y0), Stra.N, N, K)
  tau.Stra.diff.ReSEM.ST[iter] = Weighted.diff.in.means(assign.ReSEM.ST.all[iter,], Y = obs.outcome(assign.ReSEM.ST.all[iter,], Y1, Y0), Stra.N, N, K)
  tau.Stra.diff.CRSE.reg[iter]=tau.Stra.diff.CRSE[iter]-Weighted.diff.regfix(assignment=assign.CRSE.all[iter,], Y=obs.outcome(assign.CRSE.all[iter,], Y1, Y0), C, E, Stra.N, N, K)
  tau.Stra.diff.ReSEM.S.reg[iter]=tau.Stra.diff.ReSEM.S[iter]-Weighted.diff.regfix(assignment=assign.ReSEM.S.all[iter,], Y=obs.outcome(assign.ReSEM.S.all[iter,], Y1, Y0), C, E, Stra.N, N, K)
  tau.Stra.diff.ReSEM.T.reg[iter]=tau.Stra.diff.ReSEM.T[iter]-Weighted.diff.regfix(assignment=assign.ReSEM.T.all[iter,], Y=obs.outcome(assign.ReSEM.T.all[iter,], Y1, Y0), C, E, Stra.N, N, K)
  tau.Stra.diff.ReSEM.ST.reg[iter]=tau.Stra.diff.ReSEM.ST[iter]-Weighted.diff.regfix(assignment=assign.ReSEM.ST.all[iter,], Y=obs.outcome(assign.ReSEM.ST.all[iter,], Y1, Y0), C, E, Stra.N, N, K)
}


## generate CI

CI.diff.SRSE <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.CRSE.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.CRSE[iter],assign.CRSE.all[iter,], Y, X, W, C, E, design = "SRSE")
    }
  CI
  
}

CI.diff.SReSEM.S <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.S.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.S[iter],assign.ReSEM.S.all[iter,], Y, X, W, C, E, design = "S")
  }
  CI
  
}

CI.diff.SReSEM.T <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.T.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.T[iter],assign.ReSEM.T.all[iter,], Y, X, W, C, E, design = "T")
    }
  CI
  
}

CI.diff.SReSEM.ST <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.ST.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.ST[iter],assign.ReSEM.ST.all[iter,], Y, X, W, C, E, design = "ST")
    
  }
  CI
}

CI.diff.SRSE.reg <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.CRSE.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.CRSE.reg[iter],assign.CRSE.all[iter,], Y, X, W, C, E, design = "adjusted")
  }
  CI
}

CI.diff.SReSEM.S.reg <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.S.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.S.reg[iter],assign.ReSEM.S.all[iter,], Y, X, W, C, E, design = "adjusted")
  }
  CI
  
}

CI.diff.SReSEM.T.reg <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.T.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.T.reg[iter],assign.ReSEM.T.all[iter,], Y, X, W, C, E, design = "adjusted")
  }
  CI
  
}

CI.diff.SReSEM.ST.reg <- {
  CI = matrix(NA, nrow=iter.max, ncol=2)
  for (iter in 1:iter.max ) {
    Y = obs.outcome(assign.ReSEM.ST.all[iter,], Y1, Y0)
    CI[iter,] = CI.SReSEM(tau.Stra.diff.ReSEM.ST.reg[iter],assign.ReSEM.ST.all[iter,], Y, X, W, C, E, design = "adjusted")
    
  }
  CI
}


# average length of CIs
average.length.CI.diff.SRSE = mean(CI.diff.SRSE[,2] - CI.diff.SRSE[,1])
average.length.CI.diff.SReSEM.S = mean(CI.diff.SReSEM.S[,2] - CI.diff.SReSEM.S[,1])
average.length.CI.diff.SReSEM.T = mean(CI.diff.SReSEM.T[,2] - CI.diff.SReSEM.T[,1])
average.length.CI.diff.SReSEM.ST = mean(CI.diff.SReSEM.ST[,2] - CI.diff.SReSEM.ST[,1])

average.length.CI.diff.SRSE.reg = mean(CI.diff.SRSE.reg[,2] - CI.diff.SRSE.reg[,1])
average.length.CI.diff.SReSEM.S.reg = mean(CI.diff.SReSEM.S.reg[,2] - CI.diff.SReSEM.S.reg[,1])
average.length.CI.diff.SReSEM.T.reg = mean(CI.diff.SReSEM.T.reg[,2] - CI.diff.SReSEM.T.reg[,1])
average.length.CI.diff.SReSEM.ST.reg = mean(CI.diff.SReSEM.ST.reg[,2] - CI.diff.SReSEM.ST.reg[,1])


# coverage of CIs
coverage.CI.diff.SRSE = mean( mean(tau) > CI.diff.SRSE[,1] & mean(tau) < CI.diff.SRSE[,2] )
coverage.CI.diff.SReSEM.S = mean( mean(tau) > CI.diff.SReSEM.S[,1] & mean(tau) < CI.diff.SReSEM.S[,2] )
coverage.CI.diff.SReSEM.T = mean( mean(tau) > CI.diff.SReSEM.T[,1] & mean(tau) < CI.diff.SReSEM.T[,2] )
coverage.CI.diff.SReSEM.ST = mean( mean(tau) > CI.diff.SReSEM.ST[,1] & mean(tau) < CI.diff.SReSEM.ST[,2] )

coverage.CI.diff.SRSE.reg = mean( mean(tau) > CI.diff.SRSE.reg[,1] & mean(tau) < CI.diff.SRSE.reg[,2] )
coverage.CI.diff.SReSEM.S.reg = mean( mean(tau) > CI.diff.SReSEM.S.reg[,1] & mean(tau) < CI.diff.SReSEM.S.reg[,2] )
coverage.CI.diff.SReSEM.T.reg = mean( mean(tau) > CI.diff.SReSEM.T.reg[,1] & mean(tau) < CI.diff.SReSEM.T.reg[,2] )
coverage.CI.diff.SReSEM.ST.reg = mean( mean(tau) > CI.diff.SReSEM.ST.reg[,1] & mean(tau) < CI.diff.SReSEM.ST.reg[,2] )

#save 


result_C3_K_2<-list(tau.Stra.diff.CRSE=tau.Stra.diff.CRSE,tau.Stra.diff.ReSEM.S=tau.Stra.diff.ReSEM.S,tau.Stra.diff.ReSEM.T=tau.Stra.diff.ReSEM.T,tau.Stra.diff.ReSEM.ST=tau.Stra.diff.ReSEM.ST,tau.Stra.diff.CRSE.reg=tau.Stra.diff.CRSE.reg,tau.Stra.diff.ReSEM.S.reg=tau.Stra.diff.ReSEM.S.reg,tau.Stra.diff.ReSEM.T.reg=tau.Stra.diff.ReSEM.T.reg,tau.Stra.diff.ReSEM.ST.reg=tau.Stra.diff.ReSEM.ST.reg,tau=tau,n=n,average.length.CI.diff.SRSE=average.length.CI.diff.SRSE,average.length.CI.diff.SReSEM.S=average.length.CI.diff.SReSEM.S,average.length.CI.diff.SReSEM.T=average.length.CI.diff.SReSEM.T,average.length.CI.diff.SReSEM.ST=average.length.CI.diff.SReSEM.ST,average.length.CI.diff.SRSE.reg=average.length.CI.diff.SRSE.reg,average.length.CI.diff.SReSEM.S.reg=average.length.CI.diff.SReSEM.S.reg,average.length.CI.diff.SReSEM.T.reg=average.length.CI.diff.SReSEM.T.reg,average.length.CI.diff.SReSEM.ST.reg=average.length.CI.diff.SReSEM.ST.reg,coverage.CI.diff.SRSE=coverage.CI.diff.SRSE,coverage.CI.diff.SReSEM.S=coverage.CI.diff.SReSEM.S,coverage.CI.diff.SReSEM.T=coverage.CI.diff.SReSEM.T,coverage.CI.diff.SReSEM.ST=coverage.CI.diff.SReSEM.ST,coverage.CI.diff.SRSE.reg=coverage.CI.diff.SRSE.reg,coverage.CI.diff.SReSEM.S.reg=coverage.CI.diff.SReSEM.S.reg,coverage.CI.diff.SReSEM.T.reg=coverage.CI.diff.SReSEM.T.reg,coverage.CI.diff.SReSEM.ST.reg=coverage.CI.diff.SReSEM.ST.reg)
save(result_C3_K_2,file = "C3_K_2.RData")
# histogram