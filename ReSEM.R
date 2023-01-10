

## covariance matrix

rh<-c(1,0.5,0.5^2,0.5^3,
             0.5,1,0.5,0.5^2,
             0.5^2,0.5,1,0.5,
             0.5^3,0.5^2,0.5,1)
      
#rh <- c(1,0.5,0.5^2,0.5^3,0.5^4,
#        0.5,1,0.5,0.5^2,0.5^3,
 #       0.5^2,0.5,1,0.5,0.5^2,
 #       0.5^3,0.5^2,0.5,1,0.5,
 #       0.5^4,0.5^3,0.5^2,0.5,1)

## generate covariate

covariate <- function(p,n){
  
  sigma <- matrix(p,4,4)
  
  x <- mvrnorm(n, rep(0, 4), sigma)
  return(x)
}


assign.Stra.ReSEM.sampling <- function(num.assign, N, n, W, prob.S, K, Stra.N, Stra.n){
  # Generates rejective sampling indicators.
  #
  # Args:
  #   num.assign: Number of assignments to be generated.
  #   N: Size of finite population.
  #   n: Number of sampled units.
  #   W: Available covariate vector at the sampling stage.
  #   prob.S: Asymptotic acceptance probability for the sampling stage.
  #   K: K Strata
  #   Stra.N: stratum information for population
  #   Stra.n: stratum information for sample
  # Returns: 
  #   Sampling indicators generated under rejective sampling, and the corresponding Mahalanobis distances. 
  #   Each row of the sampling matrix is a N-dimensional vector that represents one sampling result, 
  #   where we use 0 for sampled units, and -1 for units that will not enter the experiment.
 if(ncol(W)>1){
   ind.assign = matrix(-1, nrow = num.assign, ncol = N)
   M.S.all = rep(NA, num.assign)
   
   a.S = qchisq(prob.S, df = ncol(W))
   Stra.W.cov=rep(list(matrix(0,ncol = ncol(W),nrow = ncol(W))),K)
   res=matrix(0,ncol(W),ncol(W))
   for(i in 1:K){
     Stra.W.cov[[i]]=cov(W[Stra.ind[[i]],])
     res=res+Stra.W.cov[[i]]*Stra.N[i]*(Stra.N[i]-Stra.n[i])/(N^2*Stra.n[i])
   }
   S2.W.inv = solve(res)
   W.bar = colMeans(W)
   
   k = 1
   while(k <= num.assign){
     # selected = sample(c(1:N), n)
     selected=c()
     for(i in 1:K){selected=c(selected, sample(Stra.ind[[i]],Stra.n[i]))
     }
     tau.hat.W = colMeans( W[selected, ] ) - W.bar
     M.S = as.numeric(  t(tau.hat.W) %*% S2.W.inv %*% tau.hat.W )
     if(M.S <= a.S){
       ind.assign[k, selected] = 0
       M.S.all[k] = M.S
       k = k+1
     }
   }
   return(list(ind.assign=ind.assign, M.S.all=M.S.all))
   
 }
  
  
  
  
  if(ncol(W)==1){
   ind.assign = matrix(-1, nrow = num.assign, ncol = N)
  M.S.all = rep(NA, num.assign)
 
  a.S = qchisq(prob.S, df = ncol(W))
  Stra.W.cov=rep(list(matrix(0,ncol = ncol(W),nrow = ncol(W))),K)
  res=matrix(0,ncol(W),ncol(W))
  for(i in 1:K){
    Stra.W.cov[[i]]=var(W[Stra.ind[[i]],])
    res=res+Stra.W.cov[[i]]*Stra.N[i]*(Stra.N[i]-Stra.n[i])/(N^2*Stra.n[i])
  }
  S2.W.inv = solve(res)
   W.bar = mean(W)
  
  k = 1
  while(k <= num.assign){
    # selected = sample(c(1:N), n)
    selected=c()
    for(i in 1:K){selected=c(selected, sample(Stra.ind[[i]],Stra.n[i]))
    }
    tau.hat.W = mean( W[selected, ] ) - W.bar
    M.S = as.numeric(  t(tau.hat.W) %*% S2.W.inv %*% tau.hat.W )
      if(M.S <= a.S){
      ind.assign[k, selected] = 0
      M.S.all[k] = M.S
      k = k+1
    }
  }
  return(list(ind.assign=ind.assign, M.S.all=M.S.all))
  }
  }







assign.Stra.ReSEM.treatment <- function(n, n1, X, prob.T,K,Stra.n,Stra.n1) {
  # Generates one treatment assignment under rerandomization.
  #
  # Args:
  #   n: Number of units to assign.
  #   n1: Number of treated units.
  #   X: Available covariate vector at the treatment assignment stage.
  #   prob.T: Asymptotic acceptance probability for the treatment assignment stage.
  #
  # Returns: 
  #   Treatment assignments under rerandomization, and the corresponding Mahalanobis distances. 
  #   The assignment vector will be a n-dimensional vector that represents one assignment, where 
  #   we use 1 for treated units, and 0 for control units.
  ind.assign = rep(0, n)
  a.T = qchisq(prob.T, df = ncol(X))
  repeat {s=0
  treated=c()
  for(i in 1:K){
    treated=c(treated,sample(c(1:Stra.n[i]),Stra.n1[i])+s)
    s=Stra.n[i]+s
  }
  
  Stra.X.cov=rep(list(matrix(0,ncol = ncol(X),nrow = ncol(X))),K)
  tmp=0
  res=matrix(0,ncol(X),ncol(X))
  for(i in 1:K){
    Stra.X.cov[[i]]=cov(X[c(1:Stra.n[i])+tmp,])
    tmp=tmp+Stra.n[i]
    res=res+Stra.X.cov[[i]]*Stra.N[i]^2*Stra.n[i]/(N^2*Stra.n1[i]*(Stra.n[i]-Stra.n1[i]))
  }
  S2.X.inv = solve(res)
  
  
    tau.hat.X = colMeans(X[treated, ]) - colMeans(X[-treated, ])
    
    M.T = as.numeric( t(tau.hat.X) %*% S2.X.inv %*% tau.hat.X )
    
    
    if (M.T <= a.T) {
      break
    }
  }
  ind.assign[treated] = 1
  return(list(ind.assign=ind.assign, M.T=M.T))
}


obs.outcome <- function(assignment, Y1, Y0){
  # Computes observed outcomes using potential outcomes and treatment assignment.
  #
  # Args:
  #   assignment: Treatment assignment vector, which contains -1, 0, 1 values. We use 1 for treated 
  #               units, 0 for control units, and -1 for units that will not enter the experiment.  
  #   Y1: Potential outcomes under treatment.
  #   Y0: Potential outcomes under control.
  #
  # Returns:
  #   Observed outcomes for all units. If a unit will not enter the experiment, its observed outcome 
  #   will be NA.
  N = length(Y1)
  Y = rep(NA, N)
  Y[assignment == 1] = Y1[assignment == 1]
  Y[assignment == 0] = Y0[assignment == 0]
  return(Y)
}







Weighted.diff.in.means <- function(assignment, Y, Stra.N, N, K){
 Stra.ind = list()
 for(i in 1:K){
   Stra.ind[[i]]=c(1:Stra.N[i])+sum(Stra.N[1:i])-Stra.N[i]
 }
  tmp = rep(NA,K)
  res=0
  for(i in 1:K){
    S=rep(NA,Stra.N[i])
    As=rep(NA,Stra.N[i])
    S=Y[Stra.ind[[i]]]  
    As=assignment[Stra.ind[[i]]]
    tmp[i]=mean(S[As==1])-mean(S[As==0])
  res=tmp%*%(Stra.N/N)
  }
  return(res)
}



Weighted.diff.regfix <- function(assignment, Y, C, E, Stra.N, N, K){
  # Calculates the regression-adjusted difference-in-means estimator for the average treatment effect.
  #
  # Args:
  #   assignment: Treatment assignment vector, which contains -1, 0, 1 values. We use 1 for treated 
  #               units, 0 for control units, and -1 for units that will not enter the experiment. 
  #   Y: Observed outcomes for all units. Will be NA for units that will not enter the experiment.
  #   C: Available covariate vector for sampled units, at the analysis stage.
  #   E: Available covariate vector for all units, at the analysis stage.
  #
  # Returns:
  #   Regression-adjusted difference-in-means estimator, with estimated adjustment coefficients.
  
   
  ##comput  delta.E and  tau.C and beta and gamma
  V_ee=matrix(rep(0,ncol(E)*ncol(E)),nrow=ncol(E),ncol=ncol(E))
  v_cc=matrix(rep(0,ncol(C)*ncol(C)),nrow=ncol(C),ncol=ncol(C))
  delta.E=matrix(rep(0,ncol(E)),nrow=1)
  tau.C=matrix(rep(0,ncol(C)),nrow=1)
  beta=matrix(rep(0,ncol(C)),ncol =1)
  gamma=matrix(rep(0,ncol(E)),ncol=1)
  
  for(i in 1:K){
    S=matrix(NA,nrow=Stra.N[i],ncol=ncol(E))
    R=matrix(NA,nrow=Stra.N[i],ncol=ncol(C))
    Z=rep(NA,Stra.N[i])
    As=rep(NA,Stra.N[i])
    S=E[Stra.ind[[i]],] 
    R=C[Stra.ind[[i]],]
    Z=Y[Stra.ind[[i]]]
    As=assignment[Stra.ind[[i]]]
    delta.E=(colMeans(S[As>=0,])-colMeans(S))*Stra.Pi[i]+delta.E
    tau.C=(colMeans(R[As==1,])-colMeans(R[As==0,]))*Stra.Pi[i]+tau.C
    
    v_cc=v_cc+Stra.weight[i]/(Stra.e1[i]*Stra.e0[i])*cov(R[As>=0,])
    V_ee=V_ee+Stra.weight[i]*(1-Stra.f[i])*cov(S)
    
    beta=beta+Stra.weight[i]*((cov(R[As==1,],Z[As==1]))/Stra.e1[i]+cov(R[As==0,],Z[As==0])/Stra.e0[i])
    gamma=gamma+Stra.weight[i]*(1-Stra.f[i])*(cov(S[As==1,],Z[As==1])-cov(S[As==0,],Z[As==0]))
  }
  beta= solve(v_cc,beta)
  gamma=solve(V_ee,gamma)
  fix= sum(beta*t(tau.C))+sum(gamma*t(delta.E))
  return (fix)
}




CI.SReSEM <- function(ht,assignment, Y, X, W, C, E, alpha=0.05, design){ 
if(ncol(W)==1){
  V.tt=0
  V.tw=matrix(0,nrow = 1,ncol = ncol(W))
  V.ww=matrix(0,nrow = ncol(W),ncol = ncol(W))
  V.tx=matrix(0,nrow = 1,ncol = ncol(X))
  V.xx=matrix(0,nrow = ncol(X),ncol = ncol(X))
  V.te=matrix(0,nrow = 1,ncol = ncol(E))
  V.ee=matrix(0,nrow = ncol(E),ncol = ncol(E))
  V.tc=matrix(0,nrow = 1,ncol = ncol(C))
  V.cc=matrix(0,nrow = ncol(C),ncol = ncol(C))
  
  
  for(k in 1:K){
    Wk=W[Stra.ind[[k]]]
    Xk=X[Stra.ind[[k]],]
    Ek=E[Stra.ind[[k]],]
    Ck=C[Stra.ind[[k]],]
    
    Yk=Y[Stra.ind[[k]]]
    As=assignment[Stra.ind[[k]]]
    
    s.1Ck = cov(Yk[As==1], Ck[As==1,])
    s.0Ck = cov(Yk[As==0], Ck[As==0,])
    s2.Ck.1 = var(Ck[As == 1, ])
  # s2.Ck.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.1))
   
    s2.Ck.0 = var(Ck[As == 0, ])
   # s2.Ck.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.0))
    
   # s2.taumidCk = t( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) ) %*% ( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) )
    
    s.1Wk=cov(Yk[As==1],Wk[As==1])
    s.0Wk=cov(Yk[As==0],Wk[As==0])
    
    s.1Xk=cov(Yk[As==1],Xk[As==1,])
    s.0Xk=cov(Yk[As==0],Xk[As==0,])
    
    s.1Ek=cov(Yk[As==1],Ek[As==1,])
    s.0Ek=cov(Yk[As==0],Ek[As==0,])
    
    
  V.tw=V.tw+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*(s.1Wk-s.0Wk)
  V.ww=V.ww+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*var(Wk)
  
  V.tx=V.tx+Stra.Pi[k]^2/Stra.pi[k]*(s.1Xk/Stra.e1[k]+s.0Xk/Stra.e0[k])
  V.xx=V.xx+Stra.Pi[k]^2/Stra.pi[k]*cov(Xk[As>=0,])/(Stra.e1[k]*Stra.e0[k])
  
  V.te=V.te+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*(s.1Ek-s.0Ek)
  V.ee=V.ee+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Ek)
  
  V.tc=V.tc+Stra.Pi[k]^2/Stra.pi[k]*(s.1Ck/Stra.e1[k]+s.0Ck/Stra.e0[k])
  V.cc=V.cc+Stra.Pi[k]^2/Stra.pi[k]*cov(Ck[As>=0,])/(Stra.e1[k]*Stra.e0[k])
  
  
  V.tt =V.tt+Stra.Pi[k]^2/Stra.pi[k]*(var(Yk[As==1])/Stra.e1[k]+var(Yk[As==0])/Stra.e0[k])
  }

  V.wt=t(V.tw)
  V.xt=t(V.tx)
  V.et=t(V.te)
  V.ct=t(V.tc)
  
  R2.S =V.tw%*%solve(V.ww,V.wt)/V.tt
    
  R2.T =V.tx%*%solve(V.xx,V.xt)/V.tt  
   
  R2.E =V.te%*%solve(V.ee,V.et)/V.tt
  
  R2.C =V.tc%*%solve(V.cc,V.ct)/V.tt
    
  
  
  
  if (design == "SRSE") {
    QR = c(qnorm(alpha/2), qnorm(1-alpha/2))
    CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "S") {
    QR = quantile( ifelse(R2.S<=1, sqrt(1-R2.S), 0) * epsilon + sqrt(R2.S) * L.J.aS,
                   probs = c(alpha/2, 1-alpha/2) ) 
    CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "T") {
    QR = quantile( ifelse(R2.T<=1, sqrt(1-R2.T), 0) * epsilon + sqrt(R2.T) * L.K.aT,
                   probs = c(alpha/2, 1-alpha/2) )
    CI =ht + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  
  if (design == "ST") {
    QR = quantile( ifelse(R2.S+R2.T<=1, sqrt(1-R2.S-R2.T), 0) * epsilon 
                   + sqrt(R2.S) * L.J.aS + sqrt(R2.T) * L.K.aT, 
                   probs = c(alpha/2, 1-alpha/2) ) 
    CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  
  if (design == "adjusted") {
    QR = ifelse(R2.E+R2.C<=1, sqrt(1-R2.E-R2.C), 0) * c(qnorm(alpha/2), qnorm(1-alpha/2))
    CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  return(CI)
}
  if(ncol(W)>1){
    V.tt=0
    V.tw=matrix(0,nrow = 1,ncol = ncol(W))
    V.ww=matrix(0,nrow = ncol(W),ncol = ncol(W))
    V.tx=matrix(0,nrow = 1,ncol = ncol(X))
    V.xx=matrix(0,nrow = ncol(X),ncol = ncol(X))
    V.te=matrix(0,nrow = 1,ncol = ncol(E))
    V.ee=matrix(0,nrow = ncol(E),ncol = ncol(E))
    V.tc=matrix(0,nrow = 1,ncol = ncol(C))
    V.cc=matrix(0,nrow = ncol(C),ncol = ncol(C))
    
    
    for(k in 1:K){
      Wk=W[Stra.ind[[k]],]
      Xk=X[Stra.ind[[k]],]
      Ek=E[Stra.ind[[k]],]
      Ck=C[Stra.ind[[k]],]
      
      Yk=Y[Stra.ind[[k]]]
      As=assignment[Stra.ind[[k]]]
      
      s.1Ck = cov(Yk[As==1], Ck[As==1,])
      s.0Ck = cov(Yk[As==0], Ck[As==0,])
      s2.Ck.1 = var(Ck[As == 1, ])
     # s2.Ck.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.1))
      
      s2.Ck.0 = var(Ck[As == 0, ])
     # s2.Ck.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.0))
      
     #s2.taumidCk = t( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) ) %*% ( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) )
      
      s.1Wk=cov(Yk[As==1],Wk[As==1,])
      s.0Wk=cov(Yk[As==0],Wk[As==0,])
      
      s.1Xk=cov(Yk[As==1],Xk[As==1,])
      s.0Xk=cov(Yk[As==0],Xk[As==0,])
      
      s.1Ek=cov(Yk[As==1],Ek[As==1,])
      s.0Ek=cov(Yk[As==0],Ek[As==0,])
      
      
      V.tw=V.tw+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*(s.1Wk-s.0Wk)
      V.ww=V.ww+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Wk)
      
      V.tx=V.tx+Stra.Pi[k]^2/Stra.pi[k]*(s.1Xk/Stra.e1[k]+s.0Xk/Stra.e0[k])
      V.xx=V.xx+Stra.Pi[k]^2/Stra.pi[k]*cov(Xk[As>=0,])/(Stra.e1[k]*Stra.e0[k])
      
      V.te=V.te+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*(s.1Ek-s.0Ek)
      V.ee=V.ee+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Ek)
      
      V.tc=V.tc+Stra.Pi[k]^2/Stra.pi[k]*(s.1Ck/Stra.e1[k]+s.0Ck/Stra.e0[k])
      V.cc=V.cc+Stra.Pi[k]^2/Stra.pi[k]*cov(Ck[As>=0,])/(Stra.e1[k]*Stra.e0[k])
      
      
      V.tt =V.tt+Stra.Pi[k]^2/Stra.pi[k]*(var(Yk[As==1])/Stra.e1[k]+var(Yk[As==0])/Stra.e0[k])
    }
    
    V.wt=t(V.tw)
    V.xt=t(V.tx)
    V.et=t(V.te)
    V.ct=t(V.tc)
    
    R2.S =V.tw%*%solve(V.ww,V.wt)/V.tt
    
    R2.T =V.tx%*%solve(V.xx,V.xt)/V.tt  
    
    R2.E =V.te%*%solve(V.ee,V.et)/V.tt
    
    R2.C =V.tc%*%solve(V.cc,V.ct)/V.tt
    
    
    
    
    if (design == "SRSE") {
      QR = c(qnorm(alpha/2), qnorm(1-alpha/2))
      CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
    }
    if (design == "S") {
      QR = quantile( ifelse(R2.S<=1, sqrt(1-R2.S), 0) * epsilon + sqrt(R2.S) * L.J.aS,
                     probs = c(alpha/2, 1-alpha/2) ) 
      CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
    }
    if (design == "T") {
      QR = quantile( ifelse(R2.T<=1, sqrt(1-R2.T), 0) * epsilon + sqrt(R2.T) * L.K.aT,
                     probs = c(alpha/2, 1-alpha/2) )
      CI =ht + sqrt(V.tt) * unname(QR) / sqrt(n)
    }
    
    if (design == "ST") {
      QR = quantile( ifelse(R2.S+R2.T<=1, sqrt(1-R2.S-R2.T), 0) * epsilon 
                     + sqrt(R2.S) * L.J.aS + sqrt(R2.T) * L.K.aT, 
                     probs = c(alpha/2, 1-alpha/2) ) 
      CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
    }
    
    if (design == "adjusted") {
      QR = ifelse(R2.E+R2.C<=1, sqrt(1-R2.E-R2.C), 0) * c(qnorm(alpha/2), qnorm(1-alpha/2))
      CI = ht + sqrt(V.tt) * unname(QR) / sqrt(n)
    }
    return(CI)
  }
  }
  






generate.constrained.Gaussian <- function(num, K, a){
  # Generates constrained Gaussian random variables.
  #
  # Args:
  #   num: number of constrained Gaussian random variables to generate.
  #   K: dimension of the Gaussian vector used to generate constrained Gaussian.
  #   a: threshold used to generate constrained Gaussian.
  #
  # Returns:
  #   Constrained Gaussian random variables.
  L.K.a = rep(NA, num)
  for (i in 1:num) {
    repeat{
      D = rnorm(K)
      if (sum(D^2) <= a) {
        break
      }
    }
    L.K.a[i] = D[1]
  }
  return(L.K.a) 
  }

Rprint <- function( Y, X, W, C, E){ 
  if(ncol(W)==1){
    V.tt=0
    V.tw=matrix(0,nrow = 1,ncol = ncol(W))
    V.ww=matrix(0,nrow = ncol(W),ncol = ncol(W))
    V.tx=matrix(0,nrow = 1,ncol = ncol(X))
    V.xx=matrix(0,nrow = ncol(X),ncol = ncol(X))
    V.te=matrix(0,nrow = 1,ncol = ncol(E))
    V.ee=matrix(0,nrow = ncol(E),ncol = ncol(E))
    V.tc=matrix(0,nrow = 1,ncol = ncol(C))
    V.cc=matrix(0,nrow = ncol(C),ncol = ncol(C))
    
    
    for(k in 1:K){
      Wk=W[Stra.ind[[k]]]
      Xk=X[Stra.ind[[k]],]
      Ek=E[Stra.ind[[k]],]
      Ck=C[Stra.ind[[k]],]
      
      Y1k=Y1[Stra.ind[[k]]]
      Y0k=Y0[Stra.ind[[k]]]
      tauk=tau[Stra.ind[[k]]]
      
        #s2.Ck.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.1))
      
       #s2.Ck.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.0))
      
      # s2.taumidCk = t( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) ) %*% ( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) )
      
  
      S.1Xk=cov(Y1k,Xk)
      S.0Xk=cov(Y0k,Xk)
      
      S.1Ck=cov(Y1k,Ck)
      S.0Ck=cov(Y0k,Ck)
      
      
      V.tw=V.tw+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(tauk,Wk)
      V.ww=V.ww+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*var(Wk)
      
      V.tx=V.tx+Stra.Pi[k]^2/Stra.pi[k]*(S.1Xk/Stra.e1[k]+S.0Xk/Stra.e0[k])
      V.xx=V.xx+Stra.Pi[k]^2/Stra.pi[k]*cov(Xk)/(Stra.e1[k]*Stra.e0[k])
      
      V.te=V.te+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(tauk,Ek)
      V.ee=V.ee+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Ek)
      
      V.tc=V.tc+Stra.Pi[k]^2/Stra.pi[k]*(S.1Ck/Stra.e1[k]+S.0Ck/Stra.e0[k])
      V.cc=V.cc+Stra.Pi[k]^2/Stra.pi[k]*cov(Ck)/(Stra.e1[k]*Stra.e0[k])
      
      
      V.tt =V.tt+Stra.Pi[k]^2/Stra.pi[k]*(var(Y1k)/Stra.e1[k]+var(Y0k)/Stra.e0[k]-Stra.f[k]*var(tauk))
    }
    
    V.wt=t(V.tw)
    V.xt=t(V.tx)
    V.et=t(V.te)
    V.ct=t(V.tc)
    
    R2.S =V.tw%*%solve(V.ww,V.wt)/V.tt
    
    R2.T =V.tx%*%solve(V.xx,V.xt)/V.tt  
    
    R2.E =V.te%*%solve(V.ee,V.et)/V.tt
    
    R2.C =V.tc%*%solve(V.cc,V.ct)/V.tt
    
    R2=list(R2.S=R2.S, R2.T=R2.T, R2.E=R2.E, R2.C=R2.C,V.tt=V.tt)
    return(R2)
  }
  if(ncol(W)>1){V.tt=0
  V.tw=matrix(0,nrow = 1,ncol = ncol(W))
  V.ww=matrix(0,nrow = ncol(W),ncol = ncol(W))
  V.tx=matrix(0,nrow = 1,ncol = ncol(X))
  V.xx=matrix(0,nrow = ncol(X),ncol = ncol(X))
  V.te=matrix(0,nrow = 1,ncol = ncol(E))
  V.ee=matrix(0,nrow = ncol(E),ncol = ncol(E))
  V.tc=matrix(0,nrow = 1,ncol = ncol(C))
  V.cc=matrix(0,nrow = ncol(C),ncol = ncol(C))
  
  
  for(k in 1:K){
    Wk=W[Stra.ind[[k]],]
    Xk=X[Stra.ind[[k]],]
    Ek=E[Stra.ind[[k]],]
    Ck=C[Stra.ind[[k]],]
    
    Y1k=Y1[Stra.ind[[k]]]
    Y0k=Y0[Stra.ind[[k]]]
    tauk=tau[Stra.ind[[k]]]
    
    #s2.Ck.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.1))
    
    #s2.Ck.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.Ck.0))
    
    # s2.taumidCk = t( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) ) %*% ( s2.Ck.1.neg.half.power %*% t(s.1Ck) - s2.Ck.0.neg.half.power %*% t(s.0Ck) )
    
    
    S.1Xk=cov(Y1k,Xk)
    S.0Xk=cov(Y0k,Xk)
    
    S.1Ck=cov(Y1k,Ck)
    S.0Ck=cov(Y0k,Ck)
    
    
    V.tw=V.tw+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(tauk,Wk)
    V.ww=V.ww+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Wk)
    
    V.tx=V.tx+Stra.Pi[k]^2/Stra.pi[k]*(S.1Xk/Stra.e1[k]+S.0Xk/Stra.e0[k])
    V.xx=V.xx+Stra.Pi[k]^2/Stra.pi[k]*cov(Xk)/(Stra.e1[k]*Stra.e0[k])
    
    V.te=V.te+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(tauk,Ek)
    V.ee=V.ee+Stra.Pi[k]^2/Stra.pi[k]*(1-Stra.f[k])*cov(Ek)
    
    V.tc=V.tc+Stra.Pi[k]^2/Stra.pi[k]*(S.1Ck/Stra.e1[k]+S.0Ck/Stra.e0[k])
    V.cc=V.cc+Stra.Pi[k]^2/Stra.pi[k]*cov(Ck)/(Stra.e1[k]*Stra.e0[k])
    
    
    V.tt =V.tt+Stra.Pi[k]^2/Stra.pi[k]*(var(Y1k)/Stra.e1[k]+var(Y0k)/Stra.e0[k]-Stra.f[k]*var(tauk))
  }
  
  V.wt=t(V.tw)
  V.xt=t(V.tx)
  V.et=t(V.te)
  V.ct=t(V.tc)
  
  R2.S =V.tw%*%solve(V.ww,V.wt)/V.tt
  
  R2.T =V.tx%*%solve(V.xx,V.xt)/V.tt  
  
  R2.E =V.te%*%solve(V.ee,V.et)/V.tt
  
  R2.C =V.tc%*%solve(V.cc,V.ct)/V.tt
  
  R2=list(R2.S=R2.S, R2.T=R2.T, R2.E=R2.E, R2.C=R2.C)
  return(R2)
  }
  }


