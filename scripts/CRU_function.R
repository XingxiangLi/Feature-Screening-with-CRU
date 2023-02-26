######################################################################################
##  This script provides R functions needed to run "CRU_ReplicationCode_Full.R" or  ##
##  "CRU_ReplicationCode_ShortRunTime.R".                                           ##
######################################################################################

## Last update: February 25 2023
## Require R >= 4.2.1

## Packages needed to load this file:
# library(mnormt)   # version 2.1.1
# library(class)    # version 7.3-20


###############################
## Functions to estimate MV  ##
###############################

## Function name: MV_component_matrix 
## Purpose: estimate the component parameters in MV utility.
## Description: a helper function for function "MV_OSFA_Value".

MV_component_matrix = function(X, Y, class){
  
  ## Initialization
  R = length(class) 
  n = length(Y) 
  p = dim(X)[2]  
  
  ## Indicator matrix to estimate component parameters 
  Yheng = matrix(Y, n, R)
  Yshu = t(matrix(class, R, n))
  Ycha = Yheng-Yshu 
  
  MY = matrix(0, n, R)
  MY[which(Ycha == 0)] = 1 
  Nr = colSums(MY) 
  Num = n*(n-1)*(n-2)  
  Res= list()
  
  ## Begin loop for each category
  for (r in 1:R){
    
    VYr = MY[, r]
    MYr = matrix(VYr, n, n)    
    Er = c()
    
    if(Nr[r] == 0){Er = matrix(0, 2, p)}else{ 
      
      for (j in 1:p){
        
        Xj = X[, j]
        MXj = shixingM(Xj)  
        diag(MXj) = 0         
        MXY = MXj*MYr  
        
        ## First component parameter
        Ar = MXY%*%t(MXY)  
        diag(Ar) = 0       
        E1 = sum(Ar)/Num   
        
        ## Second component parameter
        Br = MXY%*%t(MXj)
        diag(Br) = 0
        E2 = sum(Br)/Num   
        
        ## Combine parameters 
        res = c(E1, E2)        
        Er = cbind(Er, res)        
      }
    }
    
    Res[[r]] = Er    
  }
  
  ## Output result
  Res[[R+1]] = Nr
  Res 
}


## Function name: MV_OSFA_Value
## Purpose: distributed estimation of MV via OSFA.
## Description: this function is used to construct Tables 3-7 in Section 4.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector; the input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical response with R classes.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments. 
## Outputs:
##   Value: a vector of estimated MV.
##   TC: computational time (in seconds) needed to carry out the estimation.

MV_OSFA_Value = function(X, Y, M_num){
  
  print("MV_OSFA")
  
  ## Obtain class labels
  m = length(M_num)
  class = c()     
  for (ii in 1:m){
    DY = Y[M_num[[ii]]]
    class_ii = unique(DY)  
    class = union(class, class_ii)  
  }
  
  ## Initialization
  R = length(class)
  p = dim(X)[2]  
  TT = c() 
  Vw = c()
  Pr = rep(0, R)
  theta1 = theta2 = matrix(0, 1, p)
  MCP = lapply(rep(0, R), function(h) return(h)) 
  
  ## Begin estimation
  for (ii in 1:m){

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ] 
    DY = Y[M_num[[ii]]]     
    ptm = proc.time()
    
    ## Compute weights of component parameters
    nii = length(DY)
    w = floor(nii/3)
    Vw = c(Vw, w)
    
    ## Estimate component parameters via local U-statistics
    res_ii = MV_component_matrix(DX, DY, class)
    Pr = Pr + res_ii[[R+1]]
    
    ## Combine U-statistics
    for (r in 1:R){
      MCP[[r]] = MCP[[r]]+res_ii[[r]]*w    
    }
    
    ct = proc.time()-ptm
    TT = c(TT, ct[3])    
  } 
  
  ## Functional aggregation of the estimated component parameters
  ptm = proc.time()  
  MV = matrix(1/3, 1, p)
  n3 = sum(Vw)  
  Pr = Pr/sum(Pr)     
 
  for (r in 1:R){
    aa = MCP[[r]]/n3  
    MV = MV+aa[1, ]/Pr[r]-2*aa[2, ]
  }
  
  ## Output results
  at = proc.time()-ptm  
  TC = mean(TT)+at[3]    
  list(Value = c(MV), TC = TC)
  
}


#################################
## Functions to estimate FAIR  ##
#################################

## Function name: FAIR_index
## Purpose: non-distributed estimation of FAIR for a single feature.
## Description: this function is used to construct Figure 1 in Section 2.1. 
## Inputs:
##   X: an n*1 vector of numerical feature.
##   Y: an n*1 vector of binary response.
## Output:
##   W: the estimated value of FAIR.

FAIR_index = function(X, Y){
  
  ## Obtain class labels
  class = unique(Y)
  weizhi1 = which(Y == class[1]) 
  weizhi2 = which(Y == class[2]) 
  num1 = length(weizhi1)
  num2 = length(weizhi2)
  
  ## Compute conditional mean and variance for class 1
  X1 = X[weizhi1]
  junzhi1 = mean(X1) 
  Cha1 = X1-junzhi1 
  Var1 = sum(Cha1^2)/(num1-1) 
  
  ## Compute conditional mean and variance for class 2
  X2 = X[weizhi2]  
  junzhi2 = mean(X2) 
  Cha2 = X2-junzhi2
  Var2 = sum(Cha2^2)/(num2-1) 
  
  ## Output the estimate 
  W = abs((junzhi1-junzhi2)/(sqrt(Var1/num1+Var2/num2)))   
  W  
}


## Function name: FAIR_OSFA_Value
## Purpose: distributed estimation of FAIR via OSFA.
## Description: this function is used to construct Tables 3, 4, and 6 in Section 4.1.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of binary responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.  
## Outputs:
##   Value: a vector of estimated FAIR.
##   TC: computational time (in seconds) needed to carry out the estimation.

FAIR_OSFA_Value = function(X, Y, M_num){
  
  print("FAIR_OSFA")
  
  ## Obtain class labels
  m = length(M_num)
  class = c()    
  for (ii in 1:m){
    DY = Y[M_num[[ii]]]
    class_ii = unique(DY) 
    class = union(class, class_ii)  
  }
  
  ## Initialization
  R = length(class) 
  p = dim(X)[2]  
  Num = matrix(0, m, R) 
  He = matrix(0, R, p) 
  SHe = matrix(0, R, p) 
  TT = c()
  
  ## Begin estimation
  for (ii in 1:m) {

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ] 
    DY = Y[M_num[[ii]]] 
    
    ## Create working matrices to carry out estimation
    ptm = proc.time()
    nii=length(DY)        
    Yheng = matrix(DY, nii, R)
    Yshu = t(matrix(class, R, nii))
    Ycha = Yheng-Yshu     
    M = matrix(0, nii, R)
    M[which(Ycha == 0)] = 1
    num = colSums(M)   
    Num[ii, ] = num        
    HE_ii = matrix(0, R, p) 
    SHE_ii = matrix(0, R, p) 
    
    ## Estimate component parameters via local U-statistics
    for (jj in 1:R){

      MY = matrix(rep(M[, jj], p), nrow = nii)
      HE_ii[jj, ] = colSums(DX*MY) 
      SHE_ii[jj, ] = colSums((DX*MY)^2)
    }

    ## Combine U-statistics over data segments
    He = He+HE_ii 
    SHe = SHe+SHE_ii 
    ct = proc.time()-ptm
    TT = c(TT, ct[3])
  }
  
  ## Functional aggregation of the estimated component parameters
  ptm = proc.time()
  N12 = colSums(Num) 
  junzhi = He/matrix(rep(N12, 1, p), 2, p) 
  junzhi_s = SHe/matrix(rep(N12, 1, p), 2, p) 
  
  junzhi1 = junzhi[1, ]
  junzhi2 = junzhi[2, ]
  num1 = N12[1]
  num2 = N12[2]

  fangcha1 = (junzhi_s[1, ]-junzhi1^2)*(num1/(num1-1))
  fangcha2 = (junzhi_s[2, ]-junzhi2^2)*(num2/(num2-1))
  W = abs((junzhi1-junzhi2)/(sqrt(fangcha1/num1+fangcha2/num2)))
  at = proc.time()-ptm
  
  ## Output results
  TC = mean(TT)+at[3]
  list(Value = W, TC = TC)
  
}


#################################
## Functions to estimate PSIS  ##
#################################

## Function name: PSIS_OSFA_Value
## Purpose: distributed estimation of PSIS via OSFA. 
## Description: this function is used to construct Tables 5 and 7 in Section 4.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.   
## Outputs:
##   Value: a vector of estimated PSIS.
##   TC: computational time (in seconds) needed to carry out the estimation.

PSIS_OSFA_Value = function(X, Y, M_num){
  
  print("PSIS_OSFA")
  
  ## Obtain class labels
  m = length(M_num)
  class = c()     
  for (ii in 1:m){
    DY = Y[M_num[[ii]]]
    class_ii = unique(DY)
    class = union(class, class_ii) 
  }
  
  ## Initialization
  R = length(class) 
  p = dim(X)[2]
  X = scale(X) 
  Num = matrix(0, m, R) 
  He = matrix(0, R, p) 
  TT = c()
  
  ## Begin estimation
  for (ii in 1:m) {

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ]
    DY = Y[M_num[[ii]]]
    
    ## Create working matrices to carry out estimation
    ptm  =  proc.time()
    nii = length(DY)        
    Yheng = matrix(DY, nii, R)
    Yshu = t(matrix(class, R, nii))
    Ycha = Yheng-Yshu 
    M = matrix(0, nii, R)
    M[which(Ycha ==  0)] = 1
    num = colSums(M)   
    Num[ii, ] = num         
    HE_ii = matrix(0, R, p)
 
    ## Estimate component parameters via local U-statistics
    for (jj in 1:R){
      MY = matrix(rep(M[, jj], p), nrow = nii)
      HE_ii[jj, ] = colSums(DX*MY) 
    }

    ## Combine U-statistics over data segments
    He = He+HE_ii
    ct = proc.time()-ptm
    TT = c(TT, ct[3])    
  }

  ## Functional aggregation of the estimated component parameters  
  ptm = proc.time()
  Nm = colSums(Num) 
  junzhi = He/matrix(rep(Nm, p), R, p) 
  PairCha = c()

  for (i1 in 1:(R-1)){
    for (i2 in (i1+1):R){
      aa = abs(junzhi[i1, ]-junzhi[i2, ])
      PairCha = rbind(PairCha, aa)
    }}
  
  ## Output results
  PSIS = apply(PairCha, 2, max)  
  at = proc.time()-ptm 
  TC = mean(TT)+at[3]
  list(Value = PSIS, TC = TC)
  
}


################################
## Functions to estimate KF   ##
################################

## Function name: KF_Value
## Purpose: non-distributed estimation of KF.
## Description: a helper function for function "KF_SWA_Value".

KF_Value = function(X, Y){ 
  
  ## Obtain sample size and the number of features
  n = length(Y)     
  p = dim(X)[2]  

  ## Obtain class labels
  class = unique(Y) 
  R = length(class) 
  
  if (R ==  1){KF_vec = rep(0, p)}else{
    
    ## Create working matrices
    Yheng = matrix(Y, n, R)
    Yshu = t(matrix(class, R, n))
    Ycha = Yheng-Yshu 
    M = matrix(0, n, R)
    M[which(Ycha == 0)] = 1
    Pr = colMeans(M)     
    KF_vec = c()

    ## Begin estimation
    for (j in 1:p){

      Xj = X[, j]  
      MX = shixingM(Xj)
      FM = matrix(0, R, n)
      
      ## Compute the conditional cdf for each class
      for (ii in 1:R){      
        MY = matrix(M[, ii], n, n)
        MXY = MX*MY
        FXY = colMeans(MXY)/Pr[ii] 
        FM[ii, ] = FXY
      }
      
      ## Compute Kolmogorov-Smirnov test statistic
      KF = max(abs(FM[1, ]-FM[2, ]))
      KF_vec = c(KF_vec, KF)}}
  
  KF_vec
}

## Function name: KF_SWA_Value
## Purpose: distributed estimation of KF via SWA.
## Description: this function is used to construct Tables 3, 4, and 6 in Section 4.1.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of binary responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.  
## Outputs:
##   Value: a vector of estimated KF.
##   TC: computational time (in seconds) needed to carry out the estimation.

KF_SWA_Value = function(X, Y, M_num){
  
  print("KF_SWA")
  
  ## Initialization
  WM = c()
  ns = c()
  TT = c()
  m = length(M_num)
  
  ## Begin estimation
  for (ii in 1:m) {

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ] 
    DY = Y[M_num[[ii]]] 

    ## Compute weights for averaging
    nii = length(DY)
    ns = c(ns, nii)
    
    ## Local KF estimate
    ptm = proc.time()
    Wii = KF_Value(DX, DY)*nii
    ct = proc.time()-ptm
    TT = c(TT, ct[3])   
    WM = rbind(WM, Wii)
  }
  
  ## Output results
  W_KF = colSums(WM)/sum(ns) 
  TC = mean(TT) 
  list(Value = W_KF, TC = TC)
  
}


################################
## Functions to estimate FKF  ##
################################

## Function name: FKF_Value
## Purpose: non-distributed estimation of KF. 
## Description: a helper function for function "FKF_SWA_Value".

FKF_Value = function(X, Y){
  
  ## Obtain sample size and the number of features
  n = length(Y)     
  p = dim(X)[2]

  ## Obtain class labels
  class = unique(Y) 
  R = length(class) 
  
  if (R ==  1){FKF_vec = rep(0, p)}else{
    
    ## Create working matrices
    Yheng = matrix(Y, n, R)
    Yshu = t(matrix(class, R, n))
    Ycha = Yheng-Yshu     
    M = matrix(0, n, R)
    M[which(Ycha ==  0)] = 1
    Pr = colMeans(M)     
    FKF_vec = c()

    ## Begin estimation
    for (j in 1:p){
      Xj = X[, j]  
      MX = shixingM(Xj)
      FM = matrix(0, R, n)

      ## Compute the conditional cdf for each class
      for (ii in 1:R){
        MY = matrix(M[, ii], n, n)
        MXY = MX*MY       
        FXY = colMeans(MXY)/Pr[ii] 
        FM[ii, ] = FXY
      }
      
      ## Compute pairwise Kolmogorov-Smirnov test statistic 
      PairM = c()
      for (i1 in 1:(R-1)){
        for (i2 in (i1+1):R){
          aa = max(abs(FM[i1, ]-FM[i2, ])) 
          PairM = c(PairM, aa)
        }
      }

      ## Final estimation
      KF = max(PairM)
      FKF_vec = c(FKF_vec, KF)
    }
  }

  FKF_vec
}

## Function name: FKF_SWA_Value
## Purpose: distributed estimation of FKF via SWA.
## Description: this function is used to construct Tables 5 and 7 in Section 4.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.   
## Outputs:
##   Value: a vector of estimated FKF.
##   TC: computational time (in seconds) needed to carry out the estimation.

FKF_SWA_Value = function(X, Y, M_num){
  
  print("FKF_SWA")
  
  ## Initialization
  WM = c() 
  ns = c()
  TT = c()
  m = length(M_num)

  ## Begin estimation
  for (ii in 1:m) {

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ]
    DY = Y[M_num[[ii]]] 
    
    ## Compute weights for averaging
    ptm = proc.time()
    nii = length(DY) 
    ns = c(ns, nii)

    ## Local FKF estimate
    Wii = FKF_Value(DX, DY)*nii
    ct = proc.time()-ptm
    TT = c(TT, ct[3])    
    WM = rbind(WM, Wii)
  }
  
  ## Output results
  W_FKF = colSums(WM)/sum(ns)  
  TC = mean(TT)  
  list(Value = W_FKF, TC = TC)
  
}


################################
## Functions to estimate CRU  ##
################################

## Function name: CRU_index
## Purpose: non-distributed estimation of CRU for a single feature.  
## Description: this function is used to construct Figure 1 in Section 2.1. 
## Inputs:
##   X: an n*1 vector of numerical feature.
##   Y: an n*1 vector of categorical response.
## Output:
##   Wj: the estimated CRU between X and Y.

CRU_index = function(X, Y){
  
  ## Obtain sample size and class labels
  n = length(Y)  
  class = unique(Y) 
  R = length(class) 
  
  ## Begin estimation
  if (R ==  1){Wj = 0}else{
     
    ## Create working matrices
    Yheng = matrix(Y, n, R)
    Yshu = t(matrix(class, R, n))
    Ycha = Yheng-Yshu 
    
    ## Estimate the first component parameter
    M = matrix(0, n, R)
    M[which(Ycha ==  0)]=1
    Pr = colMeans(M) 
    
    ## Estimate the second component parameter
    MX = shixingM_CRU(X)
    MXY = t(M)%*%MX
    MXY_mean = rowSums(MXY)/(n*(n-1))
    
    ## Functional aggregation
    Wj = sum((MXY_mean-0.5*Pr)^2)
  }
  
  Wj
}


## Function name: CRU_Value
## Purpose: non-distributed estimation of CRU.
## Description: this function is used to construct Figures 1-2 and Table 2 in Sections 2.1 
##              and 4.1. It is also a helper function for function "CRU_SWA_Value".
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
## Outputs:
##   Value: a vector of estimated CRU.
##   TC: computational time (in seconds) needed to carry out the estimation.

CRU_Value = function(X, Y){
  
  ## Timer
  ptm = proc.time()
  
  ## Obtain sample size and the number of features
  n = length(Y)    
  p = dim(X)[2]
  
  ## Obtain class labels
  class = unique(Y) 
  R = length(class) 
  
  ## Begin estimation
  if (R ==  1){CRU_vec = rep(0, p)}else{
    
    ## Create working matrices
    Yheng = matrix(Y, n, R)
    Yshu = t(matrix(class, R, n))
    Ycha = Yheng-Yshu 
    
    ## Estimate the first component parameter
    M = matrix(0, n, R)
    M[which(Ycha ==  0)] = 1
    Pr = colMeans(M)     
    CRU_vec = c()

    ## Estimate the second component parameter
    for (j in 1:p){
      Xj = X[, j]  
      MX = shixingM_CRU(Xj)
      MXY = t(M)%*%MX
      MXY_mean = rowSums(MXY)/(n*(n-1))      
      Wj = sum((MXY_mean-0.5*Pr)^2)
      CRU_vec = c(CRU_vec, Wj)
    }
  }
  
  ## Output results
  ct = proc.time()-ptm
  TT = ct[3]  
  list(Value = CRU_vec, TC = TT)
}


## Function name: CRU_SWA_Value
## Purpose: distributed estimation of CRU via SWA.
## Description: this function is used to construct Figure 2 in Section 4.1.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.   
## Outputs:
##   Value: a vector of estimated CRU.
##   TC: computational time (in seconds) needed to carry out the estimation.

CRU_SWA_Value = function(X, Y, M_num){
  
  ## Initialization
  WM = c()
  ns = c()
  TT = c()
  m = length(M_num)

  ## Begin estimation
  for (ii in 1:m) {

    ## Create data segments from X and Y
    DX = X[M_num[[ii]], ] 
    DY = Y[M_num[[ii]]] 
    
    ## Compute weights for averaging
    ptm = proc.time()
    nii = length(DY)
    ns = c(ns, nii)

    ## Local CRU estimates
    Wii = CRU_Value(DX, DY)$Value*nii
    WM = rbind(WM, Wii)    
    ct = proc.time()-ptm
    TT = c(TT, ct[3])
  }
    
  ## Output results
  W_CRU = colSums(WM)/sum(ns)
  TC = mean(TT)  
  list(Value = W_CRU, TC = TC)
  
}


## Name: CRU_OSFA_Value
## Purpose: distributed estimation of CRU via OSFA. 
## Description: this function is used to construct Tables 3-7 in Section 4.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.   
## Outputs:
##   Value: a vector of estimated CRU.
##   TC: computational time (in seconds) needed to carry out the estimation.

CRU_OSFA_Value = function(X, Y, M_num){
  
  print("CRU_OSFA")
  
  ## Obtain class labels
  m = length(M_num)
  class = c()     
  for (ii in 1:m){
    DY = Y[M_num[[ii]]]
    class_ii = unique(DY)  
    class = union(class, class_ii)  
  }  

  ## Initialization
  R = length(class)
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## Begin estimation
  if(R ==  1){CRU_vec = rep(0, p); TC = 0}else{
    MCPy = 0 
    MCPj = 0 
    Sw = 0
    TT = c()
    
    for (ii in 1:m){

      ## Create data segments from X and Y
      DX = X[M_num[[ii]], ]
      DY = Y[M_num[[ii]]] 
      
      ## Create working matrices to carry out estimation
      ptm = proc.time()
      nii = length(DY)    
      Yheng = matrix(DY, nii, R)
      Yshu = t(matrix(class, R, nii))
      Ycha = Yheng-Yshu      
      M = matrix(0, nii, R)
      M[which(Ycha ==  0)] = 1
      nr = colSums(M) 

      ## Estimate component parameters via local U-statistics 
      MCPy = MCPy+matrix(rep(nr, p), R, p)/n   
      wm = floor(nii/2)
      Sw = Sw+wm      
      MCPjm = c()

      for (j in 1:p){
        Xj = DX[, j]  
        MX = shixingM_CRU(Xj)
        MXY = t(M)%*%MX   
        MXY_mean = wm*rowSums(MXY)/(nii*(nii-1))
        MCPjm = cbind(MCPjm, MXY_mean)
      }
      
      ## Combine U-statistics over data segments
      MCPj = MCPj+MCPjm
      ct = proc.time()-ptm
      TT = c(TT, ct[3])
      
    }
    
    ## Functional aggregation of the estimated component parameters  
    ptm = proc.time()
    CRU_vec = colSums((MCPj/Sw-.5*MCPy)^2)
    at = proc.time()-ptm
    TC = mean(TT)+at[3]
    
  }
  
  list(Value = CRU_vec, TC = TC)
}


#########################################
## Functions used in data generation   ##
#########################################

## Function name: dis_gen
## Purpose: generate a categorical response from a discrete distribution.
## Description: this function is used to construct Figures 1-2 and Tables 2-5 in Sections 2.1 and 4.1. 
## Inputs:
##   n: number of responses to be simulated.
##   pp: a vector of probabilities specifying the distribution of the categorical response.
## Outputs:
##   Y: the simulated categorical response vector of size n.
##   SR: a list indicating the positions of Y corresponding to the categories.

dis_gen = function(n, pp){ 

  ## Compute the cumulated probabilities to generate Y
  lp = length(pp) 
  cp = c()
  lc = lp-1
  for (i in 1:lc){
    cp = c(cp, sum(pp[1:i]))
  }
  Q = qunif(cp) 
  
  ## Generate the response Y
  x = runif(n) 
  X = x
  X[which(x < Q[1])] = 1
  if (lp > 2){
    for (ii in 2:lc){X[which(x >= Q[ii-1] & x < Q[ii])] = ii}
  }
  X[which(x >= Q[lc])] = lp 
  
  ## Obtain the positions of Y for different categories
  S = list()
  for (r in 1:lp){
    S[[r]] = which(X ==  r)
  }
  
  list(Y = X, SR = S)
}


## Function name: logistic_inde_dg
## Purpose: generate data in setup (d) of Example 3 in Section 4.1.
## Description: this function is used to construct Table 6 in Section 4.1. 
## Inputs:
##   n: number of observations to be generated.
##   p: total number of features.
##   q: number of relevant features.
##   Beta: non-zero coefficients in the model; the default value is specified in setup (d) of Example 3.
##   index: the index set of relevant features; the default value is specified in setup (d) of Example 3.
## Outputs:
##   Y: a binary response vector of length n.
##   X: the n by p feature matrix associated with Y. 
##   Beta: a vector of p regression coefficients used to generate the data.
##   B_index: the index set of relevant features used to generate the data. 

logistic_inde_dg = function(n = 3000, p = 5000, q = 8, Beta = "default", index = "default")
{
  
  ## Default index set of relevant features
  if(index ==  "default")
  {
    index = 1:q
  }
  
  ## Default non-zero coefficients
  if(Beta[1] ==  "default")
  {
    Beta = 1.2 + abs(rnorm(q, 0, 4))
    Beta = (-1)^(rbinom(q, 1, 0.5)) * Beta
  }
  
  ## Generate p coefficients
  BB = matrix(0, nrow = p, ncol = 1)
  BB[index, 1] = Beta  
  
  ## Generate feature matrix
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  
  ## Generate response
  theta = X %*% BB
  pi = exp(theta) / (1 + exp(theta))
  Y = rbinom(n, size = 1, prob = pi)
  
  ## Output results
  list(Y = Y, X = X, Beta = BB, B_index = index)
  
}


## Function name: logistic_moav_dg
## Purpose: generate data in setup (e) of Example 3 in Section 4.1.
## Description: this function is used to construct Table 6 in Section 4.1. 
## Inputs:
##   n: number of observations to be generated.
##   p: total number of features.
##   q: number of relevant features.
##   Beta: non-zero coefficients in the model; the default value is specified in setup (e) of Example 3.
##   index: the index set of relevant features; the default value is specified in setup (e) of Example 3.
## Outputs:
##   Y: a binary response vector of length n.
##   X: the n by p feature matrix associated with Y. 
##   Beta: a vector of p regression coefficients used to generate the data.
##   B_index: the index set of relevant features used to generate the data.

logistic_moav_dg = function(n = 3000, p = 5000, q = 5, Beta = "default", index = "default")
{
  
  ## Default index set of relevant features
  if(index[1] ==  "default"){
    index = seq(1,(2*q-1), 2) 
  }
  
  ## Default non-zero coefficients
  if(Beta[1] ==  "default")
  {
    Beta = 2.4 + abs(rnorm(q))
    Beta = (-1)^(rbinom(q, 1, 0.4)) * Beta
  }
  
  ## Generate p coefficients
  BB = matrix(0, nrow = p, ncol = 1)
  BB[index, 1] = Beta
  
  ## Generate feature matrix
  Z = matrix(rnorm(n*(p+2)), nrow = n, ncol = p+2)
  X = matrix(0, nrow = n, ncol = p)
  for (j in 1:p){ 
    X[, j] = Z[, j+2] + Z[, j+1] + Z[, j] 
  }
  X = X/sqrt(3)  
  
  ## Generate response
  theta = X %*% BB
  pi = exp(theta) / (1 + exp(theta))
  Y = rbinom(n, size = 1, prob = pi)
  
  ## Output results
  list(Y = Y, X = X, Beta = BB, B_index = index)
  
}


## Function name: logistic_comp_dg
## Purpose: generate data in setup (f) of Example 3 in Section 4.1.
## Description: this function is used to construct Table 6 in Section 4.1. 
## Inputs:
##   n: number of observations to be generated.
##   p: total number of features.
##   p1: number of the correlated features.
##   q: number of relevant features.
##   Beta: non-zero coefficients in the model; the default value is specified in setup (f) of Example 3.
##   index: the index set of relevant features; the default value is specified in setup (f) of Example 3.
##   rr: the correlation strength among p1 features.
## Outputs:
##   Y: a binary response vector of length n.
##   X: the n by p feature matrix associated with Y. 
##   Beta: a vector of p regression coefficients used to generate the data.
##   B_index: the index set of relevant features used to generate the data.

logistic_comp_dg = function(n = 3000, p = 5000, p1 = 1000, q = 6, Beta = "default", index = "default", rr = 0.4)
{ 
  
  ## Default index set of relevant features
  if(index[1] ==  "default")
  {
    index = 1:q  
  }
  
  ## Default non-zero coefficients
  if(Beta[1] ==  "default")
  {
    Beta = c(0.4, 0.4, 0.4, 1, 1, 2)
  }
  
  ## Generate p coefficients
  BB = matrix(0, nrow = p, ncol = 1)
  BB[index, 1] = Beta 
  
  ## Generate feature matrix
  CC = matrix(rr, ncol = p1, nrow = p1)
  diag(x = CC) = 1
  Z1 = rmnorm(n = n, mean = rep(0, p1), varcov = CC)
  Z2 = matrix(rnorm(n*(p-p1)), nrow = n, ncol = (p-p1))
  X = cbind(Z1, Z2)
  
  ## Generate response
  theta = X %*% BB 
  pi = exp(theta) / (1 + exp(theta))
  Y = rbinom(n, size = 1, prob = pi)
  
  ## Output results
  list(Y = Y, X = X, Beta = BB, B_index = index)
  
}


#########################################################
## Functions to conduct distributed feature screening  ##
#########################################################

## Function name: Dis_Scr
## Purpose: conduct distributed feature screening in Example 2 of Section 4.1.
## Description: this function is used to construct Tables 3-5 in Section 4.1.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   S: index set of relevant features.
##   paux: the number of auxiliary features to compute the screening threshold. 
##   m: the number of data segments.
## Output:
##   a list containing the summary of screening results for different utilities.

Dis_Scr = function(X, Y, S, paux, m){
  
  ## Initialization
  n = length(Y) 
  p = dim(X)[2] 
  R = length(unique(Y))
  RR_FAIR = c()
  RR_PSIS = c()
  RR_KF = c()
  RR_FKF = c()
  RR_MV = c()
  RR_CRU = c()
  
  ## Data partition
  Sn = rep((n/m), m)
  M_num = fenzu(Sn) 
  
  ## Generate auxiliary features
  paux_max = max(paux) 
  Saux = sample(1:p, paux_max)
  Xaux = X[, Saux] 
  for (ii in 1:m){ 
    wzii = M_num[[ii]] 
    nii = length(wzii)
    Xwzii = Xaux[wzii, ] 
    wzrd = sample(1:nii, nii) 
    Xaux[wzii, ] = Xwzii[wzrd, ] 
  } 
  npaux = length(paux)
  
  ## Begin feature screening
  if (R == 2){ 
    
    ## Methods for binary classification ##
    
    ## FAIR
    WT = FAIR_OSFA_Value(X, Y, M_num)
    W_FAIR = WT$Value
    WTaux = FAIR_OSFA_Value(Xaux, Y, M_num)
    W_FAIRaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_FAIRaux[1:paux[ii]])
      res = Evaluation(W_FAIR, S, thresh)
      RR_FAIR = rbind(RR_FAIR, c(res, (WT$TC+WTaux$TC)))}
    
    # KF
    WT = KF_SWA_Value(X, Y, M_num)
    W_KF = WT$Value
    WTaux = KF_SWA_Value(Xaux, Y, M_num)
    W_KFaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_KFaux[1:paux[ii]])
      res = Evaluation(W_KF, S, thresh)
      RR_KF = rbind(RR_KF, c(res, (WT$TC+WTaux$TC)))}
    
  }else{  
    
    ## Methods for multi-classification ##
    
    ## PSIS
    WT = PSIS_OSFA_Value(X, Y, M_num)
    W_PSIS = WT$Value
    WTaux = PSIS_OSFA_Value(Xaux, Y, M_num)
    W_PSISaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_PSISaux[1:paux[ii]])
      res = Evaluation(W_PSIS, S, thresh)
      RR_PSIS = rbind(RR_PSIS, c(res, (WT$TC+WTaux$TC)))}
    
    ## FKF
    WT = FKF_SWA_Value(X, Y, M_num)
    W_FKF = WT$Value
    WTaux = FKF_SWA_Value(Xaux, Y, M_num)
    W_FKFaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_FKFaux[1:paux[ii]])
      res=Evaluation(W_FKF, S, thresh)
      RR_FKF = rbind(RR_FKF, c(res, (WT$TC+WTaux$TC)))}
  }
  
  ## Methods for general classification ##
  
  ## MV
  WT = MV_OSFA_Value(X, Y, M_num)
  W_MV = WT$Value
  WTaux = MV_OSFA_Value(Xaux, Y, M_num)
  W_MVaux = WTaux$Value
  for (ii in 1:npaux){
    thresh = max(W_MVaux[1:paux[ii]])
    res = Evaluation(W_MV, S, thresh)
    RR_MV = rbind(RR_MV, c(res, (WT$TC+WTaux$TC)))}
  
  ## CRU
  WT = CRU_OSFA_Value(X, Y, M_num)
  W_CRU = WT$Value
  WTaux = CRU_OSFA_Value(Xaux, Y, M_num)
  W_CRUaux = WTaux$Value
  for (ii in 1:npaux){
    thresh = max(W_CRUaux[1:paux[ii]])
    res = Evaluation(W_CRU, S, thresh)
    RR_CRU = rbind(RR_CRU, c(res, (WT$TC+WTaux$TC)))}
  
  ## Output screening results
  list(RR_FAIR = RR_FAIR,
       RR_PSIS = RR_PSIS,
       RR_KF = RR_KF,
       RR_FKF = RR_FKF,
       RR_MV = RR_MV,
       RR_CRU = RR_CRU)  
}


## Function name: Dis_Scr_K
## Purpose: conduct distributed feature screening in Example 3 of Section 4.1.
## Description: this function is used to construct Table 6 in Section 4.1.
## Inputs:
##   X: an N by p feature matrix with each column denoting a feature and each row denoting
##      an observation vector. The input should be a 'matrix' object for numerical data.
##   Y: a vector of categorical responses.
##   S: index set of relevant features.
##   K: the number of features to be retained after screening.
##   m: the number of data segments.
## Output:
##   a list containing the summary of screening results for different utilities.

Dis_Scr_K = function(X, Y, S, K, m){
  
  n = length(Y)
  p = dim(X)[2]
  R = length(unique(Y))
  nK = length(K)
  
  RR_FAIR = c()
  RR_PSIS = c()
  RR_KF = c()
  RR_FKF = c()
  RR_MV = c()
  RR_CRU = c()
  
  ## split the data set
  Sn = rep(floor(n/m), m)
  Nr = n-sum(Sn)
  if (Nr>0){Sn[1:Nr] = Sn[1:Nr]+1}
  print(Sn)
  M_num = fenzu(Sn) 
  
  if (R ==  2){
    
    ## FAIR estimated by OSFA
    WT = FAIR_OSFA_Value(X, Y, M_num)
    W_FAIR = WT$Value
    for (ii in 1:nK){
      res = Evaluation_K(W_FAIR, S, K[ii])
      RR_FAIR = rbind(RR_FAIR, c(res, WT$TC))}
    
    ## KF estimated by SWA
    WT = KF_SWA_Value(X, Y, M_num)
    W_KF = WT$Value
    for (ii in 1:nK){
      res = Evaluation_K(W_KF, S, K[ii])
      RR_KF = rbind(RR_KF, c(res, WT$TC))}
    
  }else{
    
    ## PSIS estimated by OSFA
    WT = PSIS_OSFA_Value(X, Y, M_num)
    W_PSIS = WT$Value
    for (ii in 1:nK){
      res = Evaluation_K(W_PSIS, S, K[ii])
      RR_PSIS = rbind(RR_PSIS, c(res, WT$TC))}
    
    ## FKF estimated by SWA
    WT = FKF_SWA_Value(X, Y, M_num)
    W_FKF = WT$Value
    for (ii in 1:nK){
      res = Evaluation_K(W_FKF, S, K[ii])
      RR_FKF = rbind(RR_FKF, c(res, WT$TC))}
  }
  
  # MV estimated by OSFA
  WT = MV_OSFA_Value(X, Y, M_num)
  W_MV = WT$Value
  for (ii in 1:nK){
    res = Evaluation_K(W_MV, S, K[ii])
    RR_MV = rbind(RR_MV, c(res, WT$TC))}
  
  ## CRU estimated by OSFA
  WT = CRU_OSFA_Value(X, Y, M_num)
  W_CRU = WT$Value
  for (ii in 1:nK){
    res = Evaluation_K(W_CRU, S, K[ii])
    RR_CRU = rbind(RR_CRU, c(res, WT$TC))}
  
  list(RR_FAIR = RR_FAIR, 
       RR_PSIS = RR_PSIS,
       RR_KF = RR_KF,
       RR_FKF = RR_FKF,
       RR_MV = RR_MV,
       RR_CRU = RR_CRU)
}


##############################################
## Functions to summarize screening result  ##
##############################################

## Function name: Evaluation
## Purpose: evaluate a feature screening result with a given screening threshold 
## Description: this function is used to construct Tables 3-5 in Section 4.1; it is
##              also a helper function to function "Dis_Scr". 
## Inputs:  
##   Omega: a vector of estimated screening utilities.
##   S: index set of relevant features.
##   threshold: the screening threshold gamma, as defined in the paper. 
## Outputs: 
##   Jishu: an binary indicator for successful screening.
##   PS: positive selection rate.
##   FD: false discovery rate.
##   dn: the number of features retained after screening.
##   Fugaizhi: the rank of the weakest relevant feature.

Evaluation = function(Omega, S, threshold){
  
  ## Sort utilities
  sort_omega = sort.int(Omega, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_omega = sort_omega$ix 
  
  ## Find the rank of relevant features
  pd = order_omega %in% S 
  Fugaizhi = max(which(pd ==  1)) 
  
  ## Features retained after screening
  weizhi2 = which(Omega >= threshold) 
  
  ## Output evaluation results
  d = length(S)
  dn = length(weizhi2)  
  a = S %in% weizhi2    
  Jishu = all(a)           
  PS = sum(a)/d           
  if(dn>0){FD = (dn-sum(a))/dn}else{FD = 1}    
  c(Jishu, PS, FD, dn, Fugaizhi)
}


## Function name: Evaluation_K
## Purpose: evaluate a screening result with a fixed number of retained features 
## Description: this function is used to construct Table 6 in Section 4.1. 
## Inputs:  
##   Omega: a vector of estimated screening utilities.
##   S: index set of relevant features.
##   K: the number of features to be retained after screening.
## Outputs: 
##   Jishu: an binary indicator for successful screening.
##   PS: positive selection rate.
##   FD: false discovery rate.
##   Fugaizhi: the rank of the weakest relevant feature.

Evaluation_K = function(Omega, S, K){
  
  ## Sort utilities
  sort_omega = sort.int(Omega, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_omega = sort_omega$ix 
  
  ## Find the rank of relevant features
  pd = order_omega %in% S 
  Fugaizhi=max(which(pd ==  1)) 
  
  ## Features retained after screening
  weizhi2 = order_omega[1:K]
  
  ## Output evaluation results
  d = length(S)
  dn = K                
  a = S %in% weizhi2    
  Jishu = all(a)       
  PS = sum(a)/d           
  if(dn>0){FD = (dn-sum(a))/dn}else{FD = 1} 
  c(Jishu, PS, FD, Fugaizhi)
}


###############################
## Other helper functions    ##
###############################

## Function name: shixingM
## Purpose: generate an n*n indicator matrix to efficiently compute cdf or 
##          conditional cdf of a numerical feature.
## Description: a helper function for functions "FKF_Value", "KF_Value", and "MV_component_matrix".

shixingM = function(X){
  
  ## Initialization
  n = length(X)
  Xheng = matrix(X, n, n)
  Xshu = t(Xheng)

  ## Generate indicator matrix
  Xcha = Xheng-Xshu
  M = matrix(0, n, n)
  M[which(Xcha <= 0)] = 1  
  return(M) 
}


## Function name: shixingM_CRU
## Purpose: generate an n*n indicator matrix to efficiently estimate the component parameters in CRU.
## Description: a helper function for functions "CRU_index", "CRU_Value", and "CRU_OSFA_Value".

shixingM_CRU = function(X){
  
  ## Initialization
  n = length(X)
  Xheng = matrix(X, n, n)
  Xshu = t(Xheng)
  
  ## Generate indicator matrix
  Xcha = Xheng-Xshu  
  M = matrix(0, n, n)
  M[which(Xcha>0)] = 1  
  return(M) 
}


## Function name: fenzu
## Purpose: create a random partition of data
## Description: this function is used to construct Figure 2 and Tables 2-6 in Section 4.1 
## Input: 
##   Sn: a vector of length m indicating the local sample sizes of m data segments.
## Output: 
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments.

fenzu = function(Sn){
  
  ## Initialization
  m = length(Sn)
  M_num = lapply(1:m, function(h) return(h)) 
  N = sum(unlist(Sn)) 
  Weizhi = 1:N 
  
  ## Data partition by sampling 
  for (j in 1:(m-1)){
    num = sample(Weizhi, Sn[[j]])
    M_num[[j]] = num
    Weizhi = setdiff(Weizhi, num)
  }

  ## Output membership of observations for different data segments
  M_num[[m]] = Weizhi 
  M_num 
}


######################################################
## Functions to conduct KNN with selected features  ##
######################################################

## Function name: Benchmark_KNN_pred.
## Purpose: KNN classification based on K randomly selected features.
## Description: this function is used to construct Table 7 in Section 4.2.
## Inputs:
##   X_train: a feature matrix of trainning set.
##   X_test: a feature matrix of testing set.
##   Y_train: a vector of class labels corresponding to X_train.
##   Y_test: a vector of class labels corresponding to X_test.
##   K: the number of randomly selected features used in KNN classification; the default is 20. 
##   kk: the number of neighbors considered in KNN; the default is 200. 
## Output:
##   AR: accuracy of the classifier on X_test

Benchmark_KNN_pred = function(X_train, X_test, Y_train, Y_test, K = 20, kk = 200){
    
  ## Create trainning and testing feature matrices with K random selected features
  Ind_random = sample(1:p, K)
  X_train_K = X_train[, Ind_random]
  X_test_K = X_test[, Ind_random]
  
  ## KNN classification
  res_knn = knn(X_train_K, X_test_K, Y_train, k = kk)
  AR = mean(res_knn == Y_test)
  AR  
}


## Function name: Scr_KNN_pred
## Purpose: KNN classification based on distributed feature screening.
## Description: this function is used to construct Table 7 in Section 4.2.
## Inputs: 
##   X_train: a feature matrix of trainning set.
##   X_test: a feature matrix of testing set.
##   Y_train: a vector of class labels corresponding to X_train.
##   Y_test: a vector of class labels corresponding to X_test.
##   M_num: a list containing m vectors to indicate the membership of observations
##          belonging to each of the m data segments. 
##   K: the number of features to be retained; the default is 20. 
##   kk: the number of neighbors considered in KNN; the default is 200. 
##   method: the method to be used for distributed feature screening;
##           choose from CRU (default), PSIS, FKF, and MV. 
## Output:
##   AR: accuracy of the classifier on X_test

Scr_KNN_pred = function(X_train, X_test, Y_train, Y_test, M_num, K = 20, kk = 200, method = c("CRU", "PSIS", "FKF", "MV")){
  
  ## Estimate screening utilities 
  if(method == "PSIS"){W = PSIS_OSFA_Value(X_train, Y_train, M_num)}
  if(method == "FKF"){W = FKF_SWA_Value(X_train, Y_train, M_num)}
  if(method == "MV"){W = MV_OSFA_Value(X_train, Y_train, M_num)}
  if(method == "CRU"){W = CRU_OSFA_Value(X_train, Y_train, M_num)}
  
  ## Feature screening by retaining top K features 
  sort_w = sort.int(W$Value, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_w = sort_w$ix 
  Ind = order_w[1:K]

  ## Create trainning and testing feature matrices based on the retained features
  X_train_K = X_train[, Ind]
  X_test_K = X_test[, Ind]
  
  ## KNN classification
  res_knn = knn(X_train_K, X_test_K, Y_train, k = kk)
  AR = mean(res_knn == Y_test)
  AR
  
}

