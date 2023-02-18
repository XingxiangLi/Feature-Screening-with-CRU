#########################################################################################
# This script provides R functions needed to produce the numerical results in the paper.#
#########################################################################################
## Last update: February 3 2022
## Require R > 4.0.0

#############################################################
########################    MV    ###########################
#############################################################

# "helper" for function "MV_OSFA_Value"
MV_component_matrix<-function(X,Y,class){
  
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # class: category vector
  
  R = length(class) 
  n = length(Y)  
  p = dim(X)[2]  
  
  Yheng = matrix(Y,n,R)
  Yshu = t(matrix(class,R,n))
  Ycha = Yheng-Yshu 
  
  MY=matrix(0,n,R)
  MY[which(Ycha==0)]=1 
  Nr=colSums(MY) 
  Num = n*(n-1)*(n-2)
  
  Res= list()
  
  for (r in 1:R){
    
    VYr = MY[,r]
    MYr = matrix(VYr,n,n)
    
    Er = c()
    
    if(Nr[r]==0){Er=matrix(0,2,p)}else{ 
      
      for (j in 1:p){
        
        Xj = X[,j]
        MXj = shixingM(Xj)  
        diag(MXj) = 0 
        
        MXY = MXj*MYr  
        
        Ar = MXY%*%t(MXY)  
        diag(Ar) = 0       
        E1 = sum(Ar)/Num   
        
        Br = MXY%*%t(MXj)
        diag(Br) = 0
        E2 = sum(Br)/Num   
        
        res = c(E1,E2)
        
        Er = cbind(Er,res)
        
      }}
    
    Res[[r]]=Er
    
  }
  
  Res[[R+1]]=Nr
  
  Res
  
}

# compute the mean variance (MV) between X and Y via OSFA estimation 
# this function is used in Section 4.1 (Example 2 and 3)
MV_OSFA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("MV_OSFA")
  
  m=length(M_num)
  class=c()     
  for (ii in 1:m){
    DY=Y[M_num[[ii]]]
    class_ii=unique(DY)  
    class=union(class,class_ii)  
  }
  
  R = length(class)
  p = dim(X)[2]
  
  TT = c(); Vw = c(); Pr = rep(0,R)
  theta1 = theta2 = matrix(0,1,p)
  MCP = lapply(rep(0,R), function(h) return(h)) 
  
  for (ii in 1:m) {
    DX=X[M_num[[ii]],] 
    DY=Y[M_num[[ii]]] 
    
    ptm <- proc.time()
    
    nii = length(DY)
    w = floor(nii/3);Vw=c(Vw,w)
    
    res_ii = MV_component_matrix(DX,DY,class)
    Pr = Pr + res_ii[[R+1]]
    
    for (r in 1:R){
      MCP[[r]] = MCP[[r]]+res_ii[[r]]*w    
    }
    
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
  } 
  
  ptm <- proc.time()
  
  MV = matrix(1/3,1,p)
  n3 = sum(Vw)
  
  Pr = Pr/sum(Pr)     #theta_{y,r}
  
  for (r in 1:R){
    aa = MCP[[r]]/n3  #n3 = sum[n_l/3]
    MV = MV+aa[1,]/Pr[r]-2*aa[2,]
  }
  
  at <- proc.time()-ptm
  
  TC = mean(TT)+at[3]  
  
  list(Value=c(MV),TC=TC)
  
}


#############################################################
########################    FAIR   ##########################
#############################################################

# the non-distributed estimator of FAIR between two vectors in Section 2.1 (Figure 1)
FAIR_index<-function(X,Y){
  
  class=unique(Y);
  weizhi1=which(Y==class[1])
  weizhi2=which(Y==class[2])
  num1=length(weizhi1);
  num2=length(weizhi2);
  
  X1 = X[weizhi1]
  junzhi1=mean(X1)
  Cha1 = X1-junzhi1
  Var1 = sum(Cha1^2)/(num1-1)
  
  X2 = X[weizhi2]
  junzhi2=mean(X2)
  Cha2 = X2-junzhi2
  Var2 = sum(Cha2^2)/(num2-1)
  
  W=abs((junzhi1-junzhi2)/(sqrt(Var1/num1+Var2/num2)))
  
  W
  
}

# compute the FAIR utility between X and Y via OSFA estimation 
# this function is used in Section 4.1 (Example 2 and 3)
FAIR_OSFA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("FAIR_OSFA")
  
  m=length(M_num)
  class=c();    
  for (ii in 1:m){
    DY=Y[M_num[[ii]]]
    class_ii=unique(DY) 
    class=union(class,class_ii)  
  }
  
  R=length(class); 
  p=dim(X)[2]
  
  Num=matrix(0,m,R) 
  He=matrix(0,R,p) 
  SHe=matrix(0,R,p) 
  TT=c()
  
  for (ii in 1:m) {
    DX=X[M_num[[ii]],] 
    DY=Y[M_num[[ii]]] 
    
    ptm <- proc.time()
    nii=length(DY)        
    Yheng=matrix(DY,nii,R)
    Yshu=t(matrix(class,R,nii))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,nii,R)
    M[which(Ycha==0)]=1
    num=colSums(M)   
    Num[ii,]=num    
    
    HE_ii=matrix(0,R,p) 
    SHE_ii=matrix(0,R,p)  
    for (jj in 1:R){
      MY=matrix(rep(M[,jj],p),nrow = nii)
      HE_ii[jj,]=colSums(DX*MY) 
      SHE_ii[jj,]=colSums((DX*MY)^2)
    }
    He=He+HE_ii 
    SHe=SHe+SHE_ii 
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
  }
  
  ptm <- proc.time()
  N12=colSums(Num) 
  junzhi=He/matrix(rep(N12,1,p),2,p) 
  junzhi_s=SHe/matrix(rep(N12,1,p),2,p) 
  
  junzhi1=junzhi[1,]
  junzhi2=junzhi[2,]
  num1=N12[1]
  num2=N12[2]
  fangcha1=(junzhi_s[1,]-junzhi1^2)*(num1/(num1-1));
  fangcha2=(junzhi_s[2,]-junzhi2^2)*(num2/(num2-1));
  W=abs((junzhi1-junzhi2)/(sqrt(fangcha1/num1+fangcha2/num2)))
  at <- proc.time()-ptm
  
  TC = mean(TT)+at[3]
  list(Value=W,TC=TC)
  
}

#############################################################
######################     PSIS     #########################
#############################################################

# compute the PSIS utility between X and Y via OSFA estimation 
# this function is used in Section 4.1 (Example 2)
PSIS_OSFA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("PSIS_OSFA")
  
  m=length(M_num)
  class=c();     
  for (ii in 1:m){
    DY=Y[M_num[[ii]]]
    class_ii=unique(DY)
    class=union(class,class_ii) 
  }
  
  R=length(class); 
  p=dim(X)[2]
  
  X=scale(X) 
  Num=matrix(0,m,R) 
  He=matrix(0,R,p) 
  TT=c()
  
  for (ii in 1:m) {
    DX=X[M_num[[ii]],]
    DY=Y[M_num[[ii]]]
    
    ptm <- proc.time()
    nii=length(DY);        
    Yheng=matrix(DY,nii,R)
    Yshu=t(matrix(class,R,nii))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,nii,R)
    M[which(Ycha==0)]=1
    num=colSums(M)   
    Num[ii,]=num 
    
    HE_ii=matrix(0,R,p) 
    for (jj in 1:R){
      MY=matrix(rep(M[,jj],p),nrow = nii)
      HE_ii[jj,]=colSums(DX*MY) 
    }
    He=He+HE_ii
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
  }
  
  ptm <- proc.time()
  Nm=colSums(Num); 
  junzhi=He/matrix(rep(Nm,p),R,p) 
  
  PairCha=c()
  for (i1 in 1:(R-1)){
    for (i2 in (i1+1):R){
      aa=abs(junzhi[i1,]-junzhi[i2,])
      PairCha=rbind(PairCha,aa)
    }}
  
  PSIS = apply(PairCha,2,max)
  
  at <- proc.time()-ptm
  
  TC = mean(TT)+at[3]
  
  list(Value=PSIS,TC=TC)
  
}

#############################################################
#######################      KF     #########################
#############################################################

# "helper" for function "KF_SWA_Value"
KF_Value<-function(X,Y){
  
  n=length(Y)     
  p=dim(X)[2]
  class=unique(Y) 
  R=length(class) 
  
  if (R==1){KF_vec=rep(0,p)}else{
    
    Yheng=matrix(Y,n,R)
    Yshu=t(matrix(class,R,n))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,n,R)
    M[which(Ycha==0)]=1
    Pr=colMeans(M) 
    
    KF_vec = c()
    for (j in 1:p){
      Xj=X[,j]  
      MX=shixingM(Xj)
      FM=matrix(0,R,n)
      for(ii in 1:R){
        MY=matrix(M[,ii],n,n)
        MXY=MX*MY
        FXY=colMeans(MXY)/Pr[ii] 
        FM[ii,]=FXY
      }
      
      KF=max(abs(FM[1,]-FM[2,]))
      KF_vec = c(KF_vec,KF)}}
  
  KF_vec
}

# compute the KF utility between X and Y via SWA estimation 
# this function is used in Section 4.1 (Example 2 and 3)
KF_SWA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("KF_SWA")
  
  WM=c();ns=c();TT=c()
  m=length(M_num)
  for (ii in 1:m) {
    DX=X[M_num[[ii]],] 
    DY=Y[M_num[[ii]]] 
    nii=length(DY);ns=c(ns,nii)
    
    ptm <- proc.time()
    Wii = KF_Value(DX,DY)*nii
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
    WM=rbind(WM,Wii)}
  
  W_KF = colSums(WM)/sum(ns) 
  
  TC = mean(TT)
  
  list(Value=W_KF,TC=TC)
  
}


#############################################################
#########################      FKF     ######################
#############################################################

# "helper" for function "FKF_SWA_Value"
FKF_Value<-function(X,Y){
  
  n=length(Y)     
  p=dim(X)[2]
  class=unique(Y) 
  R=length(class) 
  
  if (R==1){FKF_vec=rep(0,p)}else{
    
    Yheng=matrix(Y,n,R)
    Yshu=t(matrix(class,R,n))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,n,R)
    M[which(Ycha==0)]=1
    Pr=colMeans(M) 
    
    FKF_vec = c()
    for (j in 1:p){
      Xj=X[,j]  
      MX=shixingM(Xj)
      FM=matrix(0,R,n)
      for(ii in 1:R){
        MY=matrix(M[,ii],n,n)
        MXY=MX*MY
        FXY=colMeans(MXY)/Pr[ii] 
        FM[ii,]=FXY
      }
      
      PairM=c()
      for (i1 in 1:(R-1)){
        for (i2 in (i1+1):R){
          aa=max(abs(FM[i1,]-FM[i2,])); 
          PairM=c(PairM,aa)
        }}
      KF = max(PairM)
      FKF_vec = c(FKF_vec,KF)
    }}
  FKF_vec
}

# compute the FKF utility between X and Y via SWA estimation 
# this function is used in Section 4.1 (Example 2)
FKF_SWA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("FKF_SWA")
  
  WM=c();ns=c();TT=c()
  m=length(M_num)
  for (ii in 1:m) {
    DX=X[M_num[[ii]],]
    DY=Y[M_num[[ii]]] 
    
    ptm <- proc.time()
    nii=length(DY);ns=c(ns,nii)
    Wii = FKF_Value(DX,DY)*nii
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
    WM=rbind(WM,Wii)
  }
  
  W_FKF = colSums(WM)/sum(ns)
  
  TC = mean(TT)
  
  list(Value=W_FKF,TC=TC)
  
}

#############################################################
######################      CRU     #########################
#############################################################

# the non-distributed estimator of CRU between two vectors in Section 2.1 (Figures 1)
CRU_index<-function(X,Y){
  
  n=length(Y)    
  class=unique(Y) 
  R=length(class) 
  
  if (R==1){Wj=0}else{
    Yheng=matrix(Y,n,R)
    Yshu=t(matrix(class,R,n))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,n,R)
    M[which(Ycha==0)]=1
    Pr=colMeans(M) 
    
    MX=shixingM_CRU(X)
    MXY = t(M)%*%MX
    MXY_mean = rowSums(MXY)/(n*(n-1))
    Wj=sum((MXY_mean-0.5*Pr)^2)}
  
    Wj
}

# the non-distributed estimator of CRU used in Section 2.1 and 4.1 (Figures 1-2, Table 2)
# the following Y is the response vector, the following X is the matrix of covariates
CRU_Value<-function(X,Y){
  
  n=length(Y)    
  p=dim(X)[2]
  class=unique(Y) 
  R=length(class) 
  
  if (R==1){CRU_vec=rep(0,p)}else{
    
    Yheng=matrix(Y,n,R)
    Yshu=t(matrix(class,R,n))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,n,R)
    M[which(Ycha==0)]=1
    Pr=colMeans(M) 
    
    CRU_vec = c()
    for (j in 1:p){
      Xj=X[,j]  
      MX=shixingM_CRU(Xj)
      MXY = t(M)%*%MX
      MXY_mean = rowSums(MXY)/(n*(n-1))
      Wj=sum((MXY_mean-0.5*Pr)^2)
      CRU_vec = c(CRU_vec,Wj)
    }}
  CRU_vec
}

# the non-distributed (centralized) computational time for CRU 
CRU_Central_Time<-function(X,Y){
  ptm <- proc.time()
  n=length(Y)    
  p=dim(X)[2]
  class=unique(Y) 
  R=length(class) 
  
  if (R==1){CRU_vec=rep(0,p)}else{
    
    Yheng=matrix(Y,n,R)
    Yshu=t(matrix(class,R,n))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,n,R)
    M[which(Ycha==0)]=1
    Pr=colMeans(M) 
    
    CRU_vec = c()
    for (j in 1:p){
      Xj=X[,j]  
      MX=shixingM_CRU(Xj)
      MXY = t(M)%*%MX
      MXY_mean = rowSums(MXY)/(n*(n-1))
      Wj=sum((MXY_mean-0.5*Pr)^2)
      CRU_vec = c(CRU_vec,Wj)
    }}
  
  ct <- proc.time()-ptm
  TT = ct[3]
  
  list(Value=CRU_vec, TC=TT)
}


# compute the CRU between X and Y via SWA estimation 
# this function is used in Section 4.1 (Examples 1-3)
CRU_SWA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  WM=c();ns=c();TT=c()
  m=length(M_num)
  for (ii in 1:m) {
    DX=X[M_num[[ii]],] 
    DY=Y[M_num[[ii]]] 
    
    ptm <- proc.time()
    nii=length(DY);ns=c(ns,nii)
    Wii = CRU_Value(DX,DY)*nii
    WM=rbind(WM,Wii)
    
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
  }
  
  W_CRU = colSums(WM)/sum(ns)
  
  TC = mean(TT)
  
  list(Value=W_CRU,TC=TC)
  
}

# compute the CRU between X and Y via OSFA estimation 
# this function is used in Section 4.1 (Examples 1-3)
CRU_OSFA_Value<-function(X,Y,M_num){
  
  # input:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # M_num is the index sets of local data segments
  # output:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  print("CRU_OSFA")
  
  m=length(M_num)
  class=c();     
  for (ii in 1:m){
    DY=Y[M_num[[ii]]]
    class_ii=unique(DY)  
    class=union(class,class_ii)  
  }
  
  R=length(class)
  n=dim(X)[1]
  p=dim(X)[2]
  
  if (R==1){CRU_vec=rep(0,p);TC = 0}else{
    MCPy=0 
    MCPj=0 
    Sw=0
    TT=c()
    
    for (ii in 1:m) {
      DX=X[M_num[[ii]],]
      DY=Y[M_num[[ii]]] 
      
      ptm <- proc.time()
      nii=length(DY)    
      Yheng=matrix(DY,nii,R)
      Yshu=t(matrix(class,R,nii))
      Ycha=Yheng-Yshu 
      
      M=matrix(0,nii,R)
      M[which(Ycha==0)]=1
      nr=colSums(M) 
      MCPy=MCPy+matrix(rep(nr,p),R,p)/n 
      
      wm = floor(nii/2)
      Sw= Sw+wm
      
      MCPjm=c()
      for (j in 1:p){
        Xj=DX[,j]  
        MX=shixingM_CRU(Xj)
        MXY = t(M)%*%MX
        MXY_mean = wm*rowSums(MXY)/(nii*(nii-1))
        MCPjm = cbind(MCPjm,MXY_mean)
      }
      
      MCPj=MCPj+MCPjm
      ct <- proc.time()-ptm
      TT=c(TT,ct[3])
      
    }
    
    ptm <- proc.time()
    CRU_vec = colSums((MCPj/Sw-.5*MCPy)^2)
    at <- proc.time()-ptm
    TC = mean(TT)+at[3]
    
  }
  
  list(Value=CRU_vec,TC=TC)
}

#############################################################
################# Example 2: data generation     ############
#############################################################

## "helper" function to generate categorical response in Sections 2.1 and 4.1 (Figures 1-2, Examples 1-2)  
dis_gen <- function(n,pp){ 
  
  # input:
  # n: sample size
  # pp: the probabilities of all categories (p_r in our paper)
  # output:
  # Y: the generated categorical response
  # TC: the index sets of all categories 
  
  lp = length(pp) 
  cp = c()
  lc = lp-1
  for (i in 1:lc){
    cp = c(cp, sum(pp[1:i]))
  }
  Q=qunif(cp) 
  
  x=runif(n) 
  X=x
  X[which(x<Q[1])]=1
  if (lp>2){
    for (ii in 2:lc){X[which(x>=Q[ii-1] & x<Q[ii])]=ii}
  }
  X[which(x>=Q[lc])]=lp 
  
  S=list()
  for (r in 1:lp){
    S[[r]]=which(X==r)
  }
  
  list(Y=X,SR=S)
}

## "helper" function to generate different local sample sizes n_1,...,n_m in Sections 4.1 (setup (c) in Examples 2)  
nlgen<-function(n,m){
  
  # input:
  # n: sample size
  # m: the number of data segments
  # output:
  # PL: local sample sizes (n_1, n_2,..., n_m)
  
  a = n/(2.5*m)
  PL=c(rep(a,m/4),rep(2*a,m/4),rep(3*a,m/4),rep(4*a,m/4))
}

## "helper" function for "FKF_Value", "KF_Value", and "MV_component_matrix"
shixingM <- function(Y){
  
  n=length(Y)
  Yheng=matrix(Y,n,n)
  Yshu=t(Yheng)
  
  Ycha=Yheng-Yshu
  
  M=matrix(0,n,n)
  M[which(Ycha<=0)]=1
  
  return(M) 
}

## "helper" function for "CRU_index", "CRU_Value", "CRU_Central_Time", "CRU_OSFA_Value"
shixingM_CRU <- function(Y){
  
  n=length(Y)
  Yheng=matrix(Y,n,n)
  Yshu=t(Yheng)
  
  Ycha=Yheng-Yshu
  
  M=matrix(0,n,n)
  M[which(Ycha>0)]=1
  
  return(M) 
}

## "helper" function to generate local index sets S1,...,S_m for distributed feature screening
fenzu<-function(Sn){
  
  # input: Sn is a list containing n_1,...,n_m
  # output: the index sets of local data segments
  
  m = length(Sn)
  M_num=lapply(1:m, function(h) return(h)) 
  N =  sum(unlist(Sn)) 
  Weizhi=1:N 
  
  for (j in 1:(m-1)){
    num=sample(Weizhi,Sn[[j]])
    M_num[[j]]=num
    Weizhi=setdiff(Weizhi,num)
  }
  M_num[[m]]=Weizhi 
  M_num 
}


# "helper" function to evaluate and summarize the result after screening
Evaluation<-function(Omega,S,threshold){
  
  # input: the screening utilities (Omega), index set of relevant features (S), screening threshold (gamma in our manuscript)  
  # output: SSR, PSR, FDR, Size, wRank, Time
  
  # sort utilities
  sort_omega<-sort.int(Omega, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_omega<-sort_omega$ix 
  
  # find the rank of relevant features
  pd = order_omega %in% S 
  Fugaizhi=max(which(pd==1)) ##wRank
  
  ## the features retained
  weizhi2=which(Omega>=threshold) 
  
  d = length(S)
  dn = length(weizhi2)  #Size
  a = S %in% weizhi2    
  Jishu = all(a)        # SSR      
  PS = sum(a)/d         # positive slection rate
  
  if(dn>0){FD=(dn-sum(a))/dn}else{FD=1}  # false discovery rate
  
  c(Jishu,PS,FD,dn,Fugaizhi)
}

# "helper" function to compute the mean SSR, PSR, FDR, Size, wRank, Time
Summarize<-function(X){
  Means = colMeans(X)
  Means 
}


## the function for distributed screening in Example 2
## methods include FAIR-OSFA, PSIS-OSFA, KF-SWA, FKF-SWA, MV-OSFA, CRU-OSFA
Dis_Scr<-function(X,Y,S,paux,m,balance){
  
  # inputs:
  # X: the matrix of all covariates (features)
  # Y: categorical response vector
  # S: index set of relevant features
  # paux: the number of auxiliary features, which are used to obtain the threshold (gamma) 
  # m: the number of data segments
  # balance==1 corresponds to the equal segmentation, balance==0 corresponds to the unequal segmentation
  # output:
  # the values of screening utilities (FAIR, PSIS, KF, FKF, MV, CRU)
  
  n=length(Y) # sample size
  p=dim(X)[2] # feature size
  R=length(unique(Y)) # the number of categories
  
  RR_FAIR = c()
  RR_PSIS = c()
  RR_KF = c()
  RR_FKF = c()
  RR_MV = c()
  RR_CRU = c()
  
  ## split the data set
  if (balance==1){
    Sn = rep((n/m),m)
    M_num = fenzu(Sn) 
  }else{
    Sn=nlgen(n,m)
    M_num = fenzu(Sn) 
  }
  
  ## generate auxiliary features
  paux_max = max(paux) 
  Saux=sample(1:p,paux_max)
  Xaux=X[,Saux] 
  for (ii in 1:m){ 
    wzii = M_num[[ii]] 
    nii = length(wzii)
    Xwzii = Xaux[wzii,] 
    wzrd = sample(1:nii,nii) 
    Xaux[wzii,] = Xwzii[wzrd,] 
  } 
  npaux = length(paux)
  
  if (R==2){ # classification for R=2
    
    # FAIR estimated by OSFA
    WT = FAIR_OSFA_Value(X,Y,M_num)
    W_FAIR = WT$Value
    WTaux = FAIR_OSFA_Value(Xaux,Y,M_num)
    W_FAIRaux = WTaux$Value
    
    for (ii in 1:npaux){
      thresh = max(W_FAIRaux[1:paux[ii]])
      res=Evaluation(W_FAIR,S,thresh)
      RR_FAIR = rbind(RR_FAIR,c(res,(WT$TC+WTaux$TC)))}
    
    # KF estimated by SWA
    WT = KF_SWA_Value(X,Y,M_num)
    W_KF = WT$Value
    WTaux = KF_SWA_Value(Xaux,Y,M_num)
    W_KFaux = WTaux$Value
    for (ii in 1:npaux){
    thresh = max(W_KFaux[1:paux[ii]])
    res=Evaluation(W_KF,S,thresh)
    RR_KF = rbind(RR_KF,c(res,(WT$TC+WTaux$TC)))}
    
  }else{  # classification for R>2
    
    # PSIS estimated by OSFA
    WT = PSIS_OSFA_Value(X,Y,M_num)
    W_PSIS = WT$Value
    WTaux = PSIS_OSFA_Value(Xaux,Y,M_num)
    W_PSISaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_PSISaux[1:paux[ii]])
      res=Evaluation(W_PSIS,S,thresh)
      RR_PSIS = rbind(RR_PSIS,c(res,(WT$TC+WTaux$TC)))}
    
    # FKF estimated by SWA
    WT = FKF_SWA_Value(X,Y,M_num)
    W_FKF = WT$Value
    WTaux = FKF_SWA_Value(Xaux,Y,M_num)
    W_FKFaux = WTaux$Value
    for (ii in 1:npaux){
      thresh = max(W_FKFaux[1:paux[ii]])
      res=Evaluation(W_FKF,S,thresh)
      RR_FKF = rbind(RR_FKF,c(res,(WT$TC+WTaux$TC)))}
  }
  
  # MV estimated by OSFA
  WT = MV_OSFA_Value(X,Y,M_num)
  W_MV = WT$Value
  WTaux = MV_OSFA_Value(Xaux,Y,M_num)
  W_MVaux = WTaux$Value
  for (ii in 1:npaux){
    thresh = max(W_MVaux[1:paux[ii]])
    res=Evaluation(W_MV,S,thresh)
    RR_MV = rbind(RR_MV,c(res,(WT$TC+WTaux$TC)))}

  # CRU estimated by OSFA
  WT = CRU_OSFA_Value(X,Y,M_num)
  W_CRU = WT$Value
  WTaux = CRU_OSFA_Value(Xaux,Y,M_num)
  W_CRUaux = WTaux$Value
  for (ii in 1:npaux){
    thresh = max(W_CRUaux[1:paux[ii]])
    res=Evaluation(W_CRU,S,thresh)
    RR_CRU = rbind(RR_CRU,c(res,(WT$TC+WTaux$TC)))}
  
  list(RR_FAIR = RR_FAIR,
       RR_PSIS = RR_PSIS,
       RR_KF = RR_KF,
       RR_FKF = RR_FKF,
       RR_MV = RR_MV,
       RR_CRU = RR_CRU)
  
}


#############################################################
######## Example 3 in Section 4.1: Logistic model ###########
#############################################################

## "helper" function for data generation of logistic model
## setup (d) independent features
logistic_inde_dg<-function(n=3000, p=10000, q=8, Beta="default", inc=0, index=-1, sig=1, tag="N")
{
  
  # inputs:
  # n: sample size
  # p: dimension of features
  # q: the number of relevant features
  # Beta: non-zero coefficients in logistic model
  # inc: value of intercept term 
  # index: index set of relevant features, index==-1 means the set is designed by this function
  # sig: the standard deviation of random error in linear model, other models can skip it
  # tag: "N" indicates linear model, "B" indicates logistic model, "P" indicates Poisson model
  
  # output:
  # Y: response
  # X: covariate matrix
  # Beta: regression coefficients
  
  if(index[1]==-1)
  {
    index <- 1:q
  }
  
  if(Beta[1]=="default")
  {
    Beta<-  1.2 + abs(rnorm(q,0,4))
    Beta<-  (-1)^(rbinom(q, 1,0.5)) * Beta
  }
  
  BB<-matrix(0, nrow=p, ncol=1)
  BB[index, 1] <- Beta  
  
  X<-matrix(rnorm(n*p), nrow=n, ncol=p)
  theta <- X %*% BB + inc
  
  if(tag=="N")
  {
    Y<- theta + rnorm(n, mean=0, sd=sig)
  }
  
  if(tag=="B")
  {
    pi <- exp(theta) / (1 + exp(theta))
    Y  <- rbinom(n, size=1, prob=pi)
  }
  
  if(tag=="P")
  {
    mmu <- exp(theta)
    Y  <- rpois(n, lambda=mmu)
  }
  
  list(Y=Y, X=X, Beta=BB, B_index=index)
  
}

## "helper" function for data generation of logistic model
## setup (e) correlated features 
logistic_moav_dg<-function(n=3000, p=10000, q=5, Beta="default", inc=0, index=-1, sig=1, tag="N")
{
  
  # inputs:
  # n: sample size
  # p: dimension of features
  # q: the number of relevant features
  # Beta: non-zero coefficients in logistic model
  # inc: value of intercept term 
  # index: index set of relevant features, index==-1 means the set is designed by this function
  # sig: the standard deviation of random error in linear model, other models can skip it
  # tag: "N" indicates linear model, "B" indicates logistic model, "P" indicates Poisson model
  
  # output:
  # Y: response
  # X: covariate matrix
  # Beta: regression coefficients
 
  if(index[1]==-1){
  index <- seq(1,(2*q-1),2)    # select relevant features randomly
  }

  if(Beta[1]=="default")
  {
    Beta<-  2.4 + abs(rnorm(q))
    Beta<-  (-1)^(rbinom(q, 1, 0.4)) * Beta
  }
  
  BB<-matrix(0, nrow=p, ncol=1)
  BB[index, 1] <- Beta
  
  Z<-matrix(rnorm(n*(p+2)), nrow=n, ncol=p+2)
  X<-matrix(0, nrow=n, ncol=p)
  
  for(j in 1:p)
  { X[,j]<- Z[,j+2] +  Z[,j+1] + Z[,j] }
  
  X<-X/sqrt(3)  
  
  theta <- X %*% BB + inc
  
  if(tag=="N")
  {
    Y<- theta + rnorm(n, mean=0, sd=sig)
  }
  
  if(tag=="B")
  {
    pi <- exp(theta) / (1 + exp(theta))
    Y  <- rbinom(n, size=1, prob=pi)
  }
  
  if(tag=="P")
  {
    mmu <- exp(theta)
    Y  <- rpois(n, lambda=mmu)
  }
  
  list(Y=Y, X=X, Beta=BB, B_index=index)
  
}
 
## "helper" function for data generation of logistic model  
## setup (f) correlated features
logistic_comp_dg<-function(n=100, p=10000, p1=1000, q=8, Beta=c(3,3,3,3,3), index=-1, inc=0, rr=0.5, sig=1, tag="N")
{ 
  
  # inputs:
  # n: sample size
  # p: dimension of features
  # p1: the number of the correlated features
  # q: the number of relevant features
  # Beta: non-zero coefficients in logistic model
  # inc: value of intercept term 
  # rr: the given degree of correlations 
  # index: index set of relevant features, index==-1 means the set is designed by this function
  # sig: the standard deviation of random error in linear model, other models can skip it
  # tag: "N" indicates linear model, "B" indicates logistic model, "P" indicates Poisson model
  
  # output:
  # Y: response
  # X: covariate matrix
  # Beta: regression coefficients
  
  if(index[1]==-1)
  {
    index <- 1:q  #select relevant features randomly
  }
  
  if(Beta[1]=="default")
  {
    Beta <- 2 + abs(rnorm(q))
    Beta <- (-1)^(rbinom(q, 1, 0.5)) * Beta
  }
  
  BB<-matrix(0, nrow=p, ncol=1)
  BB[index, 1] <- Beta 
  
  CC<-matrix(rr, ncol=p1, nrow=p1)
  diag(x=CC)<-1
  
  
  Z1<- rmnorm(n=n, mean=rep(0,p1), varcov=CC)
  Z2<-matrix(rnorm(n*(p-p1)), nrow=n, ncol=(p-p1))
  X=cbind(Z1,Z2)
  
  theta <- X %*% BB + inc
  
  if(tag=="N")
  {
    Y<- theta + rnorm(n, mean=0, sd=sig)
  }
  
  if(tag=="B")
  {
    pi <- exp(theta) / (1 + exp(theta))
    Y  <- rbinom(n, size=1, prob=pi)
  }
  
  if(tag=="P")
  {
    mmu <- exp(theta)
    Y  <- rpois(n, lambda=mmu)
  }
  
  list(Y=Y, X=X, Beta=BB, B_index=index)
  
}

# "helper" function to evaluate and summarize the result after screening, 
Evaluation_K<-function(Omega,S,K){
  
  # input: the screening utilities (Omega), index set of relevant features (S), the number of retained features (K) 
  # output: SSR, PSR, FDR, wRank
  
  # sort utilities
  sort_omega<-sort.int(Omega, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_omega<-sort_omega$ix 
  
  # find the rank of relevant features
  pd = order_omega %in% S 
  Fugaizhi=max(which(pd==1)) ##wRank
  
  ## the features retained
  weizhi2=order_omega[1:K]
  
  d = length(S)
  dn = K                #Size
  a = S %in% weizhi2    
  Jishu = all(a)        # SSR      
  PS = sum(a)/d         # positive slection rate
  
  if(dn>0){FD=(dn-sum(a))/dn}else{FD=1}  # false discovery rate
  
  c(Jishu,PS,FD,Fugaizhi)
}


## the function for distributed screening in Example 3
## methods include FAIR-OSFA, KF-SWA, MV-OSFA, CRU-OSFA
Dis_Scr_K<-function(X,Y,S,K,m,balance){
  
  # inputs:
  # X: the matrix of all covariates (features)
  # Y: categorical response 
  # K: the number of features retained, K can be a vector with multiple choices of K
  # m: the number of data segments
  # balance==1 corresponds to the equal segmentation
  # output:
  # the values of screening utilities (FAIR, KF, MV, CRU)
  
  n=length(Y)
  p=dim(X)[2]
  R=length(unique(Y))
  nK = length(K)
  
  RR_FAIR = c()
  RR_PSIS = c()
  RR_KF = c()
  RR_FKF = c()
  RR_MV = c()
  RR_CRU = c()
  
  ## split the data set
  if (balance==1){
    Sn = rep(floor(n/m),m);
    Nr = n-sum(Sn)
    if (Nr>0){Sn[1:Nr]=Sn[1:Nr]+1}
    print(Sn)
    M_num = fenzu(Sn) 
  }else{
    Sn=nlgen(n,m)
    M_num = fenzu(Sn) 
  }
  
  if (R==2){
    
    ## FAIR estimated by OSFA
    WT = FAIR_OSFA_Value(X,Y,M_num)
    W_FAIR = WT$Value
    for (ii in 1:nK){
      res=Evaluation_K(W_FAIR,S,K[ii])
      RR_FAIR = rbind(RR_FAIR,c(res,WT$TC))}
    
    ## KF estimated by SWA
    WT = KF_SWA_Value(X,Y,M_num)
    W_KF = WT$Value
    for (ii in 1:nK){
      res=Evaluation_K(W_KF,S,K[ii])
      RR_KF = rbind(RR_KF,c(res,WT$TC))}
    
  }else{

    ## PSIS estimated by OSFA
    WT = PSIS_OSFA_Value(X,Y,M_num)
    W_PSIS = WT$Value
    for (ii in 1:nK){
      res=Evaluation_K(W_PSIS,S,K[ii])
      RR_PSIS = rbind(RR_PSIS,c(res,WT$TC))}
    
    ## FKF estimated by SWA
    WT = FKF_SWA_Value(X,Y,M_num)
    W_FKF = WT$Value
    for (ii in 1:nK){
      res=Evaluation_K(W_FKF,S,K[ii])
      RR_FKF = rbind(RR_FKF,c(res,WT$TC))}
  }
  
  # MV estimated by OSFA
  WT = MV_OSFA_Value(X,Y,M_num)
  W_MV = WT$Value
  for (ii in 1:nK){
    res=Evaluation_K(W_MV,S,K[ii])
    RR_MV = rbind(RR_MV,c(res,WT$TC))}
  
  ## CRU estimated by OSFA
  WT = CRU_OSFA_Value(X,Y,M_num)
  W_CRU = WT$Value
  for (ii in 1:nK){
    res=Evaluation_K(W_CRU,S,K[ii])
    RR_CRU = rbind(RR_CRU,c(res,WT$TC))}
  
  list(RR_FAIR = RR_FAIR,
       RR_PSIS = RR_PSIS,
       RR_KF = RR_KF,
       RR_FKF = RR_FKF,
       RR_MV = RR_MV,
       RR_CRU = RR_CRU)
  
}


#####################################
######## Real data analysis #########
#####################################

# compute the mean variance (MV) between separate data sets XD and YD via OSFA estimation 
# this function is used in Section 4.2 (real data analysis)
MV_OSFA_segment<-function(XD,YD){

  # inputs:
  # XD: the sets of local X segments
  # YD: the sets of local Y segments
  # outputs:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  m=length(YD)
  class=c()     
  for (ii in 1:m){
    DY=YD[[ii]]
    class_ii=unique(DY) 
    class=union(class,class_ii)  
  }
  
  R = length(class)
  p=dim(XD[[1]])[2]
  
  TT = c(); Vw = c(); Pr = rep(0,R)
  theta1 = theta2 = matrix(0,1,p)
  MCP = lapply(rep(0,R), function(h) return(h)) 
  
  for (ii in 1:m) {
    DX=XD[[ii]] 
    DY=YD[[ii]]
    
    ptm <- proc.time()
    
    nii = length(DY)
    w = floor(nii/3);Vw=c(Vw,w)
    
    res_ii = MV_component_matrix(DX,DY,class)
    Pr = Pr + res_ii[[R+1]]
    
    for (r in 1:R){
      MCP[[r]] = MCP[[r]]+res_ii[[r]]*w  
    }
    
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
  } 
  
  ptm <- proc.time()
  
  MV = matrix(1/3,1,p)
  n3 = sum(Vw)
  
  Pr = Pr/sum(Pr) #theta_{y,r}
  
  for (r in 1:R){
    aa = MCP[[r]]/n3     #n3 = sum[n_l/3]
    MV = MV+aa[1,]/Pr[r]-2*aa[2,]
  }
  
  at <- proc.time()-ptm
  
  TC = mean(TT)+at[3]  
  
  list(Value=c(MV),TC=TC)
  
}

# compute the PSIS utility between separate data sets XD and YD via OSFA estimation 
# this function is used in Section 4.2 (real data analysis)
PSIS_OSFA_segment<-function(XD,YD){
  
  # inputs:
  # XD: the sets of local X segments. Note that each X_j in X has been standardized during the implement of PSIS.
  # YD: the sets of local Y segments
  # outputs:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  m=length(YD)
  class=c();     
  for (ii in 1:m){
    DY=YD[[ii]]
    class_ii=unique(DY)
    class=union(class,class_ii) 
  }
  
  R=length(class); 
  p=dim(XD[[1]])[2]
  
  
  Num=matrix(0,m,R) 
  He=matrix(0,R,p) 
  TT=c()
  
  for (ii in 1:m) {
    DX=XD[[ii]] 
    DY=YD[[ii]] 
    
    ptm <- proc.time()
    nii=length(DY);        
    Yheng=matrix(DY,nii,R)
    Yshu=t(matrix(class,R,nii))
    Ycha=Yheng-Yshu 
    
    M=matrix(0,nii,R)
    M[which(Ycha==0)]=1
    num=colSums(M)   
    Num[ii,]=num 
    
    HE_ii=matrix(0,R,p) 
    for (jj in 1:R){
      MY=matrix(rep(M[,jj],p),nrow = nii)
      HE_ii[jj,]=colSums(DX*MY) 
    }
    He=He+HE_ii
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
  }
  
  ptm <- proc.time()
  Nm=colSums(Num); 
  junzhi=He/matrix(rep(Nm,p),R,p) 
  
  PairCha=c()
  for (i1 in 1:(R-1)){
    for (i2 in (i1+1):R){
      aa=abs(junzhi[i1,]-junzhi[i2,])
      PairCha=rbind(PairCha,aa)
    }}
  
  PSIS = apply(PairCha,2,max)
  
  at <- proc.time()-ptm
  
  TC = mean(TT)+at[3]
  
  list(Value=PSIS,TC=TC)
  
}

# compute the FKF utility between separate data sets XD and YD via SWA estimation 
# this function is used in Section 4.2 (real data analysis)
FKF_SWA_segment<-function(XD,YD){
  
  # inputs:
  # XD: the sets of local X segments
  # YD: the sets of local Y segments
  # outputs:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  WM=c();ns=c();TT=c()
  m=length(YD)
  for (ii in 1:m) {
    DX=XD[[ii]] 
    DY=YD[[ii]] 
    nii=length(DY);ns=c(ns,nii)
    
    ptm <- proc.time()
    Wii = FKF_Value(DX,DY)*nii
    ct <- proc.time()-ptm
    TT=c(TT,ct[3])
    
    WM=rbind(WM,Wii)
  }
  
  W_FKF = colSums(WM)/sum(ns)
  
  TC = mean(TT)
  
  list(Value=W_FKF,TC=TC)
  
}

# compute the mean variance (MV) between separate data sets XD and YD via OSFA estimation 
# this function is used in Section 4.2 (real data analysis)
CRU_OSFA_segment<-function(XD,YD){
  
  # inputs:
  # XD: the sets of local X segments
  # YD: the sets of local Y segments
  # outputs:
  # Value: the values of screening utilities corresponding to all features 
  # TC: the computational time 
  
  m=length(YD)
  class=c();nl=c()     
  for (ii in 1:m){
    DY=YD[[ii]] 
    nl=c(nl,length(DY)) 
    class_ii=unique(DY)  
    class=union(class,class_ii)   
  }
  
  R=length(class)
  n=sum(nl)
  p=dim(XD[[1]])[2]
  
  TT_prepare = c()
  
  if (R==1){CRU_vec=rep(0,p);TC = 0}else{
    
    MCPy=0 
    MCPj=0 
    Sw=0
    TT=c()
    
    for (ii in 1:m) {
      DX=XD[[ii]] 
      DY=YD[[ii]] 
      
      ptm <- proc.time()
      nii=length(DY)   
      Yheng=matrix(DY,nii,R)
      Yshu=t(matrix(class,R,nii))
      Ycha=Yheng-Yshu 
      
      M=matrix(0,nii,R)
      M[which(Ycha==0)]=1
      nr=colSums(M) 
      MCPy=MCPy+matrix(rep(nr,p),R,p)/n 
      
      wm = floor(nii/2)
      Sw = Sw+wm
      
      MCPjm=c()
      for (j in 1:p){
        Xj=DX[,j]  
        MX=shixingM_CRU(Xj)
        MXY = t(M)%*%MX
        MXY_mean = wm*rowSums(MXY)/(nii*(nii-1))
        MCPjm = cbind(MCPjm,MXY_mean)
        
      }
      MCPj=MCPj+MCPjm
      ct <- proc.time()-ptm
      TT=c(TT,ct[3])
      
    }
    
    ptm <- proc.time()
    CRU_vec = colSums((MCPj/Sw-.5*MCPy)^2)
    at <- proc.time()-ptm
    TC = mean(TT)+at[3]
  }
  
  list(Value=CRU_vec,TC=TC)
  
}

## "helper" function to obtain the submodel based on screening methods
Eva_K<-function(Omega,K){
  
  # input: 
  # Omega: the screenig utilities, 
  # K: the number of features retained  
  # output: 
  # S_hat_K: the submodel based on the screening
  
  sort_omega<-sort.int(Omega, partial = NULL, na.last = NA, decreasing = TRUE, index.return = TRUE)
  order_omega<-sort_omega$ix 
  S_hat_K=order_omega[1:K]
}

## the classifier based on 20 randomly selected features (RS) 
Benchmark_KNN_pred <- function(K){
  
  # inputs: 
  # K: the number of features randomly selected 
  # outputs:
  # AR: classification accuracy
  # TC: computational time 
  
  ptm <- proc.time() 
  
  Ind_random = sample(1:p,K)
  X_train_K = X_train[,Ind_random]
  X_test_K = X_test[,Ind_random]
  
  ## knn classifier
  res_knn = knn(X_train_K, X_test_K, Y_train, k=KNN_number)
  AR = mean(res_knn == Y_test)
  at <- proc.time()-ptm
  
  TC = at[3]
  
  list(AR=AR,TC=TC)
  
}

# the classifiers based on the distributed screening methods (PSIS, FKF, MV, and CRU)  
Scr_KNN_pred <- function(K,method){
  
  # inputs: 
  # K: the number of features retained by screening method
  # method: indicate which distributed screening method is used  
  # outputs:
  # AR: classification accuracy
  # TC: computational time 
  
  ## screening ##
  if(method == "PSIS"){W = PSIS_OSFA_segment(XD_psis,YD)}
  if(method == "FKF"){W = FKF_SWA_segment(XD,YD)}
  if(method == "MV"){W = MV_OSFA_segment(XD,YD)}
  if(method == "CRU"){W = CRU_OSFA_segment(XD,YD)}
  
  ptm <- proc.time() 
  Ind = Eva_K(W$Value,K)
  print(method)
  X_train_K = X_train[,Ind]
  X_test_K = X_test[,Ind]
  
  ## knn
  res_knn = knn(X_train_K, X_test_K, Y_train, k=201)
  AR = mean(res_knn == Y_test)
  at <- proc.time()-ptm
  
  TC = W$TC+at[3]
  
  list(AR=AR,TC=TC)
  
}

