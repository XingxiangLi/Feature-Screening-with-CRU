
############################################################################################
# This script provides R code that produces the full numerical results shown in the paper. #
############################################################################################
## Last update: February 8 2023
## Require R > 4.0.0

# Packages needed to run this code 
# Package "here" is needed to create relative paths.
# Package "mnormt" is needed to generate the X follow multiple normal distribution with given mean and covariance in Example 3 (Table 6-setup (f)).
# Package "class" is needed to implement KNN classifier in our real data analysis (Table 7).

library(here)     # version 1.0.1
library(mnormt)   # version 2.1.1
library(class)    # version 7.3-20

## The source file "CRU_function.r" provides functions needed to run CRU screening.
function_path = here("scripts","CRU_function.r")
source(function_path)

###################################################
#   Code for Figure 1 in Section 2.1              #
###################################################

## The running time for Figure 1 is approximately 3 minutes on 
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter = 500
n=1000
R=2
D=c(0,0.25,0.5,0.75,1)

M_CRU_cen=c();M_FAIR_cen=c()
for (iter in 1:Num_iter){
  
  R_CRU_cen=c();R_FAIR_cen=c()
  
  for (ii in 1:length(D)){
    
    # data generation
    pp = rep(1/R,R)
    res = dis_gen(n,pp)
    Y = res$Y
    weizhi = res$SR
    
    X=rnorm(n)
    X[weizhi[[1]]]=X[weizhi[[1]]]+D[ii]
    
    # case (a) no outliers
    a_CRU = CRU_index(X,Y);R_CRU_cen=c(R_CRU_cen,a_CRU)
    a_FAIR = FAIR_index(X,Y);R_FAIR_cen=c(R_FAIR_cen,a_FAIR)
    
    # case (b) with outliers
    nn=10
    Y1=c(Y,rep(1,nn))
    X1=c(X,rnorm(nn,20,1))
    a_CRU = CRU_index(X1,Y1);R_CRU_cen=c(R_CRU_cen,a_CRU)
    a_FAIR = FAIR_index(X1,Y1);R_FAIR_cen=c(R_FAIR_cen,a_FAIR)
    
    # case (c) with outliers
    Y2=c(Y,rep(2,nn))
    X2=c(X,rnorm(nn,20,1))
    a_CRU = CRU_index(X2,Y2);R_CRU_cen=c(R_CRU_cen,a_CRU)
    a_FAIR = FAIR_index(X2,Y2);R_FAIR_cen=c(R_FAIR_cen,a_FAIR)
    
  } 
  
  M_CRU_cen=rbind(M_CRU_cen,R_CRU_cen);
  M_FAIR_cen=rbind(M_FAIR_cen,R_FAIR_cen)
  
}

FAIR_res = colMeans(M_FAIR_cen)
FAIR_case_a = FAIR_res[seq(1,15,3)]
FAIR_case_b = FAIR_res[seq(2,15,3)]
FAIR_case_c = FAIR_res[seq(3,15,3)]

CRU_res = colMeans(M_CRU_cen)
CRU_case_a = CRU_res[seq(1,15,3)]
CRU_case_b = CRU_res[seq(2,15,3)]
CRU_case_c = CRU_res[seq(3,15,3)]

## Plotting Figure 1
figure1_path = here("figures","Figure1_full.png")
png(filename=figure1_path)

par(mfcol=c(2,1)) 
method=c("case(a) no outliers","case(b) with outliers","case(c) with outliers")
color = c("red","black","blue")
a1=a2=a3=0.75

L_limit1=c(-0.5,17)
L_limit2=c(-0.001,0.01)

plot(c(1:5),FAIR_case_a,type="b", col="red", pch=8,lty=1, ylim=L_limit1,xlab="u",ylab="Estimates",
     xaxt="n",cex.lab=a1, cex.axis=a1, cex.main=a1, cex.sub=a1)
lines(c(1:5),FAIR_case_b,type="b", col="black", pch=1,lty=2)
lines(c(1:5),FAIR_case_c,type="b", col="blue", pch=2,lty=3)
legend("topleft",legend=method,bty="n", col=color,pch=c(8,1,2),lty=c(1,2,3), ncol=1, cex=a3)
axis(1,at=c(1,2,3,4,5),labels=c("0","0.25","0.5","0.75","1"),cex.axis=a2)
title("FAIR",cex.main=a1)

plot(c(1:5),CRU_case_a,type="b", col="red", pch=8,lty=1, ylim=L_limit2,xlab="u",ylab="Estimates",
     xaxt="n",cex.lab=a1, cex.axis=a1, cex.main=a1, cex.sub=a1)
lines(c(1:5),CRU_case_b,type="b", col="black", pch=1,lty=2)
lines(c(1:5),CRU_case_c,type="b", col="blue", pch=2,lty=3)
legend("topleft",legend=method,bty="n", col=color,pch=c(8,1,2),lty=c(1,2,3), ncol=1, cex=a3)
axis(1,at=c(1,2,3,4,5),labels=c("0","0.25","0.5","0.75","1"),cex.axis=a2)
title("CRU",cex.main=a1)

dev.off()

###################################################
#   Code for Figure 2 in Section 4.1 (Example 1)  #
###################################################

## The running time for Figure 2 is approximately 13 minutes on 
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

p=1000

#### setup I ####

intercept=c(0,0.25,0.5,0.75,1)

## Benchmark: non-distributed (central) estimator 
n=2400
R=2
a_CRU_cen=c(0)
for (i in 2:length(intercept)){
  pp = rep(1/R,R)
  res = dis_gen(n,pp)
  Y = res$Y
  weizhi = res$SR
  X=matrix(rnorm(n*p),n,p)
  X[weizhi[[1]],]=X[weizhi[[1]],]+intercept[i]
  a_CRU = mean(CRU_Value(X,Y))
  a_CRU_cen=c(a_CRU_cen,a_CRU)
}

## distributed estimator  
n=2400
m=240
R=2
R_CRU_SWA=c()
R_CRU_OSFA=c()

for (i in 1:length(intercept)){
  pp=rep(1/R,R)
  res = dis_gen(n,pp)
  Y = res$Y
  weizhi = res$SR
  X=matrix(rnorm(n*p),n,p)
  X[weizhi[[1]],]=X[weizhi[[1]],]+intercept[i]
  Sn = rep((n/m),m)
  M_num = fenzu(Sn)
  R_CRU_SWA=rbind(R_CRU_SWA,CRU_SWA_Value(X,Y,M_num)$Value)
  R_CRU_OSFA=rbind(R_CRU_OSFA,CRU_OSFA_Value(X,Y,M_num)$Value)
}

SetupI_SWA = rowMeans(R_CRU_SWA)
SetupI_OSFA = rowMeans(R_CRU_OSFA)
SetupI_BM = a_CRU_cen

#### setup II ####

n=2400
R=2
X=matrix(rnorm(p*n),n,p)
pp=rep(1/R,R)
res = dis_gen(n,pp)
Y = res$Y

RM_CRU_SWA=c();RM_CRU_OSFA=c();

Mset=c(20,40,80,240,480)

for (m in 1:length(Mset)){
  Sn = rep((n/Mset[m]),Mset[m])
  M_num = fenzu(Sn)
  a1=CRU_SWA_Value(X,Y,M_num)
  RM_CRU_SWA=rbind(RM_CRU_SWA,a1$Value)
  a2=CRU_OSFA_Value(X,Y,M_num)
  RM_CRU_OSFA=rbind(RM_CRU_OSFA,a2$Value)
}

SetupII_SWA = rowMeans(RM_CRU_SWA)
SetupII_OSFA = rowMeans(RM_CRU_OSFA)
SetupII_BM = rep(0,5)

## Plotting Figure 2 ##
figure2_path = here("figures","Figure2_full.png")
png(filename=figure2_path)

par(mfcol=c(2,1)) 
method=c("BM","SWA","OSFA")
color = c("red","blue","black")
a1=a2=a3=0.75

L_limit=c(-0.001,0.014)

plot(c(1:5),SetupI_BM,type="b", col="red", pch=8,lty=1, ylim=L_limit,xlab="u",ylab="Estimates",
     xaxt="n",cex.lab=a1, cex.axis=a1, cex.main=a1, cex.sub=a1)
lines(c(1:5),SetupI_SWA,type="b", col="blue", pch=2,lty=2)
lines(c(1:5),SetupI_OSFA,type="b", col="black", pch=1,lty=3)
legend("topleft",legend=method, col=color,pch=c(8,2,1),lty=c(1,2,3), ncol=3, cex=a3)
axis(1,at=c(1,2,3,4,5),labels=c("0","0.25","0.5","0.75","1"),cex.axis=a2)

plot(c(1:5),SetupII_BM,type="b", col="red", pch=8,lty=1, ylim=L_limit,xlab="m",ylab="Estimates",
     xaxt="n",cex.lab=a1, cex.axis=a1, cex.main=a1, cex.sub=a1)
lines(c(1:5),SetupII_SWA,type="b", col="blue", pch=2,lty=2)
lines(c(1:5),SetupII_OSFA,type="b", col="black", pch=1,lty=3)
legend("topleft",legend=method, col=color,pch=c(8,2,1),lty=c(1,2,3), ncol=3, cex=a3)
axis(1,at=c(1,2,3,4,5),labels=c("20","40","80","240","480"),cex.axis=a2)
dev.off()

###################################################
#   Code for Table 2 in Section 4.1 (Example 1)   #
###################################################

## Table 2 shows the mean computational time (in seconds) of m local computers
## needed to carry out p=1000 repeated CRU estimations. 
## The running time for Table 2 is approximately 5 minutes on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

p=1000
n=2400
R=2
X=matrix(rnorm(p*n),n,p)
pp=rep(1/R,R)
res = dis_gen(n,pp)
Y = res$Y

## computational time of non-distributed (central) estimator
ac=CRU_Central_Time(X,Y)
T_central=ac$TC

T_CRU_SWA=c()
T_CRU_OSFA=c()

## computational time of SWA and OSFA estimator
Mset=c(5,10,20,30,40)
for (m in 1:length(Mset)){
  
  Sn = rep((n/Mset[m]),Mset[m])
  M_num = fenzu(Sn)
  
  a=CRU_SWA_Value(X,Y,M_num)
  T_CRU_SWA=rbind(T_CRU_SWA,a$TC)
  
  aa=CRU_OSFA_Value(X,Y,M_num)
  T_CRU_OSFA=rbind(T_CRU_OSFA,aa$TC)
  
}

TT = c(T_CRU_SWA,T_CRU_OSFA)
TT = t(matrix(TT,length(Mset),2))
Table2 = cbind(rep(T_central,2),TT)

row=c("SWA","OSFA")
column=c("m=1","m=5","m=10","m=20","m=30","m=40")
dimnames(Table2)=list(row,column)
Table2 = round(Table2, 3)

table2_path = here("tables","table2_full.txt")
write.table(Table2,file=table2_path)

###############################################################
#   Code for Table 3 in Section 4.1 (Setup (a) in Example 2)  #
###############################################################

## The running time for Table 3 is approximately 118 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter=100 
p=10000
n=3000
pp1=c(0.5,0.5)
pp2=c(0.3,0.7)
S=c(1:8)
MM=c(30,50)
paux=p/10
intercept=0.35
0

RES=list()

for (ii in 1:Num_iter){
  
  print(ii)
  
  R=c()
  
  ## data generation with pp1
  res = dis_gen(n,pp1) 
  Y = res$Y
  cate_loc = res$SR
  X=matrix(rnorm((p)*n),nrow=n)
  
  for (jj in 1:length(S)){
    X[cate_loc[[1]],S[jj]]=X[cate_loc[[1]],S[jj]]+intercept
  }
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  ## data generation with pp2
  res = dis_gen(n,pp2) 
  Y = res$Y
  cate_loc = res$SR
  X=matrix(rnorm((p)*n),nrow=n)
  
  for (jj in 1:length(S)){
    X[cate_loc[[1]],S[jj]]=X[cate_loc[[1]],S[jj]]+intercept
  }
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

## Summarize the results in Table 3 
rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:16, function(h) return(matrix(0,num_iter,6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,16,6)
for (ii in 1:16){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,9),c(2,10),c(3,11),c(4,12),c(5,13),c(6,14),c(7,15),c(8,16))
Table3 = Summ[wz,]

row=c(rep("FAIR",2),rep("KF",2),rep("MV",2),rep("CRU",2))
row=rep(row,2)
column=c("SSR","PSR","FDR","Size","wRank","Time")
dimnames(Table3)=list(row,column)
Table3 = round(Table3, 3)

table3_path = here("tables","table3_full.txt")
write.table(Table3,file=table3_path)

###############################################################
#   Code for Table 4 in Section 4.1 (Setup (b) in Example 2)  #
###############################################################

## The running time for Table 4 is approximately 58 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter=100
p=10000
n=3000
paux=p/10
pp1=c(0.5,0.5)
pp2=c(0.2,0.8)
S=c(1:4)
MM=c(75,150)
intercept=0.35


RES=list()

for (ii in 1:Num_iter){
  
  print(ii)
  
  R=c()
  
  ## data generation with pp1
  res = dis_gen(n,pp1) 
  Y = res$Y
  cate_loc = res$SR
  
  X=matrix(exp(rnorm(p*n)),nrow=n)
  X[cate_loc[[1]],S[1]]=X[cate_loc[[1]],S[1]]+intercept
  X[cate_loc[[1]],S[2]]=X[cate_loc[[1]],S[2]]+intercept
  X[cate_loc[[2]],S[3]]=X[cate_loc[[2]],S[3]]+intercept
  X[cate_loc[[2]],S[4]]=X[cate_loc[[2]],S[4]]+intercept
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  ## data generation pp2
  res = dis_gen(n,pp2) 
  Y = res$Y
  cate_loc = res$SR
  
  X=matrix(exp(rnorm(p*n)),nrow=n)
  X[cate_loc[[1]],S[1]]=X[cate_loc[[1]],S[1]]+intercept
  X[cate_loc[[1]],S[2]]=X[cate_loc[[1]],S[2]]+intercept
  X[cate_loc[[2]],S[3]]=X[cate_loc[[2]],S[3]]+intercept
  X[cate_loc[[2]],S[4]]=X[cate_loc[[2]],S[4]]+intercept
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

## Summarize the results in Table 4 
rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:16, function(h) return(matrix(0,num_iter,6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,16,6)
for (ii in 1:16){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,9),c(2,10),c(3,11),c(4,12),c(5,13),c(6,14),c(7,15),c(8,16))
Table4 = Summ[wz,]

row=c(rep("FAIR",2),rep("KF",2),rep("MV",2),rep("CRU",2))
row=rep(row,2)
column=c("SSR","PSR","FDR","Size","wRank","Time")
dimnames(Table4)=list(row,column)
Table4 = round(Table4, 3)

table4_path = here("tables","table4_full.txt")
write.table(Table4,file=table4_path)


###############################################################
#   Code for Table 5 in Section 4.1 (Setup (c) in Example 2)  #
###############################################################

## The running time for Table 5 is approximately 141 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter=100
p=10000
n=3000
pp1=c(0.2,0.2,0.2,0.2,0.2)
pp2=c(0.05,0.05,0.05,0.05,0.8)
Rnum=length(pp1)
S=c(1,2,3,4,5)
paux=p/10
MM=c(75,300)
intercept=1.2

RES=list()

for (iter in 1:Num_iter){
  
  print(iter)
  
  R=c()
  
  ## data generation with pp1
  res = dis_gen(n,pp1)
  Y = res$Y
  cate_loc = res$SR
  
  X=matrix(rt((p)*n,2),nrow=n)
  for (ii in 1:Rnum){
    X[cate_loc[[ii]],S[ii]]=X[cate_loc[[ii]],S[ii]]+intercept
  }
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_PSIS)
  R=rbind(R,Res$RR_FKF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_PSIS)
  R=rbind(R,Res$RR_FKF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  ## data generation with pp2
  res = dis_gen(n,pp2) 
  Y = res$Y
  cate_loc = res$SR
  
  X=matrix(rt((p)*n,2),nrow=n)
  for (ii in 1:Rnum){
    X[cate_loc[[ii]],S[ii]]=X[cate_loc[[ii]],S[ii]]+intercept
  }
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_PSIS)
  R=rbind(R,Res$RR_FKF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr(X,Y,S,paux,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_PSIS)
  R=rbind(R,Res$RR_FKF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

## Summarize the results in Table 5 
rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:16, function(h) return(matrix(0,num_iter,6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,16,6)
for (ii in 1:16){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,9),c(2,10),c(3,11),c(4,12),c(5,13),c(6,14),c(7,15),c(8,16))
Table5 = Summ[wz,]

row=c(rep("PSIS",2),rep("FKF",2),rep("MV",2),rep("CRU",2))
row=rep(row,2)
column=c("SSR","PSR","FDR","Size","wRank","Time")
dimnames(Table5)=list(row,column)
Table5 = round(Table5, 3)

table5_path = here("tables","table5_full.txt")
write.table(Table5,file=table5_path)


###############################################################
#   Code for Setup (d) in Example 3 (Table 6 in Section 4.1)  #
###############################################################

## The running time for Setup (d) is approximately 15 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter=100
p=5000
n=2980
Num_noise = 20
q=8
MM=c(50,100)
K=50

RES = list()

for (ii in 1:Num_iter){
  
  print(ii)
  
  R=c()
  
  res = logistic_inde_dg(n, p, q, Beta="default", inc=0, index=-1, sig=1, tag="B")
  Y=res$Y
  X=res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  #add some noise
  noise_value=runif(Num_noise*p)*20
  X_noise=matrix(noise_value,Num_noise,p)
  Y_noise=rep(0,Num_noise)
  
  X = rbind(X,X_noise)
  Y = c(Y,Y_noise)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:8, function(h) return(matrix(0,num_iter,5)))

for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,8,5)
for (ii in 1:8){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,5),c(2,6),c(3,7),c(4,8))
Table6_d = Summ[wz,]

###############################################################
#   Code for Setup (e) in Example 3 (Table 6 in Section 4.1)  #
###############################################################

## The running time for Setup (e) is approximately 17 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

Num_iter=100
p=5000
n=2980
Num_noise = 20
q=5
MM=c(100,200)
K=50

RES = list()
for (ii in 1:Num_iter){
  
  print(ii)
  
  R=c()
  
  res = logistic_moav_dg(n, p, q, Beta="default", inc=0, index=-1, sig=1, tag="B")
  Y=res$Y
  X=res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  #add some noise
  noise_value=runif(Num_noise*p)*20
  X_noise=matrix(noise_value,Num_noise,p)
  Y_noise=rep(0,Num_noise)
  
  X = rbind(X,X_noise)
  Y = c(Y,Y_noise)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:8, function(h) return(matrix(0,num_iter,5)))

for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,8,5)
for (ii in 1:8){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,5),c(2,6),c(3,7),c(4,8))
Table6_e = Summ[wz,]


###############################################################
#   Code for Setup (f) in Example 3 (Table 6 in Section 4.1)  #
###############################################################

## The running time for Setup (f) is approximately 25 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.


Num_iter=100
p=5000
p1=1000
n=2980
Num_noise=20
q=6
MM=c(250,500)
K=50

RES = list()
for (ii in 1:Num_iter){
  
  print(ii)
  
  R=c()
  
  res = logistic_comp_dg(n, p, p1, q, Beta=c(0.4,0.4,0.4,1,1,2), index=-1, rr=0.4, tag="B")
  Y=res$Y
  X=res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  #add some noise
  noise_value=runif(Num_noise*p)*20
  X_noise=matrix(noise_value,Num_noise,p)
  Y_noise=rep(0,Num_noise)
  
  X = rbind(X,X_noise)
  Y = c(Y,Y_noise)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[1],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  Res=Dis_Scr_K(X,Y,S,K,m=MM[2],balance=1)
  
  R=rbind(R,Res$RR_FAIR)
  R=rbind(R,Res$RR_KF)
  R=rbind(R,Res$RR_MV)
  R=rbind(R,Res$RR_CRU)
  
  RES=c(RES,list(R))
  
}

rr1 = RES
num_iter = length(rr1)
RESf<-lapply(1:8, function(h) return(matrix(0,num_iter,5)))

for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii,]=rr1[[ii]][jj,]
  } 
}

Summ = matrix(0,8,5)
for (ii in 1:8){
  Summ[ii,] = Summarize(RESf[[ii]])
}

wz = c(c(1,5),c(2,6),c(3,7),c(4,8))
Table6_f = Summ[wz,]

# Combine the above Table6_d, Table6_e and Table6_f to be Table 6

Table6 = rbind(Table6_d,Table6_e)
Table6 = rbind(Table6,Table6_f)

row=c(rep("FAIR",2),rep("KF",2),rep("MV",2),rep("CRU",2))
row=rep(row,3)
column=c("SSR","PSR","FDR","wRank","Time")
dimnames(Table6)=list(row,column)

Table6 = round(Table6, 3)

table6_path = here("tables","table6_full.txt")
write.table(Table6,file=table6_path)

#############################################################################
#           Code for Table 7 in real data analysis (Section 4.2)            #
#############################################################################

## The running time for Table 7 is approximately 35 hours on
## the computers with Intel Core i7-8700 3.2Hz + 32 GB RAM.

## The dataset is available at https://archive.ics.uci.edu/ml/datasets/Multiple+Features and contain 6 feature matrices.
## More details of features can be found in "M. van Breukelen, R.P.W. Duin, D.M.J. Tax, and J.E. den Hartog, Handwritten  
## digit recognition by combined classifiers, Kybernetika, vol. 34, no. 4,1998, 381-386."
## Each digit is represented in terms of the following six feature sets.
## Regarding more details of data dictionary, please see README.txt.

Num_iter=100

X_fac = as.matrix(read.table(here("data","fac.txt"),header=T))
X_fou = as.matrix(read.table(here("data","fou.txt"),header=T))
X_kar = as.matrix(read.table(here("data","kar.txt"),header=T))
X_mor = as.matrix(read.table(here("data","mor.txt"),header=T))
X_pix = as.matrix(read.table(here("data","pix.txt"),header=T))
X_zer = as.matrix(read.table(here("data","zer.txt"),header=T))

# X generation
X_raw = cbind(X_fac,X_fou)
X_raw = cbind(X_raw,X_kar)
X_raw = cbind(X_raw,X_mor)
X_raw = cbind(X_raw,X_pix)
X_raw = cbind(X_raw,X_zer)

# Y generation
Y=c()
for(k in 0:9){Y = c(Y,rep(k,200))}

# some setups
balance = 1
K=20
nl=10
p1=1351 
KNN_number=51 
noise_number_MM=c(0,1000,2000,3000)

RE_benchmark = c()
RE_psis = c()
RE_fkf = c()
RE_mv = c()
RE_cru = c()  

## repetitions
for (iter in 1:Num_iter){
  
  set.seed(iter)
  
  print(iter)
  ## add some irrelevant features
  Xar = matrix(rt(p1*2000,2),2000,p1)
  X_new = cbind(X_raw,Xar)
  X_new = scale(X_new)
  p = dim(X_new)[2]
  
  # training and testing set
  test_index=c()
  for(k in 0:9){test_index = c(test_index,k*200+sort(sample(200,20)))} #n_train=1800, n_test=200
  train_index=setdiff(1:2000,test_index)
  X_train = X_new[train_index,]
  X_test = X_new[test_index,]
  Y_train = Y[train_index]
  Y_test = Y[test_index]
  
  X_train_raw = X_train
  Y_train_raw = Y_train
  
  RR_benchmark = c()
  RR_psis = c()
  RR_fkf = c()
  RR_mv = c()
  RR_cru = c()  
  
  for (noise_number in noise_number_MM){
    
    if (noise_number==0){ ## the analysis on the original data set
      print(paste("N1=",noise_number,sep=""))
      n_train = length(Y_train)
      X_train = X_train+0.0001*matrix(runif(n_train*p),n_train,p) ## differentiate the same values in some features
      X_train_psis = scale(X_train)
    }
    
    if (noise_number>0){
      print(paste("N1=",noise_number,sep=""))
      ## add some outliers into training set
      Xout = matrix(runif(noise_number*p),noise_number,p)*1000
      Yout = sample(0:9,noise_number,replace = TRUE)
      X_train = rbind(X_train,Xout)
      Y_train = c(Y_train,Yout)
      n_train = length(Y_train)
      X_train = X_train+0.0001*matrix(runif(n_train*p),n_train,p) ## differentiate the same values 
      X_train_psis = scale(X_train)
    }
    
    ## divide the data set for distributed screening
    m=n_train/nl
    if (balance==1){
      Sn = rep(floor(n_train/m),m);
      Nr = n_train-sum(Sn)
      if (Nr>0){Sn[1:Nr]=Sn[1:Nr]+1}
      M_num = fenzu(Sn) 
    }
    
    XD=c();YD=c();XD_psis=c()
    for(l in 1:m)
    {
      Xl=X_train[M_num[[l]],]
      Yl=Y_train[M_num[[l]]]
      Xl_psis=X_train_psis[M_num[[l]],]
      XD=c(XD,list(Xl))   
      YD=c(YD,list(Yl))
      XD_psis=c(XD_psis,list(Xl_psis))
    }
    
    ## RS
    res_benchmark = Benchmark_KNN_pred(K)
    RR_benchmark = c(RR_benchmark,res_benchmark$AR)
    
    ## PSIS  
    res_psis = Scr_KNN_pred(K,method="PSIS")
    RR_psis = c(RR_psis,res_psis$AR)
    
    ## FKF 
    res_fkf = Scr_KNN_pred(K,method="FKF")
    RR_fkf = c(RR_fkf,res_fkf$AR)
    
    ## MV 
    res_mv = Scr_KNN_pred(K,method="MV")
    RR_mv = c(RR_mv,res_mv$AR)
    
    ## CRU 
    res_cru = Scr_KNN_pred(K,method="CRU")
    RR_cru = c(RR_cru,res_cru$AR)
    
    X_train = X_train_raw
    Y_train = Y_train_raw 
    
  } 
  
  RE_benchmark = rbind(RE_benchmark,RR_benchmark)
  RE_psis = rbind(RE_psis,RR_psis)
  RE_fkf = rbind(RE_fkf,RR_fkf)
  RE_mv = rbind(RE_mv,RR_mv)
  RE_cru = rbind(RE_cru,RR_cru)
  
  a1=list(RE_benchmark,RE_psis,RE_fkf,RE_mv,RE_cru)
  
}

Table7=c()
for (ii in 1:5){
  res = colMeans(a1[[ii]])
  Table7 = rbind(Table7,res)
}

row=c("RS","PSIS","FKF","MV","CRU")
column=c("N_1=0","N_1=1000","N_1=2000","N_1=3000")
dimnames(Table7)=list(row,column)

Table7 = round(Table7, 3)

table7_path = here("tables","table7_full.txt")
write.table(Table7,file=table7_path)



