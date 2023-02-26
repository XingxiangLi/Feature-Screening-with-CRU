
#########################################################################
## R code to fully reproduce the numerical results shown in "Feature    #
## Screening with Conditional Rank Utility for Big-data Classification" #
#########################################################################

## Last update: February 25 2023
## Require R >= 4.2.1

## Runtime to construct each figure or table is given at the beginning of the
## relevant code sections; the reported runtimes are based on a Windows computer 
## with 3.2 GHz CPUs and 32 GB memory. 


###################
## Preliminaries ##
###################

## Packages needed to run this code:
## Package "here" - run the code with relative paths linked to the top-level directory 
## Package "mnormt" - generate random sample from a multi-normal distribution (Table 6)
## Package "class" - implement KNN classifier (Table 7)
## Package "writexl" - export Tables 2-7 in ".xlsx" format

## Uncomment the command below to install the required packages
# install.packages(c("here", "mnormt", "class", "writexl"))

## Load packages
library(here)     # version 1.0.1
library(mnormt)   # version 2.1.1
library(class)    # version 7.3-20
library(writexl)  # version 1.4.2

## Source in functions needed to run this code
source(here("scripts", "CRU_function.r"))


#######################################
## Code for Figure 1 in Section 2.1  ##
## Runtime: about 3 minutes          ##
#######################################

print("---Figure 1---")

## Initial setup
n = 1000
R = 2
D = c(0, 0.25, 0.5, 0.75, 1) 
M_CRU_cen = c() 
M_FAIR_cen = c()

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 500

## Begin repetitions
for (iter in 1:Num_iter){
  
  R_CRU_cen = c() 
  R_FAIR_cen = c()
  
  ## loop for utility estimates
  for (ii in 1:length(D)){
    
    ## Data generation   
    pp = rep(1/R, R) 
    res = dis_gen(n, pp) 
    Y = res$Y
    weizhi = res$SR
    
    X = rnorm(n) 
    X[weizhi[[1]]] = X[weizhi[[1]]]+D[ii]
    
    ## Case (a) no outliers
    a_CRU = CRU_index(X, Y)
    R_CRU_cen = c(R_CRU_cen, a_CRU)
    a_FAIR = FAIR_index(X, Y)
    R_FAIR_cen = c(R_FAIR_cen, a_FAIR)
    
    ## Case (b) with outliers
    nn = 10
    Y1 = c(Y, rep(1, nn))
    X1 = c(X, rnorm(nn, 20, 1))
    a_CRU = CRU_index(X1, Y1)
    R_CRU_cen = c(R_CRU_cen, a_CRU)
    a_FAIR = FAIR_index(X1, Y1)
    R_FAIR_cen = c(R_FAIR_cen, a_FAIR)
    
    ## Case (c) with outliers
    Y2 = c(Y, rep(2, nn))
    X2 = c(X, rnorm(nn, 20, 1))
    a_CRU = CRU_index(X2, Y2)
    R_CRU_cen = c(R_CRU_cen, a_CRU)
    a_FAIR = FAIR_index(X2, Y2)
    R_FAIR_cen = c(R_FAIR_cen, a_FAIR)    
  } 
  
  M_CRU_cen = rbind(M_CRU_cen, R_CRU_cen)
  M_FAIR_cen = rbind(M_FAIR_cen, R_FAIR_cen)
}

## FAIR estimates
FAIR_res = colMeans(M_FAIR_cen)
FAIR_case_a = FAIR_res[seq(1, 15, 3)]
FAIR_case_b = FAIR_res[seq(2, 15, 3)]
FAIR_case_c = FAIR_res[seq(3, 15, 3)]

## CRU estimates
CRU_res = colMeans(M_CRU_cen)
CRU_case_a = CRU_res[seq(1, 15, 3)]
CRU_case_b = CRU_res[seq(2, 15, 3)]
CRU_case_c = CRU_res[seq(3, 15, 3)]

## Output Figure 1 in subfolder "figures"
pdf(here("figures", "fig1_full.pdf"), width = 6, height = 5)

par(mfcol = c(2, 1))
par(mgp = c(2, 1, 0), mar = c(3.5, 3, 1.2, 1))
method = c("case(a) no outliers", "case(b) with outliers", "case(c) with outliers")
color = c("red", "black", "blue")
a1 = a2 = a3 = 0.75

L_limit1 = c(-0.5, 17)
L_limit2 = c(-0.001, 0.01)

plot(c(1:5), FAIR_case_a, type = "b", col = "red", pch = 8, lty = 1, ylim = L_limit1, xlab = expression(mu), ylab = "Estimates", 
     xaxt = "n", cex.lab = 0.8, cex.axis = a1, cex.main = a1, cex.sub = a1)
lines(c(1:5), FAIR_case_b, type = "b", col = "black", pch = 1, lty = 2)
lines(c(1:5), FAIR_case_c, type = "b", col = "blue", pch = 2, lty = 3)
legend("topleft", legend = method, bty = "n", col = color, pch = c(8, 1, 2), lty = c(1, 2, 3), ncol = 1, cex = a3)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = a2)
title("FAIR", cex.main = a1)

plot(c(1:5), CRU_case_a, type = "b", col = "red", pch = 8, lty = 1, ylim = L_limit2, xlab = expression(mu), ylab = "Estimates",
     xaxt = "n", cex.lab = 0.8, cex.axis = a1, cex.main = a1, cex.sub = a1)
lines(c(1:5), CRU_case_b, type = "b", col = "black", pch = 1, lty = 2)
lines(c(1:5), CRU_case_c, type = "b", col = "blue", pch = 2, lty = 3)
legend("topleft", legend = method, bty = "n", col = color, pch = c(8, 1, 2), lty = c(1, 2, 3), ncol = 1, cex = a3)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = a2)
title("CRU", cex.main = a1)

dev.off()


##################################################
## Code for Figure 2 in Section 4.1 (Example 1) ##
## Runtime: about 15 minutes                    ##
##################################################

print("---Figure 2---")

## Setup I ## 

## Initialization
n = 2400
R = 2
m = 240
intercept = c(0, 0.25, 0.5, 0.75, 1) 

## Vectors to store CRU estimates  
a_CRU_cen = c()
R_CRU_SWA = c()
R_CRU_OSFA = c()

## Set seed and number of features
## Reduce p for faster runtimes if needed
set.seed(1)
p = 1000

for (i in 1:length(intercept)){
  
  ## Generate response Y
  pp = rep(1/R, R)  
  res = dis_gen(n, pp) 
  Y = res$Y
  
  ## Generate feature X
  weizhi = res$SR 
  X = matrix(rnorm(n*p), n, p) 
  X[weizhi[[1]], ] = X[weizhi[[1]], ]+intercept[i] 
  
  ## Non-distributed estimation
  a_CRU = mean(CRU_Value(X, Y)$Value) 
  a_CRU_cen = c(a_CRU_cen, a_CRU)
  
  ## Distributed estimation
  Sn = rep((n/m), m)
  M_num = fenzu(Sn) 
  R_CRU_SWA = rbind(R_CRU_SWA, CRU_SWA_Value(X, Y, M_num)$Value)
  R_CRU_OSFA = rbind(R_CRU_OSFA, CRU_OSFA_Value(X, Y, M_num)$Value)
}

## Summarize results
SetupI_SWA = rowMeans(R_CRU_SWA)
SetupI_OSFA = rowMeans(R_CRU_OSFA)
SetupI_BM = a_CRU_cen


## Setup II ##

## Data generation
n = 2400
R = 2
X = matrix(rnorm(p*n), n, p)
pp = rep(1/R, R)
res = dis_gen(n, pp)
Y = res$Y

## Vectors to store distributed estimates 
RM_CRU_SWA = c()
RM_CRU_OSFA = c()

## Set number of data segments
Mset = c(20, 40, 80, 240, 480) 

## Distributed CRU estimation with Mset
for (m in 1:length(Mset)){
  Sn = rep((n/Mset[m]), Mset[m])
  M_num = fenzu(Sn) 
  a1 = CRU_SWA_Value(X, Y, M_num)
  RM_CRU_SWA = rbind(RM_CRU_SWA, a1$Value)
  a2 = CRU_OSFA_Value(X, Y, M_num)
  RM_CRU_OSFA = rbind(RM_CRU_OSFA, a2$Value)
}

## Summarize results
SetupII_SWA = rowMeans(RM_CRU_SWA)
SetupII_OSFA = rowMeans(RM_CRU_OSFA)
SetupII_BM = rep(0, 5)

## Output Figure 2 in subfolder "figures"
pdf(here("figures", "fig2_full.pdf"), width = 6, height = 5)

par(mfcol = c(2, 1)) 
par(mgp = c(2, 1, 0), mar=c(3.5, 3, 1.2, 1))
method = c("BM", "SWA", "OSFA")
color = c("red", "blue", "black")
a1 = a2 = a3 = 0.75
L_limit = c(-0.001, 0.014)

plot(c(1:5), SetupI_BM, type = "b", col = "red", pch = 8, lty = 1, ylim = L_limit, xlab = expression(mu), ylab = "Estimates",
     xaxt = "n", cex.lab = 0.8, cex.axis = a1, cex.main = a1, cex.sub = a1)
lines(c(1:5), SetupI_SWA, type = "b", col = "blue", pch = 2, lty = 2)
lines(c(1:5), SetupI_OSFA, type = "b", col = "black", pch = 1, lty = 3)
legend("topleft", legend = method, bty = "n", col = color, pch = c(8, 2, 1), lty = c(1, 2, 3), ncol = 1, cex = a3)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = a2)

plot(c(1:5), SetupII_BM, type = "b", col = "red", pch = 8, lty = 1, ylim = L_limit, xlab = expression(m), ylab = "Estimates",
     xaxt = "n", cex.lab = 0.8, cex.axis = a1, cex.main = a1, cex.sub = a1)
lines(c(1:5), SetupII_SWA, type = "b", col = "blue", pch = 2, lty = 2)
lines(c(1:5), SetupII_OSFA, type = "b", col = "black", pch = 1, lty = 3)
legend("topleft", legend = method, bty = "n", col = color, pch = c(8, 2, 1), lty = c(1, 2, 3), ncol = 1, cex = a3)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("20", "40", "80", "240", "480"), cex.axis = a2)

dev.off()


###################################################
##   Code for Table 2 in Section 4.1 (Example 1) ##
##   Runtime: about 4 minutes                    ##
###################################################

print("---Table 2---")

## Set sample size and number of categories 
n = 2400
R = 2

## Set seed and number of features
## Reduce p for faster runtimes if needed
set.seed(1)
p = 1000

## Data generation
X = matrix(rnorm(p*n), n, p)
pp = rep(1/R, R) 
res = dis_gen(n, pp) 
Y = res$Y

## Record runtime for non-distributed estimator
ac = CRU_Value(X, Y)
T_central = ac$TC

## Initialization for distributed estimation
Mset = c(5, 10, 20, 30, 40)
T_CRU_SWA = c()
T_CRU_OSFA = c()

## Record runtime for distributed estimators
for (m in 1:length(Mset)){
  
  Sn = rep((n/Mset[m]), Mset[m])
  M_num = fenzu(Sn) 
  
  t_swa = CRU_SWA_Value(X, Y, M_num)
  T_CRU_SWA = rbind(T_CRU_SWA, t_swa$TC)
  
  t_osfa = CRU_OSFA_Value(X, Y, M_num)
  T_CRU_OSFA = rbind(T_CRU_OSFA, t_osfa$TC)
  
}

## Output Table 2 in subfolder "tables"
TT = c(T_CRU_SWA, T_CRU_OSFA)
TT = t(matrix(TT, length(Mset), 2))
Table2 = cbind(rep(T_central, 2), TT)
Table2 = round(Table2, 3)
colnames(Table2) = c("m=1", "m=5", "m=10", "m=20", "m=30", "m=40")
Estimation = c("SWA", "OSFA")
Table2 = cbind(Estimation, Table2)
Table2 = as.data.frame(Table2)
write_xlsx(Table2, here("tables", "table2_full.xlsx"))


###############################################################
##  Code for Table 3 in Section 4.1 (Example 2 - Setup (a))  ##
##  Runtime: about 126 hours                                 ##
###############################################################

print("---Table 3---")

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 100

## Set number of features and observations
## Reduce p for faster runtimes if needed
p = 10000
n = 3000 

## Model setup
pp1 = c(0.5, 0.5) 
pp2 = c(0.3, 0.7)
S = c(1:8) 
MM = c(30, 50)
paux = p/10 
intercept = 0.35 
RES = list()

## Begin repetitions
for (ii in 1:Num_iter){
  
  print(ii)
  R = c()
  
  ## Balanced classes ##
  
  ## Generate response 
  res = dis_gen(n, pp1) 
  Y = res$Y
  cate_loc = res$SR
  
  ## Generate feature matrix
  X = matrix(rnorm((p)*n), nrow = n) 
  for (jj in 1:length(S)){
    X[cate_loc[[1]], S[jj]] = X[cate_loc[[1]], S[jj]]+intercept
  }
  
  ## Distributed feature screening with m = 30
  Res = Dis_Scr(X, Y, S, paux, m = MM[1])   
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 50
  Res = Dis_Scr(X, Y, S, paux, m = MM[2])  
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Unbalanced classes ##
  
  ## Generate response
  res = dis_gen(n, pp2) 
  Y = res$Y
  cate_loc = res$SR
  
  ##  Generate feature matrix 
  X = matrix(rnorm((p)*n), nrow = n)
  for (jj in 1:length(S)){
    X[cate_loc[[1]], S[jj]] = X[cate_loc[[1]], S[jj]]+intercept
  }
  
  ## Distributed feature screening with m = 30
  Res = Dis_Scr(X, Y, S, paux, m = MM[1]) 
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 50
  Res = Dis_Scr(X, Y, S, paux, m = MM[2]) 
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R))  
}

## Summarize the results
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:16, function(h) return(matrix(0, num_iter, 6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}

Summ = matrix(0, 16, 6)
for (ii in 1:16){
  Summ[ii, ] = colMeans(RESf[[ii]])
}

wz = c(c(1, 9), c(2, 10), c(3, 11), c(4, 12), c(5, 13), c(6, 14), c(7, 15), c(8, 16))
RES3 = Summ[wz, ]

## Output Table 3 in subfolder "tables"
RES3[, 1:5] = round(RES3[, 1:5], 2)
RES3[, 6] = round(RES3[, 6], 4)
RES3[seq(2, 16, 2), 6]=""
colnames(RES3) = c("SSR", "PSR", "FDR", "Size", "wRank", "Time")
m = c(MM[1], rep(" ", 7), MM[2], rep(" ", 7))
Utility = c("FAIR", " ", "KF", " ", "MV", " ", "CRU", " ")
Utility = rep(Utility, 2)
pi_r = c("equal", "unequal")
pi_r = rep(pi_r, 8)
Table3 = cbind(m, Utility, pi_r, RES3)
Table3 = as.data.frame(Table3)
write_xlsx(Table3, here("tables", "table3_full.xlsx"))


###############################################################
##  Code for Table 4 in Section 4.1 (Example 2 - Setup (b))  ##
##  Runtime: about 62 hours                                  ##
###############################################################

print("---Table 4---")

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 100

## Set number of features and observations
## Reduce p for faster runtimes if needed
p = 10000
n = 3000 

## Model setup
pp1 = c(0.5, 0.5) 
pp2 = c(0.2, 0.8)
S = c(1:4) 
MM = c(75, 150) 
paux = p/10 
intercept = 0.35 
RES = list()

## Begin repetitions
for (ii in 1:Num_iter){
  
  print(ii)  
  R = c()
  
  ## Balanced classes ##
  
  ## Generate response
  res = dis_gen(n, pp1) 
  Y = res$Y
  cate_loc = res$SR
  
  ## Generate feature matrix
  X = matrix(exp(rnorm(p*n)), nrow = n)
  X[cate_loc[[1]], S[1]] = X[cate_loc[[1]], S[1]]+intercept
  X[cate_loc[[1]], S[2]] = X[cate_loc[[1]], S[2]]+intercept
  X[cate_loc[[2]], S[3]] = X[cate_loc[[2]], S[3]]+intercept
  X[cate_loc[[2]], S[4]] = X[cate_loc[[2]], S[4]]+intercept
  
  ## Distributed feature screening with m = 75
  Res = Dis_Scr(X, Y, S, paux, m = MM[1])  
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 150
  Res = Dis_Scr(X, Y, S, paux, m = MM[2])  
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Unbalanced classes ##
  
  ## Generate response
  res = dis_gen(n, pp2) 
  Y = res$Y
  cate_loc = res$SR
  
  ## Generate feature matrix 
  X = matrix(exp(rnorm(p*n)), nrow = n)
  X[cate_loc[[1]], S[1]] = X[cate_loc[[1]], S[1]]+intercept
  X[cate_loc[[1]], S[2]] = X[cate_loc[[1]], S[2]]+intercept
  X[cate_loc[[2]], S[3]] = X[cate_loc[[2]], S[3]]+intercept
  X[cate_loc[[2]], S[4]] = X[cate_loc[[2]], S[4]]+intercept
  
  ## Distributed feature screening with m = 75
  Res = Dis_Scr(X, Y, S, paux, m = MM[1])  
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 150
  Res = Dis_Scr(X, Y, S, paux, m = MM[2]) 
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R))
  
}

## Summarize the results
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:16, function(h) return(matrix(0, num_iter, 6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}

Summ = matrix(0, 16, 6)
for (ii in 1:16){
  Summ[ii, ] = colMeans(RESf[[ii]])
}

wz = c(c(1, 9), c(2, 10), c(3, 11), c(4, 12), c(5, 13), c(6, 14), c(7, 15), c(8, 16))
RES4 = Summ[wz, ]

## Output Table 4 in subfolder "tables"
RES4[, 1:5] = round(RES4[, 1:5], 2)
RES4[, 6] = round(RES4[, 6], 4)
RES4[seq(2, 16, 2), 6]=""
colnames(RES4) = c("SSR", "PSR", "FDR", "Size", "wRank", "Time")
m = c(MM[1], rep(" ", 7), MM[2], rep(" ", 7))
Utility = c("FAIR", " ", "KF", " ", "MV", " ", "CRU", " ")
Utility = rep(Utility, 2)
pi_r = c("equal", "unequal")
pi_r = rep(pi_r, 8)
Table4 = cbind(m, Utility, pi_r, RES4)
Table4 = as.data.frame(Table4)
write_xlsx(Table4, here("tables", "table4_full.xlsx"))



##############################################################
## Code for Table 5 in Section 4.1 (Example 2 - Setup (c))  ##
## Runtime: about 157 hours                                 ##
##############################################################

print("---Table 5---")

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 100

## Set number of features and observations
## Reduce p for faster runtimes if needed
p = 10000
n = 3000 

## Model setup
pp1 = c(0.2, 0.2, 0.2, 0.2, 0.2) 
pp2 = c(0.05, 0.05, 0.05, 0.05, 0.8)
Rnum = length(pp1)
S = c(1, 2, 3, 4, 5)
MM = c(75, 300)
paux = p/10 
intercept = 1.2 
RES = list()

## Begin repetitions
for (iter in 1:Num_iter){
  
  print(iter)  
  R = c()
  
  ## Balanced classes ##
  
  ## Generate response 
  res = dis_gen(n, pp1)
  Y = res$Y
  cate_loc = res$SR
  
  ## Generate feature matrix 
  X = matrix(rt((p)*n, 2), nrow = n)
  for (ii in 1:Rnum){
    X[cate_loc[[ii]], S[ii]] = X[cate_loc[[ii]], S[ii]]+intercept
  }
  
  ## Distributed feature screening with m = 75
  Res = Dis_Scr(X, Y, S, paux, m = MM[1])  
  R = rbind(R, Res$RR_PSIS)
  R = rbind(R, Res$RR_FKF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 300
  Res = Dis_Scr(X, Y, S, paux, m = MM[2]) 
  R = rbind(R, Res$RR_PSIS)
  R = rbind(R, Res$RR_FKF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Unbalanced classes ##
  
  ## Generate response 
  res = dis_gen(n, pp2) 
  Y = res$Y
  cate_loc = res$SR
  
  ## Generate feature matrix  
  X = matrix(rt((p)*n, 2), nrow = n)
  for (ii in 1:Rnum){
    X[cate_loc[[ii]], S[ii]] = X[cate_loc[[ii]], S[ii]]+intercept
  }
  
  ## Distributed feature screening with m = 75
  Res = Dis_Scr(X, Y, S, paux, m = MM[1]) 
  R = rbind(R, Res$RR_PSIS)
  R = rbind(R, Res$RR_FKF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 300
  Res = Dis_Scr(X, Y, S, paux, m = MM[2]) 
  R = rbind(R, Res$RR_PSIS)
  R = rbind(R, Res$RR_FKF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R))
  
}

## Summarize the results 
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:16, function(h) return(matrix(0, num_iter, 6)))

for (ii in 1:num_iter){
  for (jj in 1:16){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}

Summ = matrix(0, 16, 6)
for (ii in 1:16){
  Summ[ii, ] = colMeans(RESf[[ii]])
}

wz = c(c(1, 9), c(2, 10), c(3, 11), c(4, 12), c(5, 13), c(6, 14), c(7, 15), c(8, 16))
RES5 = Summ[wz, ]

## Output Table 5 in subfolder "tables"
RES5[, 1:5] = round(RES5[, 1:5], 2)
RES5[, 6] = round(RES5[, 6], 4)
RES5[seq(2, 16, 2), 6]=""
colnames(RES5) = c("SSR", "PSR", "FDR", "Size", "wRank", "Time")
m = c(MM[1], rep(" ", 7), MM[2], rep(" ", 7))
Utility = c("PSIS", " ", "FKF", " ", "MV", " ", "CRU", " ")
Utility = rep(Utility, 2)
pi_r = c("equal", "unequal")
pi_r = rep(pi_r, 8)
Table5 = cbind(m, Utility, pi_r, RES5)
Table5 = as.data.frame(Table5)
write_xlsx(Table5, here("tables", "table5_full.xlsx"))


#########################################
## Code for Table 6 in Section 4.1     ##
## Runtime: about 48 hours             ##
#########################################

print("---Table 6---")

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 100

## Set number of features and observations
## Reduce p for faster runtimes if needed
p = 5000
n = 2980 

## Setup (d) ##

## Initialization
Num_noise = 20 
q = 8 
MM = c(50,100) 
K = 50 
RES = list()

## Begin repetitions
for (ii in 1:Num_iter){
  
  print(ii)  
  R = c()
  
  ## Data generation
  res = logistic_inde_dg(n, p, q)
  Y = res$Y
  X = res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  ## Add noisy observations
  noise_value = runif(Num_noise*p)*20
  X_noise = matrix(noise_value, Num_noise, p)
  Y_noise = rep(0, Num_noise)
  X = rbind(X, X_noise)
  Y = c(Y, Y_noise)
  
  ## Distributed feature screening with m = 50
  Res = Dis_Scr_K(X, Y, S, K, m = MM[1])  
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 100
  Res = Dis_Scr_K(X, Y, S, K, m = MM[2]) 
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R)) 
}

## Summarize the results of Setup (d)
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:8, function(h) return(matrix(0, num_iter, 5)))
for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}
Summ = matrix(0, 8, 5)
for (ii in 1:8){
  Summ[ii, ] = colMeans(RESf[[ii]])
}
wz = c(c(1, 5), c(2, 6), c(3, 7), c(4, 8))
Table6_d = Summ[wz, ]


## Setup (e) ##

## Initialization
Num_noise = 20 
q = 5
MM = c(100, 200) 
K = 50 
RES = list()

# Begin repetitions
for (ii in 1:Num_iter){
  
  print(ii) 
  R = c()
  
  ## Data generation
  res = logistic_moav_dg(n, p, q)
  Y = res$Y
  X = res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  ## Add noisy observations
  noise_value = runif(Num_noise*p)*20
  X_noise = matrix(noise_value, Num_noise, p)
  Y_noise = rep(0, Num_noise)  
  X = rbind(X, X_noise)
  Y = c(Y, Y_noise)
  
  ## Distributed feature screening with m = 100
  Res = Dis_Scr_K(X, Y, S, K, m = MM[1])
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 200
  Res = Dis_Scr_K(X, Y, S, K, m = MM[2])
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R))
  
}

## Summarize the results of Setup (e)
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:8, function(h) return(matrix(0, num_iter, 5)))

for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}
Summ = matrix(0, 8, 5)
for (ii in 1:8){
  Summ[ii, ] = colMeans(RESf[[ii]])
}
wz = c(c(1, 5), c(2, 6), c(3, 7), c(4, 8))
Table6_e = Summ[wz, ]


## Setup (f) ##

## Initialization 
Num_noise = 20
p1 = p/5
q = 6 
MM = c(250, 500) 
K = 50 
RES = list()

# Begin repetitions
for (ii in 1:Num_iter){
  
  print(ii) 
  R = c()
  
  ## Data generation
  res = logistic_comp_dg(n, p, p1, q)
  Y = res$Y
  X = res$X
  S = res$B_index
  Beta = res$Beta[S]
  
  ## Add noisy observations
  noise_value = runif(Num_noise*p)*20
  X_noise = matrix(noise_value, Num_noise, p)
  Y_noise = rep(0, Num_noise)
  X = rbind(X, X_noise)
  Y = c(Y, Y_noise)
  
  ## Distributed feature screening with m = 250
  Res = Dis_Scr_K(X, Y, S, K, m = MM[1])
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  ## Distributed feature screening with m = 500
  Res = Dis_Scr_K(X, Y, S, K, m = MM[2])
  R = rbind(R, Res$RR_FAIR)
  R = rbind(R, Res$RR_KF)
  R = rbind(R, Res$RR_MV)
  R = rbind(R, Res$RR_CRU)
  
  RES = c(RES, list(R))  
}

## Summarize the results of Setup (f)
rr1 = RES
num_iter = length(rr1)
RESf = lapply(1:8, function(h) return(matrix(0, num_iter, 5)))
for (ii in 1:num_iter){
  for (jj in 1:8){
    RESf[[jj]][ii, ] = rr1[[ii]][jj, ]
  } 
}
Summ = matrix(0, 8, 5)
for (ii in 1:8){
  Summ[ii, ] = colMeans(RESf[[ii]])
}
wz = c(c(1, 5), c(2, 6), c(3, 7), c(4, 8))
Table6_f = Summ[wz, ]


## Generate Table 6  ##

## Merge results from Setups (d) (e) (f)
RES6 = rbind(Table6_d, Table6_e, Table6_f)

## Output Table 6 in in subfolder "tables"
RES6[, 1:4] = round(RES6[, 1:4], 2)
RES6[, 5] = round(RES6[, 5], 4)
colnames(RES6) = c("SSR", "PSR", "FDR", "wRank", "Time")
Setup = c("(d)", rep(" ", 7), "(e)", rep(" ",7), "(f)", rep(" ", 7))
Utility = c("FAIR", " ", "KF", " ", "MV", " ", "CRU", " ")
Utility = rep(Utility, 3)
MM1 = c(50, 100)
MM2 = c(100, 200)
MM3 = c(250, 500)
m = c(rep(MM1, 4), rep(MM2, 4), rep(MM3, 4))
Table6 = cbind(Setup, Utility, m, RES6)
Table6 = as.data.frame(Table6)
write_xlsx(Table6, here("tables", "table6_full.xlsx"))


######################################
## Code for Table 7 in Section 4.2  ## 
## Runtime: about 40 hours          ##
######################################

print("---Table 7---")

## The dataset contains 649 features of 2000 handwritten numerals ("0"--"9")
## extracted from a collection of Dutch utility maps. 
## For each of the 2000 images, six groups of features (649 in total) are collected based on different
## processing methods. The measurements of those features are respectively stored in six separeate feature
## sets (files): "fac.txt", "fou.txt", "kar.txt", "mor.txt", "pix.txt", and "zer.txt".
## The dataset is available for download at https://archive.ics.uci.edu/ml/datasets/Multiple+Features 
## See readme.txt for more details about this dataset.

## Load data files
X_fac = as.matrix(read.table(here("data", "fac.txt"), header = T))
X_fou = as.matrix(read.table(here("data", "fou.txt"), header = T))
X_kar = as.matrix(read.table(here("data", "kar.txt"), header = T))
X_mor = as.matrix(read.table(here("data", "mor.txt"), header = T))
X_pix = as.matrix(read.table(here("data", "pix.txt"), header = T))
X_zer = as.matrix(read.table(here("data", "zer.txt"), header = T))

## Generate feature matrix and response 
X_raw = cbind(X_fac, X_fou, X_kar, X_mor, X_pix, X_zer)
Y = c()
for (k in 0:9){
  Y = c(Y, rep(k, 200))
}

## Set seed and number of repetitions
## Reduce Num_iter for faster runtimes if needed
set.seed(1)
Num_iter = 100

## Set data segment size and number of features to be retained
nl = 10 
K = 20

## Set number of noisy features and fake observations to be added
p1 = 1351 
noise_number_MM = c(0, 1000, 2000, 3000) 

## Matrices to store final output
RE_benchmark = c()
RE_psis = c()
RE_fkf = c()
RE_mv = c()
RE_cru = c()  

## Begin repetitions
for (iter in 1:Num_iter){
  
  print(iter)
  
  ## Add noisy features
  Xar = matrix(rt(p1*2000, 2), 2000, p1)
  X_new = cbind(X_raw, Xar)
  X_new = scale(X_new)
  p = dim(X_new)[2]
  
  ## Generate training and testing sets
  test_index = c()
  for (k in 0:9){
    test_index = c(test_index, k*200+sort(sample(200, 20)))
  } 
  train_index = setdiff(1:2000, test_index) 
  X_train = X_new[train_index, ] 
  X_test = X_new[test_index, ]
  Y_train = Y[train_index]
  Y_test = Y[test_index]
  X_train_raw = X_train
  Y_train_raw = Y_train
  
  ## Vectors to store prediction errors
  RR_benchmark = c()
  RR_psis = c()
  RR_fkf = c()
  RR_mv = c()
  RR_cru = c()  
  
  ## Analysis for different numbers of fake observations
  for (noise_number in noise_number_MM){
    
    ## Case with no fake observations
    if (noise_number==0){ 
      print(paste("N1=", noise_number, sep = ""))
      n_train = length(Y_train)
      X_train = X_train+0.0001*matrix(runif(n_train*p), n_train, p) 
    }
    
    ## Cases with fake observations
    if (noise_number>0){
      print(paste("N1=", noise_number, sep=""))
      Xout = matrix(runif(noise_number*p), noise_number, p)*1000
      Yout = sample(0:9, noise_number, replace = TRUE)
      X_train = rbind(X_train, Xout)
      Y_train = c(Y_train, Yout)
      n_train = length(Y_train) 
      X_train = X_train+0.0001*matrix(runif(n_train*p), n_train, p) 
    }
    
    ## Create data segments
    m = n_train/nl
    Sn = rep(floor(n_train/m), m)
    Nr = n_train-sum(Sn)
    if (Nr>0){Sn[1:Nr] = Sn[1:Nr]+1}
    M_num = fenzu(Sn)       
    
    ## KNN prediction errors with different screening methods
    res_benchmark = Benchmark_KNN_pred(X_train, X_test, Y_train, Y_test, K, kk=200)
    RR_benchmark = c(RR_benchmark, res_benchmark)    
    res_psis = Scr_KNN_pred(X_train, X_test, Y_train, Y_test, M_num, K, kk=200, method =  "PSIS")
    RR_psis = c(RR_psis, res_psis)    
    res_fkf = Scr_KNN_pred(X_train, X_test, Y_train, Y_test, M_num, K, kk=200, method = "FKF")
    RR_fkf = c(RR_fkf, res_fkf) 
    res_mv = Scr_KNN_pred(X_train, X_test, Y_train, Y_test, M_num, K, kk=200, method =  "MV")
    RR_mv = c(RR_mv, res_mv)
    res_cru = Scr_KNN_pred(X_train, X_test, Y_train, Y_test, M_num, K, kk=200, method = "CRU")
    RR_cru = c(RR_cru, res_cru)
    
    ## Restore training set
    X_train = X_train_raw
    Y_train = Y_train_raw 
  } 
  
  ## Store results for the current repetition
  RE_benchmark = rbind(RE_benchmark, RR_benchmark)
  RE_psis = rbind(RE_psis, RR_psis)
  RE_fkf = rbind(RE_fkf, RR_fkf)
  RE_mv = rbind(RE_mv, RR_mv)
  RE_cru = rbind(RE_cru, RR_cru)
  a1 = list(RE_benchmark, RE_psis, RE_fkf, RE_mv, RE_cru)
}

## Output Table 7 in subfolder "tables"
Table7 = c()
for (ii in 1:5){
  res = colMeans(a1[[ii]])
  Table7 = rbind(Table7, res)
}
colnames(Table7) = c("N_1=0", "N_1=1000", "N_1=2000", "N_1=3000")
Method = c("RS", "PSIS", "FKF", "MV", "CRU")
Table7 = cbind(Method, Table7)
Table7 = as.data.frame(Table7)
write_xlsx(Table7, here("tables", "table7_full.xlsx"))

