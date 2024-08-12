#-------------------------------------------------------------------------------
# Simulation1 of FANMI
#-------------------------------------------------------------------------------
rm(list=ls())

######################################
## Load Packages: 
######################################
library(numDeriv)
library(alabama)
library(splines)
library(Matrix)
library(fda)
library(cluster)
library(plyr)
library(splines2)
library(nlme)
library(mgcv)
library(MASS)
library(conquer)
library(tcltk)
library(ggplot2)
library(plot3D)
library(plotly)
library(akima)
library(phtt)
library(progress)
library(xtable)
library(NbClust)
library(mclust)
library(mgcv)
######################################
## End Load Packages: 
######################################

######################################
## Setup: 
######################################
N <- 300
degree <- 3
basisnumber1 <- 6
basisnumber2 <- 5
n<-100
iteration <- 500 # number of iteration 
fhat1MC<-array(NA,dim=c(iteration,N,N))
fhat2MC<-array(NA,dim=c(iteration,N,N))
fhat3MC<-array(NA,dim=c(iteration,N,N))
fhat1nMC<-array(NA,dim=c(iteration,n,n))
fhat2nMC<-array(NA,dim=c(iteration,n,n))
fhat3nMC<-array(NA,dim=c(iteration,n,n))
zeta1MC<-array(NA,dim=c(iteration,N))
zeta2MC<-array(NA,dim=c(iteration,N))
zeta3MC<-array(NA,dim=c(iteration,N))
CI1_up<-array(NA,dim=c(iteration,n,n))
CI1_low<-array(NA,dim=c(iteration,n,n))
CI2_up<-array(NA,dim=c(iteration,n,n))
CI2_low<-array(NA,dim=c(iteration,n,n))
CI3_up<-array(NA,dim=c(iteration,n,n))
CI3_low<-array(NA,dim=c(iteration,n,n))
zMC<-array(NA,dim=c(iteration,N))
######################################
## End Setup: 
######################################

#######################################
## Data Generating Propcess:
#######################################
M = 20 # that many summands for X 
u = seq(1,100,1)/100 # xgrid
K = 3 # number of additive functions
n = 100
t1 <- seq(0.001,0.99,length.out = n)
tz <- seq(0.001,0.99,length.out = n)


#f1 <- function(s,t){5*sin(s)*cos(t)}
#f2 <- function(s,t){5*sin(t)*cos(s)}
#f3 <- function(s,t){5*t*cos(s)}
#f1 <- function(s,t){sin(s)*(1+0.5*cos(t))}
#f2 <- function(s,t){(1+0.5*sin(t))*s}
#f3 <- function(s,t){cos(s)*(1+t)}
f1 <- function(s,t){2*cos(2*pi*s)*t^2}
f2 <- function(s,t){(3*(s-0.5)^2+sin(2.5*pi*s))*t}
f3 <- function(s,t){5*t*cos(s)}
f1 <- function(s,t){0.2*cos(2*pi*s)*t^2}
f2 <- function(s,t){0.1*(3*(s-0.5)^2+sin(2.5*pi*s))*t}
f3 <- function(s,t){0.5*t*cos(pi*s)}
#f1 <- function(s,t){sin(2*pi*s)*(1+cos(2*pi*t))}
#f2 <- function(s,t){(1+c[k]*sin(2.5*pi*t))*s}
#f2 <- function(s,t){cos(2*pi*s)*(1+2*t)}

#######################################
## End Data Generating Propcess:
#######################################


#######################################
######## Monte Carlo 
#######################################
pb <- tkProgressBar(title="进度",label="已完成 %", min=0, max=100, initial = 0, width = 300) 
for (it in 1:iteration){
  info<- sprintf("已完成 %d%%", round(it*100/length(1:iteration)))  
  setTkProgressBar(pb, value = it*100/length(1:iteration), title = sprintf("进度 (%s)",info),label = info) 
  ## -------------Generate Xi(t)--------------------------
  ## Generate Xi(t)
  xii = matrix(rep(0,N*length(u)),nrow=(N),ncol=length(u))
  xij = matrix(rep(0,N*M),nrow=N,ncol=M)
  for (i in 1:N){
    xi = 0
    for (j in 1:M){ 
      xij[i,j] = rnorm(1,0,1/((j-0.5)*pi))
      phij = sqrt(2)*sin(pi*u*(j-0.5))
      xi = xi+xij[i,j]*phij }
    xii[i,] = xi
  }
  ## -------------Scale the FPC scores--------------------------
  ## Scale the FPC scores
  zeta = matrix(rep(0,N*K),nrow=N,ncol=K)
  zeta[,1] = pnorm(xij[,1]/(1/((1-0.5)*pi)))
  zeta[,2] = pnorm(xij[,2]/(1/((2-0.5)*pi)))
  zeta[,3] = pnorm(xij[,3]/(1/((3-0.5)*pi)))
  ## -------------Generate Z--------------------------
  ## Generate Z
  z = matrix(runif(N),N,1)
  #z = sort(z)
  ## -------------Generate g(u_i)=xbeta--------------------------
  ## Generate g(u_i)=xbeta
  mui = matrix(rep(0,N),nrow=N,ncol=1)
  yi = matrix(rep(0,N),nrow=N,ncol=1)
  mui = f1(zeta[,1],z)+f2(zeta[,2],z)+f3(zeta[,3],z)
  meanf1 = matrix(NA,nrow=N,ncol=1)
  meanf2 = matrix(NA,nrow=N,ncol=1)
  meanf3 = matrix(NA,nrow=N,ncol=1)
  for(i in 1:N){
    meanf1[i,1] = mean(f1(zeta[,1],z[i]))
    meanf2[i,1] = mean(f2(zeta[,2],z[i]))
    meanf3[i,1] = mean(f3(zeta[,3],z[i]))
  }
  mui_c = (meanf1 + meanf2 + meanf3) + f1(zeta[,1],z)-meanf1 + f2(zeta[,2],z)-meanf2 + f3(zeta[,3],z)-meanf3

  ## ---------------------------------------
  # Normal 
  yi = mui_c + rnorm(1,0,0.1)
  ## ---------------------------------------
  # Poisson
  #p = exp(mui_c)
  #yi = rpois(N,lambda=p)
  ## ---------------------------------------
  # Binomial
  #mui_c = mui_c 
  #p = exp(mui_c)/(1+exp(mui_c))
  #yi = rbinom(N,prob = p,size=1)
  ## ---------------------------------------
  ## Estimate pca scores
  X = matrix(rep(0,N*length(u)),nrow=N,ncol=length(u))
  Sig = array(NA,dim=c(length(u),length(u)))
  m_hat = rep(NA)
  #X[,]<- xii[,]-matrix(colSums(xii[,])/nrow(xii[,]),N,length(u),byrow=T) 
  X[,] = xii[,]
  Sig[,] = t(X[,])%*%X[,]/(N*length(u))
  #suppressWarnings(m_hat <-c(as.list(OptDim(t(X[1:N,]),criteria = "ER"))[[1]]$`Optimal Dimension`))
  
  m_hat=3
  Phi_hat = array(NA,dim=c(length(u),m_hat))
  Lmbd = matrix(NA,m_hat,1)
  xi_hat = matrix(NA,N,m_hat)
  zeta_hat = matrix(NA,nrow=N,ncol=m_hat)
  
  Phi_hat[,] = eigen(Sig[,])$vectors[,1:m_hat]*(length(u)^(1/2))
  Lmbd[,1] = eigen(Sig[,])$values[1:m_hat]
  xi_hat = X%*%Phi_hat/length(u)
  sign = array(NA,m_hat)
  for(m in 1:m_hat){
    sign[m] = sign(xi_hat[1,m])*sign(xij[1,m])
    zeta_hat[,m] = pnorm(sign[m]*(xi_hat[,m])/sqrt(Lmbd[m,1]))
  }
  #######################################
  ## End Data Generating Propcess:
  #######################################
  
  #######################################
  ## Estimation 
  #######################################  
  knot1 = as.vector(quantile(zeta_hat[,1],seq(0,1,length=4)))
  knot2 = as.vector(quantile(zeta_hat[,2],seq(0,1,length=4)))
  knot3 = as.vector(quantile(zeta_hat[,3],seq(0,1,length=4)))
  
  basis_1 = create.bspline.basis(breaks = knot1,nbasis = 5,norder = 3)
  basis_2 = create.bspline.basis(breaks = knot2,nbasis = 5,norder = 3)
  basis_3 = create.bspline.basis(breaks = knot3,nbasis = 5,norder = 3)
  basis_z = create.bspline.basis(rangeval = range(z),nbasis = 4,norder = 2)
  tbval1 = eval.basis(zeta_hat[,1],basis_1)
  tbval2 = eval.basis(zeta_hat[,2],basis_2)
  tbval3 = eval.basis(zeta_hat[,3],basis_3)
  tbval_z = eval.basis(as.vector(z),basis_z)

  ## ---------------------------------------
  nbasis = dim(tbval1)[2]
  tbval1 = tbval1[,-nbasis]
  tbval2 = tbval2[,-nbasis]
  tbval3 = tbval3[,-nbasis]
  tbval_z = tbval_z[,-nbasis]
  
  basis0 = (nbasis - 1)
  basis1 = (nbasis - 1)^2
  basis2 = (nbasis - 1)^2
  basis3 = (nbasis - 1)^2
  b0 <- matrix(rep(0,N*basis0),nrow=N,ncol=basis0)
  b1 <- matrix(rep(0,N*basis1),nrow=N,ncol=basis1)
  b2 <- matrix(rep(0,N*basis2),nrow=N,ncol=basis2)
  b3 <- matrix(rep(0,N*basis3),nrow=N,ncol=basis3)
  b0 <- tbval_z
  for(i in 1:N){
      b1[i,] <- kronecker(tbval1[i,],tbval_z[i,])
      b2[i,] <- kronecker(tbval2[i,],tbval_z[i,])
      b3[i,] <- kronecker(tbval3[i,],tbval_z[i,])
  }
  
  ## ---------------------------------------
  # Normal
  f<-glm(yi[1:N]~0+b0+b1+b2+b3)
  ## ---------------------------------------
  # Poisson
  #f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "poisson")
  ## ---------------------------------------
  # Binomial
  #f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "binomial")
  ## ---------------------------------------

  fhat0=matrix(NA,nrow=N,ncol=1)
  fhat1=matrix(NA,nrow=N,ncol=N)
  fhat2=matrix(NA,nrow=N,ncol=N)
  fhat3=matrix(NA,nrow=N,ncol=N)
  fhat0 = tbval_z %*% f$coefficients[1:basis0]
  fhat1 = b1 %*% f$coefficients[(basis0+1):(basis0+basis1)]
  fhat2 = b2%*%f$coefficients[(basis0+basis1+1):(basis0+basis1+basis2)]
  fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+1):(basis0+basis1+basis2+basis3)]
  
  ## Step2 Estimation
  tbval1n = eval.basis(zeta_hat[,1],basis_1)
  tbval2n = eval.basis(zeta_hat[,2],basis_2)
  tbval3n = eval.basis(zeta_hat[,3],basis_3)
  tbval_z = eval.basis(as.vector(z),basis_z)
  basis1n = (nbasis)*(nbasis-1)
  basis2n = (nbasis)*(nbasis-1)
  basis3n = (nbasis)*(nbasis-1)
  b1n <- matrix(rep(0,N*basis1n),nrow=N,ncol=basis1n)
  b2n <- matrix(rep(0,N*basis2n),nrow=N,ncol=basis2n)
  b3n <- matrix(rep(0,N*basis3n),nrow=N,ncol=basis3n)
  fhat1n=matrix(NA,nrow=N,ncol=N)
  fhat2n=matrix(NA,nrow=N,ncol=N)
  fhat3n=matrix(NA,nrow=N,ncol=N)
  for(i in 1:N){
    b1n[i,] <- kronecker(tbval1n[i,],tbval_z[i,])
    b2n[i,] <- kronecker(tbval2n[i,],tbval_z[i,])
    b3n[i,] <- kronecker(tbval3n[i,],tbval_z[i,])
  }
  
  knot1h = as.vector(quantile(t1,seq(0,1,length=4)))
  knot2h = as.vector(quantile(t1,seq(0,1,length=4)))
  knot3h = as.vector(quantile(t1,seq(0,1,length=4)))
  basis_1h = create.bspline.basis(breaks = knot1h,nbasis = 5,norder = 3)
  basis_2h = create.bspline.basis(breaks = knot2h,nbasis = 5,norder = 3)
  basis_3h = create.bspline.basis(breaks = knot3h,nbasis = 5,norder = 3)
  basis_zh = create.bspline.basis(rangeval = range(tz),nbasis = 4,norder = 2)
  tbval1h = eval.basis(t1,basis_1h)
  tbval2h = eval.basis(t1,basis_2h)
  tbval3h = eval.basis(t1,basis_3h)
  tbval_zh= eval.basis(tz,basis_zh)
  basis1h = (nbasis)*(nbasis-1)
  basis2h = (nbasis)*(nbasis-1)
  basis3h = (nbasis)*(nbasis-1)
  b1h <- matrix(rep(0,n*basis1h),nrow=n,ncol=basis1h)
  b2h <- matrix(rep(0,n*basis2h),nrow=n,ncol=basis2h)
  b3h <- matrix(rep(0,n*basis3h),nrow=n,ncol=basis3h)
  fhat1h=matrix(NA,nrow=n,ncol=n)
  fhat2h=matrix(NA,nrow=n,ncol=n)
  fhat3h=matrix(NA,nrow=n,ncol=n)
  for(i in 1:n){
    b1h[i,] <- kronecker(tbval1h[i,],tbval_zh[i,])
    b2h[i,] <- kronecker(tbval2h[i,],tbval_zh[i,])
    b3h[i,] <- kronecker(tbval3h[i,],tbval_zh[i,])
  }
  
    
  # Normal
  yi1 = yi[1:N] - fhat0 - fhat2 - fhat3
  g1<-glm(yi1~0+b1n)
  ## ---------------------------------------
  # Poisson
  #mui1 = mui_c- fhat0c - diag(fhat2c) - diag(fhat3c)
  #p1 = exp(mui1)
  #yi1 = rpois(N,lambda=p1)
  
  #mui1 = mui_c- fhat0 - fhat2 - fhat3
  #g1<-glm(mui1~0+b1n)
  ## ---------------------------------------
  # Binomial
  #mui1 = mui_c- fhat0c - diag(fhat2c) - diag(fhat3c)
  #mui1 = mui_c- fhat0 - fhat2 - fhat3
  #g1<-glm(mui1~0+b1n)
  
  #p1 = exp(mui1)/(1+exp(mui1))
  #yi1 = rbinom(N,prob = p1,size=1)
  #g1<-glm(yi1~0+b1n,family = "binomial")
  ## ---------------------------------------
  cv_q = sqrt(2*log(16)) - (1/sqrt(2*log(16)))*(log(-0.5*log(0.95))+0.5*(log(log(16))+log(4*pi)))
  
  #for (i in 1:N) {for(j in 1:N){fhat1n[i,j] = kronecker(tbval1n[i,],tbval_z[j,])%*%g1$coefficients[1:(basis1n)]}}
  fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
  for (i in 1:n) {for(j in 1:n){
    fhat1h[i,j] = kronecker(tbval1h[i,],tbval_zh[j,])%*%g1$coefficients[1:(basis1h)]
  }}
  fhat1hc = fhat1h - t(matrix(rep(colSums(fhat1h)/n,n,each=1),n,n))
  for (i in 1:n) {for(j in 1:n){
    sgm1 <- sqrt(kronecker(tbval1h[i,],tbval_zh[j,])%*%vcov(g1)%*%kronecker(tbval1h[i,],tbval_zh[j,]))
    CI1_up[it,i,j] <-fhat1hc[i,j] + cv_q * sgm1
    CI1_low[it,i,j] <- fhat1hc[i,j] - cv_q * sgm1
  }}
  
    
  
  # Normal
  yi2 = yi[1:N] - fhat0 - fhat1n - fhat3
  g2<-glm(yi2~0+b2n)
  ## ---------------------------------------
  # Poisson
  #mui2 = mui_c - fhat0c - diag(fhat1n) - diag(fhat3c)
  #p2 = exp(mui2)
  #yi2 = rpois(N,lambda=p2)
  #g2<-glm(yi2~0+b2n,family = "poisson")
  
  #mui2 = mui_c - (fhat0) - fhat1n - fhat3
  #g2<-glm(mui2~0+b2n)
  ## ---------------------------------------
  # Binomial
  #mui2 = mui_c - (fhat0c) - diag(fhat1n) - diag(fhat3c)
  #mui2 = mui_c - (fhat0) - fhat1n - fhat3
  #g2<-glm(mui2~0+b2n)
  #p2 = exp(mui2)/(1+exp(mui2))
  #yi2 = rbinom(N,prob = p2,size=1)
  #g2<-glm(yi2~0+b2n,family = "binomial")
  ## ---------------------------------------
  fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
  for (i in 1:n) {for(j in 1:n){fhat2h[i,j] = kronecker(tbval2h[i,],tbval_zh[j,])%*%g2$coefficients[1:(basis2h)]}}
  fhat2hc = fhat2h - t(matrix(rep(colSums(fhat2h)/n,n,each=1),n,n))
  for (i in 1:n) {for(j in 1:n){
    sgm2 <- sqrt(kronecker(tbval2h[i,],tbval_zh[j,])%*%vcov(g2)%*%kronecker(tbval2h[i,],tbval_zh[j,]))
    CI2_up[it,i,j] <-fhat2hc[i,j] + cv_q * sgm2
    CI2_low[it,i,j] <- fhat2hc[i,j] - cv_q * sgm2
  }}
  #for (i in 1:N) {for(j in 1:N){fhat2n[i,j] = kronecker(tbval2n[i,],tbval_z[j,])%*%g2$coefficients[1:(basis2n)]}}
  #fhat2nc = fhat2n - t(matrix(rep(colSums(fhat2n)/N,N,each=1),N,N))
  
  # Normal
  yi3 = yi[1:N] - fhat0 - fhat1n - fhat2n
  g3<-glm(yi3~0+b3n)
  ## ---------------------------------------
  # Poisson
  #mui3 = mui_c- fhat0c - diag(fhat1n) - diag(fhat2n)
  #p3 = exp(mui3)
  #yi3 = rpois(N,lambda=p3)
  #g3<-glm(yi3~0+b3n,family = "poisson")
  #g3<-glm(mui3~0+b3n)
  ## ---------------------------------------
  # Binomial
  
  #mui3 = mui_c- fhat0 - fhat1n - fhat2n
  #g3 <- glm(mui3~0+b3n)
  
  #p3 = exp(mui3)/(1+exp(mui3))
  #yi3 = rbinom(N,prob = p3,size=1)
  #g3 <- glm(yi3~0+b3n,family = "binomial")
  ## ---------------------------------------
 
  fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
  for (i in 1:n) {for(j in 1:n){fhat3h[i,j] = kronecker(tbval3h[i,],tbval_zh[j,])%*%g3$coefficients[1:(basis3h)]}}
  fhat3hc = fhat3h - t(matrix(rep(colSums(fhat3h)/n,n,each=1),n,n))
  for (i in 1:n) {for(j in 1:n){
    sgm3 <- sqrt(kronecker(tbval3h[i,],tbval_zh[j,])%*%vcov(g3)%*%kronecker(tbval3h[i,],tbval_zh[j,]))
    CI3_up[it,i,j] <-fhat3hc[i,j] + cv_q * sgm3
    CI3_low[it,i,j] <- fhat3hc[i,j] - cv_q * sgm3
  }}
  #for (i in 1:N) {for(j in 1:N){fhat3n[i,j] = kronecker(tbval3n[i,],tbval_z[j,])%*%g3$coefficients[1:(basis3n)]}}
  #fhat3nc = fhat3n - t(matrix(rep(colSums(fhat3n)/N,N,each=1),N,N))
  ## ---------------------------------------
  
  #fhat1MC[it,,]<-fhat1c[order(zeta_hat[,1]),order(z[1:N])]
  #fhat2MC[it,,]<-fhat2c[order(zeta_hat[,2]),order(z[1:N])]
  #fhat3MC[it,,]<-fhat3c[order(zeta_hat[,3]),order(z[1:N])]
  #fhat1nMC[it,,]<-fhat1hc[order(zeta_hat[,1]),order(z[1:N])]
  #fhat2nMC[it,,]<-fhat2hc[order(zeta_hat[,2]),order(z[1:N])]
  #fhat3nMC[it,,]<-fhat3hc[order(zeta_hat[,3]),order(z[1:N])]
  fhat1nMC[it,,]<-fhat1hc[order(t1),order(tz)]
  fhat2nMC[it,,]<-fhat2hc[order(t1),order(tz)]
  fhat3nMC[it,,]<-fhat3hc[order(t1),order(tz)]
  
  zeta1MC[it,]<-sort(zeta_hat[,1])
  zeta2MC[it,]<-sort(zeta_hat[,2])
  zeta3MC[it,]<-sort(zeta_hat[,3])
  zMC[it,]<-sort(z)
}
close(pb)

fhat1Sum = fhat1MC[1,,]
fhat2Sum = fhat2MC[1,,]
fhat3Sum = fhat3MC[1,,]
fhat1nSum = fhat1nMC[1,,]
fhat2nSum = fhat2nMC[1,,]
fhat3nSum = fhat3nMC[1,,]
CI1lowSum = CI1_low[1,,]
CI2lowSum = CI2_low[1,,]
CI3lowSum = CI3_low[1,,]
CI1upSum = CI1_up[1,,]
CI2upSum = CI2_up[1,,]
CI3upSum = CI3_up[1,,]
for(i in 2:iteration){
  fhat1nSum = fhat1nSum + fhat1nMC[i,,]
  fhat2nSum = fhat2nSum + fhat2nMC[i,,]
  fhat3nSum = fhat3nSum + fhat3nMC[i,,]
  CI1lowSum = CI1lowSum + CI1_low[i,,]
  CI2lowSum = CI2lowSum + CI2_low[i,,]
  CI3lowSum = CI3lowSum + CI3_low[i,,]
  CI1upSum = CI1lowSum + CI1_up[i,,]
  CI2upSum = CI2lowSum + CI2_up[i,,]
  CI3upSum = CI3lowSum + CI3_up[i,,]
}
fhat1mean = fhat1Sum/iteration
fhat2mean = fhat2Sum/iteration
fhat3mean = fhat3Sum/iteration
fhat1nmean = fhat1nSum/iteration
fhat2nmean = fhat2nSum/iteration
fhat3nmean = fhat3nSum/iteration
CI1lowmean = CI1lowSum/iteration
CI2lowmean = CI2lowSum/iteration
CI3lowmean = CI3lowSum/iteration
CI1upmean = CI1upSum/iteration
CI2upmean = CI2upSum/iteration
CI3upmean = CI3upSum/iteration
zeta1mean = colSums(zeta1MC)/iteration
zeta2mean = colSums(zeta2MC)/iteration
zeta3mean = colSums(zeta3MC)/iteration
zmean = colSums(zMC)/iteration

func1 = matrix(NA,nrow = n,ncol=n)
func2 = matrix(NA,nrow = n,ncol=n)
func3 = matrix(NA,nrow = n,ncol=n)

for(i in 1:n){
  for(j in 1:n){
    func1[i,j] = f1(t1[i],tz[j])
    func2[i,j] = f2(t1[i],tz[j])
    func3[i,j] = f3(t1[i],tz[j])
  }
}
func1 = func1 - t(matrix(rep(colSums(func1)/n,n,each=1),n,n))
func2 = func2 - t(matrix(rep(colSums(func2)/n,n,each=1),n,n))
func3 = func3 - t(matrix(rep(colSums(func3)/n,n,each=1),n,n))

persp3D(sort(t1),sort(tz),fhat1hc[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),CI1_low[1,order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),CI1_up[1,order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),fhat2hc[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),fhat3hc[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),func1[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),func2[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")
persp3D(sort(t1),sort(tz),func3[order(t1),order(tz)],theta=45,phi=40,xlab="",ylab="",zlab="",main="")

####################################### 
## End Estimation 
####################################### 

