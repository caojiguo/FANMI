# Simulation3 of FANMI: Interaction test
#-------------------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
######################################
## Load Packages: 
######################################
library(MASS)
library(fda)
library(alabama)
library(pracma)
library(rmutil)
library(cubature)
library(parallel)
library(foreach)
library(doParallel)

######################################
## Setup: 
######################################

c = c(0,0.5,1,1.25,1.5)
N = 500
M = 20 # that many summands for X 
u = seq(1,100,1)/100 # xgrid
K = 3
n = 100
t1 <-  seq(0.001,0.99,length.out = 100)
tz <- t1
type0 = c("Normal","Poisson","Binomial")

# power
level = 0.05
pfun = array(0, dim = c(5,500,3))
powerfun = array(0, dim = c(5,500,3))
T = array(0, dim = c(5,500,3))

## -------------Generate Xi(t)--------------------------
## Generate Xi(t)
Generate_X <- function(N,M,u,K){
  xii = matrix(rep(0,N*length(u)),nrow=(N),ncol=length(u))
  xij = matrix(rep(0,N*M),nrow=N,ncol=M)
  for (i in 1:N){
    xi = 0
    for (j in 1:M){ 
      xij[i,j] = rnorm(1,0,1/((j-0.5)*pi))
      phij = sqrt(2)*sin(pi*u*(j-0.5))
      xi = xi+xij[i,j]*phij 
    }
    xii[i,] = xi
  }
  zeta = matrix(rep(0,N*K),nrow=N,ncol=K)
  zeta[,1] = pnorm(xij[,1]/(1/((1-0.5)*pi)))
  zeta[,2] = pnorm(xij[,2]/(1/((2-0.5)*pi)))
  zeta[,3] = pnorm(xij[,3]/(1/((3-0.5)*pi)))
  ## Estimate pca scores
  X = matrix(rep(0,N*length(u)),nrow=N,ncol=length(u))
  Sig = array(NA,dim=c(length(u),length(u)))
  m_hat = rep(NA)
  X[,] = xii[,]
  Sig[,] = t(X[,])%*%X[,]/(N*length(u))
  m_hat = 3
  
  ## Scale the FPC scores
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
  return(list('zeta_hat'=zeta_hat,'zeta'=zeta))
}

## -------------Generate Z-------------------------------
## Generate Z
Generate_z <- function(N){
  z = matrix(runif(N),N,1)
  return(z)
}

## -------------Generate g(u_i)=xbeta--------------------
## Generate g(u_i)=xbeta
Generate_g <- function(f1,f2,f3,zeta,z,N,type="Normal"){
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
  if(type=="Normal"){
    yi = mui_c + rnorm(1,0,0.1)
  }
  if(type=="Poisson"){
    p = exp(mui_c)
    yi = rpois(N,lambda=p)
  }
  if(type=="Binomial"){
    p = exp(mui_c)/(1+exp(mui_c))
    yi = rbinom(N,prob = p,size=1)
  }
  return(list('mui_c'=mui_c,'yi'=yi))
}

## -------------H1--------------------
GFVM = function(zeta_hat,z,yi,mui_c,N,type="Normal",t1= seq(0.001,0.99,length.out = 100),tz=t1 <- seq(0.001,0.99,length.out = 100),
                weights=NULL){
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
  
  
  if(type=="Normal"){
    f<-glm(yi[1:N]~0+b0+b1+b2+b3)
    #r = f$residuals
    fhat0 = tbval_z %*% f$coefficients[1:basis0]
    fhat1 = b1%*%f$coefficients[(basis0+1):(basis0+basis1)]
    fhat2 = b2%*%f$coefficients[(basis0+basis1+1):(basis0+basis1+basis2)]
    fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+1):(basis0+basis1+basis2+basis3)]
    
    yi1 = yi[1:N] - fhat0 - fhat2 - fhat3
    g1<-glm(yi1~0+b1n)
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    for (i in 1:n) {for(j in 1:n){fhat1h[i,j] = kronecker(tbval1h[i,],tbval_zh[j,])%*%g1$coefficients[1:(basis1h)]}}
    fhat1hc = fhat1h - t(matrix(rep(colSums(fhat1h)/n,n,each=1),n,n))
    
    yi2 = yi[1:N] - fhat0 - fhat1n - fhat3
    g2<-glm(yi2~0+b2n)
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    for (i in 1:n) {for(j in 1:n){fhat2h[i,j] = kronecker(tbval2h[i,],tbval_zh[j,])%*%g2$coefficients[1:(basis2h)]}}
    fhat2hc = fhat2h - t(matrix(rep(colSums(fhat2h)/n,n,each=1),n,n))
    
    yi3 = yi[1:N] - fhat0 - fhat1n - fhat2n
    g3<-glm(yi3~0+b3n)
    r = g3$residuals
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    for (i in 1:n) {for(j in 1:n){fhat3h[i,j] = kronecker(tbval3h[i,],tbval_zh[j,])%*%g3$coefficients[1:(basis3h)]}}
    fhat3hc = fhat3h - t(matrix(rep(colSums(fhat3h)/n,n,each=1),n,n))
  }
  
  if(type=="Poisson"){
    f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "poisson",weights = weights)
    fhat0 = tbval_z %*% f$coefficients[1:basis0]
    fhat1 = b1 %*% f$coefficients[(basis0+1):(basis0+basis1)]
    fhat2 = b2 %*% f$coefficients[(basis0+basis1+1):(basis0+basis1+basis2)]
    fhat3 = b3 %*% f$coefficients[(basis0+basis1+basis2+1):(basis0+basis1+basis2+basis3)]
    
    mui1 = mui_c- fhat0 - fhat2 - fhat3
    g1<-glm(mui1~0+b1n)
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    for (i in 1:n) {for(j in 1:n){fhat1h[i,j] = kronecker(tbval1h[i,],tbval_zh[j,])%*%g1$coefficients[1:(basis1h)]}}
    fhat1hc = fhat1h - t(matrix(rep(colSums(fhat1h)/n,n,each=1),n,n))
    
    mui2 = mui_c - fhat0 - fhat1n - fhat3
    g2<-glm(mui2~0+b2n)
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    for (i in 1:n) {for(j in 1:n){fhat2h[i,j] = kronecker(tbval2h[i,],tbval_zh[j,])%*%g2$coefficients[1:(basis2h)]}}
    fhat2hc = fhat2h - t(matrix(rep(colSums(fhat2h)/n,n,each=1),n,n))
    
    mui3 = mui_c- fhat0 - fhat1n - fhat2n
    g3 <- glm(mui3~0+b3n)
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    r = g3$residuals
    for (i in 1:n) {for(j in 1:n){fhat3h[i,j] = kronecker(tbval3h[i,],tbval_zh[j,])%*%g3$coefficients[1:(basis3h)]}}
    fhat3hc = fhat3h - t(matrix(rep(colSums(fhat3h)/n,n,each=1),n,n))
  }
  
  if(type=="Binomial"){
    f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "binomial",weights = weights)
    fhat0 = tbval_z %*% f$coefficients[1:basis0]
    fhat1 = b1%*%f$coefficients[(basis0+1):(basis0+basis1)]
    fhat2 = b2%*%f$coefficients[(basis0+basis1+1):(basis0+basis1+basis2)]
    fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+1):(basis0+basis1+basis2+basis3)]
    
    mui1 = mui_c- fhat0 - fhat2 - fhat3
    g1<-glm(mui1~0+b1n)
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    for (i in 1:n) {for(j in 1:n){fhat1h[i,j] = kronecker(tbval1h[i,],tbval_zh[j,])%*%g1$coefficients[1:(basis1h)]}}
    fhat1hc = fhat1h - t(matrix(rep(colSums(fhat1h)/n,n,each=1),n,n))
    
    mui2 = mui_c - fhat0 - fhat1n - fhat3
    g2<-glm(mui2~0+b2n)
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    for (i in 1:n) {for(j in 1:n){fhat2h[i,j] = kronecker(tbval2h[i,],tbval_zh[j,])%*%g2$coefficients[1:(basis2h)]}}
    fhat2hc = fhat2h - t(matrix(rep(colSums(fhat2h)/n,n,each=1),n,n))
    
    mui3 = mui_c- fhat0 - fhat1n - fhat2n
    g3 <- glm(mui3~0+b3n)
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    r = g3$residuals
    for (i in 1:n) {for(j in 1:n){fhat3h[i,j] = kronecker(tbval3h[i,],tbval_zh[j,])%*%g3$coefficients[1:(basis3h)]}}
    fhat3hc = fhat3h - t(matrix(rep(colSums(fhat3h)/n,n,each=1),n,n))
  }
  
  return(list("f1" = fhat1hc[order(t1),order(tz)],
              "f2" = fhat2hc[order(t1),order(tz)],
              "f3" = fhat3hc[order(t1),order(tz)],
              "zeta" = cbind(sort(zeta_hat[,1]),sort(zeta_hat[,2]),sort(zeta_hat[,3])),
              "z" = sort(z),
              "theta1" = g1$coefficients[1:(basis1n)],
              "theta2" = g2$coefficients[1:(basis2n)],
              "theta3" = g3$coefficients[1:(basis3n)],
              "f1n" = fhat1n,
              "f2n" = fhat2n,
              "f3n" = fhat3n,
              "f0" = fhat0,
              "r" = r))
}

## -------------H0--------------------
H0_GFVM = function(zeta_hat,z,yi,mui_c,N,type="Normal",t1= seq(0.001,0.99,length.out = 100),tz = seq(0.001,0.99,length.out = 100),
                   weights = NULL){
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
  tbval_z = tbval_z[,-nbasis]
  
  basis0 = (nbasis - 1)
  basis1 = (nbasis - 1)
  basis2 = (nbasis - 1)
  basis3 = (nbasis - 1)
  b0 <- tbval_z
  b1 <- tbval1
  b2 <- tbval2
  b3 <- tbval3
  
  
  tbval1n = eval.basis(zeta_hat[,1],basis_1)
  tbval2n = eval.basis(zeta_hat[,2],basis_2)
  tbval3n = eval.basis(zeta_hat[,3],basis_3)
  tbval_z = eval.basis(as.vector(z),basis_z)
  basis1n = (nbasis)
  basis2n = (nbasis)
  basis3n = (nbasis)
  b1n <- tbval1n
  b2n <- tbval2n
  b3n <- tbval3n
  fhat1n=matrix(NA,nrow=N,ncol=1)
  fhat2n=matrix(NA,nrow=N,ncol=1)
  fhat3n=matrix(NA,nrow=N,ncol=1)

  
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
  basis1h = (nbasis)
  basis2h = (nbasis)
  basis3h = (nbasis)
  b1h <- tbval1h
  b2h <- tbval2h
  b3h <- tbval3h
  fhat1h=matrix(NA,nrow=n,ncol=1)
  fhat2h=matrix(NA,nrow=n,ncol=1)
  fhat3h=matrix(NA,nrow=n,ncol=1)
  
  
  if(type=="Normal"){
    f<-glm(yi[1:N]~0+b1+b2+b3,weights = weights)
    fhat1 = b1%*%f$coefficients[1:basis1]
    fhat2 = b2%*%f$coefficients[(basis1+1):(basis1+basis2)]
    fhat3 = b3%*%f$coefficients[(basis1+basis2+1):(basis1+basis2+basis3)]
    
    yi1 = yi[1:N] - fhat2 - fhat3
    g1<-glm(yi1~0+b1n)
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    r1 = g1$residuals
    for (i in 1:n) {fhat1h[i,1] = tbval1h[i,]%*%g1$coefficients[1:(basis1h)]}
    fhat1hc = fhat1h -mean(fhat1h)
    
    yi2 = yi[1:N]- fhat1n - fhat3
    g2<-glm(yi2~0+b2n)
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    r2 = g2$residuals
    for (i in 1:n) {fhat2h[i,1] = tbval2h[i,]%*%g2$coefficients[1:(basis2h)]}
    fhat2hc = fhat2h - mean(fhat2h)
    
    yi3 = yi[1:N]- fhat1n - fhat2n
    g3<-glm(yi3~0+b3n)
    r3 = g3$residuals
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    for (i in 1:n) {fhat3h[i,1] = tbval3h[i,]%*%g3$coefficients[1:(basis3h)]}
    fhat3hc = fhat3h - mean(fhat3h)
  }
  
  if(type=="Poisson"){
    f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "poisson",weights = weights)
    #fhat0 = tbval_z %*% f$coefficients[1:basis0]
    fhat1 = b1%*%f$coefficients[1:basis1]
    fhat2 = b2%*%f$coefficients[(basis1+1):(basis1+basis2)]
    fhat3 = b3%*%f$coefficients[(basis1+basis2+1):(basis1+basis2+basis3)]
    
    mui1 = mui_c - fhat2 - fhat3
    g1<-glm(mui1~0+b1n)
    r1 = g1$residuals
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    for (i in 1:n) {fhat1h[i,1] = tbval1h[i,]%*%g1$coefficients[1:(basis1h)]}
    fhat1hc = fhat1h -mean(fhat1h)
    
    mui2 = mui_c - fhat1n - fhat3
    g2<-glm(mui2~0+b2n)
    r2 = g2$residuals
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    for (i in 1:n) {fhat2h[i,1] = tbval2h[i,]%*%g2$coefficients[1:(basis2h)]}
    fhat2hc = fhat2h - mean(fhat2h)
    
    mui3 = mui_c - fhat1n - fhat2n
    g3 <- glm(mui3~0+b3n)
    r3 = g3$residuals
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    for (i in 1:n) {fhat3h[i,1] = tbval3h[i,]%*%g3$coefficients[1:(basis3h)]}
    fhat3hc = fhat3h - mean(fhat3h)
  }
  
  if(type=="Binomial"){
    f<-glm(yi[1:N]~0+b0+b1+b2+b3,family = "binomial",weights = weights)
    #fhat0 = tbval_z %*% f$coefficients[1:basis0]
    fhat1 = b1%*%f$coefficients[1:basis1]
    fhat2 = b2%*%f$coefficients[(basis1+1):(basis1+basis2)]
    fhat3 = b3%*%f$coefficients[(basis1+basis2+1):(basis1+basis2+basis3)]
    
    mui1 = mui_c - fhat2 - fhat3
    g1<-glm(mui1~0+b1n)
    r1 = g1$residuals
    fhat1n = b1n%*%g1$coefficients[1:(basis1n)]
    for (i in 1:n) {fhat1h[i,1] = tbval1h[i,]%*%g1$coefficients[1:(basis1h)]}
    fhat1hc = fhat1h -mean(fhat1h)
    
    mui2 = mui_c - fhat1n - fhat3
    g2<-glm(mui2~0+b2n)
    r2 = g2$residuals
    fhat2n = b2n%*%g2$coefficients[1:(basis2n)]
    for (i in 1:n) {fhat2h[i,1] = tbval2h[i,]%*%g2$coefficients[1:(basis2h)]}
    fhat2hc = fhat2h - mean(fhat2h)
    
    mui3 = mui_c - fhat1n - fhat2n
    g3 <- glm(mui3~0+b3n)
    r3 = g3$residuals
    fhat3n = b3n%*%g3$coefficients[1:(basis3n)]
    for (i in 1:n) {fhat3h[i,1] = tbval3h[i,]%*%g3$coefficients[1:(basis3h)]}
    fhat3hc = fhat3h - mean(fhat3h)
  }
  
  return(list("f1" = fhat1hc[order(t1)],
              "f2" = fhat2hc[order(t1)],
              "f3" = fhat3hc[order(t1)],
              "f1n" = fhat1n,
              "f2n" = fhat2n,
              "f3n" = fhat3n,
              "r1" = r1,
              "r2" = r2,
              "r3" = r3,
              "zeta" = cbind(sort(zeta_hat[,1]),sort(zeta_hat[,2]),sort(zeta_hat[,3])),
              "z" = sort(z),
              "theta1" = g1$coefficients[1:(basis1n)],
              "theta2" = g2$coefficients[1:(basis2n)],
              "theta3" = g3$coefficients[1:(basis3n)]))
}


#######################################
## Step1 Estimate parameters and calculate T 
####################################### 

#source(paste(filepathfun,'GFVM.R',sep=''))
#source(paste(filepathfun,'H0_GFVM.R',sep=''))


for(b in 1:1){
  type = type0[b]
  for(k in 1:7){
    f1 <- function(s,t){sin(s)*(1+c[k]*cos(t))}
    f2 <- function(s,t){(1+c[k]*sin(t))*s}
    f3 <- function(s,t){cos(s)*(1+c[k]*2*t)}
    
    f1_0 <- function(s,t){sin(s)}
    f2_0 <- function(s,t){s}
    f3_0 <- function(s,t){cos(s)}
    for(l in 1:500){
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
      
      zeta = matrix(rep(0,N*K),nrow=N,ncol=K)
      zeta[,1] = pnorm(xij[,1]/(1/((1-0.5)*pi)))
      zeta[,2] = pnorm(xij[,2]/(1/((2-0.5)*pi)))
      zeta[,3] = pnorm(xij[,3]/(1/((3-0.5)*pi)))
      ## Estimate pca scores
      X = matrix(rep(0,N*length(u)),nrow=N,ncol=length(u))
      Sig = array(NA,dim=c(length(u),length(u)))
      m_hat = rep(NA)
      X[,] = xii[,]
      Sig[,] = t(X[,])%*%X[,]/(N*length(u))
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
      
      z = Generate_z(N)
      yi = Generate_g(f1,f2,f3,zeta,z,N,type = type)$yi
      mui_c = Generate_g(f1,f2,f3,zeta,z,N,type = type)$mui_c
      result <- GFVM(zeta_hat,z,yi,mui_c,N,type = type)
      resultH0 <- H0_GFVM(zeta_hat,z,yi,mui_c,N,type = type)
      
      f1hat = result$f1
      f2hat = result$f2
      f3hat = result$f3
      f1hatn = result$f1n
      f2hatn = result$f2n
      f3hatn = result$f3n
      f0hat = result$f0
      theta1 = result$theta1
      theta2 = result$theta2
      theta3 = result$theta3
      r = result$r
      
      f1hat0 = resultH0$f1
      f2hat0 = resultH0$f2
      f3hat0 = resultH0$f3
      f1hatn0 = resultH0$f1n
      f2hatn0 = resultH0$f2n
      f3hatn0 = resultH0$f3n
      f0hatn0 = resultH0$f0
      r1 = resultH0$r1
      r2 = resultH0$r2
      r3 = resultH0$r3
      theta01 = resultH0$theta1
      theta02 = resultH0$theta2
      theta03 = resultH0$theta3
      
      # calculate T
      g <- function(t){
        basis_1h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
        basis_2h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
        basis_3h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
        basis_zh = create.bspline.basis(rangeval = c(0,1),nbasis = 4,norder = 2)
        tbval1h = eval.basis(t[1],basis_1h)
        tbval2h = eval.basis(t[1],basis_2h)
        tbval3h = eval.basis(t[1],basis_3h)
        tbval_zh= eval.basis(t[2],basis_zh)
        f1hat = kronecker(tbval1h,tbval_zh)%*%theta1
        f2hat = kronecker(tbval2h,tbval_zh)%*%theta2
        f3hat = kronecker(tbval3h,tbval_zh)%*%theta3
        f1hat0 = tbval1h%*%theta01
        f2hat0 = tbval2h%*%theta02
        f3hat0 = tbval3h%*%theta03
        
        res = as.numeric(crossprod(t(cbind(f1hat,f2hat,f3hat)-cbind(f1hat0,f2hat0,f3hat0))))
        return(res)
      }
      T[k,l,b] = adaptIntegrate(g, lowerLimit = c(0, 0), upperLimit = c(1,1))[[1]]
      
      # bootstrap 
      no_cores <- detectCores(nocores) - 1
      cl2 = makeCluster(no_cores)
      registerDoParallel(cl2)
      clusterEvalQ(cl2,library(fda))
      clusterEvalQ(cl2,library(rmutil))
      clusterEvalQ(cl2,library(MASS))
      clusterEvalQ(cl2,library(pracma))
      clusterEvalQ(cl2,library(alabama))
      clusterEvalQ(cl2,library(cubature))
      
      Tb = array(0, dim = c(100,1))
      res2 <- foreach(B = 1:100, .combine = 'c')%dopar%
        {
          weights = NULL
          
          #Normal
          if(type=="Normal"){
            mui_c0 = f1hatn0 + f2hatn0 + f3hatn0 + f0hatn0
            yi0 = mui_c0 + r3*rnorm(N,0,1)
          }
          
          #Poisson
          if(type=="Poisson"){
            mui_c0 = f1hatn0 + f2hatn0 + f3hatn0 + f0hatn0
            yi0 = rpois(N,exp(mui_c0))
          }
          
          #Binomial
          if(type=="Binomial"){
            mui_c0 = f1hatn0 + f2hatn0 + f3hatn0 + f0hatn0
            yi0 = rbinom(N,exp(mui_c0)/(1+exp(mui_c0)),size=1)
          }
          
          bresult <- GFVM(zeta_hat,z,yi0,mui_c0,N,type = type,weights = weights)
          bresultH0 <- H0_GFVM(zeta_hat,z,yi0,mui_c0,N,type = type,weights = weights)
          bf1hat = bresult$f1
          bf2hat = bresult$f2
          bf3hat = bresult$f3
          btheta1 = bresult$theta1
          btheta2 = bresult$theta2
          btheta3 = bresult$theta3
          
          bf1hat0 = bresultH0$f1
          bf2hat0 = bresultH0$f2
          bf3hat0 = bresultH0$f3
          btheta01 = bresultH0$theta1
          btheta02 = bresultH0$theta2
          btheta03 = bresultH0$theta3
          
          rm(bresult)
          rm(bresultH0)
          
          gb = function(t){
            basis_1h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
            basis_2h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
            basis_3h = create.bspline.basis(rangeval = c(0,1),nbasis = 5,norder = 3)
            basis_zh = create.bspline.basis(rangeval = c(0,1),nbasis = 4,norder = 2)
            tbval1h = eval.basis(t[1],basis_1h)
            tbval2h = eval.basis(t[1],basis_2h)
            tbval3h = eval.basis(t[1],basis_3h)
            tbval_zh= eval.basis(t[2],basis_zh)
            bf1hat = kronecker(tbval1h,tbval_zh)%*%btheta1
            bf2hat = kronecker(tbval2h,tbval_zh)%*%btheta2
            bf3hat = kronecker(tbval3h,tbval_zh)%*%btheta3
            bf1hat0 = tbval1h%*%btheta01
            bf2hat0 = tbval2h%*%btheta02
            bf3hat0 = tbval3h%*%btheta03
            
            res = as.numeric(crossprod(t(cbind(bf1hat,bf2hat,bf3hat)-cbind(bf1hat0,bf2hat0,bf3hat0))))
            return(res)
          }
          Tb[B] = adaptIntegrate(gb, lowerLimit = c(0, 0), upperLimit = c(1,1))[[1]]
          rm(bf1hat,bf2hat,bf3hat,bf1hat0,bf2hat0,bf3hat0)
          res1 <- Tb[B]
        }
      stopCluster(cl2)
      rm(result,resultH0,f1hat,f2hat,f3hat,f1hat0,f2hat0,f3hat0)
      # p-value
      Tb <- res2
      pvalue = mean(Tb>=T[k,l,b])
      pfun[k,l,b] = pvalue
      # power
      cv = quantile(Tb,1-level)
      if(T[k,l,b]>cv){powerfun[k,l,b] = 1}
      print(c(b,k,l))
    }
  }
  
}

save(powerfun, file = "power_N_N300.rda")
save(pfun, file = "pvalue_N_N300.rda")
save(T,file = "T_N_N300.rda")


#i:c,j:iteration,b:type
#Normal
outputNormal <-c(mean(powerfun[1,,1]),mean(pfun[1,,1]),mean(powerfun[2,,1]),mean(pfun[2,,1]),
                 mean(powerfun[3,,1]),mean(pfun[3,,1]),mean(powerfun[4,,1]),mean(pfun[4,,1]),
                 mean(powerfun[5,,1]),mean(pfun[5,,1]))
outputNormal

#Poisson
outputPoisson <-c(mean(powerfun[1,,2]),mean(pfun[1,,2]),mean(powerfun[2,,2]),mean(pfun[2,,2]),
                 mean(powerfun[3,,2]),mean(pfun[3,,2]),mean(powerfun[4,,2]),mean(pfun[4,,2]),
                 mean(powerfun[5,,2]),mean(pfun[5,,2]))
outputPoisson

#Binomial
outputBinomial <-c(mean(powerfun[1,,3]),mean(pfun[1,,3]),mean(powerfun[2,,3]),mean(pfun[2,,3]),
                  mean(powerfun[3,,3]),mean(pfun[3,,3]),mean(powerfun[4,,3]),mean(pfun[4,,3]),
                  mean(powerfun[5,,3]),mean(pfun[5,,3]))
outputBinomial




