# Simulation2 of FANMI: goodness of fit test
#-------------------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------

library(MASS)
library(fda) # bspline
library(alabama) # which includes constrOptim.nl
library(pracma) # the sqrt of a matrix
library(rmutil) # two-dimensional integrate
library(cubature)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)

######################################
## Setup: 
######################################
N = 700
M = 20 # that many summands for X 
u = seq(1,100,1)/100 # xgrid
K = 3 # number of additive functions
n = 100
t1 <-  seq(0.001,0.99,length.out = 100)
tz <- t1
type0 = c("Normal","Poisson","Binomial")

# power
level = 0.05
pfun = array(0, dim = c(7,500,3))
powerfun = array(0, dim = c(7,500,3))
T = array(0, dim = c(7,500,3))

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
  return(list(zeta_hat=zeta_hat,zeta=zeta))
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
  meanf0 = mean(f0(z) + meanf1 + meanf2 + meanf3)
  mui_c = meanf0 + (f0(z) + meanf1 + meanf2 + meanf3 - meanf0) + f1(zeta[,1],z)-meanf1 + f2(zeta[,2],z)-meanf2 + f3(zeta[,3],z)-meanf3
  ## ---------------------------------------
  if(type=="Normal"){
    yi = mui_c + rnorm(n,0,0.1)
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

Generate_g0 <- function(f1,f2,f3,zeta,z,N,type="Normal"){
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
  meanf0 = mean(f0(z) + meanf1 + meanf2 + meanf3)
  mui_c = meanf0
  ## ---------------------------------------
  if(type=="Normal"){
    yi = mui_c + rnorm(N,0,0.1)
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

## -------------H0--------------------
H0_MT = function(zeta_hat,z,yi,mui_c,N,type="Normal",t1= seq(0.001,0.99,length.out = 100),tz = seq(0.001,0.99,length.out = 100),
                   weights = NULL){
  knot1 = as.vector(quantile(zeta_hat[,1],seq(0,1,length=4)))
  knot2 = as.vector(quantile(zeta_hat[,2],seq(0,1,length=4)))
  knot3 = as.vector(quantile(zeta_hat[,3],seq(0,1,length=4)))
  
  basis_1 = create.bspline.basis(breaks = knot1,nbasis = 5,norder = 3)
  basis_2 = create.bspline.basis(breaks = knot2,nbasis = 5,norder = 3)
  basis_3 = create.bspline.basis(breaks = knot3,nbasis = 5,norder = 3)
  basis_z = create.bspline.basis(rangeval = range(z),nbasis = 5,norder = 3)
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
  basis1 = (nbasis -1)
  basis2 = (nbasis -1)
  basis3 = (nbasis -1)
  b0 <- tbval_z
  b1 <- tbval1
  b2 <- tbval2
  b3 <- tbval3
  
  tbval1n = eval.basis(zeta_hat[,1],basis_1)
  tbval2n = eval.basis(zeta_hat[,2],basis_2)
  tbval3n = eval.basis(zeta_hat[,3],basis_3)
  tbval_zn = eval.basis(as.vector(z),basis_z)
  basis1n = (nbasis)
  basis2n = (nbasis)
  basis3n = (nbasis)
  basis0n = (nbasis)
  b1n <- tbval1n
  b2n <- tbval2n
  b3n <- tbval3n
  b0n <- tbval_zn
  fhat1n=matrix(NA,nrow=N,ncol=1)
  fhat2n=matrix(NA,nrow=N,ncol=1)
  fhat3n=matrix(NA,nrow=N,ncol=1)
  
  knot1h = as.vector(quantile(t1,seq(0,1,length=4)))
  knot2h = as.vector(quantile(t1,seq(0,1,length=4)))
  knot3h = as.vector(quantile(t1,seq(0,1,length=4)))
  basis_1h = create.bspline.basis(breaks = knot1h,nbasis = 5,norder = 3)
  basis_2h = create.bspline.basis(breaks = knot2h,nbasis = 5,norder = 3)
  basis_3h = create.bspline.basis(breaks = knot3h,nbasis = 5,norder = 3)
  basis_zh = create.bspline.basis(rangeval = range(tz),nbasis = 5,norder = 3)
  tbval1h = eval.basis(t1,basis_1h)
  tbval2h = eval.basis(t1,basis_2h)
  tbval3h = eval.basis(t1,basis_3h)
  tbval_zh= eval.basis(tz,basis_zh)
  basis1h = (nbasis)
  basis2h = (nbasis)
  basis3h = (nbasis)
  basis0h = (nbasis)
  b1h <- tbval1h
  b2h <- tbval2h
  b3h <- tbval3h
  fhat1h=matrix(NA,nrow=n,ncol=1)
  fhat2h=matrix(NA,nrow=n,ncol=1)
  fhat3h=matrix(NA,nrow=n,ncol=1)
  
  if(type=="Normal"){
    f<-glm(yi[1:N]~b0+b1+b2+b3)
    I0 = f$coefficients[1]
    fhat0 = tbval_z %*% f$coefficients[2:(basis0+1)]
    fhat1 = b1%*%f$coefficients[(basis0+2):(basis0+basis1+1)]
    fhat2 = b2%*%f$coefficients[(basis0+basis1+2):(basis0+basis1+basis2+1)]
    fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+2):(basis0+basis1+basis2+basis3+1)]
    sigma = vcov(f)
    sigma0 = sigma[2:(basis0+1),2:(basis0+1)]
    sigma1 = sigma[(basis0+2):(basis0+basis1+basis2+basis3+1),(basis0+2):(basis0+basis1+basis2+basis3+1)]
    gamma = as.matrix(f$coefficients[2:(basis0+1)])
    theta = as.matrix(f$coefficients[(basis0+2):(basis0+basis1+basis2+basis3+1)])
    
    t0 = (t(gamma)%*%solve(sigma0)%*%gamma - basis0)/sqrt(2*basis0)
    t1 = (t(theta)%*%solve(sigma1)%*%theta - (basis1+basis2+basis3))/sqrt(2*(basis1+basis2+basis3))
  }
  
  if(type=="Poisson"){
    f<-glm(yi[1:N]~b0+b1+b2+b3,family = "poisson")
    I0 = f$coefficients[1]
    fhat0 = tbval_z %*% f$coefficients[2:(basis0+1)]
    fhat1 = b1%*%f$coefficients[(basis0+2):(basis0+basis1+1)]
    fhat2 = b2%*%f$coefficients[(basis0+basis1+2):(basis0+basis1+basis2+1)]
    fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+2):(basis0+basis1+basis2+basis3+1)]
    sigma = vcov(f)
    sigma0 = sigma[2:(basis0+1),2:(basis0+1)]
    sigma1 = sigma[(basis0+2):(basis0+basis1+basis2+basis3+1),(basis0+2):(basis0+basis1+basis2+basis3+1)]
    gamma = as.matrix(f$coefficients[2:(basis0+1)])
    theta = as.matrix(f$coefficients[(basis0+2):(basis0+basis1+basis2+basis3+1)])
    
    t0 = (t(gamma)%*%solve(sigma0)%*%gamma - basis0)/sqrt(2*basis0)
    t1 = (t(theta)%*%solve(sigma1)%*%theta - (basis1+basis2+basis3))/sqrt(2*(basis1+basis2+basis3))
      }
  
  if(type=="Binomial"){
    f<-glm(yi[1:N]~b0+b1+b2+b3,family = "binomial")
    I0 = f$coefficients[1]
    fhat0 = tbval_z %*% f$coefficients[2:(basis0+1)]
    fhat1 = b1%*%f$coefficients[(basis0+2):(basis0+basis1+1)]
    fhat2 = b2%*%f$coefficients[(basis0+basis1+2):(basis0+basis1+basis2+1)]
    fhat3 = b3%*%f$coefficients[(basis0+basis1+basis2+2):(basis0+basis1+basis2+basis3+1)]
    sigma = vcov(f)
    sigma0 = sigma[2:(basis0+1),2:(basis0+1)]
    sigma1 = sigma[(basis0+2):(basis0+basis1+basis2+basis3+1),(basis0+2):(basis0+basis1+basis2+basis3+1)]
    gamma = as.matrix(f$coefficients[2:(basis0+1)])
    theta = as.matrix(f$coefficients[(basis0+2):(basis0+basis1+basis2+basis3+1)])
    
    t0 = (t(gamma)%*%solve(sigma0)%*%gamma - basis0)/sqrt(2*basis0)
    t1 = (t(theta)%*%solve(sigma1)%*%theta - (basis1+basis2+basis3))/sqrt(2*(basis1+basis2+basis3))
  }
  
  return(list("gamma" = gamma,
              "theta" = theta,
              "t0" = t0,
              "t1" = t1
              ))
}

######################################
## End 
######################################

#######################################
## Step1 Estimate parameters and calculate T 
####################################### 

for(b in 3:3){
  type = type0[b]
  for(k in 1:1){
    f0 <- function(s){0.2*sin(2*pi*s)}
    f1 <- function(s,t){0.2*2*cos(2*pi*s)*t^2}
    f2 <- function(s,t){0.2*(3*(s-0.5)^2+sin(2.5*pi*s))*t}
    f3 <- function(s,t){0.2*5*t*cos(s)}
    
    for(l in 1:500){
      GX = Generate_X(N,M,u,K)
      zeta = GX$zeta
      zeta_hat = GX$zeta_hat
      z = Generate_z(N)
      yi = Generate_g(f1,f2,f3,zeta,z,N,type = type)$yi
      mui_c = Generate_g(f1,f2,f3,zeta,z,N,type = type)$mui_c
      resultH0 <- H0_MT(zeta_hat,z,yi,mui_c,N,type = type)
      
      # calculate T
      t0 = resultH0$t0
      t1 = resultH0$t1
      T[k,l,b] = t0^2 + t1^2
      
      # p-value
      Tb <- qchisq(0.95,2)
      pfun[k,l,b] = 1 - pchisq(T[k,l,b],2)
      # power
      cv = qchisq(0.95,2)
      if(T[k,l,b]>cv){powerfun[k,l,b] = 1}
      print(c(b,k,l))
    }
  }
  
}





