##############################################################################
## Interval-censored linear quantile regression                              #
## By Taehwa Choi, Seohyeon Park, Hunyong Cho and Sangbum Choi               #
## Update: July 19, 2023                                                     #
##############################################################################
# Required packages
library(evd)
library(tidyverse)
library(quantreg)
library(survival)
library(doParallel)
library(foreach)

# Load basic functions
source("functions.R")

# Data preparation 
beta0 = c(1.5, 1, 1) # True beta
n = 200              # Sample size
tau = 0.3            # Quantile level
set.seed(230719)
x1 = runif(n, -1, 1)
x2 = rbinom(n, 1, 0.5)
xx = cbind(1, x1, x2)
eps = rgumbel(n, -1, 1)
sig = 1 + 0.5 * (1 - x1)^2
Time = exp(xx %*% beta0 + (eps - quantile(eps, tau)) * sig)
Cens = runif(n, 30, 50)
int = outer(1 : n, c(0, 0), "*")
delta = 0
prob = rep(0, n)
for (i in 1 : n) 
{
  ind = rbinom(1, 1, prob[i])
  if (ind == 1) 
  { 
    if (Time[i] <= Cens[i]) 
    {
      int[i, ] = rep(Time[i], 2) 
      delta[i] = 1
    } else 
    {
      int[i, ] = c(Cens[i], Inf)
      delta[i] = 2
    } 
  } else 
  {
    utmp = U = 0
    j=1
    while (utmp < Cens[i]) 
    {
      U[j + 1] = U[j] + runif(1, 0.1, 1)
      utmp = U[j + 1]
      j = j + 1
    }
    uu = U[2 : (j - 1)] 
    
    if (Time[i] < min(uu)) 
    { 
      int[i, ] =  c(0, min(uu))
      delta[i] = 3
    } else if (Time[i] > max(uu)) 
    { 
      int[i, ] = c(max(uu), Inf)
      delta[i] = 2
    } else 
    {
      int[i, ] = c(max(uu[Time[i] > uu]), min(uu[Time[i] < uu]))
      delta[i] = case_when(int[i, 1] == 0 ~ 3, TRUE ~ 4)
    }
  }
}
dt = data.frame(log(int), delta = delta, x1, x2, log(Time))
colnames(dt) = c("L", "R", "delta", "x1", "x2", "T")

# Fitting ICQR-KS 
icrq(dt = dt, tau = tau, type = "em")
# crq_se(dt = dt, B = 50, tau = tau, type = "em")

# Fitting ICQR-RF
icrq(dt = dt, tau = tau, type = "rf")
# crq_se(dt = dt, B = 50, tau = tau, type = "rf")