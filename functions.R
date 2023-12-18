# Kernel-smoothed nonparametric maximum likelihood
km_ic = function(dt, 
                 h, 
                 resam = rep(1, nrow(dt)))
{
  L = dt$L; R = dt$R; x1 = dt$x1; x2 = dt$x2
  x = cbind(x1, x2)
  K = dnorm(outer(x1, x1, "-"), sd = h) * dnorm(outer(x2, x2, "-"), sd = h)
  weights = (K / colSums(K)) * resam
  t = sort(unique(c(L[L > -Inf], R[R < Inf])))
  m = length(t)
  delta = ifelse(L == R, 1, 0)
  lambda = lambda_old = rep(1 / m, m)
  Rstar = ifelse(R < Inf, R, L)
  atrisk = outer(t, Rstar, "<=")
  ind_eq = outer(t, Rstar, "==")
  maxiter = 30; tol = 1e-4; error = 10; iter = 0
  while (iter < maxiter & error > tol) 
  {
    #Estep
    S1 = drop(outer(L, t, ">=") %*% lambda)
    S2 = drop(outer(R, t, ">=") %*% lambda)
    temp = 1 / (1 - exp(S1 - S2))
    temp = ifelse(S1 < S2, temp, 0)
    EW = outer(temp, lambda, "*") + t(outer(lambda, rep(1, length(L)), "*") * (1 - atrisk))
    EW = EW + t(I(delta == 1) * ind_eq)
    EW[outer(L, t, ">=")] = 0
    EW[outer(R, t, "<")] = 0
    EW[outer(R, rep(Inf, m), "==")] = 0
    #Mstep
    num = apply(crossprod(outer(Rstar, t, ">=") * EW, weights), 1, sum)
    denum = apply(crossprod(outer(Rstar, t, ">="), weights), 1, sum)
    lambda = num / pmax(denum, 1e-5)
    error = sum((lambda - lambda_old)^2)
    iter = iter + 1
    lambda_old = lambda
  }
  dt = data.frame(time = t, suv = exp(-cumsum(lambda)))
  
  return(dt)
}
d
# Weight estimation
wfunc = function(dt, 
                 tau, 
                 type, 
                 h, 
                 resam = rep(1, nrow(dt)))
{
  L = dt$L ; R = dt$R; Y = ifelse(R < Inf, R, L)
  delta = dt$delta; x1 = dt$x1; x2 = dt$x2
  if (type == "em") 
    km = km_ic(dt = dt, h = h, resam = resam)
  else if (type == "rf") {
    require(icrf)
    fit = icrf(Surv(exp(L), exp(R), type = "interval2") ~ x1 + x2, data = dt, 
               ntree = 30, nfold = 1, do.trace = FALSE, bandwidth = 0)
    survs = predict(fit)
    time = attr(survs, "time")[1 : length(colnames(survs))]
    Fstar = 0
    for (i in 1 : nrow(survs)) 
      Fstar[i] = approx(time, survs[i, ], exp(Y)[i])$y
    Fstar = pmax(1 - Fstar, 1e-5)
    km = data.frame(time = Y, suv = 1 - Fstar)
  }
  FL = 1 - approx(c(-Inf, min(km$time) - 100, km$time, 1000), c(1, 1, km$suv, 0), L)$y
  FR = 1 - approx(c(-1000, min(km$time) - 100, km$time, Inf), c(1, 1, km$suv, 0), R)$y
  wt = case_when((FL < tau & tau < FR) & delta != 1 ~ (tau - FL) / (FR - FL), 
                 FL >= tau & delta != 1 ~ 1, FR <= tau & delta != 1 ~ 0,
                 delta == 1 ~ 1, TRUE ~ 0)
  
  return(wt)
}

# Interval-censored quantile regression estimator
icrq = function(dt, 
                tau, 
                type, 
                resam = rep(1, nrow(dt)), 
                h = nrow(dt)^(-1/3))
{
  n = nrow(dt)
  delta = dt$delta; x1 = dt$x1; x2 = dt$x2; xx0 = cbind(1, x1, x2)
  L = dt$L ; R = dt$R
  Rstar = ifelse(delta == 1, L, ifelse(R < Inf, R, L))
  Mp = max(R[R < Inf]) + 100; Mm = min(L[L > -Inf]) - 100
  ww0 = wfunc(dt = dt, tau = tau, type = type, h = h, resam = resam)
  YL = ifelse(L == -Inf, Mm, L)
  YR = ifelse(R == Inf, Mp, R)
  xx = rbind(xx0, xx0)
  yy = c(YL, YR)
  ww = c(ww0, 1 - ww0)
  rr = rep(resam, 2)
  est = rq.wfit(xx, yy, tau, weights = ww * rr)$coef
  
  return(est)
}

# Bootstrap or perturbed-resampling variance
crq_se = function(dt,
                  B,
                  tau,
                  type) 
{
  mc = 2; cores = 2
  n = nrow(dt)
  if (type == "rf") 
  {
    registerDoParallel(mc, cores = cores)
    res = foreach(i = 1 : B, .combine = cbind) %dopar% icrq(dt = dt[sample(n, n, replace = TRUE), ], tau = tau, type = type)
  } else if (type == "em") 
  {
    registerDoParallel(mc, cores = cores)
    res = foreach(i = 1 : B, .combine = cbind) %dopar% icrq(dt = dt, resam = rexp(n), tau = tau, type = type)
  }
  se = apply(res, 1, sd)
  
  return(se)
}