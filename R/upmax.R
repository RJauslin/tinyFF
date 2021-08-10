
#'pik=c(0.07,0.17,0.41,0.61,0.83,0.91)
#'eps = 1e-06
#'
#'
#'
UPmaxentropy <- function (pik) 
{
  n = sum(pik)
  if (n >= 2) {
    pik2 = pik[pik != 1]
    n = sum(pik2)
    
    # piktilde 
    piktilde = UPMEpiktildefrompik(pik2)
    
    w = piktilde/(1 - piktilde)
    
    q = UPMEqfromw(w, n)
    
    s2 = UPMEsfromq(q)
    
    
    s = rep(0, times = length(pik))
    s[pik == 1] = 1
    s[pik != 1][s2 == 1] = 1
  }
  if (n == 0) 
    s = rep(0, times = length(pik))
  if (n == 1) 
    s = as.vector(rmultinom(1, 1, pik))
  s
}



UPMEpiktildefrompik <- function (pik, eps = 1e-06) 
{
  n = sum(pik)
  
  pikt = pik
  arr = 1
  while (arr > eps) {
    
    w = (pikt)/(1 - pikt)
    q = UPMEqfromw(w, n)
    pikt1 = pikt + pik - UPMEpikfromq(q)
    arr = sum(abs(pikt - pikt1))
    pikt = pikt1
  }
  pikt
}


UPMEqfromw <- function (w, n) {
  N = length(w)
  expa = array(0, c(N, n))
  for (i in 1:N){
    expa[i, 1] = sum(w[i:N])
  } 
  for (i in (N - n + 1):N){
    expa[i, N - i + 1] = exp(sum(log(w[i:N])))
  } 
  for (i in (N - 2):1){
    for (z in 2:min(N - i, n)) {
      expa[i, z] = w[i] * expa[i + 1, z - 1] + expa[i + 1, z]
    }
  } 
  q = array(0, c(N, n))
  for (i in N:1){
    q[i, 1] = w[i]/expa[i, 1]
  } 
  for (i in N:(N - n + 1)){
    q[i, N - i + 1] = 1
  }
  for (i in (N - 2):1){
    for (z in 2:min(N - i, n)){
      q[i, z] = w[i] * expa[i + 1, z - 1]/expa[i, z]
    }
  }
  return(q)
}


# selet poisson sample
UPMEsfromq <- function (q) {
  n = ncol(q)
  N = nrow(q)
  s = rep(0, times = N)
  for (k in 1:N){
    if (n != 0) {
      if (runif(1) < q[k, n]) {
      s[k] = 1
      n = n - 1
      }
    }
  }
  s
}
