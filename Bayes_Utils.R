# Given piors proportions and common coefficient of variation of the beta priors
# convert the corresponding parameters of beta priors
Beta_par_convert = function(p, tau){
  a = (1 - p - (tau^2)*p)/(tau^2)
  b = a*(1 - p)/p
  return(cbind(a, b))
}

Beta_string = function(x){
  str_tp = stringr::str_c(x, collapse = ", ")
  return(paste0("B(", str_tp, ")"))
}

Norm_string = function(x){
  str_tp = stringr::str_c(x, collapse = ", ")
  return(str_tp)
}

Norm_IA_string = function(x){
  str_tp = stringr::str_c(c("Look ", x), collapse = "")
  return(str_tp)
}

Bin_cumsum = function(sim_data, n){
  return(sum(sim_data[1:x]))
}

Beta_par_method1 = function(m, s2){
  p = m*(m*(1-m)/s2 - 1)
  q = (1 - m)*(m*(1-m)/s2 - 1)
  return(c(p, q))
}

Beta_par_method2 = function(MODE, AUC, L, U){
  f = function(b){
    res = pbeta(U, (MODE*b - 2*MODE + 1)/(1 - MODE), b) - 
      pbeta(L, (MODE*b - 2*MODE + 1)/(1 - MODE), b) - AUC
    return(res)
  }
  
  b = uniroot(f, lower = 1, upper = 200)$root
  a = (MODE*b - 2*MODE + 1)/(1 - MODE)
  return(c(a, b))
}