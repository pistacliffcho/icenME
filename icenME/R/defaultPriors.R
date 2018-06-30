alpha_dft = function(vals, info){
  ans = sum( dnorm(vals, sd = 10, log = T) )
  return(ans)
}

beta_dft = function(vals, info){
  ans = sum( dnorm(vals, sd = 2, log = T) )
  return(ans)
}

gamma_mu_dft = function(vals, info){
  ans = sum( dnorm(vals, sd = 2, log = T) )
  return(ans)
}

sqrt_S_dft = function(vals, info){
  ans = sum( dnorm(vals, sd = 2, log = T) )
  return(ans)
}