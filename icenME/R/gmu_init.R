icenME_dmvnorm = function(x, mu, sqrt_sigma){
  sigma = t(sqrt_sigma) %*% sqrt_sigma
  part1 = -0.5 * (log(det(sigma)) + length(mu) * log(pi * 2))
  part2 = -t(x - mu) %*% solve(sigma) %*% (x-mu)/2
  ans =  part1 + part2
  return(ans)
}


make_gamma_list = function(re_logDens = quick_dmvnorm, 
                           gamma_inits){
  ans = list(g_fxn = re_logDens, 
             gamma_inits = gamma_inits)
  return(ans)
}

make_gmus_list = function(gmu_prior = gamma_mu_dft, 
                          s_prior = sqrt_S_dft, 
                          gmu_inits, 
                          S_inits){
  n_gamma = length(gmu_inits)
  if(length(S_inits) != n_gamma^2){
    stop("length of s_inits should be length(g_inits)^2")
  }
  ans = list(
    g_fxn = gmu_prior, 
    s_fxn = s_prior, 
    gmu_init = gmu_inits,
    S_init = S_inits
  )
  return(ans)
}