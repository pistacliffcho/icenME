MCMC = setRefClass("MCMC", 
                   fields = c("samplers", 
                              "model"))
MCMC$methods(
  reset_accept = function(){
    samplers$alphaBeta$num_accept <<- 0
    gamma_id = names(samplers$gamma)
    for(gid in seq_along(samplers$gamma)){
      samplers$gmuS[[gid]]$num_accept <<- 0
      gnames = names(samplers$gamma[[gid]])
      these_gammas = samplers$gamma[[gid]]
      for(gn in gnames){
        these_gammas[[gn]]$num_accept <- 0
      }
    }
  }
)

MCMC$methods(
  update = function(){
    samplers$alphaBeta$update()
    num_gmuS = length(samplers$gmuS)
    for(gid in seq_len(num_gmuS)){
      these_gamma_samplers = samplers$gamma[[gid]]
      num_gamma = length(these_gamma_samplers)
      for(gn in seq_len(num_gamma)){
        these_gamma_samplers[[gn]]$update()
      }
      samplers$gmuS[[gid]]$update()
    }
  }
)

MCMC$methods(
  get_lp = function(){
    ans = fullModel_lp(model)
    return(ans)
  }
)

MCMC$methods(
  update_adapt = function(mv){
    ab_vals = cbind(mv$alpha, mv$beta)
    samplers$alphaBeta$info$chol <<- update_block_chol(ab_vals)
    
    gamma_ids = names(mv$gamma)
    for(gid in gamma_ids){
      this_gmuS_smp = samplers$gmuS[[gid]]
      these_gmuS_mv = cbind(mv$gamma_mu[[gid]], mv$sqrt_S[[gid]])
      this_gmuS_smp$info$chol = update_block_chol(these_gmuS_mv)
      gamma_names = names(samplers$gamma[[gid]])
      for(gn in gamma_names){
        these_gamma_vals = mv$gamma[[gid]][[gn]]
        samplers$gamma[[gid]][[gn]]$info$chol <<- 
          update_block_chol(these_gamma_vals)
      }
    }
  }
)

MCMC$methods(
  run = function(MC = 1000, thin = 5, 
                 burn = 1000, 
                 burn2 = 500){
  mv = modelValues(model = model, 
                   samples = MC / thin)
  for(i in 1:burn){
    update()
  }
  for(i in 1:MC){
    update()
    if(i %% thin == 0)
      mv$addRow(i/thin)
  }
  ans = mv
  return(mv)
})

#' @export
make_MCMC = function(model){
  # Building samplers
  sampler_env = new.env()
  sampler_env$alphaBeta = make_ab_sampler(model)
  
  n_gmu = n_gmus(model)
  
  gmus_list = list()
  gamma_list = list()
  
  for(i in seq_len(n_gmu)){
    this_gmus_sampler = make_gmuS_sampler(model = model, i)
    gmus_list[[i]] = this_gmus_sampler
    n_gamma = n_gamma(model, i)
    this_gamma_list = list()
    for(j in seq_len(n_gamma)){
      this_gamma_sampler = make_gamma_sampler(model = model,i, j)
      this_gamma_list[[j]] = this_gamma_sampler
    }
    gamma_list[[i]] = this_gamma_list
  }
  sampler_env$gmuS = gmus_list
  sampler_env$gamma = gamma_list
  ans = MCMC()
  ans$samplers = sampler_env
  ans$model = model
  return(ans)
}