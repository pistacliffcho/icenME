
#' @export
make_gamma_sampler = function(model, 
                              gamma1_ind, 
                              gamma2_ind){
  info = gamma_MH_info
  info$model = model
  info$sd = 1
  info$ind1 = gamma1_ind
  info$ind2 = gamma2_ind
  ans = MHsampler(info = info, 
                  propGenerator = simpleGenerator, 
                  model = model)
  return(ans)
}

gamma_MH_info = list(
  getVals = function(info){
    model = info$model
    ind1 = info$ind1
    ind2 = info$ind2
    ans = get_gamma(model, ind1, ind2)
    return(ans)
  }, 
  setVals = function(vals, info){
    model = info$model
    ind1 = info$ind1
    ind2 = info$ind2
    set_gamma(model, vals, ind1, ind2)
  }, 
  get_lp = function(vals, info){
    model = info$model
    ind1 = info$ind1
    ind2 = info$ind2
    calc_gamma(model, vals, ind1, ind2)
  }
)

#' @export
make_gmuS_sampler = function(model, gamma_id){
  info = gammaMuS_MH_info
  info$model = model
  info$sd = 1
  info$ind = gamma_id
  ans = MHsampler(info = info, 
                  propGenerator = simpleGenerator, 
                  model = model)
  return(ans)
}

# The info object should have
#   gamma_mu_node
#   sqrt_S_node
#   gamma_nodes
#   n_gamma
gammaMuS_MH_info = list(
  getVals = function(info){
    model = info$model
    ind = info$ind
    ans = get_gmus(model, ind)
    return(ans)
  }, 
  setVals = function(vals, info){
    model = info$model
    ind = info$ind
    set_gmus(model, vals, ind)
  }, 
  get_lp = function(vals, info){
    model = info$model
    ind = info$ind
    calc_gmus(model, vals, ind)
  }
)

#' @export
make_ab_sampler = function(model){
  info = alphaBeta_MH_info
  info$model = model
  info$sd = 1
  ans = MHsampler(info = info, 
                  propGenerator = simpleGenerator, 
                  model = model)
  return(ans)
}

# Info requires:
#   alpha
#   beta
#   n_alpha
#   eta
#   data_node
alphaBeta_MH_info = list(
  getVals = function(info){
    model = info$model
    get_ab(model)
  }, 
  setVals = function(vals, info){
    model = info$model
    set_ab(model, vals)
  }, 
  get_lp = function(vals, info){
    model = info$model
    calc_ab(model, vals)
  }
)

# make_sampler = function(model, info){
#   ans = MHsampler()
#   ans$getVals = info$getVals
#   ans$setVals = info$setVals
#   ans$get_lp = info$get_lp
#   ans$model = model
#   ans$info = list(sd = 0.25)
#   ans$propGenerator = simpleGenerator
#   return(ans)
# }

MHsampler = setRefClass("MHsampler", 
                        fields = c("propGenerator", 
                                   "getVals", 
                                   "setVals", 
                                   "get_lp", 
                                   "model", 
                                   "info", 
                                   "num_accept"
                                   ))
MHsampler$methods(
  initialize = function(info, propGenerator, model){
    getVals <<- info$getVals
    setVals <<- info$setVals
    get_lp <<- info$get_lp
    propGenerator <<- propGenerator
    info <<- info
    model <<- model
    num_accept <<- 0
  }
)

MHsampler$methods(
  decide = function(new_lp, 
                    old_lp){
    ratio = exp(new_lp - old_lp)
    ans = ratio > runif(1)
    return(ans)
  }
)

simpleGenerator = function(old_vals, info){
  n = length(old_vals)
  ans = old_vals + rnorm(n, sd = info$sd)
  return(ans)
}

adaptiveBlock_generator = function(old_vals, info){
  n = ncol(info$chol)
  move = rnorm(n) %*% info$chol
  move = as.numeric(move)
  ans = old_vals + move
  return(ans)
}

update_block_chol = function(sample_mat){
  cov_mat = cov(sample_mat)
  diag(cov_mat) = diag(cov_mat) + 0.01
  cov_mat = cov_mat / ncol(sample_mat)^2
  ans = chol(cov_mat)
  return(ans)
}


MHsampler$methods(
  update = function(){
    org_vals = getVals(info)
    prop_vals = propGenerator(org_vals, info)
    # Nodes will skip setVal if new_val is NULL
    new_lp = get_lp(prop_vals, info)
    old_lp = get_lp(org_vals, info)
    jump = decide(new_lp = new_lp, old_lp = old_lp)
    if(jump){
      setVals(prop_vals, info)
      num_accept <<- num_accept + 1
    }
  }
)