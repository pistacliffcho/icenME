expandREData = function(data, 
                        resp_form = NULL,
                        fixed_form = NULL, 
                        rand_forms = NULL, 
                        cenTol = 0.001){
  
  resp_mat = model.matrix(resp_form, data)
  isCen = (resp_mat[,2] - resp_mat[,1]) > cenTol
  fixed_mat = NULL
  rand_list = NULL
  
  if(!is.null(fixed_form)){
    fixed_mat = model.matrix(fixed_form, data)
  }
  if(!is.null(rand_forms)){
    all_rand_list = list()
    if(is(rand_forms, "formula")){ rand_forms = list(rand_forms) }
    for(i in seq_along(rand_forms)){
      this_form = rand_forms[[i]]
      this_id_name = as.character(this_form)[2]
      id = data[,this_id_name]
      this_mat = model.matrix(this_form, data)
      rownames(this_mat) = 1:nrow(this_mat)
      rand_list = split.data.frame(this_mat, id)
      this_rand_list = list()
      for(ii in seq_along(rand_list)){
        this_name = names(rand_list)[ii]
        this_data = rand_list[[ii]]
        this_ind = as.numeric(rownames(this_data) )
        this_list = list(data = this_data, rows = this_ind)
        this_rand_list[[this_name]] = this_list
      }
      all_rand_list[[this_id_name]] = this_rand_list
    }
  }
  ans = list(resp = resp_mat, 
             fixed_cov = fixed_mat, 
             rand_cov = all_rand_list, 
             isCen = isCen)
  return(ans)
}


# For a given set of covariates (rand or fixed)
# generates named random normal initial values
gen_ind_inits = function(xmat, sd = 0.1){
  cnames = colnames(xmat)
  ans = rnorm(length(cnames), sd = sd)
  names(ans) = cnames
  return(ans)
}

# Generate random chol(Sigma) random effects
gen_sqrt_S = function(example_mat){
  k = ncol(example_mat)
  vals = rnorm(k^2, sd = 0.5)
  ans = matrix(vals, nrow = k)
  return(ans)
}

generate_inits = function(data_list, alpha_len = 2){
  ans = list()
  ans$alpha = rnorm(alpha_len, sd = 0.5)  
  fixed_cov = data_list$fixed_cov
  if(!is.null(fixed_cov)){ ans$beta = gen_ind_inits(fixed_cov) }
  
  rand_cov = data_list$rand_cov
  if(!is.null(rand_cov)){
    gamma = list()
    gamma_mu = list()
    sqrt_S = list()
    r_id_names = names(rand_cov)
    for(i in seq_along(r_id_names)){
      this_id_name = r_id_names[i]
      r_eff_names = names(rand_cov[[this_id_name]])
      this_id_data = rand_cov[[this_id_name]]
      this_id_gamma = list()
      for(j in seq_along(r_eff_names)){
        this_name = r_eff_names[j]
        this_xmat = this_id_data[[this_name]]$data
        this_id_gamma[[this_name]] = gen_ind_inits(this_xmat)
      }
      gamma[[this_id_name]] = this_id_gamma
      gamma_mu[[this_id_name]] = rnorm(ncol(this_xmat), sd = 0.1)
      sqrt_S[[this_id_name]] = gen_sqrt_S(this_xmat)
    }
    ans$gamma = gamma
    ans$gamma_mu = gamma_mu
    ans$sqrt_S = sqrt_S 
  }
  return(ans)
}

generate_default_priors = function(init_list){
  ans = list(alpha = alpha_dft, 
       beta = beta_dft)
  gamma_mu_list = list()
  sqrt_S_list = list()
  gamma_list = list()
  for(i in seq_along(init_list$gamma)){
    this_id = names(init_list$gamma)[i]
    gamma_mu_list[[this_id]] =  gamma_mu_dft
    sqrt_S_list[[this_id]] = sqrt_S_dft
    this_gamma_list = list()
    for(j in seq_along(init_list$gamma[[i]])){
      this_name = names(init_list$gamma[[i]])[j]
      this_gamma_list[[this_name]] = icenME_dmvnorm
    }
    gamma_list[[this_id]] = this_gamma_list
  }
  ans$gamma_mu = gamma_mu_list
  ans$sqrt_S = sqrt_S_list
  ans$gamma = gamma_list
  ans$data_llk = weib_ph$computeLLK
  return(ans)
}
