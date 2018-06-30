modelValues = setRefClass("modelValues", 
                          fields = c("alphaBeta",  
                                     "gmus", 
                                     "gamma", "lp", 
                                     "model"))

modelValues$methods(
  initialize = function(model, samples = 100){
    model <<- model
    lp <<- rep(NA, samples)
    ab_vals = get_ab(model)
    n_ab = length(ab_vals)
    alphaBeta <<- matrix(NA, nrow = samples, ncol = n_ab)
    
    gamma_ids = as.character(seq_len(n_gmus(model)) )
    gmus <<- new.env()
    gamma <<- new.env()
    for(i in seq_along(gamma_ids)){
      gid = gamma_ids[i]
      gmus_vals = get_gmus(model, i)
      gmus_size = length(gmus_vals)
      
      gmus[[gid]] <<- matrix(NA, ncol = gmus_size, nrow = samples)
      
      num_gamma = n_gamma(model, i)
      gamma_names = as.character(seq_len(num_gamma))
      this_gamma_env = new.env()
      for(j in seq_along(gamma_names)){
        gid2 = gamma_names[j]
        gamma_vals = get_gamma(model, i, j)
        gamma_size = length(gamma_vals)
        this_gamma_env[[gid2]] = matrix(NA, ncol = gamma_size, 
                                        nrow = samples)
      }
      gamma[[gid]] <<- this_gamma_env
    }
  }
)

modelValues$methods(
  addRow = function(row){
    lp[row] <<- full_lp(model)
    alphaBeta[row,] <<- get_ab(model)
    num_gmus <- n_gmus(model)
    for(i in seq_len(num_gmus)){
      gmus_name = names(gmus)[i]
      gmus[[gmus_name]][row, ] <<- get_gmus(model, i)
      num_gamma = n_gamma(model, i)
      for(j in seq_len(num_gamma)){
        gamma_name = names(gamma[[gmus_name]])[j]
        gamma[[gmus_name]][[gamma_name]][row,] <<- get_gamma(model, i, j)
      }
    }
  })