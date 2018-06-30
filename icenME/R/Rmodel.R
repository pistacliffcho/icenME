#' @useDynLib icenME
#' @import Rcpp


#' @export
makeAB_list = function(aFxn = alpha_dft, 
                       bFxn = beta_dft, 
                       aInit = c(0,0), 
                       bInit){
  ans = list(
    alpha_fxn = aFxn, 
    beta_fxn = bFxn, 
    a_init = aInit, 
    b_init = bInit
  )
  return(ans)
}


#' @export
makeModel = function(initValues, 
                     initData, 
                     initPriors, 
                     respMat, LLK, isCen){
  ans = cMakeModel(initValues, 
                   initData,
                   initPriors)
  return(ans)
}

#' @export
get_ab = function(model) c_get_ab_vals(model)

#' @export
calc_ab = function(model, vals) c_setAndCalc_ab(model, vals)

#' @export
set_ab = function(model, vals) c_set_ab(model, vals)

#' @export
calc_gmus = function(model, vals, ind) c_setAndCalc_gmus(model, vals, ind)

#' @export
get_gmus = function(model, ind) c_get_gmus_vals(model, ind)

#' @export
get_gamma = function(model, ind1, ind2) c_get_gamma_vals(model, ind1, ind2)

#' @export
set_gmus = function(model, vals, ind) c_set_gmus(model, vals, ind)

#' @export 
set_gamma = function(model, vals, ind1, ind2) c_set_gamma(model, vals, ind1, ind2)

#' @export 
calc_gamma = function(model, vals, ind1, ind2) c_setAndCalc_gamma(model, vals, ind1, ind2)

#' @export
full_lp = function(model) c_full_lp(model)