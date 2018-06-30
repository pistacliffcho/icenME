# Class
baseInfo = setRefClass("baseInfo", 
                       fields = c("log_dfun", "sfun"))

# Instances
weib_base = baseInfo()
weib_base$log_dfun = function(y, alpha){
  pars = exp(alpha)
  ans = dweibull(y, shape = pars[1], scale = pars[2], 
                 log = TRUE)
  return(ans)
}
weib_base$sfun = function(y, alpha){
  pars = exp(alpha)
  ans = pweibull(y, shape = pars[1], scale = pars[2], 
                 lower.tail = FALSE)
  return(ans)
}


### LINK INFO CLASS

# Class
linkInfo = setRefClass("linkInfo", 
                       fields = c("dlink", "slink"))
# Instances
aft_info = linkInfo()
aft_info$dlink = function(y, eta, alpha, dfun, sfun){
  ans = log_dfun(y/eta, alpha) - log(eta)
  return(ans)
}
aft_info$slink = function(y, eta, alpha, sfun){
  ans = sfun(y / eta, alpha)
  return(ans)
}

ph_info = linkInfo()
ph_info$dlink = function(y, eta, alpha, dfun, sfun){
  ans = dfun(y, eta) + log(eta) + (eta - 1) * sfun(y, alpha)
  return(ans)
}
ph_info$slink = function(y, eta, alpha, sfun){
  ans = sfun(y, alpha)^eta
  return(ans)
}

## CLASS FOR COMPUTING LLK

data_llkInfo = setRefClass("data_llkInfo", 
                           fields = c("link", 
                                      "baseInfo"))

make_data_llkInfo = function(link = ph_info, 
                             base = weib_base){
  ans = data_llkInfo()
  ans$link = link
  ans$baseInfo = base
#  ans$link$sfun = ans$baseInfo
  ans$link
  return(ans)
}

weib_ph = make_data_llkInfo()

data_llkInfo$methods(
  computeLLK = function(y, eta, alpha, isCen, inds = "all"){
    ans = 0
    use_y = y
    useIsCen = isCen
    use_eta = eta
    if(!identical(inds, "all")){
      use_y = y[inds,]
      use_eta = eta[inds]
      useIsCen = isCen[inds]
    }
    unCens_y = use_y[!useIsCen, ]
    unCens_eta = use_eta[!useIsCen]
    if(length(unCens_y) > 0){
      ans = ans + link$dlink(unCens_y, 
                             unCens_eta, 
                             alpha, baseInfo$log_dfun, baseInfo$sfun)
    }
    cens_y = use_y[useIsCen,]
    cens_eta = use_eta[useIsCen]
    if(length(cens_y) > 0){
      probs_mat = link$slink(cens_y, cens_eta, 
                             alpha, baseInfo$sfun)
      probs = probs_mat[,1] - probs_mat[,2]
      ans = ans + sum(log(probs))
    }
    return(ans)
  }
)
