library(icenReg)
library(icenME)

# Testing expandREData
n = 40
simdata = simIC_weib(n)
simdata$id = rep(c('cat','dog'), n/2)
simdata$id2 = c(rep("cow", n/2), rep("hen", n/2))

resp_form = ~ cbind(l, u) + 0 
fixed_form = ~ x1 + 0
id_form = id ~ 1 + x2
id2_form = id2 ~ 1 + x2

data_list = icenME:::expandREData(simdata, 
                                      resp_form = resp_form, 
                                      fixed_form = fixed_form, 
                                      rand_forms = list(id_form, id2_form) 
)

init_vals = icenME:::generate_inits(data_list, alpha_len = 2)
init_priors = icenME:::generate_default_priors(init_vals)

my_model = makeModel(init_vals, data_list, init_priors)
get_ab(my_model)
get_gmus(my_model, 0)
get_gmus(my_model, 1)
get_gamma(my_model, 1, 2)

calc_ab(my_model, rnorm(3))
calc_gmus(my_model, rnorm(3), 0)
calc_gmus(my_model, rnorm(3), 10)
calc_gmus(my_model, rnorm(6), 2)

get_gamma(my_model, 1, 1)

full_lp(my_model)