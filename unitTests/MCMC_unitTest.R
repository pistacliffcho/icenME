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


# Make ab sampler
ab_smp = make_ab_sampler(my_model)
ab_smp$getVals(ab_smp$info)
ab_smp$get_lp(ab_smp$getVals(ab_smp$info), ab_smp$info)
ab_smp$getVals(ab_smp$info)
ab_smp$update()
ab_smp$getVals(ab_smp$info)
ab_smp$get_lp(ab_smp$getVals(ab_smp$info), ab_smp$info)

flp1 = full_lp(my_model)
prtlp1 = ab_smp$get_lp(ab_smp$getVals(ab_smp$info), ab_smp$info)
ab_smp$update()
flp2 = full_lp(my_model)
prtlp2 = ab_smp$get_lp(ab_smp$getVals(ab_smp$info), ab_smp$info)

cat("Full lp diff = ", flp2 - flp1, 
    "part lp diff = ", prtlp2 -prtlp1)


# Make gmus sampler
gmus_smp = make_gmuS_sampler(my_model, 1)
cur_vals = gmus_smp$getVals(gmus_smp$info)
gmus_smp$get_lp(cur_vals, gmus_smp$info)
gmus_smp$get_lp(cur_vals + 1, gmus_smp$info)
gmus_smp$getVals(gmus_smp$info)
gmus_smp$update()
gmus_smp$getVals(gmus_smp$info)
gmus_smp$get_lp(cur_vals, gmus_smp$info)


flp1 = full_lp(my_model)
prtlp1 = gmus_smp$get_lp(gmus_smp$getVals(gmus_smp$info), gmus_smp$info)
gmus_smp$update()
flp2 = full_lp(my_model)
prtlp2 = gmus_smp$get_lp(gmus_smp$getVals(gmus_smp$info), gmus_smp$info)

cat("Full lp diff = ", flp2 - flp1, 
    "part lp diff = ", prtlp2 -prtlp1)



# Make gamma sampler
gamma_smp = make_gamma_sampler(my_model, 1, 2)
cur_vals = gamma_smp$getVals(gamma_smp$info)
gamma_smp$get_lp(cur_vals, gamma_smp$info)
gamma_smp$get_lp(cur_vals + 1, gamma_smp$info)
gamma_smp$setVals(cur_vals + 1, gamma_smp$info)
gamma_smp$getVals(gamma_smp$info)
gamma_smp$update()
gamma_smp$getVals(gamma_smp$info)
gamma_smp$get_lp(gamma_smp$getVals(gamma_smp$info), gamma_smp$info)

flp1 = full_lp(my_model)
prtlp1 = gamma_smp$get_lp(gamma_smp$getVals(gamma_smp$info), gamma_smp$info)
gamma_smp$update()
flp2 = full_lp(my_model)
prtlp2 = gamma_smp$get_lp(gamma_smp$getVals(gamma_smp$info), gamma_smp$info)

cat("Full lp diff = ", flp2 - flp1, 
    "part lp diff = ", prtlp2 -prtlp1)

# Make mcmc

mcmc <- make_MCMC(my_model)

mv = mcmc$run(MC =10000)
plot(mv$lp, type = 'l')
plot(mv$alphaBeta[,1], type = 'l')
plot(mv$alphaBeta[,2], type = 'l')
plot(mv$alphaBeta[,3], type = 'l')
plot(mv$gmus$`1`[,1], type = 'l')
plot(mv$gmus$`2`[,2], type = 'l')
plot(mv$gamma$`1`$`1`[,1], type = 'l')
