sumstats = fit$draws(parms) %>%
  summarise_draws(mean, mcse = mcse_mean, sd, 
                  ~quantile(.x, probs = c(0.025, 0.05, .5, .95, .975)),
                  N_eff = ess_bulk, rhat)
sumstats = as.data.frame(sumstats)
row.names(sumstats) = sumstats$variable; sumstats = sumstats[,-1] 
colnames(sumstats) = c("mean", "mcse", "sd","q2.5","q5","q50","q95","q97.5","N_eff", "rhat")
mcmc = as_draws_matrix(fit$draws(variables = parms))
Nsims = nrow(mcmc)
mcmc_array = as_draws_array(fit$draws(variables = parms))
vn = colnames(mcmc); vns = row.names(sumstats)
