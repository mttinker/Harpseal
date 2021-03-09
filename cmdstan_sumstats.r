sumstats = as.data.frame(fit$summary(variables = params))
row.names(sumstats) = sumstats$variable; sumstats = sumstats[,-1] 
tmp = as.data.frame(fit$summary(variables = params, mcse = mcse_mean, ~quantile(.x, probs = c(0.025, 0.975))))
sumstats$mcse = tmp$mcse; sumstats$q2.5 = tmp$`2.5%` ; sumstats$q97.5 = tmp$`97.5%`; 
sumstats$q50 = sumstats$median; sumstats$N_eff = sumstats$ess_bulk
col_order = c("mean", "mcse", "sd","q2.5","q5","q50","q95","q97.5","N_eff", "rhat")
sumstats = sumstats[, col_order]
mcmc = as_draws_matrix(fit$draws(variables = params))
vn = colnames(mcmc); vns = row.names(sumstats)
Nsims = nrow(mcmc)
rm(tmp,col_order)