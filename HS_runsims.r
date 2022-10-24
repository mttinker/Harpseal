# Script to plot some results from model fitting
#
# USER-SET PARAMETERS---------------------------------------------------------------
reps = 10000     # number of simulation reps to use (1000 for fast, 10000 for precise)
PAmeans = c(.18,.07,.75) # future expected proportion of pups in S Gulf, N Gulf, Front
# Specify year ranges to use to parameterize future conditions for ice and climate
Y1_ice = 2000    # first year of ice cover time series to use to parameterize future sims
Y2_ice = 2020    # last year of ice cover time series to use to parameterize future sims
Y1_CI = 2000     # first year of CI index time series to use to parameterize future sims
Y2_CI = 2020     # last year of CI index time series to use to parameterize future sims
# Set AVERAGE "other" human mortality for future sims (not included in Canadian harvest)
GrnH = 50000
ArcH = 1000
Byctc = 2000
# Specify adult/pup ratios for "other" human mortality sources
Grn_A = .857  # Greenland, proportion adults
Grn_P = .143  # Greenland, proportion pups
Arc_A = .95   # Arctic, proportion adults
Arc_P = .05   # Arctic, proportion pups
Byc_A = .2    # Bycatch, proportion adults
Byc_P = .8    # Bycatch, proportion pups
#
# END USER PARAMETERS --------------------------------------------------------
#
# Load libraries
#
library(parallel)
library(foreach)
library(doParallel)
require(future) 
library(doRNG)
library(fitdistrplus)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(cowplot)
require(dplyr)
require(stats)
require(gtools)
require(betareg)
require(svDialogs)
require(rJava)
require(rChoiceDialogs)
#
if( !exists("resultsfile") ){
  # resultsfile = file.choose(new = FALSE)
  stop_quietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  file_list = list.files(path = "./Results",pattern = "HS_Results",full.names = F)
  rslt_list = grep("TAC", file_list, value = TRUE, invert = TRUE)
  rslt_list = grep("Report", rslt_list, value = TRUE, invert = TRUE)
  rdata_file = rselect.list(rslt_list, preselect = NULL, multiple = FALSE,
                            title = "Select results file" ,
                            graphics = getOption("menu.graphics")) 
  if(length(rdata_file)==0){
    dlg_message(c("No data file selected"), "ok")
    stop_quietly()
  }else{
    load(paste0("./Results/",rdata_file))
  }
}
saveSIMname = paste0("SIMS_",resultsfile)
# Function to draw from joint posteriors
init_fun <- function(rr) {list(
  sigF = mcmc[rr,vn=="sigF"],
  sigS = mcmc[rr,vn=="sigS"],
  sigH = mcmc[rr,vn=="sigH"],
  gamma_0 = mcmc[rr,vn=="gamma_0"],
  gamma_A = mcmc[rr,startsWith(vn,"gamma_A[")],
  phiF = mcmc[rr,vn=="phi[1]"],
  phiS = mcmc[rr,vn=="phi[2]"],
  beta1 = mcmc[rr,vn=="beta1"],
  beta2 = mcmc[rr,vn=="beta2"],
  thtaF = mcmc[rr,vn=="thta[1]"],
  thtaS = mcmc[rr,vn=="thta[2]"],
  psi1Gu = mcmc[rr,vn=="psi1[1]"],
  psi1Fr = mcmc[rr,vn=="psi1[2]"],
  psi2 = mcmc[rr,vn=="psi2"],
  dlta = mcmc[rr,startsWith(vn,"dlta[")],
  zeta = mcmc[rr,vn=="zeta"],
  gamma_H0_mn = mcmc[rr,vn=="gamma_H0_mn"],
  gamma_HA_mn = mcmc[rr,vn=="gamma_HA_mn"],
  gamma_H0 = mcmc[rr,startsWith(vn,"gamma_H0[")],
  gamma_HA = mcmc[rr,startsWith(vn,"gamma_HA[")]
)}
#

ncores = min(40,availableCores()-4)
cl = parallelly::makeClusterPSOCK(ncores, autoStop = F)
registerDoParallel(cl,cores=ncores)

source("HSmod_sim.r")
# Part 1: Hindcast sims------------------------------------------------------------
#
Grn_A = Grn_A*GrnH
Grn_P = Grn_P*GrnH
Arc_A = Arc_A*ArcH
Arc_P = Arc_P*ArcH
Byc_A = Byc_A*Byctc
Byc_P = Byc_P*Byctc
#
futuresim = 0
Yearst2=Year1 
Nyrs2=Nyrs
N_end = sumstats[vns=="N0",1]
IC2=IC
CI2=CI
# Create some random variables
set.seed(123)
# r_vec = matrix(0,nrow = reps, ncol = 2)
r_vec = sample(nrow(mcmc),reps,replace = T)
r_vec2 = sample(1000,reps,replace = T)
PAr = rdirichlet(1000+Nyrs2, 50*PAmeans)
epsFr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
epsSr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
ig1 = runif(reps,.3,1.5)
ig2 = runif(reps,.05,1.2)
#
sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps,r_vec=r_vec,r_vec2=r_vec2,
                 PAr=PAr,epsFr=epsFr,epsSr=epsSr,ig1=ig1,ig2=ig2,
                 PAmeans=PAmeans,futuresim=futuresim,
                 NCages=NCages,NCage1=NCage1,Nages=Nages,Nareas=Nareas,
                 ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CI=CI2,
                 omega=omega,Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
                 Byc_A=Byc_A,Byc_P=Byc_P) 
# Run sims:
rslt=HSmod_sim(init_fun,sim.data,mcmc,vn)
# Process results:
Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
N_Predict = rslt$N_Predict
P_Predict = rslt$P_Predict
P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
Np_Predict = (N_Predict + P_Predict)
dpN_HCst = data.frame(Year = Yearp2, N_pred_mean = rowMeans(Np_Predict),
                    N_pred_lo = apply(Np_Predict,1,quantile,prob=0.025),
                    N_pred_hi = apply(Np_Predict,1,quantile,prob=0.975))
Predpup = sumstats[startsWith(vns,"Pups_pred[")==T,1]; Predpup = c(Predpup,Predpup[Nyrs-1])
dpN_HCst$N_actual = sumstats[startsWith(vns,"N[")==T,1] + Predpup
titletxt = "Hindcast stochastic simulations of abundance with observed harvest levels"
subtxt = " (red line = actual model estimate of abundance) "
plt_N_hindcast = ggplot(data=dpN_HCst,aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Projected abundance") +
  geom_line(aes(y=N_actual),col="red") +
  # ggtitle(titletxt,subtitle =subtxt) + 
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(plt_N_hindcast)
SAD0 = as.data.frame(rslt$SAD)

# Part 2: K-est sims--------------------------------------------------------------
#
futuresim = 1; # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 100; Yearst2 = YearT+1 ; 
N_end = sumstats[which(vns==paste0("N[",Nyrs,"]")),1]  
reps2 = reps*5
#
set.seed(123)
# r_vec = matrix(0,nrow = reps, ncol = 2)
r_vec = sample(nrow(mcmc),reps2,replace = T)
r_vec2 = sample(1000,reps2,replace = T)
PAr = rdirichlet(1000+Nyrs2, 50*PAmeans)
epsFr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
epsSr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
ig1 = runif(reps2,.3,1.5)
ig2 = runif(reps2,.05,1.2)
#
for (sc in 1:2){
  # Future conditions: based on ice and CI indices from specified year range, 
  if (sc==1){
    #  Pre-2000 conditions
    ii = which(df.NLCI$Year>= min(df.NLCI$Year) &  df.NLCI$Year < min(max(df.NLCI$Year), Y1_CI))
    ft = fitdist(exp(df.NLCI$CI[ii]),"lnorm")
    CI2 = pmax(-1.5,pmin(1.5,log(rlnorm(1000,coef(ft)[1],coef(ft)[2]))))
    # mean(CI2)
    ii = which(df.Ice$Year>= min(df.Ice$Year) &  df.Ice$Year < min(max(df.Ice$Year), Y1_ice))
    ft = fitdist(exp(Gulf_Anom[ii]),"norm")
    icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
    ft = fitdist(exp(Front_Anom[ii]),"norm")
    icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
    IC2 = as.matrix(cbind(icg1,icg1,icg2))    
    # colMeans(IC2)
  }else{
    #  Post-2000 conditions
    ii = which(df.NLCI$Year>= max(min(df.NLCI$Year), Y1_CI) &  df.NLCI$Year <= min(max(df.NLCI$Year), Y2_CI))
    ft = fitdist(exp(df.NLCI$CI[ii]),"lnorm")
    CI2 = pmax(-1.5,pmin(1.5,log(rlnorm(1000,coef(ft)[1],coef(ft)[2]))))
    # mean(CI2)
    ii = which(df.Ice$Year>= max(min(df.Ice$Year), Y1_ice) &  df.Ice$Year <= min(max(df.Ice$Year), Y2_ice))
    ft = fitdist(exp(Gulf_Anom[ii]),"norm")
    icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
    ft = fitdist(exp(Front_Anom[ii]),"norm")
    icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
    IC2 = as.matrix(cbind(icg1,icg1,icg2))
    # colMeans(IC2)
  }
  # Update data
  sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps2,r_vec=r_vec,r_vec2=r_vec2,
                   PAr=PAr,epsFr=epsFr,epsSr=epsSr,ig1=ig1,ig2=ig2,
                   PAmeans=PAmeans,futuresim=futuresim,
                   NCages=NCages,NCage1=NCage1,Nages=Nages,Nareas=Nareas,
                   ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CI=CI2,
                   omega=omega,Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
                   Byc_A=Byc_A,Byc_P=Byc_P) 
  # Run sims:
  rsltK=HSmod_sim(init_fun,sim.data,mcmc,vn)
  # process results
  Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
  N_Predict = rsltK$N_Predict
  P_Predict = rsltK$P_Predict
  P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
  Np_Predict = (N_Predict + P_Predict)
  Nfin = apply(Np_Predict[round(Nyrs2/2):Nyrs2,],1, mean)
  # stats = as.numeric(exp(quantile(log(Nfin),probs=c(0.5,0.025,0.975))))
  Kest = mean(Nfin)
  Kest_sd = sd(Nfin)
  Kest_CI = quantile(Nfin,probs=c(0.025,0.975))
  if (sc==1){
    Ktab = data.frame(Metric = "Equilibrium abundance, pre-2000 conditions",
                      Mean = Kest/1000000, SD = Kest_sd/1000000, 
                      CI95_low=Kest_CI[1]/1000000,CI95_high=Kest_CI[2]/1000000)
    dpN_K = data.frame(Scenario = rep("Pre-2000 conditions",length(Yearp2)), Year = Yearp2, 
                       N_pred_mean = smooth.spline(apply(Np_Predict, 1, mean),spar=.1)$y /1000000,
                       N_pred_lo = smooth.spline(apply(Np_Predict,1,quantile,prob=0.1),spar=.1)$y /1000000,
                       N_pred_hi = smooth.spline(apply(Np_Predict,1,quantile,prob=0.9),spar=.1)$y /1000000 )
  }else{
    Ktab = rbind(Ktab, data.frame(Metric = "Equilibrium abundance, post-2000 conditions",
                                  Mean = Kest/1000000, SD = Kest_sd/1000000, 
                                  CI95_low=Kest_CI[1]/1000000,CI95_high=Kest_CI[2]/1000000) )
    dpN_K = rbind(dpN_K, data.frame(Scenario = rep("Post-2000 conditions", length(Yearp2)), Year = Yearp2, 
                                    N_pred_mean = smooth.spline(apply(Np_Predict, 1, mean),spar=.1)$y /1000000,
                                    N_pred_lo = smooth.spline(apply(Np_Predict,1,quantile,prob=0.1),spar=.1)$y /1000000,
                                    N_pred_hi = smooth.spline(apply(Np_Predict,1,quantile,prob=0.9),spar=.1)$y /1000000 ) )
  }
}
# Plot
titletxt1 = "A): Estimted K, pre-2000 conditions"
subtxt =  paste0("Estimated long-term equilibrium (K) = ", format(Ktab$Mean[1],digits = 3),
                 "million (CI95 = ", format(Ktab$CI95_low[1],digits = 3),
                 " - ", format(Ktab$CI95_high[1],digits = 3), ")")
ii = which(dpN_K$Scenario=="Pre-2000 conditions")
plt_N_Kest1 = ggplot(data=dpN_K[ii,],
                     aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + 
  labs(x = "Year",
       y = expression(paste("Projected abundance (millions) \n pre-2000 environmental conditions"))) +
  geom_hline(yintercept=Ktab$Mean[1], linetype="dashed", color = "red") +
  ylim(4.5,12) +
  #ggtitle(titletxt1) + 
  theme_classic() + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
# print(plt_N_Kest1)
titletxt2 = "B): Estimted K, post-2000 conditions"
subtxt =  paste0("Estimated long-term equilibrium (K) = ", format(Ktab$Mean[2],digits = 3),
                 "million (CI95 = ", format(Ktab$CI95_low[2],digits = 3),
                 " - ", format(Ktab$CI95_high[2],digits = 3), ")")
ii = which(dpN_K$Scenario=="Post-2000 conditions")
plt_N_Kest2 = ggplot(data=dpN_K[ii,],
                     aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + 
  labs(x = "Year",
       y = expression(paste("Projected abundance (millions) \n post-2000 environmental conditions"))) +
  geom_hline(yintercept=Ktab$Mean[2], linetype="dashed", color = "red") +
  ylim(4.5,12) +
  #ggtitle(titletxt2) + 
  theme_classic() + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
# print(plt_N_Kest2)
plot_grid(plt_N_Kest1,plt_N_Kest2,nrow = 2,labels=c("A","B")) #,labels = c(titletxt1,titletxt2))

# 
# Part 3: TAC-est sims-------------------------------------------------------------
#
# Set up sims for 15 years under varying harvest levels
# futuresim = 2 # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
# Nyrs2 = 15; Yearst2 = YearT+1 ; 
# N_end = sumstats[which(vns==paste0("N[",Nyrs,"]")),1]  
# #
# reps2 = reps*10
# set.seed(123)
# # r_vec = matrix(0,nrow = reps, ncol = 2)
# r_vec = sample(nrow(mcmc),reps2,replace = T)
# r_vec2 = sample(1000,reps2,replace = T)
# PAr = rdirichlet(1000+Nyrs2, 25*PAmeans)
# epsFr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
# epsSr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
# ig1 = runif(reps2,.01,1.2)
# ig2 = runif(reps2,.01,1.5)
# # source("HSmod_sim.r")
# # Future conditions: based on ice and CI indices from specified year range, 
# #  fit appropriate random sampling distributions
# ii = which(df.NLCI$Year>= max(min(df.NLCI$Year), Y1_CI) &  df.NLCI$Year <= min(max(df.NLCI$Year), Y2_CI))
# ft = fitdist(exp(df.NLCI$CI[ii]),"lnorm")
# CI2 = pmax(-1.5,pmin(1.5,log(rlnorm(1000,coef(ft)[1],coef(ft)[2]))))
# ii = which(df.Ice$Year>= max(min(df.Ice$Year), Y1_ice) &  df.Ice$Year <= min(max(df.Ice$Year), Y2_ice))
# ft = fitdist(exp(df.Ice$Gulf_Anom[ii]),"norm")
# icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
# ft = fitdist(exp(df.Ice$Front_Anom[ii]),"norm")
# icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
# IC2 = as.matrix(cbind(icg1,icg1,icg2))
# # Update data
# sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps2,r_vec=r_vec,r_vec2=r_vec2,
#                  PAr=PAr,epsFr=epsFr,epsSr=epsSr,ig1=ig1,ig2=ig2,
#                  PAmeans=PAmeans,futuresim=futuresim,
#                  NCages=NCages,NCage1=NCage1,Nages=Nages,Nareas=Nareas,
#                  ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CI=CI2,
#                  omega=omega,Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
#                  Byc_A=Byc_A,Byc_P=Byc_P) 
# # Run sims:
# rsltTAC=HSmod_sim(init_fun,sim.data,mcmc,vn)
# # process results
# N_ad = sumstats[which(startsWith(vns,"N[")),1]
# N_all = N_ad[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,1]
# N_end = N_all[Nyrs-1]
# N70 = 0.7*max(N_all)
# N50 = 0.5*max(N_all)
# N70P = N70/8000000
# N50P = N50/8000000
# 
# Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
# N_Predict = rsltTAC$N_Predict
# P_Predict = rsltTAC$P_Predict
# P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
# Np_Predict = (N_Predict + P_Predict)
# Hvp = rsltTAC$HV0_predict
# Hva = rsltTAC$HVA_predict
# Hvp_k = apply(Hvp,2,mean)/1000 
# Hva_k = apply(Hva,2,mean)/1000 
# HvT = Hvp_k + Hva_k
# ii = which(Hvp_k > 0 & Hva_k > 0 & HvT<500)
# Nmin = apply(Np_Predict,2,min); Nmin = Nmin[ii]
# HvT = Hvp_k[ii]+Hva_k[ii]
# ppA = Hva_k[ii]/HvT
# # Nmin = apply(Np_Predict,2,quantile,prob=0.05); Nmin = Nmin[ii]
# 
# # Convert to 0-1 for betareg
# NminP = Nmin / 8000000
# dat_N_Hv = data.frame(Nmin=NminP,HvT=HvT,HvT2 = HvT*HvT,ppA=ppA,ppA2=ppA*ppA,intx = HvT*ppA)
# # GLM model - convert to 0-1 distribtion and use beta regression 
# #fit = betareg(NminP ~ HvT+ppA+intx,data=dat_N_Hv,link = "logit")
# #summary(fit)
# dat_N_Hv05 = dat_N_Hv[which(dat_N_Hv$ppA>0.02 & dat_N_Hv$ppA<0.08),]
# dat_N_Hv10 = dat_N_Hv[which(dat_N_Hv$ppA>0.07 & dat_N_Hv$ppA<0.13),]
# dat_N_Hv50 = dat_N_Hv[which(dat_N_Hv$ppA>0.45 & dat_N_Hv$ppA<0.55),]
# #
# fit05 = betareg(Nmin ~ HvT+HvT2,data=dat_N_Hv05,link = "logit")
# summary(fit05)
# fit10 = betareg(Nmin ~ HvT+HvT2,data=dat_N_Hv10,link = "logit")
# summary(fit10)
# fit50 = betareg(Nmin ~ HvT+HvT2,data=dat_N_Hv50,link = "logit")
# summary(fit50)
# #
# HvTeval05 = seq(1,600)
# HvTeval10 = seq(1,600)
# HvTeval50 = seq(1,600)
# tmp = data.frame(HvT=HvTeval05, HvT2=HvTeval05^2, ppA=rep(.05,length(HvTeval05)));  
# tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
# newdat05ad = tmp; 
# tmp = data.frame(HvT=HvTeval10, HvT2=HvTeval10^2, ppA=rep(.1,length(HvTeval10))); 
# tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
# newdat10ad =   tmp;
# tmp = data.frame(HvT=HvTeval50, HvT2=HvTeval50^2, ppA=rep(.5,length(HvTeval50))); 
# tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
# newdat50ad = tmp;
# #
# tmp = predict(fit05, newdata = newdat05ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
# pred05ad = data.frame(HvT=HvTeval05,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])
# tmp = predict(fit10, newdata = newdat10ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
# pred10ad = data.frame(HvT=HvTeval10,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])
# tmp = predict(fit50, newdata = newdat50ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
# pred50ad = data.frame(HvT=HvTeval50,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])
# #
# ggplot(pred05ad,aes(x=HvT,y=Fitted)) +
#   geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
#   geom_line() +
#   labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
#        title="Harvest impacts on futue abundance, 5% adults") +
#   geom_point(data=dat_N_Hv05,
#              aes(x=HvT,y=Nmin*8000000)) +
#   geom_hline(yintercept=N70, linetype="dashed", color = "red") +
#   geom_hline(yintercept=N50, linetype="dashed", color = "purple") +
#   xlim(0,500) + ylim(0,8000000) +
#   theme_classic()
# #
# ggplot(pred10ad,aes(x=HvT,y=Fitted)) +
#   geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
#   geom_line() +
#   labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
#        title="Harvest impacts on futue abundance, 10% adults") +
#   geom_point(data=dat_N_Hv10,
#              aes(x=HvT,y=Nmin*8000000)) +
#   geom_hline(yintercept=N70, linetype="dashed", color = "red") +
#   geom_hline(yintercept=N50, linetype="dashed", color = "purple") +  
#   xlim(0,400) + ylim(0,8000000) +
#   theme_classic()
# #
# ggplot(pred50ad,aes(x=HvT,y=Fitted)) +
#   geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
#   geom_line() +
#   labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
#        title="Harvest impacts on futue abundance, 50% adults") +
#   geom_point(data=dat_N_Hv50,
#              aes(x=HvT,y=Nmin*8000000)) +
#   geom_hline(yintercept=N70, linetype="dashed", color = "red") +
#   geom_hline(yintercept=N50, linetype="dashed", color = "purple") +  
#   xlim(0,300) + ylim(0,8000000) +
#   theme_classic()
# #
# rndint = 1
# ii = which(abs(pred05ad$Lower-(N70+1000))==min(abs(pred05ad$Lower-(N70+1000))))
# TACN70pa05 = floor(HvTeval05[ii]/rndint)*rndint  
# TACN70pa05 = max(5,TACN70pa05)
# ii = which(abs(pred05ad$Lower-(N50+1000))==min(abs(pred05ad$Lower-(N50+1000))))
# TACN50pa05 = floor(HvTeval05[ii]/rndint)*rndint 
# TACN50pa05 = max(5,TACN50pa05)
# ii = which(abs(pred10ad$Lower-(N70+1000))==min(abs(pred10ad$Lower-(N70+1000))))
# TACN70pa10 = floor(HvTeval10[ii]/rndint)*rndint  
# TACN70pa10 = max(2, min(TACN70pa10,TACN70pa05-1))
# ii = which(abs(pred10ad$Lower-(N50+1000))==min(abs(pred10ad$Lower-(N50+1000))))
# TACN50pa10 = floor(HvTeval10[ii]/rndint)*rndint  
# TACN50pa10 = max(2, min(TACN50pa10,TACN50pa05-1))
# ii = which(abs(pred50ad$Lower-(N70+1000))==min(abs(pred50ad$Lower-(N70+1000))))
# TACN70pa50 = floor(HvTeval50[ii]/rndint)*rndint  
# TACN70pa50 = max(1, min(TACN70pa50,TACN70pa10-1))
# ii = which(abs(pred50ad$Lower-(N50+1000))==min(abs(pred50ad$Lower-(N50+1000))))
# TACN50pa50 = floor(HvTeval50[ii]/rndint)*rndint  
# TACN50pa50 = max(1, min(TACN50pa50,TACN50pa10-1))
# # 
# # Create table of TAC recomendations
# TACtab = data.frame(Pcnt_adlt=c(5,10,50),
#                     TAC_N70 = c(TACN70pa05,TACN70pa10,TACN70pa50),
#                     TAC_N50 =c(TACN50pa05,TACN50pa10,TACN50pa50))
# print(TACtab)
# #
# Save results ---------------------------------------------------------------
rspns = dlg_message(c("Do you wish to save sim results? (This will over-write any",
                      "existing sim results asociated with this model results file)"), "yesno")$res
if(rspns=="yes"){
  save(file=paste0("./Results/",saveSIMname), 
     Ktab,dpN_HCst,dpN_K,plt_N_hindcast,plt_N_Kest1,plt_N_Kest2)
}
#
tryCatch(stopCluster(cl), error = function(e1) e1 = print("Cluster terminted"))
tryCatch(rm(cl), error = function(e1) e1 = print("Cluster removed"),
         warning = function(e1) e1 = print("Cluster removed"))
