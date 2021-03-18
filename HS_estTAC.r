# Script to plot some results from model fitting
# - K CI, should be CI of mean, not quantiles of abundance (bootstrap?)
# - fix TAC so not multiples of 5...
#
rm(list = ls())
# USER-SET PARAMETERS:
reps = 5000     # number of simulation reps to use (1000 for fast, 10000 for precise)
PAmeans = c(.18,.07,.75) # future expected proportion of pups in S Gulf, N Gulf, Front
# Specify year ranges to use to parameterize future conditions for ice and climate
Y1_ice = 2000    # first year of ice cover time series to use to parameterize future sims
Y2_ice = 2020    # last year of ice cover time series to use to parameterize future sims
Y1_CI = 1951     # first year of CI index time series to use to parameterize future sims
Y2_CI = 2020     # last year of CI index time series to use to parameterize future sims
# Set AVERAGE "other" human mortality for future sims (not included in Canadian harvest)
GrnH = 50000
ArcH = 1000
Byctc = 2000
# Specify adult/pup breakdown of "other" human mortality
Grn_A = .857*GrnH
Grn_P = .143*GrnH
Arc_A = .95*ArcH
Arc_P = .05*ArcH
Byc_A = .2*Byctc
Byc_P = .8*Byctc
#
# END USER PARAMETERS
#
# Load results file (if not already loaded into workspace):-------------------
load(file="Results.rdata")
#
library(foreach)
library(doParallel)
library(doRNG)
library(fitdistrplus)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(stats)
require(gtools)
require(betareg)
#
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
  psi1 = mcmc[rr,vn=="psi1"],
  psi2 = mcmc[rr,vn=="psi2"],
  dlta = mcmc[rr,startsWith(vn,"dlta[")],
  zeta = mcmc[rr,vn=="zeta"],
  gamma_H0_mn = mcmc[rr,vn=="gamma_H0_mn"],
  gamma_HA_mn = mcmc[rr,vn=="gamma_HA_mn"],
  gamma_H0 = mcmc[rr,startsWith(vn,"gamma_H0[")],
  gamma_HA = mcmc[rr,startsWith(vn,"gamma_HA[")]
)}
#
cores=detectCores()
cl <- makeCluster(min(20,cores[1]-1)) 
registerDoParallel(cl)
source("HSmod_sim.r")
# Part 1: Hindcast sims------------------------------------------------------------
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
r_vec = matrix(c(),)
PAr = rdirichlet(1000+Nyrs2, 25*PAmeans)
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
  ggtitle(titletxt,subtitle =subtxt) + theme_classic()
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
PAr = rdirichlet(1000+Nyrs2, 25*PAmeans)
epsFr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
epsSr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
ig1 = runif(reps2,.3,1.5)
ig2 = runif(reps2,.05,1.2)
#
# Future conditions: based on ice and CI indices from specified year range, 
#  fit appropriate random sampling distributions
ii = which(df.NLCI$Year>= max(min(df.NLCI$Year), Y1_CI) &  df.NLCI$Year <= min(max(df.NLCI$Year), Y2_CI))
ft = fitdist(df.NLCI$CI[ii],"norm")
CI2 = rnorm(1000,coef(ft)[1],coef(ft)[2])
ii = which(df.Ice$Year>= max(min(df.Ice$Year), Y1_ice) &  df.Ice$Year <= min(max(df.Ice$Year), Y2_ice))
ft = fitdist(exp(df.Ice$Gulf_Anom[ii]),"norm")
icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
ft = fitdist(exp(df.Ice$Lab_Anom[ii]),"norm")
icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
IC2 = as.matrix(cbind(icg1,icg1,icg2))
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
Ktab = data.frame(Metric = "Equilibrium abundance",
                  Mean = Kest/1000000, SD = Kest_sd/1000000, 
                  CI95_low=Kest_CI[1]/1000000,CI95_high=Kest_CI[2]/1000000)
print(Ktab)
# Plot
dpN_K = data.frame(Year = Yearp2, 
                    N_pred_mean = smooth.spline(apply(Np_Predict, 1, mean),spar=.5)$y /1000000,
                    N_pred_lo = smooth.spline(apply(Np_Predict,1,quantile,prob=0.1),spar=.5)$y /1000000,
                    N_pred_hi = smooth.spline(apply(Np_Predict,1,quantile,prob=0.9),spar=.5)$y /1000000 )
titletxt = "Model projected abundance with zero harvest (including pups)"
subtxt =  paste0("Estimated long-term equilibrium (K) = ", format(Kest/1000000,digits = 3),
                 "million (CI95 = ", format(Kest_CI[1]/1000000,digits = 3),
                 " - ", format(Kest_CI[2]/1000000,digits = 3),")")
plt_N_Kest = ggplot(data=dpN_K,aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Projected abundance (millions)") +
  geom_hline(yintercept=Kest/1000000, linetype="dashed", color = "red") +
  ggtitle(titletxt,subtitle =subtxt) + theme_classic()
print(plt_N_Kest)
# 
# Part 3: TAC-est sims-------------------------------------------------------------
#
# Set up sims for 15 years under varying harvest levels
futuresim = 2 # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 15; Yearst2 = YearT+1 ; 
#
reps2 = reps*10
set.seed(123)
# r_vec = matrix(0,nrow = reps, ncol = 2)
r_vec = sample(nrow(mcmc),reps2,replace = T)
r_vec2 = sample(1000,reps2,replace = T)
PAr = rdirichlet(1000+Nyrs2, 25*PAmeans)
epsFr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
epsSr = matrix(rnorm(1000*Nyrs2,0,1),nrow=1000)
ig1 = runif(reps2,.1,1.5)
ig2 = runif(reps2,.01,1.2)
# source("HSmod_sim.r")
# Future conditions: based on ice and CI indices from specified year range, 
#  fit appropriate random sampling distributions
ii = which(df.NLCI$Year>= max(min(df.NLCI$Year), Y1_CI) &  df.NLCI$Year <= min(max(df.NLCI$Year), Y2_CI))
ft = fitdist(df.NLCI$CI[ii],"norm")
CI2 = rnorm(1000,coef(ft)[1],coef(ft)[2])
ii = which(df.Ice$Year>= max(min(df.Ice$Year), Y1_ice) &  df.Ice$Year <= min(max(df.Ice$Year), Y2_ice))
ft = fitdist(exp(df.Ice$Gulf_Anom[ii]),"norm")
icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
ft = fitdist(exp(df.Ice$Lab_Anom[ii]),"norm")
icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
IC2 = as.matrix(cbind(icg1,icg1,icg2))
# Update data
sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps2,r_vec=r_vec,r_vec2=r_vec2,
                 PAr=PAr,epsFr=epsFr,epsSr=epsSr,ig1=ig1,ig2=ig2,
                 PAmeans=PAmeans,futuresim=futuresim,
                 NCages=NCages,NCage1=NCage1,Nages=Nages,Nareas=Nareas,
                 ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CI=CI2,
                 omega=omega,Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
                 Byc_A=Byc_A,Byc_P=Byc_P) 
# Run sims:
rsltTAC=HSmod_sim(init_fun,sim.data,mcmc,vn)
# process results
N_ad = sumstats[which(startsWith(vns,"N[")),1]
N_all = N_ad[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,1]
N_end = N_all[Nyrs-1]
N70 = 0.7*max(N_all)
N50 = 0.5*max(N_all)
N70P = N70/8000000
N50P = N50/8000000

Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
N_Predict = rsltTAC$N_Predict
P_Predict = rsltTAC$P_Predict
P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
Np_Predict = (N_Predict + P_Predict)
Hvp = rsltTAC$HV0_predict
Hva = rsltTAC$HVA_predict
Hvp_k = apply(Hvp,2,mean)/1000 
Hva_k = apply(Hva,2,mean)/1000 
ii = which(Hvp_k>0 & Hva_k>0)
HvT = Hvp_k[ii]+Hva_k[ii]
ppA = Hva_k[ii]/HvT
# Nmin = apply(Np_Predict,2,quantile,prob=0.025); Nmin = Nmin[ii]
Nmin = apply(Np_Predict,2,min); Nmin = Nmin[ii]
# Convert to 0-1 for betareg
NminP = Nmin / 8000000
dat_N_Hv = data.frame(Nmin=NminP,HvT=HvT,HvT2 = HvT*HvT,ppA=ppA,ppA2=ppA*ppA,intx = HvT*ppA)
# GLM model - convert to 0-1 distribtion and use beta regression 
fit = betareg(Nmin~ HvT+HvT2+ppA+ppA2+intx,data=dat_N_Hv,link = "logit")
summary(fit)
#
HvTeval05 = seq(1,1200)
HvTeval10 = seq(1,1200)
HvTeval50 = seq(1,1200)
tmp = data.frame(HvT=HvTeval05, HvT2=HvTeval05^2, ppA=rep(.05,length(HvTeval05)));  
tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
newdat05ad = tmp; 
tmp = data.frame(HvT=HvTeval10, HvT2=HvTeval10^2, ppA=rep(.1,length(HvTeval10))); 
tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
newdat10ad =   tmp;
tmp = data.frame(HvT=HvTeval50, HvT2=HvTeval50^2, ppA=rep(.5,length(HvTeval50))); 
tmp$intx = tmp$HvT*tmp$ppA; tmp$ppA2 =  tmp$ppA*tmp$ppA
newdat50ad = tmp;

tmp = predict(fit, newdata = newdat05ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
pred05ad = data.frame(HvT=HvTeval05,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])

tmp = predict(fit, newdata = newdat10ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
pred10ad = data.frame(HvT=HvTeval10,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])

tmp = predict(fit, newdata = newdat50ad, type = "quantile", at = c(0.2, 0.5, 0.8)); tmp = tmp * 8000000; 
pred50ad = data.frame(HvT=HvTeval50,Fitted = tmp[,2],Lower = tmp[,1],Upper=tmp[,3])

ggplot(pred05ad,aes(x=HvT,y=Fitted)) +
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
  geom_line() +
  labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
       title="Harvest impacts on futue abundance, 5% adults") +
  geom_point(data=dat_N_Hv[which(dat_N_Hv$ppA>0.03 & dat_N_Hv$ppA<0.07),],
             aes(x=HvT,y=Nmin*8000000)) +
  xlim(0,500) + ylim(0,8000000) +
  theme_classic()

ggplot(pred10ad,aes(x=HvT,y=Fitted)) +
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
  geom_line() +
  labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
       title="Harvest impacts on futue abundance, 10% adults") +
  geom_point(data=dat_N_Hv[which(dat_N_Hv$ppA>0.08 & dat_N_Hv$ppA<0.12),],
             aes(x=HvT,y=Nmin*8000000)) +
  xlim(0,500) + ylim(0,8000000) +
  theme_classic()

ggplot(pred50ad,aes(x=HvT,y=Fitted)) +
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3) +
  geom_line() +
  labs(x="Harvest level",y="Minimum abundance, 15-yr projection",
       title="Harvest impacts on futue abundance, 50% adults") +
  geom_point(data=dat_N_Hv[which(dat_N_Hv$ppA>0.45 & dat_N_Hv$ppA<0.55),],
             aes(x=HvT,y=Nmin*8000000)) +
  xlim(0,500) + ylim(0,8000000) +
  theme_classic()

rndint = 1
ii = which(abs(pred05ad$Lower-(N70+1000))==min(abs(pred05ad$Lower-(N70+1000))))
TACN70pa05 = floor(HvTeval05[ii]/rndint)*rndint  
ii = which(abs(pred05ad$Lower-(N50+1000))==min(abs(pred05ad$Lower-(N50+1000))))
TACN50pa05 = floor(HvTeval05[ii]/rndint)*rndint 
ii = which(abs(pred10ad$Lower-(N70+1000))==min(abs(pred10ad$Lower-(N70+1000))))
TACN70pa10 = floor(HvTeval10[ii]/rndint)*rndint  
ii = which(abs(pred10ad$Lower-(N50+1000))==min(abs(pred10ad$Lower-(N50+1000))))
TACN50pa10 = floor(HvTeval10[ii]/rndint)*rndint  
ii = which(abs(pred50ad$Lower-(N70+1000))==min(abs(pred50ad$Lower-(N70+1000))))
TACN70pa50 = floor(HvTeval50[ii]/rndint)*rndint  
TACN70pa50 = min(TACN70pa50,TACN70pa10-5)
ii = which(abs(pred50ad$Lower-(N50+1000))==min(abs(pred50ad$Lower-(N50+1000))))
TACN50pa50 = floor(HvTeval50[ii]/rndint)*rndint  
TACN50pa50 = min(TACN50pa50,TACN50pa10-5)
# 
# Create table of TAC recomendations
TACtab = data.frame(Pcnt_adlt=c(5,10,50),
                    TAC_N70 = c(TACN70pa05,TACN70pa10,TACN70pa50),
                    TAC_N50 =c(TACN50pa05,TACN50pa10,TACN50pa50))
print(TACtab)

# Save results ---------------------------------------------------------------
save(file = "TAC_K_est.rdata",Kest,Kest_CI,SAD0,
     TACtab,Ktab,dpN_HCst,dpN_K,dat_N_Hv,plt_N_hindcast,plt_N_Kest)
#
stopCluster(cl)
