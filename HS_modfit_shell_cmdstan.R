# Shell script for fitting harp seal model, Cmdstan version
#   NOTE: 3 pupping areas: 1) S.GSL, 2) N.GSL, 3) Front
#       - for ice mortality, allow some pup deaths to occur before counts?
#
# User-set parameters--------------------------------------------------------
#
# Specify year range of time series for model fitting 
Year1 = 1951   #  Year1 = first year that harvest data available
YearT = 2020   #  YearT = year AFTER last available data for fitting
Nareas = 3     # Number pupping areas (assume 3: S.GSL, N.GSL, Front)  
# Prior for Initial Population size year 1 (model uses weak prior):
N0pri = 2400000 # prior guess of starting pop size, default ~ 2.2 million
# Assumed CV for harvest/bycatch totals:
CV_HV = 0.1 
# Youngest adult age class to consider for age composition fitting
NCage1 = 2; # Youngest age of adults for age composition fitting
# Suggested Prior means 
# Ice anomaly effects on pup survival: export opinion for prior
psi1 = -2    # prior for psi param1 of ice anomaly effect fxn (see figure)
psi2 = 5       # prior for psi param2 of ice anomaly effect fxn 
psi1_sd = 1    # stddev of psi1 prior mean (determines strength of prior)
psi2_sd = 1    # stddev of psi2 prior mean (determines strength of prior)
#
nburnin = 1000
nsamples = 10000
#
# End user params -----------------------------------------------------------
#
# Ice anomaly mortality fxn params (can explore graphically):
psipri1a = psi1   # prior for psi param1 mean: ice anomaly effect fxn
psipri1b = psi1_sd   # prior for psi param1, sd: ice anomaly effect fxn 
psipri2a = psi2 # prior for psi param2 mean: ice anomaly effect fxn 1.5
psipri2b = psi2_sd # prior for psi param2, sd: ice anomaly effect fxn   
source("Ice_mort_plot.r")
Ice_mort_plot(psipri1a,psipri1b,psipri2a,psipri2b)
#
# Load libraries-------------------------------------------------------------
require(parallel)
require(gtools)
#require(lattice)
#require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
require(stats)
require(bayesplot)
require(fitdistrplus)
require(readxl)
# require(loo)
# Load data ----------------------------------------------------------------
source("Importdat.r")
# Process data--------------------------------------------------------------
Yr0 = Year1-1;   
Nyrs = YearT-Yr0 
omega = -7
Nages = 36
ages = seq(1,Nages); ages2 = ages^2
agesC = ages - 10; agesC[Nages] = agesC[Nages] + 4; agesC2 = agesC^2
NCages = length(seq(NCage1,Nages))
NCobs = nrow(AgeComp)
Agects = AgeComp[,NCage1:Nages]
NPRobs = nrow(df.Rep)
NPcts = nrow(df.Pup)
NPctsA = nrow(df.Pup[!is.na(df.Pup$Sgulf_N),])
YrAGsmp = Years_ages - Yr0
NFsamp = df.Rep$N
NPrg = df.Rep$Preg
YrPRsmp = df.Rep$Year - Yr0
AgePRsmp = df.Rep$Age
sad0 = df.sad$SAD
Pups = numeric()
sdNP = numeric()
YrPct = numeric()
PAtag = numeric(length = Nyrs-1); 
NPA = matrix(0,nrow = Nyrs-1,ncol = Nareas)
aa = c(2,4,6)
for (i in 1:nrow(df.Pup)){
  YrPct[i] = df.Pup$Year[i] - Yr0
  Pups[i] = df.Pup$Total_Npup[i]
  sdNP[i] = min(df.Pup$Total_se[i],Pups[i]*.2)
  if (!is.na(df.Pup$Sgulf_N[i])){
    NPA[YrPct[i],1:Nareas] = as.numeric(df.Pup[i,aa]/sum(df.Pup[i,aa]) )
    PAtag[YrPct[i]] = 1
  }
}
PAidx = which(PAtag==1)
NPA = NPA[PAidx,]
NPActs = length(PAidx)
# 
# Set priors for beta distributions of proportion pups in areas 1 & 2
#   proportion in area 3, PA3, is then be calculated as 1 - (PA1+PA2) 
PApri = matrix(nrow = Nareas-1,ncol = 2)
ft = fitdist(NPA[,1],"beta")
PApri[1,1] = as.numeric(ft$estimate[1]); PApri[1,2] = as.numeric(ft$estimate[2]); 
ft = fitdist(NPA[,2],"beta")
PApri[2,1] = as.numeric(ft$estimate[1]); PApri[2,2] = as.numeric(ft$estimate[2]); 
# Create data inputs for observed havest, environmental CE index and ice anomalies
IC = matrix(0,nrow = Nyrs, ncol = 3)
CE = numeric(length = Nyrs)
H_0 = numeric(length = (Nyrs-1))
H_A = numeric(length = (Nyrs-1))
Q_0 = numeric(length = (Nyrs-1))
Q_A = numeric(length = (Nyrs-1))
for (y in 1:Nyrs){
  ii = which(df.Ice$Year==Yr0+y)
  if(length(ii)==1){
    IC[y,1] = df.Ice$Gulf_Anom[ii]
    IC[y,2] = df.Ice$Gulf_Anom[ii]
    IC[y,3] = df.Ice$Lab_Anom[ii]
  }
  ii = which(df.CE$Year==Yr0+y)
  if(length(ii)==1){
    CE[y] = log(df.CE$CEindex[ii])
  }
  ii = which(df.HV$YEAR==Yr0+y)
  if(length(ii)==1 & y<Nyrs){
    H_0[y] = df.HV$PUPTOT[ii]
    H_A[y] = df.HV$ADLTOT[ii]
    Q_A[y] = df.HV$ADL_prob_rec[ii]
    Q_0[y] = df.HV$PUP_prob_rec[ii]
  }
}
# Create a "dummy vector" of ice anomaly values from -1 to 1 (step size 0.05), 
#  and discretized version of IC anomaly values clasified to units of 0.05 
#  (to speed up fitting)
ICvec = seq(-1,1,by=.05)
ICcat = matrix(0,nrow = Nyrs, ncol = 3) 
for (i in 1:(3*Nyrs)){
  ICcat[i] = which(abs(IC[i] - ICvec) == min(abs(IC[i]- ICvec)))
}
#
rm(i,ii,y,aa,ft) 
#
# Set up Jags inputs --------------------------------------------------------
#
fitmodel = c("HSmodfit.stan")
#  
stan.data <- list(NPcts=NPcts,NPctsA=NPctsA,NCages=NCages,NCage1=NCage1,NPActs=NPActs,
                  NCobs=NCobs,NPRobs=NPRobs,Nyrs=Nyrs,Nages=Nages,Nareas=Nareas,
                  ages=ages,agesC=agesC,agesC2=agesC2,sad0=sad0,IC=ICcat,CE=CE,
                  Pups=Pups,NPA=NPA,sdNP=sdNP,Agects=Agects,NFsamp=NFsamp,NPrg=NPrg,
                  YrPct=YrPct,PAidx=PAidx,YrAGsmp=YrAGsmp,H_0=H_0,H_A=H_A,
                  YrPRsmp=YrPRsmp,AgePRsmp=AgePRsmp,psipri1a=psipri1a,
                  psipri1b=psipri1b,psipri2a=psipri2a,psipri2b=psipri2b,
                  omega=omega,Q_A=Q_A,Q_0=Q_0,CV_HV=CV_HV,
                  PApri=PApri,N0pri=N0pri,ICvec=ICvec) 
#
init_fun <- function() {list(sigF=runif(1, .5, .8),
                             sigH=runif(1, .5, 1),
                             sigS=runif(1, 1, 2),
                             tau=runif(1, 8, 12),
                             phi=runif(2, .1, .5),
                             zeta=runif(1, .1, .2),
                             beta1=runif(1, 2.5-.25, 2.5+.25),
                             beta2=runif(1, .25, .35),
                             psi1=runif(1, psipri1a-.25, psipri1a+.25),
                             psi2=runif(1, psipri2a-.25, psipri2a+.25),
                             dlta=runif(1, .05, .1),
                             gamma_H0_mn=runif(1, 5.8, 6.2),
                             gamma_HA_mn=runif(1, 3.5, 4),
                             thta = runif(2, 1, 2),
                             N0 = runif(1, N0pri*.95, N0pri*1.05),
                             alpha0 = runif(1, 2.5, 3),
                             alpha1 = runif(1, .05, .1),
                             alpha2 = runif(1, .005, .01),
                             nu = runif(1, .2, .3)
)}
#
params <- c("sigF","sigH","sigS","tau","phi","thta","beta1","beta2",
            "N0","alpha0","alpha1","alpha2","nu","zeta","dlta",
            "psi1","psi2","gamma_H0_mn","gamma_HA_mn","N",
            "gamma_0","S0_ld","S0_hd","gamma_A","SA_ld","SA_hd","epsF","epsS",
            "Pups_pred","gamma_H0","gamma_HA","H0_pred","HA_pred",
            "Fc1966_prdct","Fc2016_prdct","Fc8_prdct","haz_Ice",
            "nphz","log_lik","y_new1","y_new2") # 
#
cores = detectCores()
ncore = min(20,cores-1)
#
# Run STAN to fit model---------------------------------------------
#
library(cmdstanr)
library(posterior)
#
Niter = round(nsamples/ncore)
#
mod <- cmdstan_model(fitmodel)
#
suppressMessages(
  suppressWarnings ( 
    fit <- mod$sample(
      data = stan.data,
      init = init_fun,
      seed = 123,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      iter_warmup = nburnin,
      iter_sampling = Niter
    )
  )
)
# generate summary stats (sumstats, mcmc matrix)
source("cmdstan_sumstats.r")
#
# Some traceplots to inspect results
mcmc_trace(fit$draws(c("sigF","sigF","sigH")))
#
rm(fit,mod)
#
# Save results --------------------------------------------------------------
#
save.image(file=paste0("./Results/FitHSmod_Results_",
                       format(Sys.time(),  "%b%d_%H%M"),".rdata"))