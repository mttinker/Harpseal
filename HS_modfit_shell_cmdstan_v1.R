# Shell script for fitting harp seal model, Cmdstan version
#   NOTE: 3 pupping areas: 1) S.GSL, 2) N.GSL, 3) Front
#
rm(list = ls())
# User-set parameters--------------------------------------------------------
#
# Specify year range of time series for model fitting 
Year1 = 1951   #  Year1 = first year that harvest data available
YearT = 2020   #  YearT = year AFTER last available data for fitting (harvest/repro data)
Nareas = 3     # Number pupping areas (assume 3: S.GSL, N.GSL, Front)  
# Prior for Initial Population size year 1 (model uses weak prior):
N0pri = 2.5  # prior rough estimate of starting pop size, in millions
# Assumed CV for harvest/bycatch totals:
CV_HV = 0.1 
# Youngest adult age class to consider for age composition fitting
NCage1 = 8; # Note: recommend at least 4, younger ages subject to negative sampling bias
nburnin = 500
nsamples = 5000
# End user parameters--------------------------------------------------------
#
# Load libraries-------------------------------------------------------------
require(parallel)
require(gtools)
require(ggplot2)
require(dplyr)
require(reshape2)
require(stats)
require(bayesplot)
require(fitdistrplus)
require(readxl)
library(cmdstanr)
library(posterior)
rstan::rstan_options(javascript=FALSE)
# require(loo)
# Load data ----------------------------------------------------------------
source("Importdat_v1.r")
# Note: downloaded Ice data: https://iceweb1.cis.ec.gc.ca/IceGraph/page1.xhtml
# Process data--------------------------------------------------------------
Yr0 = Year1-1;   
Nyrs = YearT-Yr0 
Years = Year1:YearT
omega = -7
Nages = 36
ages = seq(1,Nages); ages2 = ages^2
NCages = length(seq(NCage1,Nages))
NCobs = nrow(AgeComp)
Agects = AgeComp
Agects = Agects[,NCage1:Nages]
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
PA = matrix(0,nrow = (Nyrs-1), ncol = 3)
PA_mean = colMeans(NPA)
PA_G_mean = PA_mean[1]/( PA_mean[1] +  PA_mean[2])
for(i in 1:(Nyrs-1)){
  ii = which(PAidx==i)
  iii = which(df.Ice$Year==Yr0+i)
  if(length(ii)>0){
    PA[i,] = NPA[ii,]
  }else{
    PA[i,] = PA_mean
  }
}
#
# Create data inputs for observed harvest, environmental CE index and ice anomalies
IC = matrix(0,nrow = Nyrs, ncol = 3)
CI = numeric(length = Nyrs)
H_0 = numeric(length = (Nyrs-1))
H_A = numeric(length = (Nyrs-1))
Q_0 = numeric(length = (Nyrs-1))
Q_A = numeric(length = (Nyrs-1))
for (y in 1:Nyrs){
  ii = which(df.Ice$Year==Yr0+y)
  if(length(ii)==1){
    IC[y,1] = Gulf_Anom[ii]
    IC[y,2] = Gulf_Anom[ii]
    IC[y,3] = Front_Anom[ii]   
  }
  ii = which(df.NLCI$Year==Yr0+y)
  if(length(ii)==1){
    CI[y] = df.NLCI$CI[ii]
  }
  ii = which(df.HV$YEAR==Yr0+y)
  if(length(ii)==1 & y<Nyrs){
    H_0[y] = df.HV$PUPTOT[ii]
    H_A[y] = df.HV$ADLTOT[ii]
    Q_A[y] = df.HV$ADL_prob_rec[ii]
    Q_0[y] = df.HV$PUP_prob_rec[ii]
  }
}
IC = IC[-Nyrs,] # 
# Create a "dummy vector" of ice anomaly values from -1 to 1 (step size 0.05), 
#  and discretized version of IC anomaly values classified to units of 0.05 
#  (to speed up fitting)
ICvec = seq(-1,1,by=.05)
#
# Newborn pup survival (accounts for abortions and deaths prior to up survey)
# - for current model fix this at 0.95
Snp=rep(0.95,Nyrs-1)
#
rm(i,ii,y,aa) 
#
# agesC =  c(0,ages) - 10; agesC[Nages+1] = agesC[Nages+1] + 4; agesC2 = agesC^2
agesC =  pmax(0,10 - c(0,ages))
agesC2 = (pmax(0,c(0,ages)-10))^2
# Set up Jags inputs --------------------------------------------------------
#  
stan.data <- list(NPcts=NPcts,NCages=NCages,NCage1=NCage1,NCobs=NCobs,
                  NPRobs=NPRobs,Nyrs=Nyrs,Nages=Nages,Nareas=Nareas,
                  ages=ages,agesC=agesC,agesC2=agesC2,sad0=sad0,IC=IC,CI=CI,
                  Pups=Pups,PA=PA,sdNP=sdNP,Agects=Agects,NFsamp=NFsamp,NPrg=NPrg,
                  YrPct=YrPct,PAidx=PAidx,YrAGsmp=YrAGsmp,H_0=H_0,H_A=H_A,
                  YrPRsmp=YrPRsmp,AgePRsmp=AgePRsmp,
                  omega=omega,Q_A=Q_A,Q_0=Q_0,CV_HV=CV_HV,
                  N0pri=N0pri,ICvec=ICvec,Snp=Snp) 
#
init_fun <- function() {list(sigF=runif(1, .4, .6),
                             sigH=runif(1, .6, .8),
                             sigS=runif(1, 1, 1.5),
                             tau10=runif(1, 10, 15),
                             phi=runif(2, c(.8,.4), c(1.2,.8)),
                             thta = runif(2, c(.5,.8), c(1,1.2)),
                             zeta=runif(1, 7, 12),
                             beta1=runif(1, 3-.25, 3+.25),
                             beta2=runif(1, .25, .35),
                             dlta=runif(2, .4, .7),
                             gamma_H0_mn=runif(1, 5.7, 6.2),
                             gamma_HA_mn=runif(1, 3.5, 4),
                             N0mil = runif(1, N0pri*.95, N0pri*1.05),
                             alpha0 = runif(1, 1.8, 1.9),
                             alpha1 = runif(1, .13, .14),
                             alpha2 = runif(1, .006, .007)
)} #  nu = runif(1, .2, .3)
#
parms <- c("sigF","sigH","sigS","tau","phi","thta","beta1","beta2",
            "N0","alpha0","alpha1","alpha2","zeta","dlta",
            "psi1","psi2","gamma_H0_mn","gamma_HA_mn","N","gamma_0",
            "S0_ld","S0_hd","SA_ld","SA_hd","gamma_A","epsF","epsS",
            "Pups_pred","gamma_H0","gamma_HA","H0_pred","HA_pred",
            "Fc1966_prdct","Fc2016_prdct","Fc8_prdct","haz_Ice",
            "Pyng_prdct","Pold_prdct","y_new1","y_new2") # 
#
cores = detectCores()
if(cores>20){
  ncore = min(20,cores-10)
}else{
  ncore = min(10,cores-3)
}
#
fitmodel = c("HSmodfit_v1.stan")
#
# Run STAN to fit model---------------------------------------------
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
      seed = 1234,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      iter_warmup = nburnin,
      iter_sampling = Niter,
      max_treedepth = 12,
      adapt_delta = 0.8
    )
  )
)
# generate summary stats (sumstats, mcmc matrix)
source("cmdstan_sumstats.r")
#
# Some traceplots to inspect results
mcmc_trace(mcmc_array,pars=("sigF"))
mcmc_trace(mcmc_array,pars=("alpha1"))
mcmc_trace(mcmc_array,pars=vn[startsWith(vn,"dlta")])
#
rm(fit,mod)
#
# Save results --------------------------------------------------------------
#
resultsfile = paste0("HS_Results_",format(Sys.time(),"%y_%b%d_%H%M"),"_agmn",NCage1,"_v1.rdata")
save.image(file=paste0("./Results/",resultsfile)) 
#
