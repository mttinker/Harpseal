# Shell script for fitting harp seal model
#   NOTE: 3 pupping areas: S.GSL, N.GSL, Front
#
# User-set parameters--------------------------------------------------------
#
# Specify year range of time series for model fitting 
Year1 = 1951   #  Year1 = first year that harvest data available
YearT = 2020   #  YearT = year AFTER last available data for fitting
Nareas = 3     # Number pupping areas (assume 3: S.GSL, N.GSL, Front)  
# Prior for Initial Population size year 1 (model uses weak prior):
N0pri = 1500000 # prior guess of starting pop size, default ~ 2 million
# Assumed CV for harvest/bycatch totals:
CV_HV = 0.1 
# Priors for demographic params
Adl_Sx = 0.94    # Adult base survival (no harvest), low density (N~ 1 million)
Juv_Sx = 0.82    # Juvenile base survival (no harvest), low density (N~ 1 million)
b0 = 2.5        # prior for logit of max adult pregancy rate (2.5 --> Fx=0.92)
# Ice anomaly effects on pup survival:
psi1 = 4        # prior for psi param1 of ice anomaly effect fxn (see figure)
psi2 = .8       # prior for psi param2 of ice anomaly effect fxn 
#
# End user params, plot ice mort prior---------------------------------------
#
# Sample ice anomaly mortality fxn params (explore graphically):
psipri1a = psi1   # prior for psi param1 mean: ice anomaly effect fxn
psipri1b = .5   # prior for psi param1 sd: ice anomaly effect fxn 
psipri2a = psi2 # prior for psi param2 mean: ice anomaly effect fxn 1.5
psipri2b = .2 # prior for psi param2 sd: ice anomaly effect fxn   .5
#source("Ice_mort_plot.r")
#Ice_mort_plot(psipri1a,psipri1b,psipri2a,psipri2b)
#
# Load libraries-------------------------------------------------------------
require(parallel)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
require(rstan)
require(stats)
require(bayesplot)
require(fitdistrplus)
require(readxl)
# require(loo)
# Set options for STAN
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# Load data ----------------------------------------------------------------
source("Importdat.r")
# Process data--------------------------------------------------------------
Yr0 = Year1-1;   
Nyrs = YearT-Yr0 
gamm0 = -7
Nages = 8
ages = seq(1,Nages); ages2 = ages^2
NFage1 = 3; # Youngest age of adults sampled for preg status
NFages = length(seq(NFage1,8))
NFobs = nrow(AgeComp)
NPRobs = nrow(df.Rep)
NPcts = nrow(df.Pup)
NPctsA = nrow(df.Pup[!is.na(df.Pup$Sgulf_N),])
YrAGsmp = Years_ages - Yr0
NFsamp = df.Rep$N
Nprg = df.Rep$Preg
YrPRsmp = df.Rep$Year - Yr0
AgePRsmp = df.Rep$Age
sad0 = df.sad$SAD
NP = numeric()
sdNP = numeric()
YrPct = numeric()
PAtag = numeric(length = Nyrs-1); 
NPA = matrix(0,nrow = Nyrs-1,ncol = Nareas)
aa = c(2,4,6)
for (i in 1:nrow(df.Pup)){
  YrPct[i] = df.Pup$Year[i] - Yr0
  NP[i] = df.Pup$Total_Npup[i]
  sdNP[i] = min(df.Pup$Total_se[i],NP[i]*.2)
  if (!is.na(df.Pup$Sgulf_N[i])){
    NPA[YrPct[i],1:Nareas] = as.numeric(df.Pup[i,aa]/sum(df.Pup[i,aa]) )
    PAtag[YrPct[i]] = 1
  }
}
# Set priors for beta distributions of proportion pups in areas 1 & 2
#  * proportion in area 3, PA3, is then be calculated as 1 - (PA1+PA2) *
PApri = matrix(nrow = Nareas-1,ncol = 2)
ft = fitdist(NPA[PAtag==1,1],"beta")
PApri[1,1] = as.numeric(ft$estimate[1]); PApri[1,2] = as.numeric(ft$estimate[2]); 
ft = fitdist(NPA[PAtag==1,2],"beta")
PApri[2,1] = as.numeric(ft$estimate[1]); PApri[2,2] = as.numeric(ft$estimate[2]); 
IC = matrix(0,nrow = Nyrs, ncol = 3)
CE = numeric(length = Nyrs)
HVp = numeric(length = Nyrs)
HVa = numeric(length = Nyrs)
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
  if(length(ii)==1){
    HVp[y] = df.HV$PUPTOT[ii]
    HVa[y] = df.HV$ADLTOT[ii]
  }
}
# Calculate adult and juvenile base log hazards at low densities, 
#  based on user-provided annual survival rate estimates
thta = 2.5
gmvals = seq(2,7,by=0.01)
SAvals = exp(-1*(exp(gamm0 + gmvals + (.033*.22*1)^thta) + exp(gamm0)))
SJvals = exp(-1*(exp(gamm0 + gmvals + (.22*1)^thta) + exp(gamm0) + exp(gamm0+.5)))
Adloghz = gmvals[which(abs(SAvals-Adl_Sx)==min(abs(SAvals-Adl_Sx)))]
Jvloghz = gmvals[which(abs(SJvals-Juv_Sx)==min(abs(SJvals-Juv_Sx)))]
# scaled DD effects (relative to juv) for 1st 8 adult age classes:
DDadlt = c(.25,.1,0,0,0,0,0,0) 
b0pri = b0
rm(i,ii,y,aa,ft) 
#
# Set up Jags inputs --------------------------------------------------------
#
fitmodel = c("HSmodfit.stan")
#  
stan.data <- list(NPcts=NPcts,NPctsA=NPctsA,NFages=NFages,NFage1=NFage1,
                  NFobs=NFobs,NPRobs=NPRobs,Nyrs=Nyrs,Nages=Nages,Nareas=Nareas,
                  ages=ages,ages2=ages2,sad0=sad0,IC=IC,CE=CE,HVp=HVp,HVa=HVa,
                  NP=NP,NPA=NPA,sdNP=sdNP,AgeF=AgeComp,NFsamp=NFsamp,
                  Nprg=Nprg,YrPct=YrPct,PAtag=PAtag,YrAGsmp=YrAGsmp,
                  YrPRsmp=YrPRsmp,AgePRsmp=AgePRsmp,DDadlt=DDadlt,psipri1a=psipri1a,
                  psipri1b=psipri1b,psipri2a=psipri2a,psipri2b=psipri2b,
                  b0pri=b0pri,psi1=psi1,psi2=psi2,Adloghz=Adloghz,Jvloghz=Jvloghz,
                  CV_HV=CV_HV,gamm0=gamm0,PApri=PApri,N0pri=N0pri) # ,thta=thta
#
init_fun <- function() {list(sigF=runif(1, .8, 1),
                             sigH=runif(1, .8, 1),
                             phiJ=runif(1, .22, .28),
                             phiF=runif(1, .26, .3),
                             b0=runif(1, b0pri-.2, b0pri+.2),
                             b1=runif(1, .18, .2),
                             psi1=runif(1, psipri1a-1, psipri1a+1),
                             psi2=runif(1, psipri2a-.5, psipri2a+.5),
                             dlta=runif(1, .02, .05),
                             gammHp_mn=runif(1, 5.5, 6.2),
                             gammHa_mn=runif(1, 3, 3.5),
                             thta = runif(1, 1, 2)
                             )}
#
# For testing inits-----------------------------------------------------
#
# source("HSmod_test.r")
# Year = seq(Year1,YearT-1)
# Yearp = seq(Year1,YearT)
# rslt=HSmod_test(init_fun,stan.data)
# N_Predict = rslt$N_Predict
# P_Predict = rslt$P_Predict
# PrPredict = rslt$PrPredict
# PrAgPred = rslt$PrAgPred
# Agepredict1 = rslt$Agepredict1 # Age dist for sample year 17 = 1967 (Agecomp row5)
# Agepredict2 = rslt$Agepredict2 # Age dist for sample year 69 = 2019 (Agecomp row48)
# HVp_predict = rslt$HVp_predict
# HVa_predict = rslt$HVa_predict
# ggplot(data.frame(Year=Yearp,Npred = rowMeans(N_Predict)),aes(x=Year,y=Npred)) +
#   geom_line()
# ggplot(data.frame(Year=Year,HvPpred = rowMeans(HVp_predict)),aes(x=Year,y=HvPpred)) +
#   geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=PUPTOT))
# ggplot(data.frame(Year=Year,HvApred = rowMeans(HVa_predict)),aes(x=Year,y=HvApred)) +
#   geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=ADLTOT))
# ggplot(data.frame(Year=Year,Ppred = rowMeans(P_Predict)),aes(x=Year,y=Ppred)) +
#   geom_line() + geom_point(data=df.Pup,aes(x=Year,y=Total_Npup))
# 
# ggplot(data.frame(Age=ages[3:8], Frequency= rowMeans(Agepredict1)),aes(x=Age,y=Frequency)) +
#   geom_line() + geom_point(data=data.frame(Age=ages[3:8], Frequency= AgeComp[5,]/sum(AgeComp[5,])),aes(x=Age,y=Frequency))
# ggplot(data.frame(Age=ages[3:8], Frequency= rowMeans(Agepredict2)),aes(x=Age,y=Frequency)) +
#   geom_line() + geom_point(data=data.frame(Age=ages[3:8], Frequency= AgeComp[48,]/sum(AgeComp[48,])),aes(x=Age,y=Frequency))
# ggplot(data.frame(Year=Year,Prt8yrPred = rowMeans(PrPredict)),aes(x=Year,y=Prt8yrPred)) +
#   geom_line() + geom_point(data=df.Rep[df.Rep$Age==8,],aes(x=Year,y=Prob))
# ggplot(data.frame(Age=ages, Pregrate= rowMeans(PrAgPred)),aes(x=Age,y=Pregrate)) +
#   geom_line() + geom_point(data=df.Rep[df.Rep$Year<=1970,],aes(x=Age,y=Prob))
# 
# # sad = rslt$SAD
# rm(rslt,N_Predict,P_Predict,PrPredict,PrAgPred,Agepredict1,
#    Agepredict2,HVp_predict,HVa_predict)
#
# Run JAGS to fit model---------------------------------------------
params <- c("sigF","sigH","phiJ","phiF","b0","b1","N0","dlta",
            "psi1","psi2","gammHp_mn","gammHa_mn","thta",
            "N","PredPup","gammHp","gammHa","HVp_pred","HVa_pred",
            "Fc1966_prdct","Fc2016_prdct","Fc8_prdct","epsF",
            "Tstat","Tstat_new","ppp","log_lik") # 
#
nsamples <- 500
nburnin <- 500
cores = detectCores()
ncore = min(20,cores-1)
#cl <- makeCluster(ncore)
nc <- ncore

out <- stan(
  file = fitmodel,     # Stan program
  data = stan.data,    # named list of data
  pars = params,       # list of params to monitor
  init= init_fun,       # initial values    "random"            
  chains = nc,         # number of Markov chains
  warmup = nburnin,        # number of warmup iterations per chain
  iter = nburnin+nsamples, # total number of iterations per chain
  cores = nc,          # number of cores 
  refresh = 100,         # show progress every 'refresh' iterations
  control = list(adapt_delta = 0.99, max_treedepth = 15) # increase to help find optimal vals
)
#
# Calclate Sumstats -------------------------------------------------
#
mcmc <- as.matrix(out)
vn = colnames(mcmc)
Nsims = nrow(mcmc)
sumstats = summary(out)$summary
vns = row.names(sumstats)

traceplot(out, pars=c("sigF","sigH"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("b1"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("b0"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("phiF","phiJ"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("gammHp_mn"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("gammHa_mn"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("dlta"), inc_warmup = F, nrow = 2)
traceplot(out, pars=c("psi1","psi2"), inc_warmup = F, nrow = 2)

save.image(file=paste0("./Results/FitHSmod",fitoption,"_Results_",
                       format(Sys.time(), "%b%d"),".rdata"))
