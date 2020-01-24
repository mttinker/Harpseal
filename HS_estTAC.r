# Script to plot some results from model fitting
# Load results file (if not already loaded into workspace):
loadfile = file.choose(new = FALSE)
load(file=loadfile)
#
library(foreach)
library(doParallel)
library(fitdistrplus)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(stats)
require(bayesplot)
#require(loo)
# Set up sims for 15 years under varying harvest levels
futuresim = 2 # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 15; Yearst2 = 2020 ; reps = 10000
PAmeans = c(.18,.07,.75) # future proportion pups in S Gulf, N Gulf, Front
# Set other human mortality not included in Canadian harvest
GrnH = 50000
Grn_A = .857*GrnH
Grn_P = .143*GrnH
ArcH = 1000
Arc_A = .95*ArcH
Arc_P = .05*ArcH
Byctc = 2000
Byc_A = .2*Byctc
Byc_P = .8*Byctc
# Future conditions: based on ice and CE indices from after year YY, 
#  fit appropriate random sampling distributions
YY = 2000
ii = which(df.CE$Year>=YY)
ft = fitdist(log(df.CE$CEindex[ii])+1,"gamma")
CE2 = rgamma(1000,coef(ft)[1],coef(ft)[2])-1
ii = which(df.Ice$Year>=YY)
ft = fitdist(exp(df.Ice$Gulf_Anom[ii]),"norm")
icg1 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
ft = fitdist(exp(df.Ice$Lab_Anom[ii]),"norm")
icg2 = log(pmax(0.368,pmin(2.715,rnorm(1000,coef(ft)[1],coef(ft)[2]))))
IC2 = as.matrix(cbind(icg1,icg1,icg2))
#
N_end = sumstats[which(vns==paste0("N[",Nyrs,"]")),1]
N_all = sumstats[which(startsWith(vns,"N[")),1]
N70 = 0.7*max(N_all)
N50 = 0.5*max(N_all)
sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps,
                  PAmeans=PAmeans,futuresim=futuresim,
                  NFages=NFages,NFage1=NFage1,Nages=Nages,Nareas=Nareas,
                  ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CE=CE2,DDadlt=DDadlt,
                  Adloghz=Adloghz,Jvloghz=Jvloghz,gamm0=gamm0,
                  Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
                  Byc_A=Byc_A,Byc_P=Byc_P) 
init_fun <- function() {list(sigF=rnorm(1,sumstats[vns=="sigF",1],sumstats[vns=="sigF",2]),
                             sigH=rnorm(1,sumstats[vns=="sigH",1],sumstats[vns=="sigH",2]),
                             phiJ=rnorm(1,sumstats[vns=="phiJ",1],sumstats[vns=="phiJ",2]),
                             phiF=rnorm(1,sumstats[vns=="phiF",1],sumstats[vns=="phiF",2]),
                             b0=rnorm(1,sumstats[vns=="b0",1],sumstats[vns=="b0",2]),
                             b1=rnorm(1,sumstats[vns=="b1",1],sumstats[vns=="b1",2]),
                             thta=rnorm(1,sumstats[vns=="thta",1],sumstats[vns=="thta",2]),
                             psi1=rnorm(1,sumstats[vns=="psi1",1],sumstats[vns=="psi1",2]),
                             psi2=rnorm(1,sumstats[vns=="psi2",1],sumstats[vns=="psi2",2]),
                             dlta=rnorm(1,sumstats[vns=="dlta",1],sumstats[vns=="dlta",2]),
                             gammHp_mn=rnorm(1,sumstats[vns=="gammHp_mn",1],sumstats[vns=="gammHp_mn",2]),
                             gammHa_mn=rnorm(1,sumstats[vns=="gammHa_mn",1],sumstats[vns=="gammHa_mn",2])
)}
#
# Set up parallel backend to use many processors --------------
cores=detectCores()
cl <- makeCluster(min(20,cores[1]-1)) 
registerDoParallel(cl)
source("HSmod_sim.r")
rslt=HSmod_sim(init_fun,sim.data,sumstats,vns)
N_Predict = rslt$N_Predict
Hvp = rslt$HVp_predict
Hva = rslt$HVa_predict
Hvp_k = apply(Hvp,2,mean)/1000 
Hva_k = apply(Hva,2,mean)/1000 
ii = which(Hvp_k>0 & Hva_k>0)
HvT = Hvp_k[ii]+Hva_k[ii]
ppA = Hva_k[ii]/HvT
Nmin = apply(N_Predict,2,min); Nmin = Nmin[ii]
dat = data.frame(Nmin=Nmin,HvT=HvT,ppA=ppA)
fit1 = lm(Nmin~HvT+ppA, data=dat)
# summary(fit1)
# predtest = predict(fit1,newdata = dat)
# plot(predtest,Nmin)
# abline(coef = c(0,1))
HvTeval = seq(51,1000)
newdat05ad = data.frame(HvT=HvTeval,ppA=rep(.05,950))
newdat10ad = data.frame(HvT=HvTeval,ppA=rep(.1,950))
newdat50ad = data.frame(HvT=HvTeval,ppA=rep(.5,950))
pred05ad = predict(fit1,newdata = newdat05ad,se.fit=T,interval = "prediction",level = 0.8)$fit
pred10ad = predict(fit1,newdata = newdat10ad,se.fit=T,interval = "prediction",level = 0.8)$fit
pred50ad = predict(fit1,newdata = newdat50ad,se.fit=T,interval = "prediction",level = 0.8)$fit
ii = which(abs(pred05ad[,2]-(N70+1000))==min(abs(pred05ad[,2]-(N70+1000))))
TACN70pa05 = floor(HvTeval[ii]/5)*5  
ii = which(abs(pred05ad[,2]-(N50+1000))==min(abs(pred05ad[,2]-(N50+1000))))
TACN50pa05 = floor(HvTeval[ii]/5)*5 
ii = which(abs(pred10ad[,2]-(N70+1000))==min(abs(pred10ad[,2]-(N70+1000))))
TACN70pa10 = floor(HvTeval[ii]/5)*5  
ii = which(abs(pred10ad[,2]-(N50+1000))==min(abs(pred10ad[,2]-(N50+1000))))
TACN50pa10 = floor(HvTeval[ii]/5)*5  
ii = which(abs(pred50ad[,2]-(N70+1000))==min(abs(pred50ad[,2]-(N70+1000))))
TACN70pa50 = floor(HvTeval[ii]/5)*5  
ii = which(abs(pred50ad[,2]-(N50+1000))==min(abs(pred50ad[,2]-(N50+1000))))
TACN50pa50 = floor(HvTeval[ii]/5)*5  
# 

# Create table of TAC recomendations
TACtab = data.frame(Pcnt_adlt=c(5,10,50),
                    TAC_N70 = c(TACN70pa05,TACN70pa10,TACN70pa50),
                    TAC_N50 =c(TACN50pa05,TACN50pa10,TACN50pa50))
print(TACtab)
stopCluster(cl)