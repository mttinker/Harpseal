# Script to plot some results from model fitting
# Load results file (if not already loaded into workspace):
load("Results.rdata")
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
require(bayesplot)

futuresim = 1; # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 50; Yearst2 = 2020 ; reps = 20000
PAmeans = c(.18,.07,.75) # future proportion pups in S Gulf, N Gulf, Front
# Future conditions: based on ice and CE indices from after year YY, 
#  fit appropriate random sampling distributions
N_end = sumstats[which(vns==paste0("N[",Nyrs,"]")),1]  
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
sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps,
                 PAmeans=PAmeans,futuresim=futuresim,
                 NFages=NFages,NFage1=NFage1,Nages=Nages,Nareas=Nareas,
                 ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CE=CE2,DDadlt=DDadlt,
                 Adloghz=Adloghz,Jvloghz=Jvloghz,gamm0=gamm0) 
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
cores=detectCores()
cl <- makeCluster(min(20,cores[1]-1)) 
registerDoParallel(cl)
source("HSmod_sim.r")
rsltK=HSmod_sim(init_fun,sim.data,sumstats,vns)
N_Predict = rsltK$N_Predict
P_Predict = rsltK$P_Predict
P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
Np_Predict = 1.1*(N_Predict + P_Predict)
Nfin = colMeans(Np_Predict[(Nyrs2-20):Nyrs2,])
Kest = mean(Nfin)
Kest_sd = sd(Nfin)
Kest_CI = quantile(Nfin, prob=c(0.025,0.975))
Ktab = data.frame(Metric = "Equilibrium abundance",
           Mean = Kest/1000000, SD = Kest_sd/1000000, 
           CI95_low=Kest_CI[1]/1000000,CI95_high=Kest_CI[2]/1000000)
print(Ktab)
# Set up sims for 15 years under varying harvest levels
futuresim = 2 # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 15; Yearst2 = 2020 ; 
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
# Future ice/env conditions: same as for K sime abive
#  fit appropriate random sampling distributions
#
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
#
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
dat = data.frame(Nmin=Nmin,HvT=HvT,ppA=ppA,intx = HvT*ppA)
# Linear model (possibly with polynomial terms)
fit1 = lm(Nmin~poly(HvT,3)+poly(ppA,3)+poly(intx,3), data=dat)
# summary(fit1)
# ii = which(dat$ppA>0.09 & dat$ppA<0.11)
# newdat1 = dat[ii,]
# newdat1 <- newdat1[with(newdat1,order(newdat1$HvT)) ,]
# predtest1 = predict(fit1,newdata = newdat1,se.fit=T,interval = "prediction",level = 0.8)$fit
# plot(predtest1[,1],newdat1$Nmin)
# abline(coef = c(0,1))
# newdat2 = data.frame(HvT=seq(301,700),ppA=rep(.1,400))
# newdat2$intx = newdat2$HvT*newdat2$ppA
# predtest2 = predict(fit1,newdata = newdat2,se.fit=T,interval = "prediction",level = 0.8)$fit
# plot(newdat1$HvT,newdat1$Nmin)
# lines(newdat2$HvT,predtest2[,1],col="red")
# lines(newdat2$HvT,predtest2[,2],col="red")
# lines(newdat2$HvT,predtest2[,3],col="red")
#
HvTeval05 = seq(201,900)
HvTeval10 = seq(151,850)
HvTeval50 = seq(101,800)
tmp = data.frame(HvT=HvTeval05,ppA=rep(.05,700)); tmp$intx = tmp$HvT*tmp$ppA
newdat05ad = tmp; 
tmp = data.frame(HvT=HvTeval10,ppA=rep(.1,700)); tmp$intx = tmp$HvT*tmp$ppA
newdat10ad = tmp;
tmp = data.frame(HvT=HvTeval50,ppA=rep(.5,700)); tmp$intx = tmp$HvT*tmp$ppA
newdat50ad = tmp;
pred05ad = predict(fit1,newdata = newdat05ad,se.fit=T,interval = "prediction",level = 0.9)$fit
pred10ad = predict(fit1,newdata = newdat10ad,se.fit=T,interval = "prediction",level = 0.9)$fit
pred50ad = predict(fit1,newdata = newdat50ad,se.fit=T,interval = "prediction",level = 0.9)$fit
ii = which(abs(pred05ad[,2]-(N70+1000))==min(abs(pred05ad[,2]-(N70+1000))))
TACN70pa05 = floor(HvTeval05[ii]/5)*5  
ii = which(abs(pred05ad[,2]-(N50+1000))==min(abs(pred05ad[,2]-(N50+1000))))
TACN50pa05 = floor(HvTeval05[ii]/5)*5 
ii = which(abs(pred10ad[,2]-(N70+1000))==min(abs(pred10ad[,2]-(N70+1000))))
TACN70pa10 = floor(HvTeval10[ii]/5)*5  
ii = which(abs(pred10ad[,2]-(N50+1000))==min(abs(pred10ad[,2]-(N50+1000))))
TACN50pa10 = floor(HvTeval10[ii]/5)*5  
ii = which(abs(pred50ad[,2]-(N70+1000))==min(abs(pred50ad[,2]-(N70+1000))))
TACN70pa50 = floor(HvTeval50[ii]/5)*5  
ii = which(abs(pred50ad[,2]-(N50+1000))==min(abs(pred50ad[,2]-(N50+1000))))
TACN50pa50 = floor(HvTeval50[ii]/5)*5  
# 
# Create table of TAC recomendations
TACtab = data.frame(Pcnt_adlt=c(5,10,50),
                    TAC_N70 = c(TACN70pa05,TACN70pa10,TACN70pa50),
                    TAC_N50 =c(TACN50pa05,TACN50pa10,TACN50pa50))
print(TACtab)
save(file = "TAC_K_est.rdata",TACtab,Ktab)
stopCluster(cl)
