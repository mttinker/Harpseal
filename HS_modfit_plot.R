# Script to plot some results from model fitting
# Load results file (if not already loaded into workspace):
load(file="FitHSmod_Results.rdata")
#
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
require(loo)
Year = seq(Year1,Year1+Nyrs-2)
Yearp = seq(Year1,Year1+Nyrs-1)
# Pop trends ----------------------------------------------------------------
dp1 = data.frame(Year=Yearp,Nexp=sumstats[startsWith(vns,"N[")==T,1],
                 N_lo = sumstats[startsWith(vns,"N[")==T,4],
                 N_hi = sumstats[startsWith(vns,"N[")==T,8])

dp1$N2exp = c(dp1$Nexp[1:(Nyrs-1)] + sumstats[startsWith(vns,"PredPup[")==T,1],NA)
dp1$N2_lo = c(dp1$N_lo[1:(Nyrs-1)] + sumstats[startsWith(vns,"PredPup[")==T,4],NA)
dp1$N2_hi = c(dp1$N_hi[1:(Nyrs-1)] + sumstats[startsWith(vns,"PredPup[")==T,8],NA)

pl1 = ggplot(data=dp1[1:(Nyrs-1),],aes(x=Year,y=Nexp)) +
      geom_ribbon(aes(ymin=N_lo,ymax=N_hi),alpha=0.3) +
      geom_line() + labs(x = "Year",y="Estimated abundance w/o pups") +
      ggtitle("Model estimated abundance (excluding pups)") + theme_bw()
print(pl1)
pl1b = ggplot(data=dp1[1:(Nyrs-1),],aes(x=Year,y=N2exp)) +
  geom_ribbon(aes(ymin=N2_lo,ymax=N2_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Estimated abundance with pups") +
  ggtitle("Model estimated abundance (including pups)") + theme_bw() 
print(pl1b)
#
# Pup counts ----------------------------------------------------------------
dp2 = data.frame(Year=Year,Pexp=sumstats[startsWith(vns,"PredPup[")==T,1],
                 P_lo = sumstats[startsWith(vns,"PredPup[")==T,4],
                 P_hi = sumstats[startsWith(vns,"PredPup[")==T,8])
dp2$Obs = numeric(length = nrow(dp2))
dp2$ObsSE = numeric(length = nrow(dp2))
for (i in 1:nrow(df.Pup)){
  ii = which(dp2$Year==df.Pup$Year[i])
  dp2$Obs[ii] = df.Pup$Total_Npup[i]
  dp2$ObsSE[ii] = df.Pup$Total_se[i]
}
dp2$Obs[dp2$ObsSE==0] = NA
dp2$ObsSE[dp2$ObsSE==0] = NA
pl2 = ggplot(data=dp2,aes(x=Year,y=Pexp)) +
  geom_ribbon(aes(ymin=P_lo,ymax=P_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Pup counts (total)") +
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-1.96*ObsSE, ymax = Obs+1.96*ObsSE)) +
  ggtitle("Model estimated vs observed pup counts") + theme_classic()
print(pl2)
#
# Preg rate by age----------------------------------------------------------
PR1966 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][16]
PR1966per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][1:32],na.rm = T)
crct1 = PR1966per/PR1966
PR2016 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][66]
PR2016per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][50:69],na.rm = T)
crct2 = PR2016per/PR2016

dp3 = data.frame(Age=ages, Year = (rep(1966,Nages)),
                 PRexp=sumstats[startsWith(vns,"Fc1966_prdct[")==T,1]*crct1,
                 PR_lo = sumstats[startsWith(vns,"Fc1966_prdct[")==T,4]*crct1,
                 PR_hi = sumstats[startsWith(vns,"Fc1966_prdct[")==T,8]*crct1)
dp3 = rbind(dp3, data.frame(Age=ages, Year = (rep(2016,Nages)),
                 PRexp=sumstats[startsWith(vns,"Fc2016_prdct[")==T,1]*crct2,
                 PR_lo = sumstats[startsWith(vns,"Fc2016_prdct[")==T,4]*crct2,
                 PR_hi = sumstats[startsWith(vns,"Fc2016_prdct[")==T,8]*crct2))
dp3$Obs = numeric(length = nrow(dp3))
dp3$ObsSE = numeric(length = nrow(dp3))
for (i in 1:nrow(dp3)){
  ii = which(df.Rep$Age==dp3$Age[i] & df.Rep$Year>(dp3$Year[i]-20) & 
               df.Rep$Year<(dp3$Year[i]+20))
  dp3$Obs[i] = mean(df.Rep$Prob[ii])
  dp3$ObsSE[i] = sd(df.Rep$Prob[ii])/sqrt(length(df.Rep$Prob[ii]))
}
dp3$Year = as.factor(dp3$Year)
pl3 = ggplot(data=dp3,aes(x=Age,y=PRexp,group=Year,color = Year, fill=Year)) +
  geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi),alpha=0.3) + ylim(0,1) +
  geom_line() + labs(x = "Age",y="Pregnancy rate") +
  geom_point(aes(y=Obs),size=2) +
  geom_errorbar(aes(ymin = Obs-1.96*ObsSE, 
                      ymax = Obs+1.96*ObsSE),width=.2) +
  ggtitle("Model estimated vs observed pregancy rates by age (early vs late time series)") + theme_classic()
print(pl3)
#
# Preg rate 8+ over time----------------------------------------------------
dp4 = data.frame(Year=Year,PRexp=sumstats[startsWith(vns,"Fc8_prdct[")==T,1],
                 PR_lo = sumstats[startsWith(vns,"Fc8_prdct[")==T,4],
                 PR_hi = sumstats[startsWith(vns,"Fc8_prdct[")==T,8])
dp4$Obs = numeric(length = nrow(dp4)) 
dp4$ObsSE = numeric(length = nrow(dp4))
for (i in 1:nrow(df.Rep)){
  if (df.Rep$Age[i] == 8 & df.Rep$N[i]>=10){
    ii = which(dp4$Year==df.Rep$Year[i])
    dp4$Obs[ii] = df.Rep$Prob[i]
    dp4$ObsSE[ii] = sqrt((df.Rep$Prob[i]*(1-df.Rep$Prob[i]))/df.Rep$N[i])
  }
}
dp4$Obs[dp4$ObsSE==0] = NA
dp4$ObsSE[dp4$ObsSE==0] = NA
pl4 = ggplot(data=dp4,aes(x=Year,y=PRexp)) +
  geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Preganancy rate (Age 8+)") +
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-1.96*ObsSE, 
                                             ymax = Obs+1.96*ObsSE),color="darkgrey") +
  ggtitle("Model estimated vs observed pregnancy rate over time") + theme_classic()
print(pl4)
#
# Harvest mortality over time ---------------------------------------------
dp5 = data.frame(Year=Year, Group = rep("Beaters",length(Year)),
                 HVexp=sumstats[startsWith(vns,"HVp_pred[")==T,1],
                 HV_lo = sumstats[startsWith(vns,"HVp_pred[")==T,4],
                 HV_hi = sumstats[startsWith(vns,"HVp_pred[")==T,8])
dp5 = rbind(dp5, data.frame(Year=Year, Group = rep("Adult",length(Year)),
                 HVexp=sumstats[startsWith(vns,"HVa_pred[")==T,1],
                 HV_lo = sumstats[startsWith(vns,"HVa_pred[")==T,4],
                 HV_hi = sumstats[startsWith(vns,"HVa_pred[")==T,8]))
dp5$Obs = numeric(length = nrow(dp5))
for (i in 1:nrow(dp5)){
  ii = which(df.HV$YEAR==dp5$Year[i])
  if(dp5$Group[i]=="Beaters"){
    dp5$Obs[i] = df.HV$PUPTOT[ii]
  }else{
    dp5$Obs[i] = df.HV$ADLTOT[ii]
  }
}
dp5$Group = as.factor(dp5$Group)
pl5 = ggplot(data=dp5,aes(x=Year,y=HVexp,group=Group,color = Group, fill=Group)) +
  geom_ribbon(aes(ymin=HV_lo,ymax=HV_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Total harvest + bycatch") +
  geom_point(aes(y=Obs),size=2) +
  ggtitle("Model estimated vs observed harvest/bycatch mortality, by age class") + theme_classic()
print(pl5)
#
# Ice anomaly effect on pup survival---------------------------------------
psi1mn = sumstats[startsWith(vns,"psi1")==T,1]  
psi1sd = sumstats[startsWith(vns,"psi1")==T,3]  
ps2mn = sumstats[startsWith(vns,"psi2")==T,1]  
ps2sd = sumstats[startsWith(vns,"psi2")==T,3]  
source("Ice_mort_plot.r")
plIce = Ice_mort_plot(psi1mn,psi1sd,ps2mn,ps2sd)
print(plIce)
#
# Fecundity vs density ----------------------------------------------------
NN = seq(.2,5,by=0.1)
iir = sample(Nsims,1000)
b0_r = mcmc[iir,vn=="b0"]
phiF_r = mcmc[iir,vn=="phiF"]
thta_r = mcmc[iir,vn=="thta"]
FC = matrix(0,nrow = 1000,ncol=length(NN))
FC_mean = numeric(length = 40)
FC_lo = numeric(length = 40)
FC_hi = numeric(length = 40)
for(r in 1:1000){
  FC[r,] = inv.logit(b0_r[r] - phiF_r[r]*NN^thta_r[r])
}
FC_mean = colMeans(FC)
FC_lo = apply(FC, 2, quantile, prob=0.055)
FC_hi = apply(FC, 2, quantile, prob=0.95)
dp6 = data.frame(Nmil = NN, Fecundity = FC_mean,
                 FC_lo = FC_lo, FC_hi = FC_hi)
pl6 = ggplot(dp6, aes(x=Nmil,y=Fecundity)) +
  geom_ribbon(aes(ymin=FC_lo,ymax=FC_hi),alpha=0.3) +
  geom_line() + labs(x = "Abundance (millions), w/o pups",y="Pregnancy rate (8+)") +
  ggtitle("Density-dependent variation in adult fecundity") + theme_classic()
print(pl6)
#
# Juvenile survival vs density --------------------------------------------
NN = seq(.2,10,by=0.1)
aJr = Jvloghz 
phiJ_r = mcmc[iir,vn=="phiJ"]
SJ = matrix(0,nrow = 1000,ncol=length(NN))
for(r in 1:1000){
  haz_J = exp(gamm0 + aJr + ( phiJ_r[r] *NN)^thta_r[r])  
  SJ[r,] = exp(-1 * (haz_J + exp(gamm0+.5) + exp(gamm0)))
}
SJ_mean = colMeans(SJ)
SJ_lo = apply(SJ, 2, quantile, prob=0.025)
SJ_hi = apply(SJ, 2, quantile, prob=0.975)
dp7 = data.frame(Nmil = NN, S_j = SJ_mean,
                 S_lo = SJ_lo, S_hi = SJ_hi)
pl7 = ggplot(dp7, aes(x=Nmil,y=SJ_mean)) +
  geom_ribbon(aes(ymin=SJ_lo,ymax=S_hi),alpha=0.3) +
  geom_line() + labs(x = "Abundance (millions), w/o pups",y="Juvenile survival rate") +
  ggtitle("Density-dependent variation in juvenile survival") + theme_classic()
print(pl7)
#
# Evaluate model sims-----------------------------------------------------
# Simulate future data with or without harvest mort
futuresim = 1; # 0 = current level of harvest, 1 = no harvest (for estimating K)
Nyrs2 = 50; Yearst2 = 2020 ; reps = 500
PAmeans = c(.18,.07,.75) # future proportion pups in S Gulf, N Gulf, Front
# Future conditions: sample ice and CE indices from after year YY 
YY = 1969
ii = which(df.CE$Year>=YY)
CE2 = log(df.CE$CEindex[sample(ii,1000,replace=T)])
ii = which(df.Ice$Year>=YY)
IC2 = as.matrix(cbind(df.Ice$Gulf_Anom[sample(ii,1000,replace=T)],
                      df.Ice$Gulf_Anom[sample(ii,1000,replace=T)],
                      df.Ice$Lab_Anom[sample(ii,1000,replace=T)]))
#
N_end = sumstats[which(vns==paste0("N[",Nyrs,"]")),1]
stan.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps,
                  PAmeans=PAmeans,futuresim=futuresim,
                  NFages=NFages,NFage1=NFage1,Nages=Nages,Nareas=Nareas,
                  ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CE=CE2,DDadlt=DDadlt,
                  b0pri=b0pri,Adloghz=Adloghz,Jvloghz=Jvloghz,gamm0=gamm0) #thta=thta
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
source("HSmod_sim.r")
Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
rslt=HSmod_sim(init_fun,stan.data)
N_Predict = rslt$N_Predict
P_Predict = rslt$P_Predict
P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
Np_Predict = N_Predict + P_Predict
Nfin = colMeans(Np_Predict[(Nyrs2-10):Nyrs2,])
Kest = mean(Nfin)
Kest_sd = sd(Nfin)
Kest_CI = quantile(Nfin, prob=c(0.025,0.975))
dpNprd = data.frame(Year = Yearp2, N_pred_mean = rowMeans(Np_Predict),
                    N_pred_lo = apply(Np_Predict,1,quantile,prob=0.025),
                    N_pred_hi = apply(Np_Predict,1,quantile,prob=0.975))
if (futuresim == 0){
  titletxt = "Model projected abundance with current harvest levels"
  subtxt = " "
}else if (futuresim == 1){
  titletxt = "Model projected abundance with zero harvest (including pups)"
  subtxt =  paste0("Estimated K = ", format(Kest/1000000,digits = 4),
                            "million (", format(Kest_CI[1]/1000000,digits = 4),
                            " - ", format(Kest_CI[2]/1000000,digits = 4),")")
}
plNprd = ggplot(data=dpNprd,aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Projected abundance") +
  ggtitle(titletxt,subtitle =subtxt) + theme_classic()
print(plNprd)

