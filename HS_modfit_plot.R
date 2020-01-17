# Script to plot some results from model fitting
load(file="./Results/FitHSmod3_Results_Jan16.rdata")
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
pl1 = ggplot(data=dp1,aes(x=Year,y=Nexp)) +
      geom_ribbon(aes(ymin=N_lo,ymax=N_hi),alpha=0.3) +
      geom_line() + labs(x = "Year",y="Estimated abundance") +
      ggtitle("Model estimated abundance") + theme_classic()
print(pl1)
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
pl2 = ggplot(data=dp2,aes(x=Year,y=Pexp)) +
  geom_ribbon(aes(ymin=P_lo,ymax=P_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Pup counts (total)") +
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-ObsSE, ymax = Obs+ObsSE)) +
  ggtitle("Model estimated vs observed pup counts") + theme_classic()
print(pl2)
# Preg rate by age----------------------------------------------------------
dp3 = data.frame(Age=ages, Year = (rep(1966,Nages)),
                 PRexp=sumstats[startsWith(vns,"Fc1966_prdct[")==T,1],
                 PR_lo = sumstats[startsWith(vns,"Fc1966_prdct[")==T,4],
                 PR_hi = sumstats[startsWith(vns,"Fc1966_prdct[")==T,8])
dp3 = rbind(dp3, data.frame(Age=ages, Year = (rep(2016,Nages)),
                 PRexp=sumstats[startsWith(vns,"Fc2016_prdct[")==T,1],
                 PR_lo = sumstats[startsWith(vns,"Fc2016_prdct[")==T,4],
                 PR_hi = sumstats[startsWith(vns,"Fc2016_prdct[")==T,8]))
dp3$Obs = numeric(length = nrow(dp3))
for (i in 1:nrow(dp3)){
  ii = which(df.Rep$Age==dp3$Age[i] & df.Rep$Year==dp3$Year[i] & df.Rep$Year==dp3$Year[i])
  dp3$Obs[i] = median(df.Rep$Prob[ii])
}
dp3$Year = as.factor(dp3$Year)
pl3 = ggplot(data=dp3,aes(x=Age,y=PRexp,group=Year,color = Year, fill=Year)) +
  geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi),alpha=0.3) + ylim(0,1) +
  geom_line() + labs(x = "Age",y="Pregnancy rate") +
  geom_point(aes(y=Obs),size=2) +
  ggtitle("Model estimated vs observed pregancy rate, by age") + theme_classic()
print(pl3)
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
  geom_line() + labs(x = "Year",y="Preganancy rate (AGe 8+)") +
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-ObsSE, ymax = Obs+ObsSE)) +
  ggtitle("Model estimated vs observed pregnancy rate over time") + theme_classic()
print(pl4)
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
  ggtitle("Model estimated vs observed harvest mortality, by age class") + theme_classic()
print(pl5)

# Ice mortality -----------------------------------------------------------
psi1mn = sumstats[startsWith(vns,"psi1")==T,1]  
psi1sd = sumstats[startsWith(vns,"psi1")==T,2]  
ps2mn = sumstats[startsWith(vns,"psi2")==T,1]  
ps2sd = sumstats[startsWith(vns,"psi2")==T,2]  
source("Ice_mort_plot.r")
Ice_mort_plot(psipri1a,psipri1b,psipri2a,psipri2b)

# Juvenile survival vs density --------------------------------------------
NN = seq(.5,12,by=0.5)
phiJ_mn = sumstats[which(startsWith(vns,"phiJ")),1]
phiJ_sd = sumstats[which(startsWith(vns,"phiJ")),3]
phiJr = rnorm(1000,phiJ_mn,phiJ_sd)
#aJ_mn = sumstats[which(startsWith(vns,"aJ")),1]
#aJ_sd = sumstats[which(startsWith(vns,"aJ")),3]
#aJr = rnorm(1000,aJ_mn,aJ_sd)
#thta_mn = sumstats[which(startsWith(vns,"thta")),1]
#thta_sd = sumstats[which(startsWith(vns,"thta")),3]
#thta_r = rnorm(1000,thta_mn,thta_sd)
aJr = Jvhzpri # replace with vector from fitted estimate
thta_r = thta # replace with vector from fitted estimate
SJ = matrix(0,nrow = 1000,ncol=length(NN))
for(r in 1:1000){
  haz_J = exp(gamm0 + aJr + ( phiJr[r] *NN)^thta_r)  
  SJ[r,] = exp(-1 * (haz_J + exp(gamm0+1) + exp(gamm0)))
}
SJ_mean = colMeans(SJ)
SJ_lo = apply(SJ, 2, quantile, prob=0.025)
SJ_hi = apply(SJ, 2, quantile, prob=0.975)
dp6 = data.frame(Nmil = NN, S_j = SJ_mean,
                 S_lo = SJ_lo, S_hi = SJ_hi)
pl6 = ggplot(dp6, aes(x=Nmil,y=SJ_mean)) +
  geom_ribbon(aes(ymin=SJ_lo,ymax=S_hi),alpha=0.3) +
  geom_line() + labs(x = "Abundance (millions)",y="Juvenile survival rate") +
  ggtitle("Density-dependent variation in juvenile survival") + theme_classic()
print(pl6)
#
# Evaluate model sims-----------------------------------------------------
# Simulate "new data sets" with or without harvest mort
futuresim = 1; Nyrs2 = 50; Yearst2 = 2020
IC2 = matrix(0,nrow = Nyrs2,ncol=3); CE2 = rep(0,Nyrs2)
#
stan.data <- list(Nyrs=Nyrs2,N0pri=7000000,futuresim=futuresim,reps=500,
                  NFages=NFages,NFage1=NFage1,Nages=Nages,Nareas=Nareas,
                  ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CE=CE2,DDadlt=DDadlt,
                  b0=b0,Adhzpri=Adhzpri,Jvhzpri=Jvhzpri,gamm0=gamm0,PApri=PApri) # thta=thta
init_fun <- function() {list(sigF=rnorm(1,sumstats[vns=="sigF",1],sumstats[vns=="sigF",2]),
                             sigH=rnorm(1,sumstats[vns=="sigH",1],sumstats[vns=="sigH",2]),
                             phiJ=rnorm(1,sumstats[vns=="phiJ",1],sumstats[vns=="phiJ",2]),
                             phiF=rnorm(1,sumstats[vns=="phiF",1],sumstats[vns=="phiF",2]),
                             b1=rnorm(1,sumstats[vns=="b1",1],sumstats[vns=="b1",2]),
                             psi1=rnorm(1,sumstats[vns=="psi1",1],sumstats[vns=="psi1",2]),
                             psi2=rnorm(1,sumstats[vns=="psi2",1],sumstats[vns=="psi2",2]),
                             dlta=rnorm(1,sumstats[vns=="dlta",1],sumstats[vns=="dlta",2]),
                             gammHp_mn = 0, gammHa_mn = 0
                             # ADD AJ AND THTA
                             #gammHp_mn=rnorm(1,sumstats[vns=="gammHp_mn",1],sumstats[vns=="gammHp_mn",2]),
                             #gammHa_mn=rnorm(1,sumstats[vns=="gammHa_mn",1],sumstats[vns=="gammHa_mn",2])
)}
#
source("HSmod_test.r")
Year = seq(Yearst2,Yearst2+Nyrs2-2)
Yearp = seq(Yearst2,Yearst2+Nyrs2-1)
rslt=HSmod_test(init_fun,stan.data)
N_Predict = rslt$N_Predict
P_Predict = rslt$P_Predict
PrPredict = rslt$PrPredict
PrAgPred = rslt$PrAgPred
Agepredict1 = rslt$Agepredict1 # Age dist for sample year 17 = 1967 (Agecomp row5)
Agepredict2 = rslt$Agepredict2 # Age dist for sample year 69 = 2019 (Agecomp row48)
HVp_predict = rslt$HVp_predict
HVa_predict = rslt$HVa_predict
Nfin = N_Predict[nrow(N_Predict),]
Kest = mean(Nfin)
Kest_sd = sd(Nfin)
Kest_CI = quantile(Nfin, prob=c(0.025,0.975))

ggplot(data.frame(Year=Yearp,Npred = rowMeans(N_Predict)),aes(x=Year,y=Npred)) +
  geom_line()
ggplot(data.frame(Year=Year,Ppred = rowMeans(P_Predict)),aes(x=Year,y=Ppred)) +
  geom_line() + geom_point(data=df.Pup,aes(x=Year,y=Total_Npup))
ggplot(data.frame(Age=ages[3:8], Frequency= rowMeans(Agepredict1)),aes(x=Age,y=Frequency)) +
  geom_line() + geom_point(data=data.frame(Age=ages[3:8], Frequency= AgeComp[5,]/sum(AgeComp[5,])),aes(x=Age,y=Frequency))
ggplot(data.frame(Age=ages[3:8], Frequency= rowMeans(Agepredict2)),aes(x=Age,y=Frequency)) +
  geom_line() + geom_point(data=data.frame(Age=ages[3:8], Frequency= AgeComp[48,]/sum(AgeComp[48,])),aes(x=Age,y=Frequency))
ggplot(data.frame(Year=Year,Prt8yrPred = rowMeans(PrPredict)),aes(x=Year,y=Prt8yrPred)) +
  geom_line() + geom_point(data=df.Rep[df.Rep$Age==8,],aes(x=Year,y=Prob))
ggplot(data.frame(Age=ages, Pregrate= rowMeans(PrAgPred)),aes(x=Age,y=Pregrate)) +
  geom_line() + geom_point(data=df.Rep[df.Rep$Year<=1970,],aes(x=Age,y=Prob))
ggplot(data.frame(Year=Year,HvPpred = rowMeans(HVp_predict)),aes(x=Year,y=HvPpred)) +
  geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=PUPTOT))
ggplot(data.frame(Year=Year,HvApred = rowMeans(HVa_predict)),aes(x=Year,y=HvApred)) +
  geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=ADLTOT))


