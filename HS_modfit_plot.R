# Script to plot some results from model fitting
# Load results file (if not already loaded into workspace):---
if( !exists("sumstats") ){
  loadfile = file.choose(new = FALSE)
  load(file=loadfile)
  #
}
require(readxl)
require(stats)
require(gtools)
require(fitdistrplus)
require(dplyr)
#require(lattice)
#require(coda)
require(ggplot2)
require(bayesplot)
library(foreach)
library(doParallel)
library(doRNG)
Year = seq(Year1,Year1+Nyrs-2)
Yearp = seq(Year1,Year1+Nyrs-1)
# Diagnostics -----------------------------------------------
#
# Posterior predictive check for pup counts
y = stan.data$Pups
rr = sample(nrow(mcmc),100,replace = T)
yrep = mcmc[rr, startsWith(vn, "y_new1[")==T]
ppc_dens_overlay(y, yrep, adjust=1.2) +
  labs(x = "Pup count",y="Relative frequency") +
  ggtitle("Posterior predictive check, pup counts",
          subtitle=" Distribution of observed (y) vs. out-of-sample (y_rep) predictions")

# Posterior predictive checks for proportion of females pregnant
ii = which(stan.data$NFsamp>30)
y = stan.data$NPrg[ii]
z = stan.data$NFsamp[ii]
yrep = mcmc[, startsWith(vn, "y_new2[")==T]
yrep = yrep[,ii]
prop_preg <- function(x) mean(x/z)
ppc_stat(y, yrep, stat = "prop_preg") +
  labs(x = "Mean proportion females pregnant",y="Relative frequency") +
  ggtitle("Posterior predictive check, reproductive data",
          subtitle=" Distribution of observed (y) vs. out-of-sample (y_rep) predictions")
#
# Pop trends ----------------------------------------------------------------
dp1 = data.frame(Year=Yearp,Nexp=sumstats[startsWith(vns,"N[")==T,1],
                 N_lo = sumstats[startsWith(vns,"N[")==T,4],
                 N_hi = sumstats[startsWith(vns,"N[")==T,8])
dp1$N2exp = c(dp1$Nexp[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,1],NA)
dp1$N2_lo = c(dp1$N_lo[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,4],NA)
dp1$N2_hi = c(dp1$N_hi[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,8],NA)
tmp = read_excel("./data/OldMod_N.xlsx")
dp1$N3exp = tmp$N3exp
dp1$N3_lo = tmp$N3_lo
dp1$N3_hi = tmp$N3_hi
tmp = read_excel("./data/OldMod_N2.xlsx")
dp1$N4exp = tmp$N3exp
dp1$N4_lo = tmp$N3_lo
dp1$N4_hi = tmp$N3_hi
dp1 = dp1[-Nyrs,]
dp1 = dp1[-1,]
pl1 = ggplot(data=dp1,aes(x=Year,y=Nexp)) +
      geom_ribbon(aes(ymin=N_lo,ymax=N_hi),alpha=0.3) +
      geom_line() + labs(x = "Year",y="Estimated abundance w/o pups") +
      ggtitle("Model estimated abundance (excluding pups)") + theme_bw()
print(pl1)
pl1b = ggplot(data=dp1,aes(x=Year,y=N2exp)) +
  geom_ribbon(aes(ymin=N2_lo,ymax=N2_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Estimated abundance with pups") +
  geom_line(aes(y=N3exp),linetype=2,colour = "blue",size=1.1) +
  geom_line(aes(y=N3_lo),linetype=3,colour = "blue",size=1.1) +
  geom_line(aes(y=N3_hi),linetype=3,colour = "blue",size=1.1) +  
  geom_line(aes(y=N4exp),linetype=2,colour = "red",size=1.1) +
  geom_line(aes(y=N4_lo),linetype=3,colour = "red",size=1.1) +
  geom_line(aes(y=N4_hi),linetype=3,colour = "red",size=1.1) +  
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  scale_y_continuous(breaks = seq(1000000,9000000,by=1000000)) +
  ggtitle("Model estimated abundance (with pups)",
    subtitle = "Deterministic est. for comparison: ‘2017 survey cv 50%’ (blue) and ‘2017 survey exclude 2014 +rpd data’ (red)") + 
  theme_bw() 
print(pl1b)
#
# Pup counts ----------------------------------------------------------------
dp2 = data.frame(Year=Year,Pexp = sumstats[startsWith(vns,"Pups_pred[")==T,1],
                 P_lo = sumstats[startsWith(vns,"Pups_pred[")==T,4],
                 P_hi = sumstats[startsWith(vns,"Pups_pred[")==T,8])
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
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  scale_y_continuous(breaks = seq(0,2100000,by=100000)) +  
  ggtitle("Model estimated vs observed pup counts")  + theme_bw()
print(pl2)
#
# Preg rate by age----------------------------------------------------------
PR1966 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][16]
PR1966per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][1:32],na.rm = T)
crct1 = PR1966per/PR1966
PR2016 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][66]
PR2016per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][50:69],na.rm = T)
crct2 = PR2016per/PR2016

dp3 = data.frame(Age=ages, Year = (rep(1966,Nages)),Period=rep("1951-1985",Nages),
                 PRexp=sumstats[startsWith(vns,"Fc1966_prdct[")==T,1]*crct1,
                 PR_lo = sumstats[startsWith(vns,"Fc1966_prdct[")==T,4]*crct1,
                 PR_hi = sumstats[startsWith(vns,"Fc1966_prdct[")==T,8]*crct1)
dp3 = rbind(dp3, data.frame(Age=ages, Year = (rep(2016,Nages)),
                 Period=rep("1990-2019",Nages),           
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
dp3$Year = as.factor(dp3$Year); dp3$Period = as.factor(dp3$Period)
dp3 = dp3[which(dp3$Age<11),]
pl3 = ggplot(data=dp3,aes(x=Age,y=PRexp,group=Period,color = Period, fill=Period)) +
  geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi),alpha=0.3) + ylim(0,1) +
  geom_line() + labs(x = "Age",y="Pregnancy rate") +
  geom_point(aes(y=Obs),size=2) +
  geom_errorbar(aes(ymin = Obs-1.96*ObsSE, 
                      ymax = Obs+1.96*ObsSE),width=.2) +
  ggtitle("Model estimated vs observed pregancy rates by age",
          subtitle = "Low density (1951-1985) vs. high density (1990-2019) population") + 
  # scale_color_discrete(palette="Harmonic") +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  theme_classic() 
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
                                            ymax = Obs+1.96*ObsSE),
                                            color="darkgrey") +
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  ggtitle("Model estimated vs observed pregnancy rate over time",
          subtitle = " (points with error bars show sampled proportions and std. errors)") +
  theme_classic()
print(pl4)
#
# Harvest mortality over time ---------------------------------------------
dp5 = data.frame(Year=Year, Group = rep("Beaters",length(Year)),
                 HVexp=sumstats[startsWith(vns,"H0_pred[")==T,1],
                 HV_lo = sumstats[startsWith(vns,"H0_pred[")==T,4],
                 HV_hi = sumstats[startsWith(vns,"H0_pred[")==T,8])
dp5 = rbind(dp5, data.frame(Year=Year, Group = rep("Adult",length(Year)),
                 HVexp=sumstats[startsWith(vns,"HA_pred[")==T,1],
                 HV_lo = sumstats[startsWith(vns,"HA_pred[")==T,4],
                 HV_hi = sumstats[startsWith(vns,"HA_pred[")==T,8]))
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
  ggtitle("Model estimated vs observed harvest/bycatch mortality, by age class") + 
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  theme_classic()
print(pl5)
#
# Fecundity stochasticity -------------------------------------------------------------
b0 = sumstats[which(vns=="beta1"),1]
dp7 = data.frame(Year=Yearp[-length(Yearp)],epsF=sumstats[startsWith(vns,"epsF[")==T,1],
                  epsF_lo = sumstats[startsWith(vns,"epsF[")==T,4],
                  epsF_hi = sumstats[startsWith(vns,"epsF[")==T,8])
dp7$Fdev = inv.logit(b0 + dp7$epsF) - inv.logit(b0)
dp7$Fdev_lo = inv.logit(b0 + dp7$epsF_lo) - inv.logit(b0)
dp7$Fdev_hi = inv.logit(b0 + dp7$epsF_hi) - inv.logit(b0)
pl7 = ggplot(data=dp7,aes(x=Year,y=Fdev)) +
  geom_ribbon(aes(ymin=Fdev_lo,ymax=Fdev_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Deviation from expected F (8+)") +
  ggtitle("Stochastic variation in fecundity") + theme_bw()
print(pl7)
#
# Juv survival stochasticity -------------------------------------------------------------
gamma_0 = sumstats[which(vns=="gamma_0"),1]
hazImn = sumstats[which(startsWith(vns,"haz_Ice[")),1][21]  
phiS = sumstats[which(vns=="phi[2]"),1]
thtaS = sumstats[which(startsWith(vns,"thta[2]")),1]
Nml =  sumstats[which(startsWith(vns,"N[")),1]; Nml = Nml[-Nyrs]/1000000
hzoth = hazImn + exp(omega)
dp8 = data.frame(Year=Yearp[-length(Yearp)],epsS=sumstats[startsWith(vns,"epsS[")==T,1],
                  epsS_lo = sumstats[startsWith(vns,"epsS[")==T,4],
                  epsS_hi = sumstats[startsWith(vns,"epsS[")==T,8])
dp8$S0_stoc = exp(-1 * (hzoth + exp(omega + gamma_0 + (phiS*Nml)^thtaS + dp8$epsS)))
dp8$S0_lo =  exp(-1 * (hzoth + exp(omega + gamma_0 + (phiS*Nml)^thtaS + dp8$epsS_lo)))
dp8$S0_hi =  exp(-1 * (hzoth + exp(omega + gamma_0 + (phiS*Nml)^thtaS + dp8$epsS_hi)))
pl8 = ggplot(data=dp8,aes(x=Year,y=S0_stoc)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Juvenile survival") +
  ggtitle("Stoachastic and D-D variation in juvenile survival") + theme_bw()
print(pl8)
#
pl8b = ggplot(data=dp8,aes(x=Year,y=epsS)) +
  geom_ribbon(aes(ymin=epsS_lo,ymax=epsS_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Juvenile log hazards, deviation from mean") +
  ggtitle("Stoachastic variation in juvenile mortality") + theme_bw()
print(pl8b)

# Age-specific survival -----------------------------------------------------------
tmp = data.frame(Density = "Low", Age = ages)
tmp$Survival = sumstats[which(startsWith(vns,"SA_ld[")),1]
tmp$Survival_lo = sumstats[which(startsWith(vns,"SA_ld[")),4]
tmp$Survival_hi = sumstats[which(startsWith(vns,"SA_ld[")),8]
dp9 = rbind(data.frame(Density = "Low", Age = 0,
                        Survival = sumstats[which(vns =="S0_ld"),1],
                        Survival_lo = sumstats[which(vns =="S0_ld"),4],
                        Survival_hi = sumstats[which(vns =="S0_ld"),8]), tmp)
tmp = data.frame(Density = "High", Age = ages)
tmp$Survival = sumstats[which(startsWith(vns,"SA_hd[")),1]
tmp$Survival_lo = sumstats[which(startsWith(vns,"SA_hd[")),4]
tmp$Survival_hi = sumstats[which(startsWith(vns,"SA_hd[")),8]
dp9 = rbind(dp9, data.frame(Density = "High", Age = 0,
                        Survival = sumstats[which(vns =="S0_hd"),1],
                        Survival_lo = sumstats[which(vns =="S0_hd"),4],
                        Survival_hi = sumstats[which(vns =="S0_hd"),8]), tmp)
dp9 = dp9[-which(dp9$Age==Nages),]
dp9$Density = factor(dp9$Density, levels = c("Low","High"))
pl9 = ggplot(data=dp9,aes(x=Age,y=Survival,group=Density,fill=Density)) +
  geom_ribbon(aes(ymin=Survival_lo,ymax=Survival_hi),alpha=0.2) +
  scale_fill_manual(name = "Abundance",labels = c("2 million","6 million"),
                     values = c("blue","red")) +  
  geom_line(aes(color=Density),show.legend = FALSE) + 
  scale_color_manual(values = c("blue","red")) +
  labs(x = "Age (0 = juvenile)",y="Annual survival rate") +
  ylim(c(floor(10*min(dp9$Survival_lo))/10 ,1)) +
  ggtitle("Age-specific variation in survival at low and high population density") + theme_bw()
print(pl9)
#
# Ice anomaly effect on pup survival---------------------------------------
dp10 = data.frame(Ice_Anomaly = ICvec)
dp10$Haz = sumstats[which(startsWith(vns,"haz_Ice[1")),1]
dp10$Haz_lo = sumstats[which(startsWith(vns,"haz_Ice[1")),4]
dp10$Haz_hi = sumstats[which(startsWith(vns,"haz_Ice[1")),8]
dp10$SvIce = exp(-1 * (dp10$Haz))
dp10$SvIce_lo = exp(-1 * (dp10$Haz_hi))
dp10$SvIce_hi = exp(-1 * (dp10$Haz_lo))
pl10 = ggplot(dp10,aes(x=Ice_Anomaly,y=SvIce)) +
  geom_ribbon(aes(ymin=SvIce_lo,ymax=SvIce_hi),alpha=0.3) +
  geom_line() + geom_vline(xintercept = 0) +
  labs(x="Ice Anomaly (deviation from 1969-2000 mean cover)",y="Pup survival from ice hazards") +
  ggtitle("Effect of ice cover on pup survival (Gulf)") +
  theme_classic()
print(pl10)
#
# Fecundity vs density ----------------------------------------------------
NN = seq(.1,6.6,by=0.1)
lngNN = length(NN)
iir = sample(Nsims,1000)
b0_r = mcmc[iir,vn=="beta1"]
phiF_r = mcmc[iir,vn=="phi[1]"]
thtaF_r = mcmc[iir,vn=="thta[1]"]
FC = matrix(0,nrow = 1000,ncol=length(NN))
FC_mean = numeric(length = lngNN)
FC_lo = numeric(length = lngNN)
FC_hi = numeric(length = lngNN)
for(r in 1:1000){
  FC[r,] = inv.logit(b0_r[r] - (phiF_r[r]*NN)^thtaF_r[r])
}
FC_mean = colMeans(FC)
FC_lo = apply(FC, 2, quantile, prob=0.055)
FC_hi = apply(FC, 2, quantile, prob=0.95)
dp11 = data.frame(Nmil = NN*1.2, Fecundity = FC_mean,
                 FC_lo = FC_lo, FC_hi = FC_hi)
pl11 = ggplot(dp11, aes(x=Nmil,y=Fecundity)) +
  geom_ribbon(aes(ymin=FC_lo,ymax=FC_hi),alpha=0.3) +
  geom_line() + labs(x = "Abundance (millions)",y="Pregnancy rate (8+)") +
  ggtitle("Density-dependent variation in adult fecundity") + theme_classic()
print(pl11)
#
# Juvenile survival vs density --------------------------------------------
NN = seq(.1,6.6,by=0.1)
gamma0_r = mcmc[iir,vn=="gamma_0"]
phiS_r = mcmc[iir,vn=="phi[2]"]
thtaS_r = mcmc[iir,vn=="thta[2]"]
hazImn = sumstats[which(startsWith(vns,"haz_Ice[1")),1][21]  
S0 = matrix(0,nrow = 1000,ncol=length(NN))
for(r in 1:1000){
  haz_J = exp(omega + gamma0_r[r] + ( phiS_r[r] *NN)^thtaS_r[r])  
  S0[r,] = exp(-1 * (haz_J + hazImn + exp(omega)))
}
S0_mean = colMeans(S0)
S0_lo = apply(S0, 2, quantile, prob=0.05)
S0_hi = apply(S0, 2, quantile, prob=0.95)
dp12 = data.frame(Nmil = NN*1.2, S0_mean = S0_mean,
                 S0_lo = S0_lo, S0_hi = S0_hi)
pl12 = ggplot(dp12, aes(x=Nmil,y=S0_mean)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  xlim(0,8) +
  geom_line() + labs(x = "Abundance (millions)",y="Juvenile survival rate") +
  ggtitle("Density-dependent variation in juvenile survival") + theme_classic()
print(pl12)
#
# Evaluate model sims-----------------------------------------------------
# Simulate future data with or without harvest mort
futuresim = 0; # 0 = past harvest, 1 = no harvest, 2 = evaluate range of harvests
Nyrs2 = 100; Yearst2 = 2021 ; reps = 5000
PAmeans = c(.18,.07,.75) # future proportion pups in S Gulf, N Gulf, Front
set.seed(123)
r_vec = sample(nrow(mcmc),reps,replace = T)
# Future conditions: based on ice and CE indices from after year YY, 
#  fit appropriate random sampling distributions
if(futuresim==0){
  Yearst2=Year1 
  Nyrs2=Nyrs
  N_end = sumstats[vns=="N0",1]
  IC2=IC
  CE2=CE
}else{
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
}
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
sim.data <- list(Nyrs=Nyrs2,N0pri=N_end,reps=reps,r_vec=r_vec,
                  PAmeans=PAmeans,futuresim=futuresim,
                  NCages=NCages,NCage1=NCage1,Nages=Nages,Nareas=Nareas,
                  ages=ages,ages2=ages2,sad0=sad0,IC=IC2,CE=CE2,
                  omega=omega,Grn_A=Grn_A,Grn_P=Grn_P,Arc_A=Arc_A,Arc_P=Arc_P,
                  Byc_A=Byc_A,Byc_P=Byc_P) 
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
  dlta = mcmc[rr,vn=="dlta"],
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
rslt=HSmod_sim(init_fun,sim.data,mcmc,vn)
Yearp2 = seq(Yearst2,Yearst2+Nyrs2-1)
N_Predict = rslt$N_Predict
P_Predict = rslt$P_Predict
P_Predict = rbind(P_Predict,colMeans(P_Predict[(Nyrs2-4):(Nyrs2-1),]))
Np_Predict = (N_Predict + P_Predict)
Nfin = colMeans(Np_Predict[(Nyrs2-20):Nyrs2,])
Kest = mean(Nfin)
Kest_sd = sd(Nfin)
Kest_CI = quantile(Nfin, prob=c(0.025,0.975))
dpNprd = data.frame(Year = Yearp2, N_pred_mean = rowMeans(Np_Predict),
                    N_pred_lo = apply(Np_Predict,1,quantile,prob=0.025),
                    N_pred_hi = apply(Np_Predict,1,quantile,prob=0.975))
if (futuresim == 0){
  titletxt = "Hindcast stochastic simulations of abundance with observed harvest levels"
  subtxt = " (red line = actual model estimate of abundance) "
}else if (futuresim == 1){
  titletxt = "Model projected abundance with zero harvest (including pups)"
  subtxt =  paste0("Estimated long-term equilibrium (K) = ", format(Kest/1000000,digits = 3),
                            "million (CI95 = ", format(Kest_CI[1]/1000000,digits = 3),
                            " - ", format(Kest_CI[2]/1000000,digits = 3),")")
}
plNprd = ggplot(data=dpNprd,aes(x=Year,y=N_pred_mean)) +
  geom_ribbon(aes(ymin=N_pred_lo,ymax=N_pred_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Projected abundance") +
  ggtitle(titletxt,subtitle =subtxt) + theme_classic()
if (futuresim == 0){
  Predpup = sumstats[startsWith(vns,"Pups_pred[")==T,1]; Predpup = c(Predpup,Predpup[Nyrs-1])
  dpNprd$N_actual = sumstats[startsWith(vns,"N[")==T,1] + Predpup
  plNprd = plNprd + geom_line(data = dpNprd,aes(y=N_actual),col="red")
}
print(plNprd)
#
# SAD = data.frame(Age = ages, SAD = rslt$SAD) 
# write.csv(SAD,file="./data/SAD0_36.csv",row.names = F)
#
stopCluster(cl)
