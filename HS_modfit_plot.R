# Script to plot some results from model fitting
# Load results file (if not already loaded into workspace):---
require(readxl)
require(stats)
require(gtools)
require(fitdistrplus)
require(dplyr)
require(ggplot2)
require(bayesplot)
library(foreach)
library(doParallel)
library(doRNG)
library(gridExtra)
require(svDialogs)
require(rJava)
require(rChoiceDialogs)
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
if( !exists("sumstats") ){
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
                 N_sd = sumstats[startsWith(vns,"N[")==T,3],
                 N_lo = sumstats[startsWith(vns,"N[")==T,4],
                 N_hi = sumstats[startsWith(vns,"N[")==T,8])
dp1$N2exp = c(dp1$Nexp[1:(Nyrs-1)] + sumstats[startsWith(vns,"Pups_pred[")==T,1],NA)
dp1$N2_sd = c(sqrt( (dp1$N_sd[1:(Nyrs-1)])^2 + 
                sumstats[startsWith(vns,"Pups_pred[")==T,3]^2),NA)
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
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  scale_y_continuous(breaks = seq(1000000,9000000,by=1000000)) +
  ggtitle("Model estimated abundance (with pups)",
    subtitle = "Deterministic est. for comparison: ‘2017 survey cv 50%’ (blue)") + 
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
dp10 = data.frame(Ice_Anomaly = ICvec,Area = "Gulf")
dp10$Haz = sumstats[which(startsWith(vns,"haz_Ice[1")),1]
dp10$Haz_lo = sumstats[which(startsWith(vns,"haz_Ice[1")),4]
dp10$Haz_hi = sumstats[which(startsWith(vns,"haz_Ice[1")),8]
dp10$SvIce = exp(-1 * (dp10$Haz))
dp10$SvIce_lo = exp(-1 * (dp10$Haz_hi))
dp10$SvIce_hi = exp(-1 * (dp10$Haz_lo))
tmp = data.frame(Ice_Anomaly = ICvec,Area = "Front")
tmp$Haz = sumstats[which(startsWith(vns,"haz_Ice[2")),1]
tmp$Haz_lo = sumstats[which(startsWith(vns,"haz_Ice[2")),4]
tmp$Haz_hi = sumstats[which(startsWith(vns,"haz_Ice[2")),8]
tmp$SvIce = exp(-1 * (tmp$Haz))
tmp$SvIce_lo = exp(-1 * (tmp$Haz_hi))
tmp$SvIce_hi = exp(-1 * (tmp$Haz_lo)) 
dp10 = rbind(dp10,tmp); dp10$Area = factor(dp10$Area, levels=c("Gulf","Front"))
pl10 = ggplot(dp10,aes(x=Ice_Anomaly,y=SvIce,fill=Area,color=Area)) +
  geom_ribbon(aes(ymin=SvIce_lo,ymax=SvIce_hi),alpha=0.25) +
  geom_line() +
  geom_vline(xintercept = 0) +
  labs(x="Ice Anomaly (proportional deviation from 1969-2000 mean cover)",y="Pup survival from ice hazards") +
  ggtitle("Effect of ice cover on pup survival") +
  theme_classic()
print(pl10)
#
# Effects of Climate Index---------------------------------------
NN = 3.75
CIvals = seq(-1.5,1.5,by = .1)
NCIvals = length(CIvals)
iir = sample(Nsims,1000)
b0_r = as.numeric(mcmc[iir,vn=="beta1"])
phiF_r = as.numeric(mcmc[iir,vn=="phi[1]"])
thtaF_r = as.numeric(mcmc[iir,vn=="thta[1]"])
phiS_r = as.numeric(mcmc[iir,vn=="phi[2]"])
thtaS_r = as.numeric(mcmc[iir,vn=="thta[2]"])
gamma0_r = as.numeric(mcmc[iir,vn=="gamma_0"])
dlta1_r =  as.numeric(mcmc[iir,vn=="dlta[1]"])
dlta2_r =  as.numeric(mcmc[iir,vn=="dlta[2]"])
FC = matrix(0,nrow = 1000,ncol=NCIvals)
S0 = matrix(0,nrow = 1000,ncol=NCIvals)
gamD_F = (mean(phiF_r)*NN)^mean(thtaF_r) 
gamD_S = (mean(phiS_r)*NN)^mean(thtaS_r)
for(r in 1:1000){
  FC[r,] = inv.logit(b0_r[r] - (phiF_r[r]*NN)^thtaF_r[r] - dlta1_r[r]*CIvals)
  hazJ = exp(omega + gamma0_r[r] + (phiS_r[r]*NN)^thtaS_r[r] + dlta2_r[r] * CIvals ) 
  S0[r,] = exp(-hazJ + 2*exp(omega)) 
}
FC_mean = colMeans(FC)
FC_lo = apply(FC, 2, quantile, prob=0.05)
FC_hi = apply(FC, 2, quantile, prob=0.95)
S0_mean = colMeans(S0)
S0_lo = apply(S0, 2, quantile, prob=0.05)
S0_hi = apply(S0, 2, quantile, prob=0.95)
dp11 = data.frame(CI = CIvals, Fecundity = FC_mean,
                  FC_lo = FC_lo, FC_hi = FC_hi,
                  Survival = S0_mean,
                  S0_lo = S0_lo, S0_hi = S0_hi)
pl11a = ggplot(dp11, aes(x=CI,y=Fecundity)) +
  geom_ribbon(aes(ymin=FC_lo,ymax=FC_hi),alpha=0.3) +
  geom_line() + labs(x = "Climate index (CI)",y="Pregnancy rate (8+)") +
  ggtitle("A) Effect of climate index on adult fecundity") + theme_classic()
#print(pl11a)

pl11b = ggplot(dp11, aes(x=CI,y=Survival)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  geom_line() + labs(x = "Climate index (CI)",y="Annual survival, ageclass 0") +
  ggtitle("B) Effect of climate index on juvenile survival") + theme_classic()
# print(pl11b)

grid.arrange(pl11a,pl11b)
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
FC_lo = apply(FC, 2, quantile, prob=0.05)
FC_hi = apply(FC, 2, quantile, prob=0.95)
dp12 = data.frame(Nmil = NN*1.2, Fecundity = FC_mean,
                 FC_lo = FC_lo, FC_hi = FC_hi)
pl12 = ggplot(dp12, aes(x=Nmil,y=Fecundity)) +
  geom_ribbon(aes(ymin=FC_lo,ymax=FC_hi),alpha=0.3) +
  geom_line() + labs(x = "Abundance (millions)",y="Pregnancy rate (8+)") +
  ggtitle("Density-dependent variation in adult fecundity") + theme_classic()
print(pl12)
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
dp13 = data.frame(Nmil = NN*1.2, S0_mean = S0_mean,
                 S0_lo = S0_lo, S0_hi = S0_hi)
pl13 = ggplot(dp13, aes(x=Nmil,y=S0_mean)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  xlim(0,8) +
  geom_line() + labs(x = "Abundance (millions)",y="Juvenile survival rate") +
  ggtitle("Density-dependent variation in juvenile survival") + theme_classic()
print(pl13)
#
