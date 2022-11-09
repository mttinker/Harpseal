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
require(cowplot)
require(svDialogs)
require(tcltk)
#require(rJava)
#require(rChoiceDialogs)
require(kableExtra)
# Load data if needed --------------------------------------
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
if( !exists("sumstats") ){
  file_list = list.files(path = "./Results",pattern = "HS_Results",full.names = F)
  rslt_list = grep("TAC", file_list, value = TRUE, invert = TRUE)
  rslt_list = grep("Report", rslt_list, value = TRUE, invert = TRUE)
  rdata_file = tk_select.list(rslt_list, preselect = NULL, multiple = FALSE,
                              title = "Select results file" ) 
  if(length(rdata_file)==0){
    dlg_message(c("No data file selected"), "ok")
    stop_quietly()
  }else{
    load(paste0("./Results/",rdata_file))
  }
}
Year = seq(Year1,Year1+Nyrs-2)
Yearp = seq(Year1,Year1+Nyrs-1)
Nyrsm1 = Nyrs - 1
df.prev = read_excel("./data/Prev_model_project.xlsx")
plotpri = 0
# Diagnostics -----------------------------------------------
#
# Posterior predictive check for pup counts
y = stan.data$Pups
rr = sample(nrow(mcmc),100,replace = T)
yrep = mcmc[rr, startsWith(vn, "y_new1[")==T]
plt0a = ppc_dens_overlay(y, yrep, adjust=1.2) +
  labs(x = "Pup count",y="Relative frequency") 
  # ggtitle("Posterior predictive check, pup counts",
  #         subtitle=" Distribution of observed (y) vs. out-of-sample (y_rep) predictions")

# Posterior predictive checks for proportion of females pregnant
ii = which(stan.data$NFsamp>30)
y = stan.data$NPrg[ii]
z = stan.data$NFsamp[ii]
yrep = mcmc[, startsWith(vn, "y_new2[")==T]
yrep = yrep[,ii]
prop_preg <- function(x) mean(x/z)
plt0b = ppc_stat(y, yrep, stat = "prop_preg") +
  labs(x = "Mean proportion females pregnant",y="Relative frequency") 
  # ggtitle("Posterior predictive check, reproductive data",
  #         subtitle=" Distribution of observed (y) vs. out-of-sample (y_rep) predictions")
plot_grid(plt0a,plt0b,nrow = 2,labels = c("A","B"))
#
# Pop trends ----------------------------------------------------------------
dp1 = data.frame(Year=Yearp,Nexp=sumstats[startsWith(vns,"N[")==T,1],
                 N_sd = sumstats[startsWith(vns,"N[")==T,3],
                 N_lo = sumstats[startsWith(vns,"N[")==T,4],
                 N_hi = sumstats[startsWith(vns,"N[")==T,8])
dp1 = dp1[1:Nyrsm1,]
Ntot_reps = mcmc[,which(startsWith(vn,"N["))]
Pup_reps = mcmc[,which(startsWith(vn,"Pups_pred["))]
Ntot_reps = Ntot_reps[,1:Nyrsm1] + Pup_reps
dp1$NT_exp = apply(Ntot_reps,2,mean)
dp1$NT_sd = apply(Ntot_reps,2,sd)
dp1$NT_lo = apply(Ntot_reps,2,quantile,prob=0.025)
dp1$NT_hi = apply(Ntot_reps,2,quantile,prob=0.975)
dp1 = dp1[-1,]
dp1$Nprev_exp = df.prev$N_exp[1:(nrow(df.prev)-1)]
dp1$Nprev_lo = df.prev$N_exp_lo[1:(nrow(df.prev)-1)]
dp1$Nprev_hi = df.prev$N_exp_hi[1:(nrow(df.prev)-1)]
pl1 = ggplot(data=dp1,aes(x=Year,y=NT_exp)) +
      geom_ribbon(aes(ymin=NT_lo,ymax=NT_hi),alpha=0.3) +
      geom_line() + labs(x = "Year",y="Estimated abundance") +
      ggtitle("Model estimated abundance (including pups)") + theme_bw()
# print(pl1)
pl1b = ggplot(data=dp1,aes(x=Year)) +
  geom_ribbon(aes(ymin=Nprev_lo,ymax=Nprev_hi,fill="existing"),alpha=0.2) +
  geom_ribbon(aes(ymin=NT_lo,ymax=NT_hi,fill="revised"),alpha=0.3) +
  geom_line(aes(y=Nprev_exp,colour="existing"),linetype=2,size=1) +
  geom_line(aes(y=NT_exp,color="revised")) + 
  labs(x = "Year",y="Estimated abundance") +
  #geom_line(aes(y=Nprev_lo),linetype=3,colour = "blue",size=1.1) +
  #geom_line(aes(y=Nprev_hi),linetype=3,colour = "blue",size=1.1) +  
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  scale_y_continuous(breaks = seq(1000000,9000000,by=2000000)) +
  scale_color_manual(name = "Legend",
                     values = c("existing" = "blue", "revised" = "darkred"),
                     labels = c("existing model", "revised model"),
                     aesthetics = c("colour", "fill")) +
#  ggtitle("Model estimated abundance (including pups)",
#    subtitle = "Current model (red) vs. previous model (blue)") + 
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(legend.position = c(0.2, 0.8),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(pl1b)
#
write.excel <- function(x,row.names=T,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}
write.excel(dp1)

# Pup counts ----------------------------------------------------------------
dp2 = data.frame(Year=Year,Pexp = sumstats[startsWith(vns,"Pups_pred[")==T,1],
                 P_sd = sumstats[startsWith(vns,"Pups_pred[")==T,3],
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

dp2$Pprev_exp[2:Nyrsm1] = df.prev$P_exp[1:(nrow(df.prev)-1)]
dp2$Pprev_lo[2:Nyrsm1] = df.prev$P_exp_lo[1:(nrow(df.prev)-1)]
dp2$Pprev_hi[2:Nyrsm1] = df.prev$P_exp_hi[1:(nrow(df.prev)-1)]

pl2 = ggplot(data=dp2,aes(x=Year)) +
  geom_ribbon(aes(ymin=Pprev_lo,ymax=Pprev_hi,fill="existing"),alpha=0.3) +
  geom_ribbon(aes(ymin=P_lo,ymax=P_hi,fill="revised"),alpha=0.3) +  
  geom_line(aes(y=Pprev_exp,colour = "existing"),linetype=1) +  
  geom_line(aes(y = Pexp,color="revised")) + 
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-1.96*ObsSE, ymax = Obs+1.96*ObsSE)) +
  scale_x_continuous(breaks = seq(Year1-1,YearT,by=5)) +
  scale_y_continuous(breaks = seq(0,2100000,by=200000)) +  
  labs(x = "Year",y="Pup counts (total)") +
  scale_color_manual(name = "Legend",
                     values = c("existing" = "blue", "revised" = "darkred"),
                     labels = c("existing model", "revised model"),
                     aesthetics = c("colour", "fill")) +
#  ggtitle("Model estimated vs observed pup counts",
#          subtitle = "Current model (red) vs. previous model (blue)") +  
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(legend.position = c(0.2, 0.8),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(pl2)
#
write.excel <- function(x,row.names=T,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}
write.excel(dp2)
#
# Preg rate by age----------------------------------------------------------
# PR1966 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][16]
# PR1966per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][1:32],na.rm = T)
# crct1 = PR1966per/PR1966
# PR2016 = sumstats[startsWith(vns,"Fc8_prdct[")==T,1][66]
# PR2016per = mean(sumstats[startsWith(vns,"Fc8_prdct[")==T,1][50:69],na.rm = T)
# crct2 = PR2016per/PR2016
# 
# dp3 = data.frame(Age=ages, Year = (rep(1966,Nages)),Period=rep("1951-1985",Nages),
#                  PRexp=sumstats[startsWith(vns,"Fc1966_prdct[")==T,1]*crct1,
#                  PR_lo = sumstats[startsWith(vns,"Fc1966_prdct[")==T,4]*crct1,
#                  PR_hi = sumstats[startsWith(vns,"Fc1966_prdct[")==T,8]*crct1)
# dp3 = rbind(dp3, data.frame(Age=ages, Year = (rep(2016,Nages)),
#                  Period=rep("1990-2019",Nages),           
#                  PRexp=sumstats[startsWith(vns,"Fc2016_prdct[")==T,1]*crct2,
#                  PR_lo = sumstats[startsWith(vns,"Fc2016_prdct[")==T,4]*crct2,
#                  PR_hi = sumstats[startsWith(vns,"Fc2016_prdct[")==T,8]*crct2))
# dp3$Obs = numeric(length = nrow(dp3))
# dp3$ObsSE = numeric(length = nrow(dp3))
# for (i in 1:nrow(dp3)){
#   ii = which(df.Rep$Age==dp3$Age[i] & df.Rep$Year>(dp3$Year[i]-20) & 
#                df.Rep$Year<(dp3$Year[i]+20))
#   dp3$Obs[i] = mean(df.Rep$Prob[ii])
#   dp3$ObsSE[i] = sd(df.Rep$Prob[ii])/sqrt(length(df.Rep$Prob[ii]))
# }
# dp3$Year = as.factor(dp3$Year); dp3$Period = as.factor(dp3$Period)
# dp3 = dp3[which(dp3$Age<11),]
# pl3 = ggplot(data=dp3,aes(x=Age,y=PRexp,group=Period,color = Period, fill=Period)) +
#   geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi),alpha=0.3) + ylim(0,1) +
#   geom_line() + labs(x = "Age",y="Pregnancy rate") +
#   geom_point(aes(y=Obs),size=2) +
#   geom_errorbar(aes(ymin = Obs-1.96*ObsSE, 
#                       ymax = Obs+1.96*ObsSE),width=.2) +
#   # ggtitle("Model estimated vs observed pregancy rates by age",
#   #         subtitle = "Low density (1951-1985) vs. high density (1990-2019) population") + 
#   # scale_color_discrete(palette="Harmonic") +
#   scale_color_brewer(palette="Dark2") +
#   scale_fill_brewer(palette="Dark2") +
#   theme_classic() 
# print(pl3)
#
# Preg rate 8+ over time----------------------------------------------------
dp4 = data.frame(Year=Year,PRexp=sumstats[startsWith(vns,"Fc8_prdct[")==T,1],
                 PR_lo = sumstats[startsWith(vns,"Fc8_prdct[")==T,4],
                 PR_hi = sumstats[startsWith(vns,"Fc8_prdct[")==T,8])

dp4$PRprev_exp[2:Nyrsm1] = df.prev$Prg_smth[1:(nrow(df.prev)-1)]
dp4$PRprev_lo[2:Nyrsm1] = df.prev$Prg_smth_lo[1:(nrow(df.prev)-1)]
dp4$PRprev_hi[2:Nyrsm1] = df.prev$Prg_smth_hi[1:(nrow(df.prev)-1)]
dp4$Obs = rep(NA,nrow(dp4)) 
dp4$ObsSE = rep(NA,nrow(dp4)) 
dp4$Obs_N = rep(NA,nrow(dp4)) 
for (i in 1:nrow(dp4)){
  ii = which(df.prev$Year==dp4$Year[i] & df.prev$N_fem_smp>9)
  if(length(ii)==1){
    dp4$Obs[i] = df.prev$Prgrt_obs[ii]
    dp4$Obs_N[i] = df.prev$N_fem_smp[ii]
    dp4$ObsSE[i] = df.prev$Prgrt_obs_se[ii]
  }
}
pl4 = ggplot(data=dp4[dp4$Year>1959,],aes(x=Year)) +
  geom_ribbon(aes(ymin=PRprev_lo,ymax=PRprev_hi,fill="existing"),alpha=0.3) +
  geom_ribbon(aes(ymin=PR_lo,ymax=PR_hi,fill="revised"),alpha=0.3) +
  geom_line(aes(y=PRprev_exp,color="existing")) +
  geom_line(aes(y=PRexp,color="revised")) + 
  geom_point(aes(y=Obs)) + geom_errorbar(aes(ymin = Obs-1.96*ObsSE, 
                                            ymax = pmin(1,Obs+1.96*ObsSE)),
                                            color="darkgrey") +
  labs(x = "Year",y="Preganancy rate (Age 8+)") +
  scale_color_manual(name = "Legend",
                     values = c("existing" = "blue", "revised" = "darkred"),
                     labels = c("existing model", "revised model"),
                     aesthetics = c("colour", "fill")) +
#  ggtitle("Model estimated vs observed pregnancy rate over time",
#          subtitle = "Current model (red) vs. previous model (blue)") +
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(legend.position = c(0.2, 0.2),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(pl4)
#
# Age composition ---------------------------------------------------------
dp14 = data.frame(Year=Year,
                  Py_exp = sumstats[startsWith(vns,"Pyng_prdct[")==T,1],
                  Py_lo = sumstats[startsWith(vns,"Pyng_prdct[")==T,4],
                  Py_hi = sumstats[startsWith(vns,"Pyng_prdct[")==T,8],
                  Po_exp = sumstats[startsWith(vns,"Pold_prdct[")==T,1],
                  Po_lo = sumstats[startsWith(vns,"Pold_prdct[")==T,4],
                  Po_hi = sumstats[startsWith(vns,"Pold_prdct[")==T,8]
                  )
dp14$Py_prev_exp = numeric(length = nrow(dp14)) * NA
dp14$Py_prev_exp_lo = numeric(length = nrow(dp14)) * NA
dp14$Py_prev_exp_hi = numeric(length = nrow(dp14)) * NA
dp14$Po_prev_exp = numeric(length = nrow(dp14)) * NA
dp14$Po_prev_exp_lo = numeric(length = nrow(dp14)) * NA
dp14$Po_prev_exp_hi = numeric(length = nrow(dp14)) * NA
dp14$Py_Obs = numeric(length = nrow(dp14))
dp14$Py_Obs_lo = numeric(length = nrow(dp14))
dp14$Py_Obs_hi = numeric(length = nrow(dp14))
dp14$Py_Obs_SE = numeric(length = nrow(dp14))
z = 1.96
for (i in 1:nrow(dp14)){
  ii = which(row.names(AgeComp)==as.character(dp14$Year[i]))
  if(length(ii)==1){
    n = sum(AgeComp[ii,8:36])
    p = sum(AgeComp[ii,8:14])/sum(AgeComp[ii,8:36])
    dp14$Py_Obs[i] = p 
    dp14$Py_Obs_SE[i] = sqrt(p*(1-p)*(1/n))     
    # Wilson score interval with continuity correction for binomial CI
    dp14$Py_Obs_lo[i] = (2*n*p + z^2 - (z*sqrt(z^2 - 1/n + 4*n*p*(1-p) + (4*p - 2)) + 1))/(2*(n+z^2))
    dp14$Py_Obs_hi[i] = (2*n*p + z^2 + (z*sqrt(z^2 - 1/n + 4*n*p*(1-p) - (4*p - 2)) + 1))/(2*(n+z^2))
  }
}
dp14$Py_Obs[dp14$Py_Obs==0] = NA
dp14$Py_Obs_SE[dp14$Py_Obs==0 | dp14$Py_Obs_SE==0] = NA
dp14$Py_Obs_lo[dp14$Py_Obs==0 | dp14$Py_Obs_SE==0] = NA
dp14$Py_Obs_hi[dp14$Py_Obs==0 | dp14$Py_Obs_SE==0] = NA
dp14$Py_prev_exp[2:Nyrsm1] = df.prev$Ppn_yng[1:(nrow(df.prev)-1)]
dp14$Py_prev_exp_lo[2:Nyrsm1] = df.prev$Ppn_yng_lo[1:(nrow(df.prev)-1)]
dp14$Py_prev_exp_hi[2:Nyrsm1] = df.prev$Ppn_yng_hi[1:(nrow(df.prev)-1)]
dp14$Po_prev_exp[2:Nyrsm1] = df.prev$Ppn_8plus[1:(nrow(df.prev)-1)]
dp14$Po_prev_exp_lo[2:Nyrsm1] = df.prev$Ppn_8plus_lo[1:(nrow(df.prev)-1)]
dp14$Po_prev_exp_hi[2:Nyrsm1] = df.prev$Ppn_8plus_hi[1:(nrow(df.prev)-1)]
#
pl14 = ggplot(data=dp14[dp14$Year>1978,],aes(x=Year)) +
  geom_ribbon(aes(ymin=Py_prev_exp_lo,ymax=Py_prev_exp_hi,fill="Py_prev_exp"),alpha=0.3) +
  geom_ribbon(aes(ymin=Py_lo,ymax=Py_hi,fill="Py_exp"),alpha=0.3) +
  geom_line(aes(y=Py_prev_exp,color = "Py_prev_exp")) +  
  geom_line(aes(y = Py_exp,color = "Py_exp")) + 
  geom_point(aes(y=Py_Obs),size = 2) +
  geom_errorbar(aes(ymin = Py_Obs_lo, ymax = Py_Obs_hi)) +
  labs(x = "Year",y="Proportion adults < 15years") +
  #ggtitle("Age structure, model estiated vs. observed",
  #        subtitle = "Current model (red) vs. previous model (blue)") +  
  scale_color_manual(name = "Legend",
                     values = c("Py_prev_exp" = "blue", "Py_exp" = "darkred"),
                     labels = c("existing model", "revised model"),
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(1975,2020,by=5)) +
  theme_bw() + theme(legend.position = c(0.85, 0.9),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(pl14)
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
  geom_line() + labs(x = "Year",y="Total removals") +
  geom_point(aes(y=Obs),size=2) +
  # ggtitle("Model estimated vs observed harvest/bycatch mortality, by age class") + 
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
  
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
  geom_line() + labs(x = "",y="Fecundity (8+), deviations") +
  # ggtitle("Stochastic variation in fecundity") + 
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
# print(pl7)
#
# Juv survival stochasticity -------------------------------------------------------------
gamma_0 = sumstats[which(vns=="gamma_0"),1]
# hazImn = sumstats[which(startsWith(vns,"haz_Ice[")),1][21]  
phiS = sumstats[which(vns=="phi[2]"),1]
dlta = sumstats[which(vns=="dlta[2]"),1]
thtaS = sumstats[which(startsWith(vns,"thta[2]")),1]
Nml =  sumstats[which(startsWith(vns,"N[")),1]; Nml = Nml[-Nyrs]/1000000
# hzoth = hazImn + exp(omega)
dp8 = data.frame(Year=Yearp[-length(Yearp)],epsS=sumstats[startsWith(vns,"epsS[")==T,1],
                  epsS_lo = sumstats[startsWith(vns,"epsS[")==T,4],
                  epsS_hi = sumstats[startsWith(vns,"epsS[")==T,8])
dp8$S0_stoc = exp(-1 * (exp(omega + gamma_0 + dlta*CI[1:Nyrsm1] + 
                              (phiS*Nml)^thtaS + dp8$epsS)))
dp8$S0_lo =  exp(-1 * (exp(omega + gamma_0 + dlta*CI[1:Nyrsm1] +
                             (phiS*Nml)^thtaS + dp8$epsS_lo)))
dp8$S0_hi =  exp(-1 * (exp(omega + gamma_0 + dlta*CI[1:Nyrsm1] +
                             (phiS*Nml)^thtaS + dp8$epsS_hi)))
#
pl8 = ggplot(data=dp8,aes(x=Year,y=epsS)) +
  geom_ribbon(aes(ymin=epsS_lo,ymax=epsS_hi),alpha=0.3) +
  geom_line() + labs(x = "",y="Juv. log-haz., deviations") +
  # ggtitle("Stochastic variation in juvenile mortality") + 
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
# print(pl8b)
pl8b = ggplot(data=dp8,aes(x=Year,y=S0_stoc)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  geom_line() + labs(x = "Year",y="Juv. survival rate") +
  # ggtitle("Stoachastic and D-D variation in juvenile survival") + 
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank())
print(pl8b)

# plot_grid(pl7,pl8,pl8b,nrow = 3,labels = c("A","B","C"))

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
  # ggtitle("Age-specific variation in survival at low and high population density") 
  theme_bw()
print(pl9)
#
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
  geom_line() + labs(x = "Climate index (NCI)",y="Pregnancy rate (8+)") +
  # ggtitle("A) Effect of climate index on adult fecundity") + 
  theme_classic()
#print(pl11a)

pl11b = ggplot(dp11, aes(x=CI,y=Survival)) +
  geom_ribbon(aes(ymin=S0_lo,ymax=S0_hi),alpha=0.3) +
  geom_line() + labs(x = "Climate index (NCI)",y="Juvenile survival rate") +
  # ggtitle("B) Effect of climate index on juvenile survival") + 
  theme_classic()
# print(pl11b)

plot_grid(pl11a,pl11b,nrow = 2,labels = c("A","B"))


# Fecundity vs density ----------------------------------------------------
NN = seq(.1,5.8,by=0.1)
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
  # ggtitle("Density-dependent variation in adult fecundity") + 
  xlim(0, 7) +
  theme_classic()
# print(pl12)
#
# Juvenile survival vs density --------------------------------------------
NN = seq(.1,5.8,by=0.1)
gamma0_r = mcmc[iir,vn=="gamma_0"]
phiS_r = mcmc[iir,vn=="phi[2]"]
thtaS_r = mcmc[iir,vn=="thta[2]"]
S0 = matrix(0,nrow = 1000,ncol=length(NN))
for(r in 1:1000){
  haz_J = exp(omega + gamma0_r[r] + ( phiS_r[r] *NN)^thtaS_r[r])  
  S0[r,] = exp(-1 * (haz_J + exp(omega)))
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
  # ggtitle("Density-dependent variation in juvenile survival") 
  xlim(0, 7) +
  theme_classic()
# print(pl13)
#
plot_grid(pl12,pl13,nrow = 2,labels = c("A","B"))

# Juvenile mortality partitioned over time -----------------------------------
omega = -7
alpha0 = sumstats[which(vn=="alpha0"),1]
alpha1 = sumstats[which(vn=="alpha1"),1]
alpha2 = sumstats[which(vn=="alpha2"),1]
thta2 = sumstats[which(vn=="thta[2]"),1]
phi2 = sumstats[which(vn=="phi[2]"),1]
dlta2 = sumstats[which(vn=="dlta[2]"),1]
psi1 = sumstats[startsWith(vns,"psi[1]")==T,1]
psi2 = sumstats[startsWith(vns,"psi[2]")==T,1]
gamma_H0 = sumstats[startsWith(vns,"gamma_H0[")==T,1]
Nexp=sumstats[startsWith(vns,"N[")==T,1]; Nexp = Nexp[-1] 
epsS = sumstats[startsWith(vns,"epsS[")==T,1]
CIadj = CI[-1]
WIadj = WI[-Nyrs]
WI_CI = (WIadj - min(WIadj)) * (CIadj - min(CIadj)) 
WI_CI = WI_CI - mean(WI_CI)
gamma_0 = alpha0 + alpha1 * agesC[1] + alpha2 * agesC2[1] 
gamma_D = exp(thta2 * log(phi2*(Nexp/1000000))) 
haz_base = exp(omega + gamma_0 + gamma_D + dlta2 * CIadj + 
                 psi1 * WIadj + psi2 * WI_CI + epsS ) 
haz_CI = haz_base - exp(omega + gamma_0 + gamma_D + epsS ) 
haz_CI = pmax(.001,haz_CI)
haz_DD = pmax(0, haz_base - haz_CI)
haz_DD_exp = exp(omega + gamma_0 + gamma_D ) 
haz_H0 = exp(omega + gamma_H0) 
haz_tot = haz_DD + haz_CI + haz_H0
Mort = 1 - exp(-haz_tot)
#Mort_DDexp = Mort * haz_DD_exp/haz_tot
Mort_DDexp = 1 - exp(-haz_DD_exp)
Mort_DDexp = Mort_DDexp * (max(Mort)/max(Mort_DDexp)) * .7

dp14 = data.frame(Year=Year,
                  Density_stochastic = Mort * haz_DD/haz_tot,
                  Climate_ice_effects = Mort * haz_CI/haz_tot,
                  Removals = Mort * haz_H0/haz_tot)
dp14 = tidyr::pivot_longer(dp14,cols = !starts_with("Year"),
                           names_to = "hazard_type",values_to = "Mortality")
dp14$hazard_type = factor(dp14$hazard_type,
                          levels = c("Removals","Climate_ice_effects","Density_stochastic"))

pl14 = ggplot(dp14, aes(x=Year)) + 
  geom_area(aes(y=Mortality, fill=hazard_type), alpha=0.5) + 
  geom_line(data=data.frame(Year=Year, Mortality = Mort_DDexp, 
                            hazard_type="Baseline + D-D, no stochasticity"),
            aes(x=Year,y=Mortality,color=hazard_type),size=1.1,linetype="dashed") +  
  scale_fill_brewer(type="qual",palette = "Set1", name="Hazard type",
                    labels=c("Removals","Climate & ice effects","Baseline + D-D + stochasticity")) +
  scale_color_manual(values = c("darkgreen"),name=NULL) + 
  labs(x="Year",y="Mortality rate") +
  scale_x_continuous(breaks = seq(1950,2020,by=5)) +
  theme_bw() + theme(legend.position = "right",
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank()) +
  guides(color=guide_legend(title="",order = 2),fill=guide_legend(title="Hazard Type",order = 1))

print(pl14)
#
# Table ----------------------------------------------------------------------
ix = c(which(vns=="N0"),which(vns=="beta1"),which(vns=="beta2"),
       which(vns=="phi[1]"),which(vns=="thta[1]"),which(vns=="dlta[1]"),
       which(vns=="sigF"),which(vns=="alpha0"),which(vns=="alpha1"),
       which(vns=="alpha2"),which(vns=="phi[2]"),which(vns=="thta[2]"),
       which(vns=="dlta[2]"),which(vns=="psi[1]"),which(vns=="psi[2]"),
       which(vns=="sigS"),which(vns=="zeta"),which(vns=="gamma_H0_mn"),
       which(vns=="gamma_HA_mn"),which(vns=="sigH"),which(vns=="tau") )

statsum = sumstats[ix,c(1,3,5,7,10)]

row.names(statsum) = c("N0(millions)","Beta_1","Beta_2","phi_F","theta_F",
                       "delta_F","sigma_F","alpha_0","alpha_1","alpha_2",
                       "phi_S","theta_S","delta_S","psi_1","psi_2",
                       "sigma_S","zeta","Gamma_H0_bar","Gamma_HA_bar",
                       "sigma_H","tau")
statsum[1,1:4] = (statsum[1,1:4])/1000000

kable(statsum,digits=3,
      caption = "<strong>Statistics for fitted parameters</strong>",
      escape = FALSE,
      format = "html") %>% kable_styling(bootstrap_options = "basic")

write.excel <- function(x,row.names=T,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}
write.excel(statsum)

# Plot priors vs posteriors ---------------------------------
if(plotpri == 1){
  ixn = row.names(statsum)
  Npri = length(ix)
  df.pri = read_excel("./data/Priordefs_v2.xlsx")
  Paramname = character()
  Priorvals = numeric()
  Priorprob = numeric() 
  Postdens = numeric()  
  reps3 = 500
  for(i in 1:Npri){
    vals = seq(df.pri$lowerbnd[i],df.pri$upperbnd[i],length.out=reps3)+.0000001
    vals = vals * (1/df.pri$rescale[i])
    dst = df.pri$PriorDist[i]
    if(dst=="normal"){
      probs = dnorm(vals,df.pri$par1[i],df.pri$par2[i])
    }else if(dst=="gamma"){
      probs = dgamma(vals,df.pri$par1[i],df.pri$par2[i])  
    }else if(dst=="cauchy"){
      probs = dcauchy(vals,df.pri$par1[i],df.pri$par2[i])  
    }
    vals = vals * df.pri$rescale[i]
    probs = probs*(1/max(probs))
    posts = mcmc[,ix[i]]; if(vn[ix[i]]=="N0"){posts = posts/1000000 }
    ft = density(posts)
    postprobs = approx(ft$x,ft$y,xout=vals); postprobs = postprobs$y
    postprobs[is.na(postprobs)] = 0
    postprobs = postprobs*(1/max(postprobs))
    Paramname = c(Paramname,rep(ixn[i],reps3))
    Priorvals = c(Priorvals,vals)
    Priorprob = c(Priorprob,probs)
    Postdens = c(Postdens,postprobs)
  }
  dp_pri = data.frame(Paramname=Paramname,
                      Priorvals=Priorvals,
                      Priorprob=Priorprob,
                      Postdens=Postdens)
  dp_pri$Paramname = factor(dp_pri$Paramname, levels = ixn)
  plt_pri = ggplot(dp_pri,aes(x=Priorvals)) +
    geom_area(aes(y=Priorprob),fill="lightgrey") +
    geom_line(aes(y=Priorprob),color="darkgrey") +
    geom_area(aes(y=Postdens),fill="black",alpha=0.5) +
    geom_line(aes(y=Priorprob),color="white",linetype="dashed") +
    labs(x = "Value",y="Probability density") +
    theme_classic() +
    facet_wrap(~ Paramname,nrow = 5,ncol = 5,scales = "free")
  print(plt_pri)
}