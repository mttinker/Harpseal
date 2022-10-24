# Import data
require(readxl)
require(stats)
df.sad = read_excel("./data/SAD0_36.xlsx")
# CI index
df.NLCI = read_excel("./data/NL_climate_index.xlsx")
# Ice Anomalies
df.Ice = read_excel("./data/IceAnom_new.xlsx")
#
# Removals (harvest)
df.HV = read_excel("./data/Raw-removal-1952.xlsx")
df.HV$Grnpup = round(df.HV$greenland*df.HV$PpupGrnlnd)
df.HV$Grn1plus = round(df.HV$greenland - df.HV$Grnpup)
df.HV$Arctpup = round(0.034*df.HV$arctic)
df.HV$Arct1plus = round(df.HV$arctic - df.HV$Arctpup)
df.HV$PUPTOT = df.HV$canpup + df.HV$bypup + df.HV$Grnpup + df.HV$Arctpup
df.HV$ADLTOT = df.HV$can1plus + df.HV$by1plus + df.HV$Grn1plus + df.HV$Arct1plus
# Estimate 1951 as mean of 1952 and 1953
tmp = df.HV[1:2,]; tmp[1:2,1] = 1951
tmp = colMeans(tmp)
# Account for Struck and Loss (SnL)
df.HV = rbind(tmp,df.HV)
df.HV$canpup_wsnl = df.HV$canpup * (1/.99); 
df.HV$canpup_wsnl[which(df.HV$YEAR>1983)] = df.HV$canpup[which(df.HV$YEAR>1983)]* (1/.95);
df.HV$can1plus_wsnl = df.HV$can1plus * (1/.5)
df.HV$Grnpup_wsnl = df.HV$Grnpup * (1/.5)
df.HV$Grn1plus_wsnl = df.HV$Grn1plus * (1/.5)
df.HV$Arctpup_wsnl = df.HV$Arctpup * (1/.5)
df.HV$Arct1plus_wsnl = df.HV$Arct1plus * (1/.5)
df.HV$PUPTOTwSNL = df.HV$canpup_wsnl + df.HV$bypup + df.HV$Grnpup_wsnl + df.HV$Arctpup_wsnl
df.HV$ADLTOTwSNL = df.HV$can1plus_wsnl + df.HV$by1plus + df.HV$Grn1plus_wsnl + df.HV$Arct1plus_wsnl
df.HV$PUP_prob_rec = df.HV$PUPTOT / df.HV$PUPTOTwSNL
df.HV$ADL_prob_rec = df.HV$ADLTOT / df.HV$ADLTOTwSNL
# Pup Production
df.Pup = read_excel("./data/PupProd.xlsx")
# Abort_rate
df.Ab = read_excel("./data/Abort_rate.xlsx")
# Reproduction data (Pregnancy rates)
df.Rep = read_excel("./data/Reprodat.xlsx")
# Create matrix of age structure data for select years
# Filter pregnancy data 
ii = which(df.Rep$N==0 | is.na(df.Rep$Preg))
df.Rep = df.Rep[-ii,]
rm(ii,tmp)

# Age samples for age composition of adults
df.Age = read_excel("./data/Ages_sampled.xlsx")
Years_ages = as.numeric(colnames(df.Age[,-1]))
AgeComp = t(df.Age[,-1])

# # Examine preg rate vs age relationship (use year as proxy for DD):
# fit = nls(Prob ~ exp(b0 - b1*(8 - Age)^b2 - b3*(Year-1950))/(1+exp(b0 - b1*(8 - Age)^b2 - b3*(Year-1950))),
#           data = df.Rep,start = list(b0=2.5,b1=.1,b2=2,b3=.03), weights=N)
# summary(fit)
# newdata = data.frame(Year = rep(1960,51), Age = seq(3,8,by=0.1))
# newdata$Pred = predict(fit, newdata)
# plot(df.Rep$Age[df.Rep$Year<1970],df.Rep$Prob[df.Rep$Year<1970])
# lines(newdata$Pred ~ newdata$Age)
#
# Calculate smoothed pregnancy rate (as per previous model)
#
# require(locfit)
# require(ggplot2)
# require(gtools)
# ii = which(df.Rep$Year>1960 & df.Rep$Age==8 & df.Rep$N>4)
# dat_Prg = data.frame(Years = df.Rep$Year[ii] - 1960,
#                  Nsamp = df.Rep$N[ii],
#                  Npreg = df.Rep$Preg[ii])
# dat_Prg$Ppn_obs= dat_Prg$Npreg/dat_Prg$Nsamp
# dat_Prg$Ppn_se= sqrt((dat_Prg$Ppn_obs * (1 - dat_Prg$Ppn_obs) )/ dat_Prg$Nsamp)
# fit = locfit(Npreg ~ lp(Years,deg=2,nn=0.95),weights=Nsamp,data=dat_Prg,family="binomial")
# newdat = seq(1,2019-1960)
# pred_dat = predict(fit,newdat,se.fit=T)
# Smthpreg = data.frame(Year = seq(1961,2019),
#                       Preg_smth = pred_dat$fit,
#                       Preg_smth_SE = pred_dat$se.fit)
# tmp = preplot(fit,newdat,band = "local")
# tmp$CI_lo = inv.logit(tmp$fit - 1.96*tmp$se.fit) 
# tmp$CI_hi = inv.logit(tmp$fit + 1.96*tmp$se.fit) 
# Smthpreg$Preg_smth_lo = tmp$CI_lo
# Smthpreg$Preg_smth_hi = tmp$CI_hi
# 
# ggplot(Smthpreg,aes(x=Year,y=Preg_smth)) +
#   geom_ribbon(aes(ymin=Preg_smth_lo,ymax=Preg_smth_hi),alpha=0.2) +
#   geom_line() +
#   geom_point(data=dat_Prg,aes(x=Years+1960,y=Ppn_obs)) +
#   geom_errorbar(data=dat_Prg,aes(x=Years+1960,y=Ppn_obs,
#                                  ymin=Ppn_obs-Ppn_se,ymax=Ppn_obs+Ppn_se)) +
#   theme_classic()
# 

