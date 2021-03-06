# Import data
require(readxl)
df.sad = read_excel("./data/SAD0_36.xlsx")
# CI index
df.NLCI = read_excel("./data/NL_climate_index.xlsx")
# Ice Anomalies
df.Ice = read_excel("./data/IceAnom.xlsx")
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
