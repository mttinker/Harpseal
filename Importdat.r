# Import data
require(readxl)
df.sad = read.csv("./data/SAD0.csv", header=T)
# CE index
tmp = as.matrix(read.table("./data/CElindex.csv", sep=";", header=F))
df.CE = as.data.frame(t(tmp),row.names = F); colnames(df.CE) <- c("Year","CEindex")  
# Ice Anomalies
df.Ice = read_excel("./data/IceAnom.xlsx")
# Create "ice anomaly index" IC by subtracting mean cover for years before 2000,
# then multiply by a consutant such that index varies between -1 and 1)
mnIceGulf =  mean(df.Ice$Gulf_conc[df.Ice$Year<2000])
mnIceLabsea =  mean(df.Ice$Labsea_conc[df.Ice$Year<2000])
df.Ice$Gulf_Anom = pmin(1,pmax(-1,(df.Ice$Gulf_conc - mnIceGulf)*3))
df.Ice$Lab_Anom = pmin(1,pmax(-1,(df.Ice$Labsea_conc - mnIceLabsea)*4))
# Removals (harvest)
tmp = as.matrix(read.table("./data/raw-removal-1952.csv", sep=";", header=T)) 
df.HV = as.data.frame(tmp,row.names = F)
tmp = as.matrix(read.table("./data/prop-green-pup-1952.csv", header=F))
df.HV$PpupGrnlnd = tmp[1:nrow(df.HV)]
df.HV$Grnpup = round(df.HV$greenland*df.HV$PpupGrnlnd)
df.HV$Grn1plus = round(df.HV$greenland - df.HV$Grnpup)
df.HV$Arctpup = round(0.034*df.HV$arctic)
df.HV$Arct1plus = round(df.HV$arctic - df.HV$Arctpup)
df.HV$PUPTOT = df.HV$canpup + df.HV$bypup + df.HV$Grnpup + df.HV$Arctpup
df.HV$ADLTOT = df.HV$can1plus + df.HV$by1plus + df.HV$Grn1plus + df.HV$Arct1plus
# Estimate 1951 as mean of 2952 and 1953
tmp = df.HV[1:2,]; tmp[1:2,1] = 1951
tmp = colMeans(tmp)
df.HV = rbind(tmp,df.HV)
# Pup Production
df.Pup = read_excel("./data/PupProd.xlsx")
# Abort_rate
df.Ab = read_excel("./data/Abort_rate.xlsx")
# Reproduction data (Pregnancy rates)
df.Rep = read_excel("./data/Reprodat.xlsx")
# Create matrix of age structure data (for 3 - 8+ year olds) for select years
Years_ages = unique(df.Rep$Year)
AgeComp = matrix(0,nrow = length(Years_ages), ncol = 6)
for (y in 1:length(Years_ages)){
  AgeComp[y,1] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==3)]
  AgeComp[y,2] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==4)]
  AgeComp[y,3] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==5)]
  AgeComp[y,4] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==6)]
  AgeComp[y,5] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==7)]
  AgeComp[y,6] = df.Rep$N[which(df.Rep$Year==Years_ages[y] & df.Rep$Age==8)]
}
# Filter pregnancy data 
ii = which(df.Rep$N==0 | is.na(df.Rep$Preg))
df.Rep = df.Rep[-ii,]
rm(ii,tmp,y)
# # Examine preg rate vs age relationship (use year as proxy for DD):
# fit = nls(Prob ~ exp(b0 - b1*(8 - Age)^b2 - b3*(Year-1950))/(1+exp(b0 - b1*(8 - Age)^b2 - b3*(Year-1950))),
#           data = df.Rep,start = list(b0=2.5,b1=.1,b2=2,b3=.03), weights=N)
# summary(fit)
# newdata = data.frame(Year = rep(1960,51), Age = seq(3,8,by=0.1))
# newdata$Pred = predict(fit, newdata)
# plot(df.Rep$Age[df.Rep$Year<1970],df.Rep$Prob[df.Rep$Year<1970])
# lines(newdata$Pred ~ newdata$Age)
