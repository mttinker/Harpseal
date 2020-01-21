# For testing inits-----------------------------------------------------
#
source("HSmod_sim.r")
Year = seq(Year1,YearT-1)
Yearp = seq(Year1,YearT)
rslt=HSmod_sim(init_fun,stan.data)
N_Predict = rslt$N_Predict
P_Predict = rslt$P_Predict
PrPredict = rslt$PrPredict
PrAgPred = rslt$PrAgPred
Agepredict1 = rslt$Agepredict1 # Age dist for sample year 17 = 1967 (Agecomp row5)
Agepredict2 = rslt$Agepredict2 # Age dist for sample year 69 = 2019 (Agecomp row48)
HVp_predict = rslt$HVp_predict
HVa_predict = rslt$HVa_predict
ggplot(data.frame(Year=Yearp,Npred = rowMeans(N_Predict)),aes(x=Year,y=Npred)) +
  geom_line()
ggplot(data.frame(Year=Year,HvPpred = rowMeans(HVp_predict)),aes(x=Year,y=HvPpred)) +
  geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=PUPTOT))
ggplot(data.frame(Year=Year,HvApred = rowMeans(HVa_predict)),aes(x=Year,y=HvApred)) +
  geom_line() + geom_point(data=df.HV,aes(x=YEAR,y=ADLTOT))
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

# sad = rslt$SAD
rm(rslt,N_Predict,P_Predict,PrPredict,PrAgPred,Agepredict1,
   Agepredict2,HVp_predict,HVa_predict)