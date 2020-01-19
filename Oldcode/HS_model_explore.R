# HS model exploratory analysis 
require(gtools)
require(ggplot2)
source("HS_matrixD.r")
# Survival params: log hazard rates for adults and juvs
a0 = -5   # Base log hazards (arbitrary constant, gives Survival of 0.99)
aA = 2.95  # Adult log hazards: strong prior
aA1 = .056
aA2 = .0015
aJ = 4 # Juv log hazards at N~0: weak prior
ages = seq(1,20); ages2 = ages^2
gammA = a0 + aA - aA1*ages + aA2*ages2; S = exp(-exp(gammA)) ; 
SJ = exp(-exp(a0+aJ))
S = c(SJ,S)
AgeY = c(0,ages); plot(AgeY,S,type="l",ylim = c(0,1))
print(paste0("Adult survival = ", format(max(S),digits = 3)))
print(paste0("Maximum Juvenile survival = ", format(SJ,digits = 3)))
alpha = c(aJ,aA, aA1, aA2)
phiJ = 0.07  # juv surviv D-D param: weak prior
# Fecunduty
ages = seq(1,20)
ages2 = ages^2
b0 = 2.5    # Logit of max adult pregancy rate (0.9): strong prior
b1 = 0.1    # Effect of age on preg rate: weak prior
beta = c(b0,b1) 
phiF = 0.195;  # fecundity D-D param: weak prior
theta = 2.4; # theta (from theta-logistic growth fxn, fixed param)
N = 2*1000000
Fi = numeric(length = 20)
Fi[1:2] = 0
Fi[3:8] = inv.logit(b0 - b1*(8-ages[3:8])^2 - (phiF*N/1000000)^theta)
Fi[9:20] = Fi[8]
plt1 = ggplot(data.frame(Age=ages[3:8],Fecundity=Fi[3:8]),aes(x=Age,y=Fecundity)) +
  geom_line() + ylim(c(0,1)) + ggtitle("Age effects on Fecundity") + theme_classic()
print(plt1)
phi = c(phiF,phiJ)  
# Psi params: % ice cover effects on pup survival 
psi1 = 4; psi2 = 1.5 # Psi has moderately strong priors (sd =0.5) based on expert opinion
psi = c(psi1,psi2)
Ice_cover = seq(-1,1,by=0.02)
S_ice = inv.logit(-psi1 + 7*(Ice_cover+1)^psi2)
#plot(Ice_cover,S_ice,type="l",col="black")
S_ice_r = matrix(0,nrow = 1000,ncol = length(Ice_cover))
for(r in 1:1000){
  psi1r = rnorm(1,psi1,1); psi2r = rnorm(1,psi2,0.5)  
  S_ice_r[r,1:length(Ice_cover)] = inv.logit(-psi1r +7*(Ice_cover+1)^psi2r)
  #lines(Ice_cover,S_ice_r[r,1:100],col="red")
}
#lines(Ice_cover,S_ice,col="black",lw=2)
dfplt = data.frame(Ice_cover=Ice_cover,
                   Mean_S_ice = S_ice,
                   Lo = apply(S_ice_r,2,quantile,0.025),
                   Hi = apply(S_ice_r,2,quantile,0.975))
plt2 = ggplot(dfplt,aes(x=Ice_cover,y=Mean_S_ice)) +
  geom_ribbon(aes(ymin=Lo,ymax=Hi),alpha=0.3) +
  geom_line() +
  labs(x="Ice Index (% Cover Anomaly)",y="Relative pup survival") +
  ggtitle("Effect of ice cover on pup survival") +
  theme_classic()
print(plt2)
# Effect of envronmental index (Et), positive values = more resources, higher fecund 
delta = 1 # Weak prior
#
# Set some arbitary values of variables
IC = c(0,0,0) # Ice cover index, by Area (S.GSL, N.GSL, Front)
CE = 0;  # Environmental index (average = 0, positive = more resources)
PA = c(0.18,0.07,0.75) # mean proprortion of pupping by area
#
# **NOTE** : use dirichlet dist (in stan) to generate random PA values each year,
#  For example: 
PApriors = 50*PA
PAvals = rdirichlet(10000, PApriors)
colMeans(PAvals)
min(PAvals[,2])
max(PAvals[,2])
apply(PAvals,2,sd)
boxplot(PAvals)
#
# Calculate matrix ----------------------------------------------
# Sim to estimate initial "Stable Age Dist" with harvest 
sad1 = seq(0.12,0.01,by = -.006); sad1 = sad1/sum(sad1)
# Use 1951 harvest estimates for beaters and all age 1+ (combined)
#  and distributed the adult harvest according to stable age dist
HV = c(207756,round(100425*sad1))
# Iterate until adult stage distribution with harvest stabilizes
for (t in 1:20){
  Nt = 2400000  # 1951 pop size for matrix calcs (9000000)
  rslt = HS_matrixD(alpha,beta,delta,phi,psi,theta,PA,
                    Nt,CE,IC,HV)
  sad1 = as.numeric(rslt$sad[2:20])/sum(as.numeric(rslt$sad[2:20]))
  HV = c(207756,round(100425*sad1))
}
SAD0 = rslt$sad
SAD0 = data.frame(Age = ages,SAD = SAD0)
write.csv(SAD0,file="./data/SAD0.csv",row.names = F)
# Sims to estimate K with no harvest
Nt = 6500000  # Current pop size for matrix calcs (9000000)
HV = rep(0,20) # Vector of age-specific Harvest values (0)
# HV = c(400000, 20000, 3000, 1000, 500, rep(50,10),rep(0,5))
# HV = .2*HV
rslt = HS_matrixD(alpha,beta,delta,phi,psi,theta,PA,
                       Nt,CE,IC,HV)
print(paste0("Lambda at N of ",format(Nt/1000000,digits = 3),
             " million = ",format(rslt$lambda,digits = 4)))
print(paste0("Pups at N of ",format(Nt/1000000,digits = 3),
             " million: S GSL = ",format(rslt$pups[1],digits = 4),
              ", N GSL = ",format(rslt$pups[2],digits = 4),
             ", Front = ",format(rslt$pups[3],digits = 4)))
#
# Estimate K: should do this without harvest and expected future conditions
if (sum(HV) == 0){
  Nsize = seq(1000000,10000000,length.out = 20)
  Grate = numeric(length = 20)
  Fecnd =  numeric(length = 20)
  for (i in 1:20){
    rslti = HS_matrixD(alpha,beta,delta,phi,psi,theta,PA,
                       Nsize[i],CE,IC,HV)
    Grate[i] = rslti$lambda
    Fecnd[i] = rslti$Fecund[8]
  }
  plt3 = ggplot(data.frame(N=Nsize,Fecundity=Fecnd),aes(x=N,y=Fecundity)) +
    geom_line() + ggtitle("Density-dependent effects on Fecundity") +
    labs(x="Population size (N)",y="Adult (8+) Fecundity") + theme_classic()
  print(plt3)
  intrp = splinefun(Nsize,Grate,method = "fmm")
  Nvals = seq(5000000,10000000,by=1000)
  Gvals = intrp(Nvals)
  ii = which(abs(Gvals-1)==min(abs(Gvals-1)))
  K = Nvals[ii[1]]
  plt4 = ggplot(data.frame(N=Nsize,Lambda=Grate),aes(x=N,y=Lambda)) +
    geom_line() + ggtitle("Density-dependent population growth") +
    geom_hline(yintercept = 1, colour = "red") +
    geom_vline(xintercept = K,colour = "green") +
    labs(x="Population size (N)",y="Lambda (annual rate of growth)") + theme_classic()  
  print(plt4)
  print(paste0("K for these conditions estimated as ",format(K/1000000,digits = 4),
               " million"))
}
