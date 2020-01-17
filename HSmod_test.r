HSmod_test <- function(init_fun,stan.data) {
  # create some counter variables: repeating year 1 20x to get sad correct
  futuresim = 0 # (note: if this is a "futer sim", loaded data will update this)
  thta = 1 # (NOTE if thta provided as fixed user input, it will replace this)
  reps=100
  for(i in 1:length(stan.data)){
    ##first extract the object value
    tempobj=stan.data[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(stan.data)[[i]],"= tempobj")))
  }
  aA = Adhzpri
  aJ = Jvhzpri   # Loop through random iterations of model
  iy = c(rep(1,20),seq(2,Nyrs-1))
  iz = numeric(length = length(iy)); iz[1:19] = 1
  # Create some arrays for saving results
  PrPredict = array(dim=c(Nyrs-1,reps))
  PrAgPred = array(dim=c(Nages,reps))
  N_Predict = array(dim=c(Nyrs,reps))
  P_Predict = array(dim=c(Nyrs-1,reps))
  Agepredict1 = array(dim=c(NFages,reps))
  Agepredict2 = array(dim=c(NFages,reps))
  HVp_predict = array(dim=c(Nyrs-1,reps))
  HVa_predict = array(dim=c(Nyrs-1,reps))
  sadmean = array(dim = c(Nages,reps))
  for(r in 1:reps){
  pars = init_fun()
  for(i in 1:length(pars)){
    ##first extract the object value
    tempobj=pars[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(pars)[[i]],"= tempobj")))
  }
  # Initialize variables
  N = numeric(length = Nyrs)
  Nml = numeric(length = Nyrs)
  n = array(dim = c(Nages,Nyrs))
  Fc = array(dim = c(Nages,Nyrs-1))
  Sice = array(dim = c(Nareas,Nyrs-1))
  SJ = numeric(length = Nyrs-1)
  S = array(dim = c(Nages,Nyrs-1))
  HVp_pred = numeric(length = Nyrs-1)
  HVa_pred = numeric(length = Nyrs-1)
  FmAgeVec = array(dim = c(NFages,Nyrs-1))
  PredPup = numeric(length = Nyrs-1)
  PredPupA = array(dim = c(Nareas,Nyrs-1))
  # Params, Random assignments
  epsF = rnorm(Nyrs-1,0,sigF)
  # epsJ = rnorm(Nyrs-1,0,sigJ)
  psi1 = rnorm(1,psipri1a,psipri1b)
  psi2 = rnorm(1,psipri2a,psipri2b)
  # Add some of the "observed" harvest patterns of major deviations
  gammHp = rnorm(Nyrs-1,gammHp_mn,sigH)
  gammHa = rnorm(Nyrs-1,gammHa_mn,sigH)
  if (futuresim==0){
    gammHp[1:4] = gammHp[1:4]*1.15
    gammHp[5:21] = gammHp[5:21]*1.3
    gammHp[33:45] = gammHp[33:45]*0.8
    gammHp[46:58] = gammHp[46:58]*1
    gammHp[59:69] = gammHp[59:69]*0.8
    gammHa[1:4] = gammHa[1:4]*1.15
    gammHa[5:21] = gammHa[5:21]*1.3
    gammHa[33:45] = gammHa[33:45]*0.8
    gammHa[46:58] = gammHa[46:58]*1
    gammHa[59:69] = gammHa[59:69]*0.8
  }
  PA = array(dim=c(3,Nyrs-1))
  for(i in 1:(Nyrs-1)){
    PA[1:3,i] = rdirichlet(1, PApri)
  }
  N[1] = N0pri 
  Nml[1] = N[1]/1000000 
  # Loop through years to calculate demographic transitions and population dynamics 
  for (ix in 1:length(iy)){
    i = iy[ix]
    if (i == 1){
      n[,i] = sad0 * N[i]
    }
    #  Declare some temporary variables for this year:
    haz_J = numeric()
    haz_A = numeric(length = Nages)
    gammIce = numeric(length = Nareas)
    haz_I_A = numeric(length = Nareas)
    haz_I = numeric()
    haz_HVp = numeric()
    haz_HVa = numeric()
    prp_HV_p = numeric()
    prp_HV_a = numeric(length = Nages)
    nt1 = numeric(length = Nages)
    #  Fecundity (compared with observed pregnancy rate: adjust for abortions?)
    Fc[1:2,i] = rep(0,2) 
    Fc[3:Nages,i] = inv.logit(b0 - b1 * (Nages-ages[3:Nages])^2 - (phiF*Nml[i])^thta + dlta*CE[i]*(Nml[i]) + epsF[i])
    #  Juvenile competing hazards and net realized survival
    haz_J = exp(gamm0 + aJ + (phiJ*Nml[i])^thta) ; #  + epsJ[i]
    for(a in 1:Nareas){ 
      gammIce[a] = psi1 * (1/exp(IC[i,a]))^psi2 - 1 
      Sice[a,i] = exp(-1 * exp(gamm0 + gammIce[a])) 
      haz_I_A[a] = exp(gamm0 + gammIce[a]) 
    } 
    haz_I = sum(PA[,i] * haz_I_A)
    haz_HVp = exp(gamm0 + gammHp[i]) 
    SJ[i] = exp(-1 * (haz_J + haz_I + haz_HVp)) 
    #  Adult competing hazards and net realized survival
    for (a in 1:Nages){
      haz_A[a] = exp(gamm0 + aA + (DDadlt[a]*phiJ*Nml[i])^thta) 
    }
    haz_HVa = exp(gamm0 + gammHa[i]) 
    S[,i] = exp(-1 * (haz_A + haz_HVa))     
    #  Calculations for proportional mortality from harvest
    prp_HV_p = haz_HVp / (haz_J + haz_I + haz_HVp) 
    HVp_pred[i] = sum((n[,i] * Fc[,i]) * 0.5) * (1 - SJ[i]) * prp_HV_p 
    prp_HV_a = haz_HVa / (haz_A + haz_HVa) 
    HVa_pred[i] = sum((n[,i] * (1 - S[,i])) * prp_HV_a) 
    #  Annual demographic transitions
    nt1[1] = sum(0.5 * n[,i] *  ((Fc[,i] * SJ[i]))) 
    nt1[2:(Nages-1)] = n[1:(Nages-2),i] * S[1:(Nages-2),i]
    nt1[Nages] = n[Nages-1,i] * S[Nages-1,i] + n[Nages,i] * S[Nages,i] ;    
    if (iz[ix] == 1){
      sad0 = nt1/sum(nt1)
    }else{
      n[,i+1] = nt1  ;
      N[i+1] = sum(n[,i+1]) ;
      Nml[i+1] = N[i+1]/1000000 ;    
      #  Female predicted age vector (ages 3 - 8+), for year i
      FmAgeVec[,i] = n[NFage1:Nages,i] / (0.00001 + sum( n[NFage1:Nages,i])) ;
      #  Predicted Pups to be surveyed (allowing for early ice mortality), 
      #    by area and total, for year i 
      for (a in 1:Nareas){
        PredPupA[a,i] = sum((n[,i] * Fc[,i]) * 0.5 * PA[a,i] * sqrt(Sice[a,i])) ;
      }
      PredPup[i] = sum(PredPupA[,i])
    }
  }
  N_Predict[,r] = N  
  P_Predict[,r] = PredPup  
  PrPredict[,r] = Fc[8,]
  PrAgPred[,r] = Fc[1:8,1]
  # Age dist for sample years 17 = 1967, 69 = 2019
  if(futuresim==0){
    Agepredict1[,r] = FmAgeVec[,17]
    Agepredict2[,r] = FmAgeVec[,69]    
  }else{
    Agepredict1[,r] = FmAgeVec[,1]
    Agepredict2[,r] = FmAgeVec[,Nyrs-1]     
  }
  HVp_predict[,r] = HVp_pred
  HVa_predict[,r] = HVa_pred
  sadmean[,r] = sad0
  }
  sadmn = rowMeans(sadmean)
  results = list(N_Predict=N_Predict,P_Predict=P_Predict,PrPredict=PrPredict,
                 Agepredict1=Agepredict1,Agepredict2=Agepredict2,PrAgPred=PrAgPred,
                 HVp_predict=HVp_predict,HVa_predict=HVa_predict,SAD=sadmn)
  return(results)  
}
