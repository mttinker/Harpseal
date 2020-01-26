HSmod_sim <- function(init_fun,stan.data,sumstats1,vns1) {
  # Function to evaluate harp seal model by simulations
  #
  # Create function to save multiple results from parallel processes
  multiResultClass <- function(result1=NULL,result2=NULL,result3=NULL,
                               result4=NULL,result5=NULL,result6=NULL,
                               result7=NULL,result8=NULL,result9=NULL){
  me <- list(
      result1 = result1,
      result2 = result2,
      result3 = result3,
      result4 = result4,
      result5 = result5,
      result6 = result6,
      result7 = result7,
      result8 = result8,
      result9 = result9
    )
    ## Set the name for the class
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  # SOme default values that may get overwritten by input files:
  futuresim = 0 # (note: if this is a "futer sim", loaded data will update this)
  thta = 2 # (NOTE if thta provided as fixed user input, it will replace this)
  PAmeans = c(.18,.07,.75)
  reps=100
  # Extract variables from input data list"
  for(i in 1:length(stan.data)){
    ##first extract the object value
    tempobj=stan.data[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(stan.data)[[i]],"= tempobj")))
  }
  # Process variables and set up arrays for simulations
  gammA = Adloghz
  gammJ = Jvloghz   
  if(nrow(IC)==1000){
    ICr = IC
    CEr = CE    
  }
  if (futuresim==2){
    # random multipliers for harvest hazards,
    # used to increase range of values considered
    set.seed(123)
    ig1 = runif(reps,.5,1.5)
    set.seed(321)
    ig2 = runif(reps,.8,1.5)
  }
  iy = c(rep(1,20),seq(2,Nyrs-1))
  iz = numeric(length = length(iy)); iz[1:19] = 1
  # Create some arrays for saving results
  N_Predict = array(dim=c(Nyrs,1))
  P_Predict = array(dim=c(Nyrs-1,1))
  PrPredict = array(dim=c(Nyrs-1,1))
  PrAgPred = array(dim=c(Nages,1))  
  Agepredict1 = array(dim=c(NFages,1))
  Agepredict2 = array(dim=c(NFages,1))
  HVp_predict = array(dim=c(Nyrs-1,1))
  HVa_predict = array(dim=c(Nyrs-1,1))
  sadmean = array(dim = c(Nages,1))
  #
  # Loop through random iterations of model
  #  using Parallel processing to speed things up
  set.seed(123)
  simresults = foreach(r=1:reps) %dorng% {
    sumstats = sumstats1
    vns=vns1
    pars = init_fun()
    for(i in 1:length(pars)){
      ##first extract the object value
      tempobj=pars[[i]]
      ##now create a new variable with the original name of the list item
      eval(parse(text=paste(names(pars)[[i]],"= tempobj")))
    }
    # Initialize variables
    if(nrow(IC)==1000){
      IC = ICr[sample(1000,Nyrs,replace = T),]
      CE = CEr[sample(1000,Nyrs,replace = T)]
    }
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
    # Harvest mort: depends on type of scenario
    if(futuresim==0){
      # Level of harvest as observed
      gammHp = rnorm(Nyrs-1,gammHp_mn,sigH)
      gammHa = rnorm(Nyrs-1,gammHa_mn,sigH)
      # Add some of the "observed" harvest patterns of major deviations
      gammHp[1:4] = gammHp[1:4]*1.1
      gammHp[5:21] = gammHp[5:21]*1.15
      gammHp[33:45] = gammHp[33:45]*0.8
      gammHp[46:58] = gammHp[46:58]*1
      gammHp[59:69] = gammHp[59:69]*0.8
      gammHa[1:4] = gammHa[1:4]*1.1
      gammHa[5:21] = gammHa[5:21]*1.15
      gammHa[33:45] = gammHa[33:45]*0.8
      gammHa[46:58] = gammHa[46:58]*1
      gammHa[59:69] = gammHa[59:69]*0.8
    }else if(futuresim==1){
      # No harvest hazards (for estimating K)
      gammHp = rep(0,Nyrs-1)
      gammHa = rep(0,Nyrs-1)
    }else if(futuresim==2){
      # Range of harvest levels, for finding TAC criteria
      # (then harvest rate remains constant for years within a sim)
      gammHp = rep(gammHp_mn*ig1[r],Nyrs-1)
      gammHa = rep(gammHa_mn*ig2[r],Nyrs-1)
    }
    PA = array(dim=c(3,Nyrs-1))
    for(i in 1:(Nyrs-1)){
      PA[1:3,i] = gtools::rdirichlet(1, 20*PAmeans)
    }
    N[1] = N0pri 
    Nml[1] = N[1]/1000000 
    # Loop through years to calculate demographic transitions and population dynamics 
    for (ix in 1:length(iy)){
      npsv = 1-(gtools::inv.logit(rnorm(1,0,1))/10) 
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
      Fc[3:Nages,i] = gtools::inv.logit(b0 - b1 * (Nages-ages[3:Nages])^2 - (phiF*Nml[i])^thta + dlta*CE[i]*(Nml[i]) + epsF[i])
      #  Juvenile competing hazards and net realized survival
      haz_J = exp(gamm0 + gammJ + (phiJ*Nml[i])^thta) ; #  + epsJ[i]
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
        haz_A[a] = exp(gamm0 + gammA + (DDadlt[a]*phiJ*Nml[i])^thta) 
      }
      haz_HVa = exp(gamm0 + gammHa[i]) 
      S[,i] = exp(-1 * (haz_A + haz_HVa))     
      #  Calculations for proportional mortality from harvest
      prp_HV_p = haz_HVp / (haz_J + haz_I + haz_HVp) 
      HVp_pred[i] = sum((n[,i] * Fc[,i]) * 0.5 * npsv) * (1 - SJ[i]) * prp_HV_p 
      if (futuresim==2){
        HVp_pred[i] = HVp_pred[i] - (Grn_P+Arc_P+Byc_P)
      }
      prp_HV_a = haz_HVa / (haz_A + haz_HVa) 
      HVa_pred[i] = sum((n[,i] * (1 - S[,i])) * prp_HV_a) 
      if (futuresim==2){
        HVa_pred[i] = HVa_pred[i] - (Grn_A+Arc_A+Byc_A)
      }    
      #  Annual demographic transitions
      nt1[1] = sum((n[,i] * Fc[,i]) * 0.5 * npsv) * (SJ[i])
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
          PredPupA[a,i] = sum((n[,i] * Fc[,i]) * 0.5 * npsv * PA[a,i] ) ;
        }
        PredPup[i] = sum(PredPupA[,i])
      }
    }
    N_Predict[,1] = N  
    P_Predict[,1] = PredPup  
    PrPredict[,1] = Fc[8,]
    PrAgPred[,1] = Fc[1:8,1]
    # Age dist for sample years 17 = 1967, 69 = 2019
    if(futuresim==0 & Nyrs>69){
      Agepredict1[,1] = FmAgeVec[,17]
      Agepredict2[,1] = FmAgeVec[,69]    
    }else{
      Agepredict1[,1] = FmAgeVec[,1]
      Agepredict2[,1] = FmAgeVec[,Nyrs-1]     
    }
    HVp_predict[,1] = HVp_pred
    HVa_predict[,1] = HVa_pred
    sadmean[,1] = sad0
    result <- multiResultClass()
    result$result1 <- N_Predict
    result$result2 <- P_Predict
    result$result3 <- PrPredict
    result$result4 <- PrAgPred
    result$result5 <- Agepredict1
    result$result6 <- Agepredict2
    result$result7 <- HVp_predict
    result$result8 <- HVa_predict
    result$result9 <- sadmean
    return(result)
  }
  #
  N_Predict = array(dim=c(Nyrs,reps))
  P_Predict = array(dim=c(Nyrs-1,reps))  
  PrPredict = array(dim=c(Nyrs-1,reps))
  PrAgPred = array(dim=c(Nages,reps))
  Agepredict1 = array(dim=c(NFages,reps))
  Agepredict2 = array(dim=c(NFages,reps))
  HVp_predict = array(dim=c(Nyrs-1,reps))
  HVa_predict = array(dim=c(Nyrs-1,reps))
  sadmean = array(dim = c(Nages,reps))
  for(r in 1:reps){
    N_Predict[,r] = simresults[[r]]$result1[,1]
    P_Predict[,r] = simresults[[r]]$result2[,1]  
    PrPredict[,r] = simresults[[r]]$result3[,1]  
    PrAgPred[,r] = simresults[[r]]$result4[,1]
    Agepredict1[,r] = simresults[[r]]$result5[,1]  
    Agepredict2[,r] = simresults[[r]]$result6[,1] 
    HVp_predict[,r] = simresults[[r]]$result7[,1]
    HVa_predict[,r] = simresults[[r]]$result8[,1]  
    sadmean[,r] = simresults[[r]]$result9[,1]   
  }
  sadmn = rowMeans(sadmean)
  results = list(N_Predict=N_Predict,P_Predict=P_Predict,PrPredict=PrPredict,
                 PrAgPred=PrAgPred,Agepredict1=Agepredict1,Agepredict2=Agepredict2,
                 HVp_predict=HVp_predict,HVa_predict=HVa_predict,SAD=sadmn)
  return(results)  
}
