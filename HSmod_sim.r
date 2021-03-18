HSmod_sim <- function(init_fun,sim.data,mcmc1,vn1) {
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
  # TO DELETE: Some default values that may get overwritten by input files:
  #  futuresim = 0 # (note: if this is a "future sim", loaded data will update this)
  #  thta = 2 # (NOTE if thta provided as fixed user input, it will replace this)
  #  PAmeans = c(.18,.07,.75)
  #  reps=100
  #
  # Extract variables from input data list"
  for(i in 1:length(sim.data)){
    ##first extract the object value
    tempobj=sim.data[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(sim.data)[[i]],"= tempobj")))
  }
  # set.seed(123)
  # r_vec2 = sample(1000,reps,replace = T)
  # PAr = gtools::rdirichlet(1000+Nyrs, 25*PAmeans)
  # epsFr = matrix(rnorm(1000*Nyrs,0,1),nrow=1000)
  # epsSr = matrix(rnorm(1000*Nyrs,0,1),nrow=1000)
  # ig1 = runif(reps,.3,1.5)
  # ig2 = runif(reps,.05,1.2)
  # Process variables and set up arrays for simulations
  if(nrow(IC)==1000){
    ICr = rbind(IC,IC[1:Nyrs,])  
    CIr = c(CI,CI[1:Nyrs])    
  }else{
    ICr = IC
    CIr = CI
  }
  if (futuresim==2){
    Mort_othP = Grn_P*(1/.5) + Arc_P*(1/.5) + Byc_P
    Mort_othA = Grn_A*(1/.5) + Arc_A*(1/.5) + Byc_A
  }
  iy = c(rep(1,20),seq(2,Nyrs-1))
  iz = numeric(length = length(iy)); iz[1:19] = 1
  npsv = 0.95
  psiadj = c(0,0,-1)
  # Create some arrays for saving results
  N_Predict = array(dim=c(Nyrs,1))
  P_Predict = array(dim=c(Nyrs-1,1))
  PrPredict = array(dim=c(Nyrs-1,1))
  PrAgPred = array(dim=c(Nages,1))  
  Agepredict1 = array(dim=c(NCages,1))
  Agepredict2 = array(dim=c(NCages,1))
  HV0_predict = array(dim=c(Nyrs-1,1))
  HVA_predict = array(dim=c(Nyrs-1,1))
  sadmean = array(dim = c(Nages,1))
  #
  # Loop through random iterations of model
  #  using Parallel processing to speed things up: try dopar
  simresults = foreach(r = 1:reps) %dopar% {
    mcmc = mcmc1
    vn = vn1
    rr = r_vec[r]
    rrr = r_vec2[r]
    # Extract parameters from joint posterior
    pars = init_fun(rr)
    for(i in 1:length(pars)){
      ##first extract the object value
      tempobj=as.numeric(pars[[i]])
      ##now create a new variable with the original name of the list item
      eval(parse(text=paste(names(pars)[[i]],"= tempobj")))
    }
    # Initialize variables
    if(nrow(ICr)>=1000){
      IC = ICr[rrr:(rrr+Nyrs),]
      CI = CIr[rrr:(rrr+Nyrs)]
    }
    N = numeric(length = Nyrs)
    Nml = numeric(length = Nyrs)
    n = array(dim = c(Nages,Nyrs))
    Fc = array(dim = c(Nages,Nyrs-1))
    S0 = numeric(length = Nyrs-1)
    SA = array(dim = c(Nages,Nyrs-1))
    HV0_pred = numeric(length = Nyrs-1)
    HVA_pred = numeric(length = Nyrs-1)
    FmAgeVec = array(dim = c(NCages,Nyrs-1))
    PredPup = numeric(length = Nyrs-1)
    PredPupA = array(dim = c(Nareas,Nyrs-1))
    # Harvest mort: depends on type of scenario
    if(futuresim==0){
      # Increased stochasticity to mimic effects of autocorrelation
      epsF = epsFr[rrr,]*as.numeric(sigF)*2 + .5 
      epsS = epsSr[rrr,]*as.numeric(sigS)*2 - .5
      # Level of harvest as observed
      gamma_H0 = gamma_H0
      gamma_HA = gamma_HA
    }else if(futuresim==1){
      # Increased stochasticity to mimic effects of autocorrelation
      epsF = epsFr[rrr,]*as.numeric(sigF)*2 + .5 
      epsS = epsSr[rrr,]*as.numeric(sigS)*2 - .5 
      # No harvest hazards (for estimating K)
      gamma_H0 = rep(0,Nyrs-1)
      gamma_HA = rep(0,Nyrs-1)
    }else if(futuresim==2){
      # Increased stochasticity to mimic effects of autocorrelation
      epsF = epsFr[rrr,]*as.numeric(sigF)
      epsS = epsSr[rrr,]*as.numeric(sigS)
      # Range of harvest levels, for finding TAC criteria
      # (then harvest rate remains constant for years within a sim)
      gamma_H0 = rep(gamma_H0_mn*ig1[r],Nyrs-1)
      gamma_HA = rep(gamma_HA_mn*ig2[r],Nyrs-1)
    }
    PA = PAr[rrr:(rrr+Nyrs),]
    N[1] = as.numeric(N0pri) 
    Nml[1] = N[1]/1000000 
    gamma_D_scale = (1/(ages + 0.5))^as.numeric(zeta)
    # Loop through years to calculate demographic transitions and population dynamics 
    for (ix in 1:length(iy)){
      # npsv = 1-(gtools::inv.logit(rnorm(1,0,1))/10) 
      i = iy[ix]
      if (i == 1){
        n[,i] = sad0 * N[i]
      }
      #  Declare some temporary variables for this year:
      haz_0 = numeric()
      haz_A = numeric(length = Nages)
      gamma_I = numeric(length = Nareas)
      haz_I_A = numeric(length = Nareas)
      haz_I = numeric()
      haz_H0 = numeric()
      haz_HA = numeric()
      prp_HV_0 = numeric()
      prp_HV_A = numeric(length = Nages)
      nt1 = numeric(length = Nages)
      #  Fecundity (compared with observed pregnancy rate: adjust for abortions?)
      Fc[1:3,i] = rep(0,3) 
      Fc[4:8,i] = gtools::inv.logit(beta1 - beta2 * (8-ages[4:8])^2 - (phiF*Nml[i])^thtaF - dlta[1]*CI[i]*(Nml[i]) + epsF[i])
      Fc[9:Nages,i] = rep(Fc[8,i],(Nages-8)) 
      #  Juvenile competing hazards and net realized survival
      gamma_D = (phiS*Nml[i])^thtaS 
      haz_0 = exp(omega + gamma_0 + gamma_D + dlta[2] * CI[i+1] + epsS[i]) ; #  + epsJ[i]
      for(a in 1:Nareas){ 
        gamma_I[a] = 8 * exp((psi1+psiadj[a]) - (IC[i,a] * psi2)) / (1 + exp((psi1+psiadj[a]) - (IC[i,a] * psi2))) 
        haz_I_A[a] = exp(omega + gamma_I[a]) 
      } 
      haz_I = sum(PA[i,] * haz_I_A)
      haz_H0 = exp(omega + gamma_H0[i]) 
      S0[i] = exp(-1 * (haz_0 + haz_I + haz_H0)) 
      #  Adult competing hazards and net realized survival
      haz_A = exp(omega + gamma_A + gamma_D * gamma_D_scale)
      haz_HA = exp(omega + gamma_HA[i]) 
      SA[,i] = exp(-1 * (haz_A + haz_HA))     
      #  Calculations for proportional mortality from harvest
      # NOTE: ACCOUNT FOR SNL FOR TAC CALCS? (futuresim2)
      prp_HV_0 = haz_H0 / (haz_0 + haz_I + haz_H0) 
      HV0_pred[i] = sum((n[,i] * Fc[,i]) * 0.5 * npsv) * (1 - S0[i]) * prp_HV_0 
      if (futuresim==0){
        HV0_pred[i] =  HV0_pred[i] * .9
      }      
      if (futuresim==2){
        HV0_can = HV0_pred[i] - Mort_othP # [r,i]
        HV0_pred[i] = HV0_can * .95
      }
      prp_HV_A = haz_HA / (haz_A + haz_HA) 
      HVA_pred[i] = sum((n[,i] * (1 - SA[,i])) * prp_HV_A) 
      if (futuresim==0){
        HVA_pred[i] =  HVA_pred[i] * .5
      }
      if (futuresim==2){
        HVA_can = HVA_pred[i] - Mort_othA # [r,i]
        HVA_pred[i] = HVA_can * .5
      }    
      #  Annual demographic transitions
      nt1[1] = sum((n[,i] * Fc[,i]) * 0.5 * npsv) * (S0[i])
      nt1[2:(Nages-1)] = n[1:(Nages-2),i] * SA[1:(Nages-2),i]
      nt1[Nages] = n[Nages-1,i] * SA[Nages-1,i] + n[Nages,i] * SA[Nages,i] ;    
      if (iz[ix] == 1){
        sad0 = nt1/sum(nt1)
      }else{
        n[,i+1] = nt1  ;
        N[i+1] = sum(n[,i+1]) ;
        Nml[i+1] = N[i+1]/1000000 ;    
        #  Female predicted age vector (ages 3 - 8+), for year i
        FmAgeVec[,i] = n[NCage1:Nages,i] / (0.00001 + sum( n[NCage1:Nages,i])) ;
        #  Predicted Pups to be surveyed (allowing for early ice mortality), 
        #    by area and total, for year i 
        for (a in 1:Nareas){
          PredPupA[a,i] = sum((n[,i] * Fc[,i]) * 0.5 * npsv * PA[i,a] ) ;
        }
        PredPup[i] = sum(PredPupA[,i])
      }
    }
    N_Predict[,1] = N  
    P_Predict[,1] = PredPup  
    PrPredict[,1] = Fc[8,]
    PrAgPred[,1] = Fc[1:Nages,1]
    # Age dist for sample years 25 = 1975, 69 = 2019
    if(futuresim==0 & Nyrs>69){
      Agepredict1[,1] = FmAgeVec[,25]
      Agepredict2[,1] = FmAgeVec[,69]    
    }else{
      Agepredict1[,1] = FmAgeVec[,1]
      Agepredict2[,1] = FmAgeVec[,Nyrs-1]     
    }
    HV0_predict[,1] = HV0_pred
    HVA_predict[,1] = HVA_pred
    sadmean[,1] = sad0
    result <- multiResultClass()
    result$result1 <- N_Predict
    result$result2 <- P_Predict
    result$result3 <- PrPredict
    result$result4 <- PrAgPred
    result$result5 <- Agepredict1
    result$result6 <- Agepredict2
    result$result7 <- HV0_predict
    result$result8 <- HVA_predict
    result$result9 <- sadmean
    return(result)
  }
  #
  N_Predict = array(dim=c(Nyrs,reps))
  P_Predict = array(dim=c(Nyrs-1,reps))  
  PrPredict = array(dim=c(Nyrs-1,reps))
  PrAgPred = array(dim=c(Nages,reps))
  Agepredict1 = array(dim=c(NCages,reps))
  Agepredict2 = array(dim=c(NCages,reps))
  HV0_predict = array(dim=c(Nyrs-1,reps))
  HVA_predict = array(dim=c(Nyrs-1,reps))
  sadmean = array(dim = c(Nages,reps))
  for(r in 1:reps){
    N_Predict[,r] = simresults[[r]]$result1[,1]
    P_Predict[,r] = simresults[[r]]$result2[,1]  
    PrPredict[,r] = simresults[[r]]$result3[,1]  
    PrAgPred[,r] = simresults[[r]]$result4[,1]
    Agepredict1[,r] = simresults[[r]]$result5[,1]  
    Agepredict2[,r] = simresults[[r]]$result6[,1] 
    HV0_predict[,r] = simresults[[r]]$result7[,1]
    HVA_predict[,r] = simresults[[r]]$result8[,1]  
    sadmean[,r] = simresults[[r]]$result9[,1]   
  }
  sadmn = rowMeans(sadmean)
  results = list(N_Predict=N_Predict,P_Predict=P_Predict,PrPredict=PrPredict,
                 PrAgPred=PrAgPred,Agepredict1=Agepredict1,Agepredict2=Agepredict2,
                 HV0_predict=HV0_predict,HVA_predict=HVA_predict,SAD=sadmn)
  return(results)  
}
