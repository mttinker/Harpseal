// HSmodfit
// This Stan program describes a harp seal population model 
//  [More explanation to go here] 
//
// Data inputs for model
data {
  int<lower=1> NPcts;             // N counts of pups, total (across 3 areas)
  int<lower=1> NFages;            // N female age classes (3 - 8+) for age composition obs 
  int<lower=1> NFage1;            // first age class for female age composition (default 3)
  int<lower=1> NFobs;             // N observations of age composition vectors
  int<lower=1> NPRobs;            // N obs of age-specifc female preg status (binomial)
  int<lower=1> Nyrs;              // N years for model to run (year 1 = 1951, year T = 2020?)
  int<lower=1> Nages;             // N age classes (default is 1:20)
  int<lower=1> Nareas;            // N breeding areas (default 3: Sglf, Nhglf, Front)
  vector[Nages] ages;             // vector of age class values 
  vector[Nages] ages2;            // vector of age class values squared
  simplex[Nages] sad0;            // vector of initial stable age distribution values
  real IC[Nyrs,Nareas];           // vector of annual ice anomaly values, by breeding area
  real CE[Nyrs];                  // array of annual environmental index values
  int<lower=0> AgeF[NFobs,NFages];// adult female samples, counts by age class 
  int<lower=0> NFsamp[NPRobs];    // Number females sampled for preg status (by year, age)  
  int<lower=0> Nprg[NPRobs];      // Number pregnant females per sample (by year, age)
  real<lower=0> HVp[Nyrs];        // Combined harvest/bycatch values, beaters (age=1,YoY)
  real<lower=0> HVa[Nyrs];        // Combined harvest/bycatch values, age>1 
  real<lower=0> NP[NPcts];        // pup counts, total
  real<lower=0> sdNP[NPcts];      // s.d. associated with total pup counts
  int<lower=1> YrPct[NPcts];      // Year of each pup count (total counts)
  int<lower=0> PAtag[Nyrs-1];     // switch variable: 1 on years w. area-specific counts   
  matrix[Nyrs-1,Nareas] NPA;      // observed pup proportions by area & year 
  int<lower=1> YrAGsmp[NFobs];    // Year of each female age composition sample 
  int<lower=1> YrPRsmp[NPRobs];   // Year of each female pregnancy status sample 
  int<lower=1> AgePRsmp[NPRobs];  // Age of each female pregnancy status sample
  real N0pri ;                    // prior estimate of starting abundance
  real gamm0;                     // Nuiscence param: base log hazards
  real DDadlt[Nages] ;            // age-sepcific scaling factor for adult D-D 
  real<lower=0> b0pri;            // Fecundity: Logit of max adult pregancy rate
  real<lower=0> psipri1a ;        // Ice anomaly effect on pup survival, fxn param 1 mn
  real<lower=0> psipri1b ;        // Ice anomaly effect on pup survival, fxn param 2 sd
  real<lower=0> psipri2a ;        // Ice anomaly effect on pup survival, fxn param 1 mn
  real<lower=0> psipri2b ;        // Ice anomaly effect on pup survival, fxn param 2 sd
  matrix[Nareas-1,2] PApri;       // beta prior estimates for proportional pups by area
  real CV_HV ;                    // CV associated with harvest counts
  real Adloghz ;                  // adult log hazard rate prior (FIXED value)
  real Jvloghz ;                  // juvenile log hazard rate prior
}
transformed data {
  real gammA;
  real gammJ;
  gammA = Adloghz;
  gammJ = Jvloghz ;    
}
// The parameters to be estimated 
parameters {
  // Annual proportion of pups in areas 1 and 2, predicted
  vector<lower=0,upper=0.45>[Nyrs-1] PA1 ; 
  vector<lower=0,upper=0.45>[Nyrs-1] PA2 ; 
  real<lower=0,upper=4> thta;     // theta parameter: controls "shape" of DD 
  real<lower=0> phiJ;             // Juvenile survival D-D effects 
  real<lower=0> phiF;             // Fecundity (preg rate) D-D effects
  real<lower=0,upper=4> b0;       // Fecundity: Logit of max adult pregancy rate
  real<lower=0,upper=0.5> b1;     // Fecundity: age effect (reduced for younger)
  real<lower=0> psi1;             // Ice anomaly effect on pup survival, fxn param 1 
  real<lower=0> psi2;             // Ice anomaly effect on pup survival, fxn param 2 
  real<lower=100000> N0;          // Initial Abundance, year 1
  real dlta;                      // Effect of environmental conditions on fecundity
  real<lower=0> sigF;             // Environmental stocasticity, var in pregnancy rates
  real epsF[Nyrs-1];              // Stochastic effects on fecundity, by year
  real nphz[Nyrs-1];              // Newborn pup hazards (e.g. abortion, abandonment)
  real<lower=0> sigH;             // Variance in Harvest log hazard rate
  real<lower=0> gammHp_mn ;       // Mean Harvest log hazard rate, Juvenile
  real<lower=0> gammHa_mn;        // Mean Harvest log hazard rate, Adult
  real gammHp[Nyrs-1];            // Annual harvest log hazard rate, Juvenile (beaters)
  real gammHa[Nyrs-1];            // Annual harvest log hazard rate, Adults (age 1+)
}
// Additional transformed parameters, including calculations for model dynamics
transformed parameters {
  matrix[Nyrs-1,Nareas] PA;             // predicted pup proportions by area & year 
  real<lower=0> N[Nyrs] ;               // Population abundance by year
  real<lower=0> Nml[Nyrs] ;             // Population abundance by year in millions
  vector<lower=0>[Nages] n[Nyrs] ;      // Population vector, by year
  vector<lower=0>[Nages] Fc[Nyrs-1] ;    // Fecundity vector, by year
  vector<lower=0>[Nareas] Sice[Nyrs-1] ;// Juv survival from ice anomalies (cause specific)
  real<lower=0> SJ[Nyrs-1] ;            // Juv Survival from all competing hazards 
  vector<lower=0>[Nages] S[Nyrs-1];     // Adult Survival from all competing hazards
  real<lower=0> HVp_pred[Nyrs-1] ;      // predicted harvest numbers by year, pups
  real<lower=0> HVa_pred[Nyrs-1] ;      // predicted harvest numbers by year, adults (age1+)
  vector<lower=0>[NFages] FmAgeVec[Nyrs-1] ; // Age comp vector for adult females ages 3 - 8+
  real<lower=0> PredPup[Nyrs-1] ;       // Predicted pups available for counting, by year
  PA[,1] = PA1;
  PA[,2] = PA2;
  PA[,3] = rep_vector(1,Nyrs-1) - (PA1 + PA2);
  N[1] = N0 ;
  Nml[1] = N[1]/1000000 ;
  n[1] = sad0 * N[1];
  // Loop through years to calculate demographic transitions and population dynamics 
  for (i in 1:(Nyrs-1)){
    // Declare some temporary variables for this year:
    real haz_J ;                  // Juv baseline hazards
    vector[Nages] haz_A ;         // Adult baseline hazards    
    row_vector[Nareas] gammIce ;  // Area specific log hazards from ice
    row_vector[Nareas] haz_I_A ;  // Area specific hazards from ice
    real haz_I ;                  // Weighted mean ice hazards (juvenile)
    real haz_HVp ;                // Harvest hazards, juvenile  
    real haz_HVa ;                // Harvest hazards, adult     
    real prp_HV_p ;               // Proportion of mortality comprised of harvest, Juve 
    vector[Nages] prp_HV_a ;      // Proportion of mortality comprised of harvest, Adult 
    vector[Nages] nt1 ;           // pop vector (temporary)
    row_vector[Nareas] PAi ;      // annual proportional pup allocation to area
    // row_vector[Nareas] PredPupA ; // annual number pups counted by area
    real npsv ;                   // Early pup survival (pre-survey)
    // Use observed pup allocation for years available, otherwise use estimated
    if (PAtag[i]==1){
      PAi =  NPA[i,1:Nareas] ;
    }else{
      PAi =  PA[i,1:Nareas] ;
    }
    npsv = 1-(inv_logit(nphz[i])/10) ;
    // Fecundity (compared with observed pregnancy rate: adjust for abortions?)
    Fc[i][1:2] = rep_vector(0,2) ;
    Fc[i][3:Nages] = inv_logit(b0 - b1 * square(Nages-ages[3:Nages]) - pow(phiF*Nml[i],thta) 
                + dlta*CE[i]*(Nml[i]) + epsF[i]) ;
    // Juvenile competing hazards and net realized survival
    haz_J = exp(gamm0 + gammJ + pow(phiJ*Nml[i], thta)) ; // + epsJ[i]
    for(a in 1:Nareas){ 
       gammIce[a] = psi1 * pow((1/exp(IC[i,a])), psi2) - 1 ;
       Sice[i][a] = exp(-1 * exp(gamm0 + gammIce[a])) ;
       haz_I_A[a] = exp(gamm0 + gammIce[a]) ;
    } 
    haz_I = sum(PAi .* haz_I_A);
    haz_HVp = exp(gamm0 + gammHp[i]) ;
    SJ[i] = exp(-1 * (haz_J + haz_I + haz_HVp)) ;
    // Adult competing hazards and net realized survival
    for (a in 1:Nages){
      haz_A[a] = exp(gamm0 + gammA + pow(DDadlt[a]*phiJ*Nml[i],thta) ) ;
    }
    haz_HVa = exp(gamm0 + gammHa[i]) ;
    S[i] = exp(-1 * (haz_A + haz_HVa)) ;    
    // Calculations for proportional mortality from harvest
    prp_HV_p = haz_HVp / (haz_J + haz_I + haz_HVp) ;
    HVp_pred[i] = sum((n[i] .* Fc[i])) * 0.5 * npsv * (1 - SJ[i]) * prp_HV_p ;
    prp_HV_a = haz_HVa ./ (haz_A + haz_HVa) ;
    HVa_pred[i] = sum((n[i] .* (1 - S[i])) .* prp_HV_a) ;
    // Annual demographic transitions
    nt1[1] = sum((n[i] .* Fc[i])) * 0.5 * npsv * SJ[i] ;
    nt1[2:(Nages-1)] = n[i][1:(Nages-2)] .* S[i][1:(Nages-2)] ;
    nt1[Nages] = n[i][Nages-1] * S[i][Nages-1] + n[i][Nages] * S[i][Nages] ;    
    n[i+1] = nt1  ;
    N[i+1] = sum(n[i+1]) ;
    Nml[i+1] = N[i+1]/1000000 ;    
    // Female predicted age vector (ages 3 - 8+), for year i
    FmAgeVec[i] = n[i][NFage1:Nages] / (0.00001 + sum(n[i][NFage1:Nages])) ;
    // Predicted Pups to be surveyed (allowing for some early season ice mortality), 
    //     by area and total, for year i 
    // for (a in 1:Nareas){
    //   PredPupA[a] = sum((n[i] .* Fc[i])) * 0.5 * npsv * PAi[a];//  * sqrt(Sice[i][a])
    // }
    PredPup[i] = sum((n[i] .* Fc[i])) * 0.5 * npsv ;
  }
}
// The model parameters, estimated by fitting to data
model {
  // Observed nodes:
  // Harvest estimates
  for (i in 1:(Nyrs-1)){
    HVp[i] ~ normal(HVp_pred[i],HVp_pred[i]*CV_HV);
    HVa[i] ~ normal(HVa_pred[i],HVa_pred[i]*CV_HV);
  }
  // Annual pup counts, total 
  for (i in 1:NPcts){
    NP[i] ~ normal(PredPup[YrPct[i]], sdNP[i]) ;
  }
  // Female age distributions (multinomial dist)
  for(i in 1:NFobs){
    AgeF[i,] ~ multinomial(FmAgeVec[YrAGsmp[i]]) ;
  }
  // Female pregancy status (binomial dist)
  for(i in 1:NPRobs){
    Nprg[i] ~ binomial(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
  //
  //Prior distributions for model parameters:
  // Initial Population size (stochastic node, weak prior):
  N0 ~ gamma( square(N0pri)/square(N0pri*.25) , N0pri/square(N0pri*.25));  
  // Random effects: fecundity and newborn pup survival (env. stochasticity) 
  epsF ~ normal(0,sigF) ;
  nphz ~ normal(0,1) ;
  // Annual harvest log hazard rates
  gammHp ~ normal(gammHp_mn,sigH);
  gammHa ~ normal(gammHa_mn,sigH);
  // Estimated annual proportion of pups for areas 1 and 2 
  //  for years not observed: Area 3 calculated as 1-(PA1+PA2) 
  PA1 ~ beta(PApri[1,1],PApri[1,2]) ;
  PA2 ~ beta(PApri[2,1],PApri[2,2]) ;
  // key parameters of interest
  thta ~ gamma(3,2) ;
  phiJ ~ cauchy(0,0.5) ;
  phiF ~ cauchy(0,0.5) ;
  b0 ~ normal(b0pri,.5);
  b1 ~ normal(0.2,0.1) ;
  psi1 ~ normal(psipri1a,psipri1b) ;
  psi2 ~ normal(psipri2a,psipri2b) ;
  dlta ~ cauchy(0,.25) ;
  sigF ~ cauchy(0,1.5) ;
  sigH ~ cauchy(0,1.5) ;
  gammHp_mn ~ normal(5.9,2) ;
  gammHa_mn ~ normal(3.3,1.5) ;
}
// Derived params (e.g. goodness of fit stats)
generated quantities {
  real PAmean[Nareas] ;  /// mean PA values
  real Fc8_prdct[Nyrs-1] ;  // Predicted fecundity rate for 8+ over years
  vector[Nages] Fc1966_prdct ;  // Predicted fecundity rate for all ages, 1966
  vector[Nages] Fc2016_prdct ;  // Predicted fecundity rate for all ages, 2016
  real log_lik[NPcts+NPRobs] ;    // Log liklihood of obs. data (for LooIC)
  real y_new1[NPcts]  ;    // New observations
  real y_new2[NPRobs]  ;    // New observations  
  real P_resid1[NPcts]  ;   // Pearson residuals from observed data
  real P_resid_new1[NPcts]; // Pearson residuals from new data
  real P_resid2[NPRobs]  ;   // Pearson residuals from observed data
  real P_resid_new2[NPRobs]; // Pearson residuals from new data
  real Tstat1 ;      // Test stat for PPC (sum of squared pearson resids)
  real Tstat_new1 ;  // New data test stat (sum of squared pearson resids)
  real Tstat2 ;      // Test stat for PPC (sum of squared pearson resids)
  real Tstat_new2 ;  // New data test stat (sum of squared pearson resids)
  real Tstat ;      // Test stat for PPC (sum of squared pearson resids)
  real Tstat_new ;  // New data test stat (sum of squared pearson resids)
  real ppp1 ;        // Posterior predictive check Bayesian P-value
  real ppp2 ;        // Posterior predictive check Bayesian P-value
  real ppp ;        // Posterior predictive check Bayesian P-value
  int<lower=0> c;
  for (i in 1:(Nyrs-1)){
    Fc8_prdct[i] = Fc[i][8] ;
  }
  Fc1966_prdct = Fc[16][1:Nages];
  Fc2016_prdct = Fc[66][1:Nages];
  c = 0 ;
  for (i in 1:NPcts){
    c = c+1 ;
    log_lik[c] = normal_lpdf(NP[i] | PredPup[YrPct[i]], sdNP[i]) ;
    y_new1[i] = normal_rng(PredPup[YrPct[i]], sdNP[i]) ;
    P_resid1[i] = (NP[i] - PredPup[YrPct[i]]) / sdNP[i];
    P_resid_new1[i] = (y_new1[i] - PredPup[YrPct[i]]) / sdNP[i];
  }
  for(i in 1:NPRobs){
    c = c+1 ;
    log_lik[c] = binomial_lpmf( Nprg[i] | NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    y_new2[i] =  binomial_rng(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    P_resid2[i] = (Nprg[i] - NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]) / sqrt(NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]*(1 - Fc[YrPRsmp[i]][AgePRsmp[i]])) ;
    P_resid_new2[i] = (y_new2[i] - NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]) / sqrt(NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]*(1 - Fc[YrPRsmp[i]][AgePRsmp[i]])) ;
  }
  Tstat1 = sum(square(P_resid1))/NPcts ;
  Tstat_new1 = sum(square(P_resid_new1))/NPcts ;
  Tstat2 = sum(square(P_resid2))/NPRobs ;
  Tstat_new2 = sum(square(P_resid_new2))/NPRobs ;  
  Tstat=Tstat1+Tstat2 ;
  Tstat_new=Tstat_new1+Tstat_new2 ;
  ppp1 = Tstat1 > Tstat_new1 ? 1 : 0;
  ppp2 = Tstat2 > Tstat_new2 ? 1 : 0;  
  ppp = Tstat > Tstat_new ? 1 : 0;
}