// HSmodfit
// This Stan program describes a harp seal population model 
//  [More explanation here] ** NOTE: make harvest an observed node! 
//
// Data inputs for model
data {
  int<lower=1> NPcts;             // N counts of pups, total (across 3 areas)
  int<lower=1> NPctsA;            // N counts of pups by area (Sglf, Nhglf, Front)
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
  real<lower=0> HVp[Nyrs];        // Combined harvest/bycatch values, beaters (age=1,YoY)
  real<lower=0> HVa[Nyrs];        // Combined harvest/bycatch values, age>1 
  real<lower=0> NP[NPcts];        // pup counts, total
  real<lower=0> NPA[NPctsA];      // pup counts by area
  real<lower=0> sdNP[NPcts];      // s.d. associated with total pup counts
  real<lower=0> sdNPA[NPctsA];    // s.d. associated with pup counts by area
  int<lower=0> AgeF[NFobs,NFages];// adult female samples, counts by age class 
  int<lower=0> NFsamp[NPRobs];    // Number females sampled for preg status (by year, age)  
  int<lower=0> Nprg[NPRobs];      // Number pregnant females per sample (by year, age)
  int<lower=1> YrPct[NPcts];      // Year of each pup count (total counts)
  int<lower=1> YrPctA[NPctsA];    // Year of each pup count (area-specific)
  int<lower=1> AreaPctA[NPctsA];  // Area of each area-specific pup count 
  int<lower=1> YrAGsmp[NFobs];    // Year of each female age composition sample 
  int<lower=1> YrPRsmp[NPRobs];   // Year of each female pregnancy status sample 
  int<lower=1> AgePRsmp[NPRobs];  // Age of each female pregnancy status sample
  real a0;                        // Nuiscence param: base hazards
  real thta;                      // theta param for theta-logistic density dependence
  vector<lower=0>[Nareas] PApri;  // prior estimate, proportional pup allocation by area
  real CV_HV ;                    // CV associated with harvest counts
  real N0pri;                     // prior for mean estimated abundance at year1 (1951)
  real N0cvpri;                   // cv for prior for mean estimated abundance at year1 
  real aA0pri;                    // prior for alpha param aA0: adult hazard ratio (intcp)
  real betapri;                   // prior for beta param (b0), logit max pregnancy rate 
  real psipri1a;                  // prior for psi param1 mean: ice anomaly effect fxn 
  real psipri1b;                  // prior for psi param1 sd: ice anomaly effect fxn 
  real psipri2a;                  // prior for psi param2 mean: ice anomaly effect fxn 
  real psipri2b;                  // prior for psi param2 sd: ice anomaly effect fxn 
  simplex[NFages] FageAv ;        // back-up probs for adult fem age distribution 
}
// The parameters to be estimated 
parameters {
  simplex[Nareas] PA[Nyrs-1] ;      // Annual proportional pup distribution across areas
  real<lower=0> aJ;               // Juvenile survival hazard ratio
  real<lower=0> aA0;              // Adult survival hazard ratio, intercept
  real<lower=0> aA1;              // Adult survival hazard ratio, age effect linear
  real<lower=0> aA2;              // Adult survival hazard ratio, age effect quadratic
  real<lower=0> phiJ;             // Juvenile survival, D-D param
  real<lower=0> phiF;             // Fecundity (preg rate), D-D param
  real<lower=0> b0;               // Fecundity: Logit of max adult pregancy rate
  real<lower=0> b1;               // Fecundity: Logit of age effect
  real<lower=0,upper=4> psi1;     // Ice anomaly effect on pup survival, fxn param 1 
  real<lower=0,upper=2> psi2;     // Ice anomaly effect on pup survival, fxn param 2 
  real<lower=500000> N0;          // Initial Abundance, year 1
  real dlta;                      // Effect of environmental conditions on fecundity
  real<lower=0> sigJ;             // Environmental stocasticity, var in juvenile survival
  real<lower=0> sigF;             // Environmental stocasticity, var in pregnancy rates
  real epsF[Nyrs-1];              // Stochastic effects on fecundity, by year
  real epsJ[Nyrs-1];              // Stochastic effects on juv survival, by year
  real<lower=0> sigH[2];          // Variance in Harvest hazard rate, 1=Juv, 2=Adult
  real<lower=0> mdnH[2];          // Median Harvest hazard rate, 1=Juv, 2=Adult
  real<lower=0> zta1[Nyrs-1];     // Annual harvest hazard rate, Juvenile (beaters)
  real<lower=0> zta2[Nyrs-1];     // Annual harvest hazard rate, Adults (age 1+)
  # real<lower=0> sigHC[2];         // Estimation error in Harvest #s, 1=beaters, 2=adults
}
// Additional transformed parameters
transformed parameters {
  real<lower=0> N[Nyrs] ;               // Population abundance by year
  real<lower=0> Nml[Nyrs] ;             // Population abundance by year
  vector<lower=0>[Nages] n[Nyrs] ;      // Population vector, by year
  vector<lower=0>[Nages] F[Nyrs-1] ;    // Fecundity vector, by year
  real gammJ[Nyrs-1] ;                  // Juv summed hazards, by year
  vector<lower=0>[Nareas] Sice[Nyrs-1] ;// Juv Survival mod factor due to ice anomalies
  real<lower=0> SJ[Nyrs-1] ;            // Juv Survival from hazards plus ice, pre-harvest
  vector[Nages] gammA[Nyrs-1];          // Adult sum hazards, by year
  vector<lower=0>[Nages] S[Nyrs-1];     // Adult Survival, all effects pre harvest
  vector<lower=0>[Nages] R[Nyrs-1]; // Recruitment rate, pre harvest, by year
  simplex[NFages] FmAgeVec[Nyrs-1] ;    // Age comp vector for adult females 3-8+
  real<lower=0> PredPup[Nyrs-1] ;       // Predicted pups available for counting, by year
  vector<lower=0>[Nareas] PredPupA[Nyrs-1];// Predicted pups by area, by year
  N[1] = N0 ;
  Nml[1] = N[1]/1000000 ;
  n[1] = sad0 * N[1];
  // Loop through years to calculate demographic transitions and population dynamics 
  for (i in 1:(Nyrs-1)){
    // Declare some temporary variables:
    //vector[Nages-1] ageDstH;
    //vector[Nages] HV;
    vector[NFages] Fage ;
    vector[Nages] nt1; 
    // matrix[Nages,Nages] M;
    F[i][1:2] = rep_vector(0,2) ;
    F[i][3:8] = inv_logit(b0 - b1 * square(8-ages[3:8]) - pow(phiF*Nml[i],thta) 
                + dlta*CE[i]*(Nml[i]/2) + epsF[i]) ;
    F[i][9:20] = rep_vector(F[i][8],12) ;
    gammJ[i] = a0 + aJ + zta1[i] + pow(phiJ*Nml[i], thta) + epsJ[i] ;
    for(a in 1:Nareas){
      Sice[i][a] = inv_logit(-psi1 + 7 * pow((IC[i,a] + 1), psi2)) ;
    }
    SJ[i] = exp(-exp(gammJ[i])) * sum(PA[i] .* Sice[i]) ;
    R[i] = 0.5 * (F[i] * SJ[i]) ;
    gammA[i] = a0 + aA0 + zta1[i] - aA1*ages + aA2*ages2 ;
    S[i] = exp(-exp(gammA[i])) ;  
    //
    // ** REPLACE MATRIX BY SIMPLER ALGEBRAIC CALCS OF TRANSITIONS
    // M[1] = R[i] ;
    // M[2:Nages,1:(Nages-1)] = diag_matrix(S[i][1:(Nages-1)]) ;
    // M[2:(Nages-1),Nages] = rep_vector(0,(Nages-2)) ;
    // M[Nages,Nages] = S[i][Nages] ;
    // ageDstH = n[i][2:Nages] / (0.000001+sum(n[i][2:Nages])) ;
    // HV[1] = HVp[i];
    // HV[2:Nages] = HVa[i] * ageDstH;
    // Matrix multiplication for demographic transitions, then add harvest mortality:
    // nt1 = (M * n[i]) - HV;
    // Next loop prevents negative "n" values ("overharvest""): ideally we can remove
    // for (a in 1:Nages){
    //   n[i+1][a] = fmax(nt1[a],0) ;
    // }
    // n[i+1] = nt1 ; // Alternative if above error-catch loop appears unnecessary
    nt1[1] = sum(n[i] .* R[i]) ;
    nt1[2:(Nages-1)] = n[i][1:(Nages-2)] .* S[i][1:(Nages-2)] ;
    nt1[Nages] = n[i][Nages-1] * S[i][Nages-1] + n[i][Nages] * S[i][Nages] ;
    // ageDstH = nt1[2:Nages] / (0.000001+sum(nt1[2:Nages])) ;
    // HV[1] = HVp[i];
    // HV[2:Nages] = HVa[i] * ageDstH;
    n[i+1] = nt1  ;
    N[i+1] = sum(n[i+1]) ;
    Nml[i+1] = N[i+1]/1000000 ;
    // Female predicted age vector (ages 3 - 8+), for year i
    Fage[1:(NFages-1)] = n[i][NFage1:(NFage1+NFages-2)] / (0.0000001 + sum(n[i][NFage1:Nages])) ;
    Fage[NFages] = sum(n[i][(NFage1+NFages-1):Nages]) / (0.0000001 + sum(n[i][NFage1:Nages])) ;
    if (sum(Fage)<=0) {
      FmAgeVec[i] = FageAv ;
    }else{
      FmAgeVec[i] = Fage / sum(Fage) ;
    }
    // Predicted Pups, by area and total for year i
    for (a in 1:Nareas){
      PredPupA[i][a] = sum((n[i] .* F[i]) * 0.5 * PA[i][a] * sqrt(Sice[i][a])) ;
    }
    PredPup[i] = sum(PredPupA[i]) ;
  }
}
// The model, parameters estimated by fitting to data
model {
  // Observed nodes:
  // HARVEST NUMBERS: FIGURE OUT HOW TO CALCULATE 
  
  // Annual pup counts, total (gamma dist)
  for (i in 1:NPcts){
    NP[i] ~ normal(PredPup[YrPct[i]], sdNP[i]) ;
  }
  // Annual pup counts, by area (gamma dist)
  for (i in 1:NPctsA){
    NPA[i] ~ normal(PredPupA[YrPctA[i]][AreaPctA[i]], sdNPA[i]) ;
  }  
  // Female age distributions (multinomial dist)
  for(i in 1:NFobs){
    AgeF[i,] ~ multinomial(FmAgeVec[YrAGsmp[i]]) ;
  }
  // Female pregancy status (binomial dist)
  for(i in 1:NPRobs){
    Nprg[i] ~ binomial(NFsamp[i], F[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
  // Random effects: fecundity and Juv survival (env. stochasticity) 
  epsF ~ normal(0,sigF) ;
  epsJ ~ normal(0,sigJ) ;
  // Initial Population size (stochastic node, weak prior):
  // N0 ~ gamma(N0pri^2/(N0pri*N0cvpri)^2, N0pri/(N0pri*N0cvpri)^2) ;
  N0 ~ normal(N0pri, N0pri*N0cvpri);
  // Proportional distribution of pupping across 3 areas for each year
  //   (stochastic node, weak dirichlet prior)
  for (y in 1:(Nyrs-1)){
    PA[y] ~ dirichlet(50*PApri); 
    zta1[y] ~ lognormal(log(mdnH[1]),sigH[1]);
    zta2[y] ~ lognormal(log(mdnH[2]),sigH[2]);
  }  
  // Prior distributions for model parameters
  aJ ~ cauchy(0,2.5) ; // normal(1.25*aA0pri,1) ;
  aA0 ~ cauchy(0,2.5) ;// normal(aA0pri,0.5) ;
  aA1 ~ cauchy(0,0.25) ;
  aA2 ~ cauchy(0,0.025) ;
  phiJ ~ cauchy(0,0.25) ;
  phiF ~ cauchy(0,0.25) ;
  b0 ~ cauchy(0,2.5) ;
  b1 ~ cauchy(0,0.25) ;
  psi1 ~ cauchy(0,2.5) ;
  psi2 ~ cauchy(0,2.5) ;
  dlta ~ cauchy(0,.25) ;
  sigJ ~ cauchy(0,2.5) ;
  sigF ~ cauchy(0,2.5) ;
  sigH[1] ~ cauchy(0,0.5) ;
  sigH[2] ~ cauchy(0,0.5) ;
  mdnH[1] ~ cauchy(0,2.5) ;
  mdnH[2] ~ cauchy(0,2.5) ;
  // sigA[1] ~ cauchy(0,2.5) ;
  // sigA[2] ~ cauchy(0,2.5) ;
  // sigA[3] ~ cauchy(0,2.5) ;
}
// Derived params (e.g. goodness of fit stats)
generated quantities {
  real log_lik[NPcts+NPctsA+NPRobs] ;    // Log liklihood of obs. data (for LooIC)
  real y_new[NPcts+NPctsA+NPRobs]  ;
  real P_resid[NPcts+NPctsA+NPRobs]  ;   // Pearson residuals from observed data
  real P_resid_new[NPcts+NPctsA+NPRobs]; // Pearson residuals from new data
  real Tstat ;      // Test stat for PPC (sum of squared pearson resids)
  real Tstat_new ;  // New data test stat (sum of squared pearson resids)
  real ppp ;        // Posterior predictive check Bayesian P-value
  real PAmean[Nareas] ;  /// mean PA values
  int<lower=0> c;
  for (a in 1:Nareas){
    PAmean[a] = mean(PA[1:(Nyrs-1)][a]) ;
  }
  c = 0 ;
  for (i in 1:NPcts){ 
    c = c+1 ;
    log_lik[c] = normal_lpdf(NP[i] | PredPup[YrPct[i]], sdNP[i]) ;
    y_new[c] = normal_rng(PredPup[YrPct[i]], sdNP[i]) ;
    P_resid[c] = (NP[i] - PredPup[YrPct[i]]) / sdNP[i];
    P_resid_new[c] = (y_new[c] - PredPup[YrPct[i]]) / sdNP[i];    
  }
  for (i in 1:NPctsA){
    c = c+1 ;
    log_lik[c] = normal_lpdf(NPA[i] | PredPupA[YrPctA[i]][AreaPctA[i]], sdNPA[i]) ;
    y_new[c] = normal_rng(PredPupA[YrPctA[i]][AreaPctA[i]], sdNPA[i]) ;
    P_resid[c] = (NPA[i] - PredPupA[YrPctA[i]][AreaPctA[i]]) / sdNPA[i];
    P_resid_new[c] = (y_new[c] - PredPupA[YrPctA[i]][AreaPctA[i]]) / sdNPA[i];
  }
  for(i in 1:NPRobs){
    c = c+1 ;
    log_lik[c] = binomial_lpmf( Nprg[i] | NFsamp[i], F[YrPRsmp[i]][AgePRsmp[i]]) ;
    y_new[c] =  binomial_rng(NFsamp[i], F[YrPRsmp[i]][AgePRsmp[i]]) ;
    P_resid[c] = (Nprg[i] - NFsamp[i]*F[YrPRsmp[i]][AgePRsmp[i]]) / sqrt(NFsamp[i]*F[YrPRsmp[i]][AgePRsmp[i]]*(1 - F[YrPRsmp[i]][AgePRsmp[i]])) ;
    P_resid_new[c] = (y_new[c] - NFsamp[i]*F[YrPRsmp[i]][AgePRsmp[i]]) / sqrt(NFsamp[i]*F[YrPRsmp[i]][AgePRsmp[i]]*(1 - F[YrPRsmp[i]][AgePRsmp[i]])) ;
  }
  Tstat = sum(square(P_resid));
  Tstat_new = sum(square(P_resid_new)) ;
  ppp = Tstat > Tstat_new ? 1 : 0;
}
