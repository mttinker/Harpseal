// HSmodfit
// This Stan code executes a harp seal population model 
//  [Refer to "HSmodel_summary" for details] 
//
// Section 1. Data inputs for model
data {
  int<lower=1> NPcts;             // N counts of pups, total (across 3 areas)
  int<lower=1> NPActs;            // N counts of pups by area 
  int<lower=1> NCages;            // N age classes used for age composition comparisons 
  int<lower=1> NCage1;            // first age class for obs. age composition (5? 7?)
  int<lower=1> NCobs;             // N observations of age composition vectors
  int<lower=1> NPRobs;            // N obs of age-specifc female pregnancy status (binomial)
  int<lower=1> Nyrs;              // N years for model to run (year 1 = 1951, year T = 2020?)
  int<lower=1> Nages;             // N age classes (default is 1:36)
  int<lower=1> Nareas;            // N breeding areas (default 3: Sglf, Nhglf, Front)
  vector[Nages] ages;             // vector of age class values 
  vector[Nages] agesC;            // vector of re-centered age class values 
  vector[Nages] agesC2;           // vector of re-centered age class values squared
  simplex[Nages] sad0;            // vector of initial stable age distribution values
  int<lower=1> IC[Nyrs,Nareas];   // Annual ice anomaly values by area (discretized, 41 steps of 0.05 units)
  real CE[Nyrs];                  // array of annual environmental index values
  int<lower=0> Agects[NCobs,NCages];// matrix of adult counts by age class (cols) & year (rows)
  int<lower=0> NFsamp[NPRobs];    // Number females sampled for pregnancy status (grouped by year, age)  
  int<lower=0> NPrg[NPRobs];      // Number pregnant females per sample (grouped by year, age)
  real<lower=0> H_0[Nyrs-1];      // Combined harvest/bycatch values, beaters (age=1,YoY)
  real<lower=0> H_A[Nyrs-1];      // Combined harvest/bycatch values, age>1 
  real<lower=0> Pups[NPcts];      // pup counts, total
  real<lower=0> sdNP[NPcts];      // s.d. associated with total pup counts
  int<lower=1> YrPct[NPcts];      // Year of each pup count (total counts)
  int<lower=0> PAidx[NPActs];     // index variable: years w. area-specific counts  
  matrix[NPActs,Nareas] NPA;      // observed pup proportions by area & year 
  int<lower=1> YrAGsmp[NCobs];    // Year of each female age composition sample 
  int<lower=1> YrPRsmp[NPRobs];   // Year of each female pregnancy status sample 
  int<lower=1> AgePRsmp[NPRobs];  // Age of each female pregnancy status sample
  real N0pri ;                    // mean for vague prior estimate of starting abundance
  real omega;                     // Nuiscance param: base log hazards
  real psipri1a ;                 // Ice anomaly effect on pup survival, fxn param 1 mn
  real<lower=0> psipri1b ;        // Ice anomaly effect on pup survival, fxn param 1 sd
  real<lower=0> psipri2a ;        // Ice anomaly effect on pup survival, fxn param 2 mn
  real<lower=0> psipri2b ;        // Ice anomaly effect on pup survival, fxn param 2 sd
  matrix[Nareas-1,2] PApri;       // beta prior estimates for proportional pups by area
  real CV_HV ;                    // CV associated with harvest counts
  vector[Nyrs-1] Q_0 ;            // probability that a harvested pup is observed (and not SnL) 
  vector[Nyrs-1] Q_A;             // probability that a harvested adult is observed (and not SnL) 
  vector[41] ICvec ;              // dummy vector of ice anomaly values, -1 to 1
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0> alpha0 ;           // base adult log haz (for 10 yr old)
  real<lower=0> alpha1 ;           // adult log haz age modifier par 1
  real<lower=0> alpha2 ;           // adult log haz age modifier par 2
  real<lower=0> nu ;               // log haz ratio for juv relative to sub-adlt   
  simplex[NCages] pie[NCobs];      // probs of age counts (for multinomial)
  vector<lower=0,upper=0.45>[Nyrs-1] PA1 ; // proportion pups in area 1
  vector<lower=0,upper=0.45>[Nyrs-1] PA2 ; // proportion pups in area 2
  real<lower=0> thta[2];           // theta params: controls "shape" of DD (for Fecundity and survival)
  real<lower=0> phi[2];            // D-D effects for fecundity (1) and survival (2)
  real<lower=0> beta1;             // Fecundity: Logit of max adult pregancy rate
  real<lower=0> beta2;            // Fecundity: age effect (reduced for younger)
  real psi1;                       // Ice anomaly effect on pup survival, fxn param 1 
  real<lower=0> psi2;              // Ice anomaly effect on pup survival, fxn param 2 
  real<lower=.001> N0mil;          // Initial Abundance, year 1 of time series (i millions)
  real dlta;                       // Effect of environmental conditions on fecundity
  real<lower=0> sigF;              // Environmental stocasticity, var in pregnancy rates
  real<lower=0> sigS;              // Environmental stocasticity, var in juv survival rates 
  real<lower=0> sigH;              // Variance in Harvest log hazard rate 
  real<lower=0> tau;               // precision param for dirichlet-multinomial age dist
  real epsF[Nyrs-1];               // Stochastic effects on fecundity, by year
  real epsS[Nyrs-1];               // Stochastic effects on juv survival, by year
  real nphz[Nyrs-1];               // Stochastic Newborn pup hazards (e.g. abortion, abandonment)
  real<lower=0> gamma_H0_mn ;      // Mean Harvest log hazard rate, Juvenile
  real<lower=0> gamma_HA_mn;       // Mean Harvest log hazard rate, Adult
  real gamma_H0[Nyrs-1];           // Annual harvest log hazard rate, Juvenile (beaters)
  real gamma_HA[Nyrs-1];           // Annual harvest log hazard rate, Adults (age 1+)
  real<lower=0,upper=1> zeta ;     // proportional D-D strength for sub-adults relative to juveniles
}
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  real N0 ;                             // initial abundance (year 1)
  vector[Nages] gamma_A ;               // log hazard ratio for adults
  real gamma_0 ;                        // log hazard ratio for juveniles
  matrix[Nyrs-1,Nareas] PA;             // predicted pup proportions by area & year 
  real<lower=0> N[Nyrs] ;               // Population abundance by year
  vector<lower=0>[Nages] n[Nyrs] ;      // Population vector, by year
  vector<lower=0>[Nages] Fc[Nyrs-1] ;   // Fecundity vector, by year
  real<lower=0> S0[Nyrs-1] ;            // Juv Survival from all competing hazards 
  vector<lower=0>[Nyrs-1] Snp ;          // Newborn pup survival (from abortion, abandonment) 
  vector<lower=0>[Nages] SA[Nyrs-1];     // Adult Survival from all competing hazards
  real<lower=0> Pups_pred[Nyrs-1] ;       // Predicted pups available for counting, by year  
  vector<lower=0>[Nyrs-1] H0_pred_ALL ;// predicted harvest numbers by year, pups (including SNL)
  vector<lower=0>[Nyrs-1] HA_pred_ALL ;// predicted harvest numbers by year, adults (age1+) (including SNL)
  vector<lower=0>[Nyrs-1] H0_pred ;    // predicted harvest numbers by year without SnL, pups
  vector<lower=0>[Nyrs-1] HA_pred ;    // predicted harvest numbers by year without SnL, adults
  vector[41] haz_Ice ; // Hazards associated with ice anomalies
  N0 = N0mil * 1000000 ;
  // Assign proportion of pups in each breeding area
  PA[,1] = PA1;
  PA[,2] = PA2;
  PA[,3] = rep_vector(1,Nyrs-1) - (PA1 + PA2);
  // overwrite PA values with observed proportions for years available
  PA[PAidx] = NPA ;
  // Newborn pup survival for each year (stochastic)
  Snp = 1 - (inv_logit(to_vector(nphz)) * 0.1) ;
  // Initize pop vector  
  N[1] = N0 ;                     // Initialize population, year 1
  n[1] = sad0 * N[1];             // Initialize population vector, year 1
  // Calculate adult and juvenile log hazards (baseline). 
  gamma_A = alpha0 - alpha1 * agesC + alpha2 * agesC2 ;
  gamma_0 = alpha0 - alpha1 * agesC[1] + alpha2 * agesC2[1] + nu ;
  // Ice hazard function
  haz_Ice = exp(omega + 8 * exp(psi1 - (ICvec * psi2)) ./ (1 + exp(psi1 - (ICvec * psi2)))) ;
  // Loop through years to calculate demographic transitions and population dynamics 
  for (i in 1:(Nyrs-1)){
    // Declare some temporary variables for this year:
    real Nml ;                    // Current population abundance in millions
    real haz_J ;                  // Juv baseline hazards
    vector[Nages] haz_A ;         // Adult baseline hazards    
    real haz_IC ;                  // Weighted mean ice hazards (juvenile)
    real haz_H0 ;                // Harvest hazards, juvenile  
    real haz_HA ;                // Harvest hazards, adult     
    real prp_HV_0 ;               // Proportion of mortality comprised of harvest, Juve 
    vector[Nages] prp_HV_A ;      // Proportion of mortality comprised of harvest, Adult
    real gamma_D ;               // density-dependent effect for juveniles/subadults
    // current N in millions of anials, for D-D calcs
    Nml = N[i]/1000000 ;
    // Annual Fecundity, with D-D, envir. effects and stochasticity 
    Fc[i][1:3] = rep_vector(0,3) ;
    Fc[i][4:8] = inv_logit(beta1 - beta2 * square(8-ages[4:8]) - pow(phi[1]*Nml,thta[1]) + dlta*CE[i]*Nml + epsF[i]) ;
    Fc[i][9:Nages] = rep_vector( Fc[i][8], (Nages - 8)) ; // ** IF >8 age classes
    // Juvenile competing hazards and net survival
    gamma_D = pow(phi[2]*Nml, thta[2]) ;
    haz_J = exp(omega + gamma_0 + gamma_D + epsS[i] ) ; // normal juv hazards
    // weighted avg ice anom haz (wt = prop. pups in each area)
    haz_IC = PA[i][1] * haz_Ice[IC[i,1]] + PA[i][2] * haz_Ice[IC[i,2]] + PA[i][3] * haz_Ice[IC[i,3]] ;
    haz_H0 = exp(omega + gamma_H0[i]) ; // hazards from bycatch/human harvest of pups/beaters
    S0[i] = exp(-1 * (haz_J + haz_IC + haz_H0)) ;  // combine hazards for Juv survival
    // Adult competing hazards and net survival (allow weak DD effects for sub-adults)
    haz_A = exp(omega + gamma_A + gamma_D * zeta * (1 ./ ages) ) ; // reg adult hazards 
    haz_HA = exp(omega + gamma_HA[i]) ; // harvest/bycatch hazards for adults
    SA[i] = exp(-1 * (haz_A + haz_HA)) ;    
    // Calculate proportional mortality from harvest from pups and adults
    prp_HV_0 = haz_H0 / (haz_J + haz_IC + haz_H0) ;
    prp_HV_A = haz_HA ./ (haz_A + haz_HA) ;
    // Calculate predicted total harvest/bycatch numbers (ALL means including SnL)
    H0_pred_ALL[i] = sum((n[i] .* Fc[i])) * 0.5 * Snp[i] * (1 - S0[i]) * prp_HV_0 ;
    HA_pred_ALL[i] = sum((n[i] .* (1 - SA[i])) .* prp_HV_A) ;
    // Calculate predicted harvest/bycatch numbers, excluding SnL
    H0_pred[i] = H0_pred_ALL[i] * Q_0[i] ;
    HA_pred[i] = HA_pred_ALL[i] * Q_A[i] ;
    // Predicted # Pups to be surveyed (accounts for abortions/early pup survival, Snp), 
    // OPTION: could also subtract some fraction of Ice motality that occurs prior to pup survey: 
    // IceMort = sum((n[i] .* Fc[i])) * 0.5 * Snp[i] * (1 - S0[i]) * (haz_I / (haz_J + haz_I + haz_HVp))
    Pups_pred[i] = sum((n[i] .* Fc[i])) * 0.5 * Snp[i] ; // optionally subtract proportion of IceMort
    // Update population vector and total abundance for next year
    n[i+1][1] = sum((n[i] .* Fc[i])) * 0.5 * Snp[i] * S0[i] ;
    n[i+1][2:(Nages-1)] = n[i][1:(Nages-2)] .* SA[i][1:(Nages-2)] ;
    n[i+1][Nages] = n[i][Nages-1] * SA[i][Nages-1] + n[i][Nages] * SA[i][Nages] ; 
    N[i+1] = sum(n[i+1]) ;    
  }
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  //  Harvest estimates (change to Gamma dist?)
  H_0 ~ normal(H0_pred, H0_pred * CV_HV) ; 
  H_A ~ normal(HA_pred, HA_pred * CV_HV) ;
  //  Annual pup counts (change to Gamma dist?)
  Pups ~ normal(Pups_pred[YrPct], sdNP) ;
  //  Ageclass distributions by year (multinomial dist)
  for(i in 1:NCobs){
    // NOTE: Use dirichlet-multinomial to handle error/variance in age counts 
    pie[i] ~ dirichlet(10 * tau * (n[YrAGsmp[i]][NCage1:Nages] / sum(n[YrAGsmp[i]][NCage1:Nages])));
    Agects[i,] ~ multinomial(pie[i]) ;
  }
  //  Pregancy status of females by age/year (binomial dist)
  for(i in 1:NPRobs){
    NPrg[i] ~ binomial(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
  // B) Prior distributions for model parameters:
  // Hierarchical random effects: NOTE: maybe make annual stochastic effects autoregressive (AR[1])?  
  epsF ~ normal(0,sigF) ;          // Annual deviations from mean fecundity
  epsS ~ normal(0,sigS) ;          // Annual deviations from mean juv hazards
  nphz ~ normal(0,1) ;             // Annual pre-survey newborn hazards
  gamma_H0 ~ normal(gamma_H0_mn,sigH); // Annual pup harvest log hazards
  gamma_HA ~ normal(gamma_HA_mn,sigH); // Annual adult harvest log hazards
  PA1 ~ beta(PApri[1,1],PApri[1,2]) ; // Annual proportion of pups area 1
  PA2 ~ beta(PApri[2,1],PApri[2,2]) ; // Annual proportion of pups area 2
  // Base parameters:
  N0mil ~ normal(0,2.5) ;
  thta ~ normal(1,1) ;
  phi ~ cauchy(0,2.5) ;
  beta1 ~ cauchy(0,2.5) ;
  beta2 ~ cauchy(0,2.5) ;
  psi1 ~ normal(psipri1a,psipri1b) ;
  psi2 ~ normal(psipri2a,psipri2b) ;
  dlta ~ cauchy(0,2.5) ;
  sigF ~ cauchy(0,2.5) ;
  sigH ~ cauchy(0,2.5) ;
  sigS ~ cauchy(0,2.5) ;
  tau ~ cauchy(0,10) ;
  gamma_H0_mn ~ cauchy(0,2.5) ;
  gamma_HA_mn ~ cauchy(0,2.5) ;
  alpha0 ~ cauchy(0,2.5) ;
  alpha1 ~ cauchy(0,2.5) ;
  alpha2 ~ cauchy(0,2.5) ;
  nu ~ cauchy(0,2.5) ;
  zeta ~ beta(1,5) ;
}
// Section 5. Derived parameters and statistics 
generated quantities {
  real PAmean[Nareas] ;  /// mean PA values
  real Fc8_prdct[Nyrs-1] ;  // Predicted fecundity rate for 8+ over years
  vector[Nages] Fc1966_prdct ;  // Predicted fecundity rate for all ages, 1966
  vector[Nages] Fc2016_prdct ;  // Predicted fecundity rate for all ages, 2016
  real S0_ld ;            // Predicted first year survival, low pop density (2 mil)
  vector[Nages] SA_ld ;  // Predicted survival rates for adult ages, low pop density (2 mil)
  real S0_hd ;            // Predicted first year survival, high pop density (6 mil)
  vector[Nages] SA_hd ;  // Predicted survival rates for adult ages, high pop density (6 mil)
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
  SA_ld = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*2,thta[2]) * zeta * (1 ./ ages)) +  exp(omega))) ;
  S0_ld =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*2, thta[2])) + haz_Ice[21] + exp(omega) )) ;
  SA_hd = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*6,thta[2]) * zeta * (1 ./ ages)) +  exp(omega))) ;
  S0_hd =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*6, thta[2])) + haz_Ice[21] + exp(omega) )) ;
  c = 0 ;
  for (i in 1:NPcts){
    c = c+1 ;
    log_lik[c] = normal_lpdf(Pups[i] | Pups_pred[YrPct[i]], sdNP[i]) ;
    y_new1[i] = normal_rng(Pups_pred[YrPct[i]], sdNP[i]) ;
    P_resid1[i] = (Pups[i] - Pups_pred[YrPct[i]]) / sdNP[i];
    P_resid_new1[i] = (y_new1[i] - Pups_pred[YrPct[i]]) / sdNP[i];
  }
  for(i in 1:NPRobs){
    c = c+1 ;
    log_lik[c] = binomial_lpmf( NPrg[i] | NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    y_new2[i] =  binomial_rng(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    P_resid2[i] = (NPrg[i] - NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]) / sqrt(NFsamp[i]*Fc[YrPRsmp[i]][AgePRsmp[i]]*(1 - Fc[YrPRsmp[i]][AgePRsmp[i]])) ;
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
