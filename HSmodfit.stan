// HSmodfit
// This Stan code executes a harp seal population model 
//  [Refer to "HSmodel_summary" for details] 
//
// Section 1. Data inputs for model
data {
  int<lower=1> NPcts;             // N counts of pups, total (across 3 areas)
  int<lower=1> NCages;            // N age classes used for age composition comparisons 
  int<lower=1> NCage1;            // first age class considered for obs. age composition
  int<lower=1> NCobs;             // N observations of age composition vectors
  int<lower=1> NPRobs;            // N obs of age-specifc female pregnancy status (binomial)
  int<lower=1> Nyrs;              // N years for model to run (year 1 = 1951, year T = 2020?)
  int<lower=1> Nages;             // N age classes (default is 1:36)
  int<lower=1> Nareas;            // N breeding areas (default 3: Sglf, Nhglf, Front)
  vector[Nages] ages;             // vector of age class values 
  vector[Nages+1] agesC;          // vector of re-centered age class values 
  vector[Nages+1] agesC2;         // vector of re-centered age class values squared
  simplex[Nages] sad0;            // vector of initial stable age distribution values
  int<lower=1> IC[Nyrs,Nareas];   // Annual ice anomaly values by area (discretized, 41 steps of 0.05 units)
  real CI[Nyrs];                // array of annual environmental climate index values
  int<lower=0> Agects[NCobs,NCages];// matrix of adult counts by age class (cols) & year (rows)
  int<lower=0> NFsamp[NPRobs];    // Number females sampled for pregnancy status (grouped by year, age)  
  int<lower=0> NPrg[NPRobs];      // Number pregnant females per sample (grouped by year, age)
  real<lower=0> H_0[Nyrs-1];      // Combined harvest/bycatch values, beaters (age=1,YoY)
  real<lower=0> H_A[Nyrs-1];      // Combined harvest/bycatch values, age>1 
  real<lower=0> Pups[NPcts];      // pup counts, total
  real<lower=0> sdNP[NPcts];      // s.d. associated with total pup counts
  int<lower=1> YrPct[NPcts];      // Year of each pup count (total counts)
  matrix[Nyrs-1,Nareas] PA;       // pup proportions by area & year (estimated for years with no area counts)
  int<lower=1> YrAGsmp[NCobs];    // Year of each female age composition sample 
  int<lower=1> YrPRsmp[NPRobs];   // Year of each female pregnancy status sample 
  int<lower=1> AgePRsmp[NPRobs];  // Age of each female pregnancy status sample
  real N0pri ;                    // mean for vague prior estimate of starting abundance
  real omega;                     // Nuiscance param: base log hazards
  real psipri1 ;                  // Ice anomaly effect on pup survival, fxn param 1 mn
  real<lower=0> psipri2 ;         // Ice anomaly effect on pup survival, fxn param 2 mn
  real CV_HV ;                    // CV associated with harvest counts
  vector[Nyrs-1] Q_0 ;            // probability that a harvested pup is observed (and not SnL) 
  vector[Nyrs-1] Q_A;             // probability that a harvested adult is observed (and not SnL) 
  vector[41] ICvec ;              // dummy vector of ice anomaly values, -1 to 1
  vector[Nyrs-1] Snp ;            // average newborn pup survival by year
  vector[Nages] age_ct_adj ;      // correction factor for age vector (bias against younger adults)
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0> alpha0 ;                  // base adult log haz (for 10 yr old)
  real<lower=0> alpha1 ;                  // adult log haz age modifier par 1
  real<lower=0> alpha2 ;                  // adult log haz age modifier par 2
  simplex[NCages] pie[NCobs];             // probs of age counts (for multinomial)
  real<lower=0,upper=4> thta[2];          // theta params: controls "shape" of DD (for Fecundity and survival)
  real<lower=0> phi[2];                   // D-D effects for fecundity (1) and survival (2)
  real<lower=0,upper=10> beta1;           // Fecundity: Logit of max adult pregancy rate
  real<lower=0,upper=1> beta2;            // Fecundity: age effect (reduced for younger)
  real psi1 ;                             // Ice anomaly effect on pup survival, fxn param 1 
  real<lower=0> psi2;                     // Ice anomaly effect on pup survival, fxn param 2 
  real<lower=.001> N0mil;                 // Initial Abundance, year 1 of time series (i millions)
  real dlta[2];                              // Effect of environmental conditions on fecundity
  real<lower=0> sigF;                     // Environmental stocasticity, var in pregnancy rates
  real<lower=0> sigS;                     // Environmental stocasticity, var in juv survival rates 
  real<lower=0> sigH;                     // Variance in Harvest log hazard rate 
  real<lower=0> tau;                      // precision param for dirichlet-multinomial age dist
  real epsFs[Nyrs-1];                     // Stochastic effects on fecundity, by year (stdzd)
  real epsSs[Nyrs-1];                     // Stochastic effects on juv survival, by year (stdzd)
  real<lower=0> gamma_H0_mn ;             // Mean Harvest log hazard rate, Juvenile
  real<lower=0> gamma_HA_mn;              // Mean Harvest log hazard rate, Adult
  real gamma_H0s[Nyrs-1];                 // Annual harvest log hazard rate, Juvenile (beaters) (stdzd)
  real gamma_HAs[Nyrs-1];                 // Annual harvest log hazard rate, Adults (age 1+) (stdzd)
  real<lower=0> zeta ;                    // adjusts proportional D-D strength for adults relative to juveniles
}
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  vector[Nyrs-1] epsF;                  // Stochastic effects on fecundity, by year 
  vector[Nyrs-1] epsS;                  // Stochastic effects on juv survival, by year 
  vector[Nyrs-1] gamma_H0;              // Annual harvest log hazard rate, Juvenile (beaters) 
  vector[Nyrs-1] gamma_HA;              // Annual harvest log hazard rate, Adults (age 1+) 
  vector[Nages] gamma_A ;               // log hazard ratio for adults
  real gamma_0 ;                        // log hazard ratio for juveniles
  real<lower=0> N[Nyrs] ;               // Population abundance by year
  vector<lower=0>[Nages] n[Nyrs] ;      // Population vector, by year
  vector<lower=0>[Nages] Fc[Nyrs-1] ;   // Fecundity vector, by year
  real<lower=0> S0[Nyrs-1] ;            // Juv Survival from all competing hazards 
  vector<lower=0>[Nages] SA[Nyrs-1];    // Adult Survival from all competing hazards
  real<lower=0> Pups_pred[Nyrs-1] ;     // Predicted pups available for counting, by year  
  vector<lower=0>[NCages] Avc[Nyrs-1] ; // predicted age vector
  vector<lower=0>[Nyrs-1] H0_pred_ALL ; // predicted harvest numbers by year, pups (including SNL)
  vector<lower=0>[Nyrs-1] HA_pred_ALL ; // predicted harvest numbers by year, adults (age1+) (including SNL)
  vector<lower=0>[Nyrs-1] H0_pred ;     // predicted harvest numbers by year without SnL, pups
  vector<lower=0>[Nyrs-1] HA_pred ;     // predicted harvest numbers by year without SnL, adults
  vector[41] haz_Ice[2] ;               // Hazards associated with ice anomalies, Gulf and front
  vector[Nages] gamma_D_scale ;         // Scaling factor for adult density dependent log hazards
  // Re-scale and re-center the random effect variables from standardized normal 
  epsF = to_vector(epsFs) * sigF ;
  epsS = to_vector(epsSs) * sigS ;
  gamma_H0 = to_vector(gamma_H0s) * sigH + gamma_H0_mn;
  gamma_HA = to_vector(gamma_HAs) * sigH + gamma_HA_mn;
  // Initize pop vector  
  N[1] = N0mil * 1000000 ;                     // Initialize population, year 1
  n[1] = sad0 * N[1];             // Initialize population vector, year 1
  // Calculate adult and juvenile log hazards (baseline). 
  gamma_A = alpha0 - alpha1 * agesC[2:(Nages+1)] + alpha2 * agesC2[2:(Nages+1)] ;
  gamma_0 = alpha0 - alpha1 * agesC[1] + alpha2 * agesC2[1]  ;
  // Density depenedent scaling factor (by age) for adults:
  for(a in 1:Nages){
    gamma_D_scale[a] = pow((1 / (a + 0.5)), zeta) ;
  }  
  // Ice hazard functions for Gulf and Front
  haz_Ice[1] = exp(omega + 8 * exp(psi1 - (ICvec * psi2)) ./ (1 + exp(psi1 - (ICvec * psi2)))) ;
  haz_Ice[2] = exp(omega + 8 * exp((psi1-1) - (ICvec * psi2)) ./ (1 + exp((psi1-1) - (ICvec * psi2)))) ;
  // Loop through years to calculate demographic transitions and population dynamics 
  for (i in 1:(Nyrs-1)){
    // Declare some temporary variables for this year:
    real Nml ;                    // Current population abundance in millions
    real haz_J ;                  // Juv baseline hazards
    vector[Nages] haz_A ;         // Adult baseline hazards    
    real haz_IC ;                 // Weighted mean ice hazards (juvenile)
    real haz_H0 ;                 // Harvest hazards, juvenile  
    real haz_HA ;                 // Harvest hazards, adult     
    real prp_HV_0 ;               // Proportion of mortality comprised of harvest, Juve 
    vector[Nages] prp_HV_A ;      // Proportion of mortality comprised of harvest, Adult
    real gamma_D ;                // density-dependent effect for juveniles/subadults
    vector[Nages] Age_ppn ;       // proportion of individuals in each age class, bias-adjusted 
    // current N in millions of animals, for D-D calcs
    Nml = N[i]/1000000 ;
    // Annual Fecundity, with D-D, envir. effects and stochasticity 
    Fc[i][1:3] = rep_vector(0,3) ;
    Fc[i][4:8] = inv_logit(beta1 - beta2 * square(8-ages[4:8]) - pow(phi[1]*Nml,thta[1]) - dlta[1]*CI[i] + epsF[i]) ;
    Fc[i][9:Nages] = rep_vector( Fc[i][8], (Nages - 8)) ; // ** IF >8 age classes
    // Juvenile competing hazards and net survival
    gamma_D = pow(phi[2]*Nml, thta[2]) ;
    haz_J = exp(omega + gamma_0 + gamma_D + dlta[2] * CI[i+1] +  epsS[i] ) ; //  juv hazards
    // weighted avg ice anom haz (wt = prop. pups in each area)
    haz_IC = PA[i][1] * haz_Ice[1][IC[i,1]] + PA[i][2] * haz_Ice[1][IC[i,2]] + PA[i][3] * haz_Ice[2][IC[i,3]] ;
    haz_H0 = exp(omega + gamma_H0[i]) ; // hazards from bycatch/human harvest of pups/beaters
    S0[i] = exp(-1 * (haz_J + haz_IC + haz_H0)) ;  // combine hazards for Juv survival
    // Adult competing hazards and net survival (allow weak DD effects for sub-adults)
    haz_A = exp(omega + gamma_A + gamma_D_scale * gamma_D ) ; // reg adult hazards 
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
    // Calculate age comp vector, adjusting for bias against observing younger animals 
    Age_ppn = ( (n[i] .* age_ct_adj) / sum(n[i]  .* age_ct_adj)) ;
    Avc[i] = Age_ppn[NCage1:Nages] ;    
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
    pie[i] ~ dirichlet(10 * tau * Avc[YrAGsmp[i]]);
    Agects[i,] ~ multinomial(pie[i]) ;
  }
  //  Pregancy status of females by age/year (beta-binomial dist)
  for(i in 1:NPRobs){
    NPrg[i] ~ binomial(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
  // B) Prior distributions for model parameters:
  // Hierarchical random effects: NOTE: maybe make annual stochastic effects autoregressive (AR[1])?  
  epsFs ~ normal(0,1) ;               // Annual deviations from mean fecundity (stdzd)
  epsSs ~ normal(0,1) ;               // Annual deviations from mean juv hazards (stdzd)
  gamma_H0s ~ normal(0,1);            // Annual pup harvest log hazards (stdzd)
  gamma_HAs ~ normal(0,1);            // Annual adult harvest log hazards (stdzd)
  // Base parameter priors:
  N0mil ~ normal(2.5,.5) ; 
  thta ~ gamma(5,4) ;
  phi ~ normal(0,1) ;
  beta1 ~ normal(0,5) ;
  beta2 ~ normal(0,1) ;
  psi1 ~ normal(-1,1) ;
  psi2 ~ normal(5,1) ;
  dlta ~ normal(0,1) ;
  sigF ~ normal(0,1) ;
  sigH ~ normal(0,1) ;
  sigS ~ normal(0,1) ;
  tau ~ cauchy(0,5) ;
  gamma_H0_mn ~ normal(6,1) ;
  gamma_HA_mn ~ normal(4,1) ;
  alpha0 ~ normal(3,1) ;
  alpha1 ~ normal(0,.5) ;
  alpha2 ~ normal(0,0.1) ;
  zeta ~ normal(4,2) ;
}
// Section 5. Derived parameters and statistics 
generated quantities {
  real N0 ;                           // initial abundance (year 1)
  real PAmean[Nareas] ;               // mean PA values
  real Fc8_prdct[Nyrs-1] ;            // Predicted fecundity rate for 8+ over years
  vector[Nages] Fc1966_prdct ;        // Predicted fecundity rate for all ages, 1966
  vector[Nages] Fc2016_prdct ;        // Predicted fecundity rate for all ages, 2016
  real S0_ld ;                        // Predicted first year survival, low pop density (2 mil)
  vector[Nages] SA_ld ;               // Predicted survival rates for adult ages, low pop density (2 mil)
  real S0_hd ;                        // Predicted first year survival, high pop density (6 mil)
  vector[Nages] SA_hd ;               // Predicted survival rates for adult ages, high pop density (6 mil)
  real log_lik[NPcts+NPRobs] ;        // Log liklihood of obs. data (for LooIC)
  real y_new1[NPcts]  ;               // New observations, pup counts
  real y_new2[NPRobs]  ;              // New observations, preg rates 
  int<lower=0> c;
  N0 = N[2] * (N[2]/N[3]) ;  
  for (i in 1:(Nyrs-1)){
    Fc8_prdct[i] = Fc[i][8] ;
  }
  Fc1966_prdct = Fc[16][1:Nages];
  Fc2016_prdct = Fc[66][1:Nages];
  S0_ld =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*1.67, thta[2])) + haz_Ice[1][21] + exp(omega) )) ;
  S0_hd =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*5, thta[2])) + haz_Ice[1][21] + exp(omega) )) ;
  SA_ld = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*1.67,thta[2]) * gamma_D_scale) +  exp(omega))) ;
  SA_hd = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*5,thta[2]) * gamma_D_scale) +  exp(omega))) ;
  c = 0 ;
  for (i in 1:NPcts){
    c = c+1 ;
    log_lik[c] = normal_lpdf(Pups[i] | Pups_pred[YrPct[i]], sdNP[i]) ;
    y_new1[i] = normal_rng(Pups_pred[YrPct[i]], sdNP[i]) ;
  }
  for(i in 1:NPRobs){
    c = c+1 ;
    log_lik[c] = binomial_lpmf( NPrg[i] | NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    y_new2[i] =  binomial_rng(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
}
