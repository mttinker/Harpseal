// HSmodfit
// This Stan code executes a harp seal population model 
//  [Refer to "HSmodel_summary" for details] 
//
functions {
  // user function to create projection matrix
  matrix makemat(int Nages, vector SA, vector F, real S0, real Sn) {
     matrix[Nages,Nages] M; 
     vector[Nages] R ; 
     R = F * (0.5 * S0 * Sn) ; 
     M = rep_matrix(0,Nages,Nages) ;
     M[1] = to_row_vector(R) ;
     M[2:Nages,1:(Nages-1)] = add_diag(M[2:Nages,1:(Nages-1)], SA[1:(Nages-1)]) ;
     M[Nages,Nages] = SA[Nages] ;
     return M;
  }
}
// Section 1. Data inputs for model
data {
  int<lower=1> NPcts;             // N counts of pups, total (sum 3 areas: Sglf, Nhglf, Front)
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
  matrix[(Nyrs-1),Nareas] IC;     // Annual ice anomaly values by area 
  array[Nyrs] real CI;            // array of annual environmental climate index values
  array[NCobs,NCages] int<lower=0> Agects;// matrix of adult counts by age class (cols) & year (rows)
  array[NPRobs] int<lower=0> NFsamp; // Number females sampled for pregnancy status (grouped by year, age)  
  array[NPRobs] int<lower=0> NPrg;// Number pregnant females per sample (grouped by year, age)
  array[Nyrs-1] real<lower=0> H_0;// Combined harvest/bycatch values, beaters (age=1,YoY)
  array[Nyrs-1] real<lower=0> H_A;// Combined harvest/bycatch values, age>1 
  array[NPcts] real<lower=0> Pups;// pup counts, total
  array[NPcts] real<lower=0> sdNP;// s.d. associated with total pup counts
  array[NPcts] int<lower=1> YrPct;// Year of each pup count (total counts)
  matrix[Nyrs-1,Nareas] PA;       // pup proportions by area & year (estimated for years with no area counts)
  array[NCobs] int<lower=1> YrAGsmp;  // Year of each female age composition sample 
  array[NPRobs] int<lower=1> YrPRsmp; // Year of each female pregnancy status sample 
  array[NPRobs] int<lower=1> AgePRsmp;// Age of each female pregnancy status sample
  real N0pri ;                    // mean for vague prior estimate of starting abundance
  real omega;                     // Nuiscance param: base log hazards
  real CV_HV ;                    // CV associated with harvest counts
  vector[Nyrs-1] Q_0 ;            // probability that a harvested pup is observed (and not SnL) 
  vector[Nyrs-1] Q_A;             // probability that a harvested adult is observed (and not SnL) 
  vector[41] ICvec ;              // dummy vector of ice anomaly values, -1 to 1
  vector[Nyrs-1] Snp ;            // average newborn pup survival by year
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0,upper=25> alpha0 ;         // base adult log haz (for 10 yr old)
  real<lower=0,upper=5> alpha1 ;          // adult log haz age modifier par 1
  real<lower=0,upper=1> alpha2 ;          // adult log haz age modifier par 2
  array[NCobs] simplex[NCages] pie;       // probs of age counts (for multinomial)
  array[2] real<lower=0.001,upper=5> thta;// theta params: controls "shape" of DD (for Fecundity and survival)
  array[2] real<lower=0,upper=10> phi;    // D-D effects for fecundity (1) and survival (2)
  real<lower=0,upper=25> beta1;           // Fecundity: Logit of max adult pregancy rate
  real<lower=0,upper=5> beta2;            // Fecundity: age effect (reduced for younger)
  array[2] real<lower=-5,upper=5> psi1 ;  // Ice anomaly effect on pup survival, fxn param 1 
  real<lower=0,upper=10> psi2;            // Ice anomaly effect on pup survival, fxn param 2 
  real<lower=.001,upper=10> N0mil;        // Initial Abundance, year 1 of time series (i millions)
  array[2] real<lower=-10,upper=10> dlta; // Effect of environmental conditions on fecundity
  real<lower=0,upper=10> sigF;            // Environmental stocasticity, var in pregnancy rates
  real<lower=0,upper=10> sigS;            // Environmental stocasticity, var in juv survival rates 
  real<lower=0,upper=10> sigH;            // Variance in Harvest log hazard rate 
  real<lower=0,upper=50> tau10;           // precision param for dirichlet-multinomial age dist
  vector[Nyrs-1] epsFs;                   // Stochastic effects on fecundity, by year (stdzd)
  vector[Nyrs-1] epsSs;                   // Stochastic effects on juv survival, by year (stdzd)
  real<lower=0> gamma_H0_mn ;             // Mean Harvest log hazard rate, Juvenile
  real<lower=0> gamma_HA_mn ;             // Mean Harvest log hazard rate, Adult
  vector[Nyrs-1] gamma_H0s ;              // Annual harvest log hazard rate, Juvenile (beaters) (stdzd)
  vector[Nyrs-1] gamma_HAs ;              // Annual harvest log hazard rate, Adults (age 1+) (stdzd)
  real<lower=0,upper=25> zeta ;           // adjusts proportional D-D strength for adults relative to juveniles
}
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  vector[Nyrs-1] epsF;                    // Stochastic effects on fecundity, by year 
  vector[Nyrs-1] epsS;                    // Stochastic effects on juv survival, by year 
  vector[Nyrs-1] gamma_H0;                // Annual harvest log hazard rate, Juvenile (beaters) 
  vector[Nyrs-1] gamma_HA;                // Annual harvest log hazard rate, Adults (age 1+) 
  vector[Nages] gamma_A ;                 // log hazard ratio for adults
  real gamma_0 ;                          // log hazard ratio for juveniles
  array[Nyrs] real<lower=0> N ;           // Population abundance by year
  array[Nyrs] vector<lower=0>[Nages] n ;  // Population vector, by year
  array[Nyrs-1] vector<lower=0>[Nages] Fc;// Fecundity vector, by year
  array[Nyrs-1] real<lower=0> S0 ;        // Juv Survival from all competing hazards 
  array[Nyrs-1] vector<lower=0>[Nages] SA;// Adult Survival from all competing hazards
  array[Nyrs-1] real<lower=0> Pups_pred ; // Predicted pups available for counting, by year  
  array[Nyrs-1] vector<lower=0>[NCages] Avc;// predicted age vector
  vector<lower=0>[Nyrs-1] H0_pred_ALL ;   // predicted harvest numbers by year, pups (including SNL)
  vector<lower=0>[Nyrs-1] HA_pred_ALL ;   // predicted harvest numbers by year, adults (age1+) (including SNL)
  vector<lower=0>[Nyrs-1] H0_pred ;       // predicted harvest numbers by year without SnL, pups
  vector<lower=0>[Nyrs-1] HA_pred ;       // predicted harvest numbers by year without SnL, adults
  vector[Nages] gamma_D_scale ;           // Scaling factor for adult density dependent log hazards
  real tau ;
  tau = tau10 * 10 ;
  // Re-scale and re-center the random effect variables from standardized normal 
  epsF = epsFs * sigF ;
  epsS = epsSs * sigS ;
  gamma_H0 = gamma_H0s * sigH + gamma_H0_mn;
  gamma_HA = gamma_HAs * sigH + gamma_HA_mn;
  // Initize pop vector  
  N[1] = N0mil * 1000000 ;                     // Initialize population, year 1
  n[1] = sad0 * N[1];             // Initialize population vector, year 1
  // Calculate adult and juvenile log hazards (baseline). 
  gamma_A = alpha0 + alpha1 * agesC[2:(Nages+1)] + alpha2 * agesC2[2:(Nages+1)] ;
  gamma_0 = alpha0 + alpha1 * agesC[1] + alpha2 * agesC2[1]  ;
  // Density depenedent scaling factor (by age) for adults:
  gamma_D_scale = exp( zeta * log(1 ./ (ages + 0.5))) ;
  // Ice hazard functions for Gulf and Front
  // haz_Ice[1] = exp(omega + 8 * exp(psi1[1] - (ICvec * psi2)) ./ (1 + exp(psi1[1] - (ICvec * psi2)))) ;
  // haz_Ice[2] = exp(omega + 8 * exp((psi1[2]) - (ICvec * psi2)) ./ (1 + exp((psi1[2]) - (ICvec * psi2)))) ;
  // Loop through years to calculate demographic transitions and population dynamics 
  for (i in 1:(Nyrs-1)){
    // Declare some temporary variables for this year:
    real haz_J ;                    // Juv baseline hazards
    vector[Nages] haz_A ;           // Adult baseline hazards    
    real haz_IC ;                   // Area-weighted mean ice hazards (juvenile)
    real haz_H0 ;                   // Harvest hazards, juvenile  
    real haz_HA ;                   // Harvest hazards, adult     
    real prp_HV_0 ;                 // Proportion of mortality comprised of harvest, Juve 
    vector[Nages] prp_HV_A ;        // Proportion of mortality comprised of harvest, Adult
    real gamma_D ;                  // density-dependent effect for juveniles/subadults
    real ddF ;                      // density-dependent effect for fecundity
    matrix[Nages,Nages] Mt;         // matrix, year t
    ddF = exp(thta[1] * log(phi[1]*(N[i]/1000000))) ;
    gamma_D = exp(thta[2] * log(phi[2]*(N[i]/1000000))) ;
    // Annual Fecundity, with D-D, envir. effects and stochasticity 
    Fc[i][1:3] = rep_vector(0,3) ;
    Fc[i][4:8] = inv_logit(beta1 - beta2 * square(8-ages[4:8]) - ddF - dlta[1]*CI[i] + epsF[i]) ;
    Fc[i][9:Nages] = rep_vector( Fc[i][8], (Nages - 8)) ; // ** IF >8 age classes
    // Juvenile competing hazards and net survival
    haz_J = exp(omega + gamma_0 + gamma_D + dlta[2] * CI[i+1] +  epsS[i] ) ; //  juv hazards
    // weighted avg ice anom haz (wt = prop. pups in each area)
    haz_IC = PA[i][1] * exp(omega + 8 * inv_logit(psi1[1] - psi2*IC[i,1])) +
             PA[i][2] * exp(omega + 8 * inv_logit(psi1[1] - psi2*IC[i,2])) +
             PA[i][3] * exp(omega + 8 * inv_logit(psi1[2] - psi2*IC[i,3])) ;
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
    // Calculate age comp vector for animals older than minimum cut-off age
    Avc[i] = n[i][NCage1:Nages] / sum(n[i][NCage1:Nages] ) ;
    // Update population vector and total abundance for next year
    Mt = makemat(Nages,SA[i],Fc[i],S0[i],Snp[i]) ;
    n[i+1] = Mt * n[i] ;
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
    pie[i] ~ dirichlet(tau * Avc[YrAGsmp[i]]);
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
  N0mil ~ gamma(1,.5) ; 
  thta ~ gamma(5,4) ;
  phi ~ cauchy(0,1) ;
  beta1 ~ cauchy(0,2.5) ;
  beta2 ~ cauchy(0,0.1) ;
  psi1 ~ normal(0,1) ;
  psi2 ~ normal(4,2) ;
  dlta ~ cauchy(0,1) ;
  sigF ~ cauchy(0,.5) ;
  sigH ~ cauchy(0,.5) ;
  sigS ~ cauchy(0,.5) ;
  tau10 ~ cauchy(0,5) ;
  gamma_H0_mn ~ normal(6,1) ;
  gamma_HA_mn ~ normal(4,1) ;
  alpha0 ~ cauchy(0,1) ;
  alpha1 ~ cauchy(0,0.1) ;
  alpha2 ~ cauchy(0,0.1) ;
  zeta ~ cauchy(0,2.5) ;
}
// Section 5. Derived parameters and statistics 
generated quantities {
  real N0 ;                           // initial abundance (year 1)
  array[Nareas] real PAmean ;         // mean PA values
  array[Nyrs-1] real Fc8_prdct ;      // Predicted fecundity rate for 8+ across years
  vector[Nages] Fc1966_prdct ;        // Predicted fecundity rate for all ages, 1966
  vector[Nages] Fc2016_prdct ;        // Predicted fecundity rate for all ages, 2016
  array[Nyrs-1] real Pyng_prdct ;     // Predicted ratio of 8-14 to 8-36 yr olds
  array[Nyrs-1] real Pold_prdct ;     // Predicted ratio of 8-36 yr olds vs all animals
  real S0_ld ;                        // Predicted first year survival, low pop density (2 mil)
  vector[Nages] SA_ld ;               // Predicted survival rates for adult ages, low pop density (2 mil)
  real S0_hd ;                        // Predicted first year survival, high pop density (6 mil)
  vector[Nages] SA_hd ;               // Predicted survival rates for adult ages, high pop density (6 mil)
  array[NPcts] real y_new1  ;         // New observations, pup counts
  array[NPRobs] real y_new2  ;        // New observations, preg rates 
  array[2] vector[41] haz_Ice ;       // Hazards associated with ice anomalies, Gulf and front
  int<lower=0> c;
  //
  haz_Ice[1] = exp(omega + 8 * exp(psi1[1] - (ICvec * psi2)) ./ (1 + exp(psi1[1] - (ICvec * psi2)))) ;
  haz_Ice[2] = exp(omega + 8 * exp((psi1[2]) - (ICvec * psi2)) ./ (1 + exp((psi1[2]) - (ICvec * psi2)))) ;
  N0 = N[2] * (N[2]/N[3]) ;  
  for (i in 1:(Nyrs-1)){
    Fc8_prdct[i] = Fc[i][8] ;
    Pyng_prdct[i] = sum(n[i][8:14]) / sum(n[i][8:Nages] ) ;
    Pold_prdct[i] = sum(n[i][8:36]) / (sum(n[i][1:Nages] ) + Pups_pred[i]);
  }
  Fc1966_prdct = Fc[16][1:Nages];
  Fc2016_prdct = Fc[66][1:Nages];
  S0_ld =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*2*.82, thta[2])) + haz_Ice[1][21] + exp(omega) )) ;
  S0_hd =  exp(-1 * (exp(omega + gamma_0 + pow(phi[2]*6*.82, thta[2])) + haz_Ice[1][21] + exp(omega) )) ;
  SA_ld = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*2*.82,thta[2]) * gamma_D_scale) +  exp(omega))) ;
  SA_hd = exp(-1 * (exp(omega + gamma_A + pow(phi[2]*6*.82,thta[2]) * gamma_D_scale) +  exp(omega))) ;
  c = 0 ;
  for (i in 1:NPcts){
    c = c+1 ;
    // log_lik[c] = normal_lpdf(Pups[i] | Pups_pred[YrPct[i]], sdNP[i]) ;
    y_new1[i] = normal_rng(Pups_pred[YrPct[i]], sdNP[i]) ;
  }
  for(i in 1:NPRobs){
    c = c+1 ;
    // log_lik[c] = binomial_lpmf( NPrg[i] | NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
    y_new2[i] =  binomial_rng(NFsamp[i], Fc[YrPRsmp[i]][AgePRsmp[i]]) ;
  }
}
