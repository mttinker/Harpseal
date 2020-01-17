HS_matrixD <- function(alpha,beta,delta,phi,psi,theta,PA,
                       Nt,CE,IC,HV) {
  # Function to create a harp seal population matrix 
  # Pre-breeding survey, so first age class in matric = ~1 year-old juveniles
  # (but # pups calculated and tracked separately)
  # NOTE: stochastic version would include non-0 process error terms (eps) for Fi and S0
  # HV = age-structured vector of harvest/catch estimates
  # Ct = ice anomalies for each of 3 areas
  # PA = proportion of pups in each of 3 areas
  # Et = envitonmental consitions in year t (goes with param delta)
  ages = seq(1,20); ages2 = ages^2
  Abund = Nt/1000000
  epsF = 0
  epsJ = 0
  PupDetect = .99 # mean probability a pup is detected (~1)
  a0 = -5         # base log hazard (nuiscence constant for haz fxn, fixed)
  aJ = alpha[1]   # Juvenile hazards: weak prior  
  aA = alpha[2]   # Adult hazards intcpt: strong prior
  aA1 = alpha[3]   # Adult hazards 1st order effect: strong prior
  aA2 = alpha[4]   # Adult hazards 2nd order effect: strong prior
  b0 = beta[1]    # Max pregancy rate: strong prior
  b1 = beta[2]    # Effect of age on preg rate: weak prior (assuming good data on age-specific pregs)
  phiF = phi[1]   # fecundity D-D param: weak prior
  phiJ = phi[2]   # juv surviv D-D param: weak prior
  psi1 = psi[1]   # Ice cover effects on pup survival, p1 (moderate prior)
  psi2 = psi[2]   # Ice cover effects on pup survival, p2 (moderate prior)
  # Calculate demographic rates * ADD ABORTIONS?
  Fi = numeric(length = 20) # Fecundity (~late term pregancy rates) (row vector in Stan)
  Fi[1:2] = 0
  Fi[3:8] = inv.logit(b0 - b1*(8-ages[3:8])^2 - (phiF*Abund)^theta + delta*CE + epsF)
  Fi[9:20] = Fi[8]
  SiceA = inv.logit(-psi1 + 7*(IC + 1)^psi2) # Pup survival from ice anomolies for 3 areas 
  gamm0 = a0 + aJ + (phiJ*Abund)^theta      # Juvenile log hazards, with D-D
  S0 = exp(-exp(gamm0 + epsJ))*sum(PA*SiceA)       # Juv survival, with ice mort
  R = Fi*S0*0.5 # Recruitment to first year class (1yr old juveniles), assume 50% sex ratio
  gammA = a0 + aA - aA1*ages[1:19] + aA2*ages2[1:19]   # Adult log hazards
  S = exp(-exp(gammA))  # Adult survival (to age 19): allow to vary by year?
  S[20] =  0.75*S[19] # Note: assume reduced survival for age 20+ 
  # Construct projection matrix  
  M = diag(S[1:19])
  M = cbind(M,rep(0,19)) 
  M[19,20] = S[20]  
  M = rbind(R,M)
  # NOTE: assuming we can use age vector of sampled females as observed multinomial variable
  # so that agevec_Yt ~ multinomial( (nt/sum(nt)) )
  # we could allow adult "S" to be quadratic fxn of age, instead of constant 
  #  (since the observed age distribution constrains age-specific survival rates)
  # Thus, for Stan version:
  #    gammA = a0 + aA - (a1*ages) + (a2*ages2) 
  #    S = exp(-exp(gammA)) 
  # where ages is a data-declared vector 1:20, ages2 = ages^2
  # and fitted params aA and a1 and a2 are all strong priors, truncated >0  
  # for example, try: a0 = -5; aA = 2.5; a1 = .2; a2 = .012
  # Also, to construct annual projection matrix M, construct all but 1st row just once:
  #    M1 = diag_matrix(S[1:19])
  #    M2 = append_col(M1, rep_vector(0,19))  
  #    M2[19,20] = 0.75*S[19]
  # and then for each year, where R is a row vector calculated as above:
  #    M = append_row(R, M2)
  lam=Re(eigen(M)$values[1])          # Calculate lambda for matrix
  W=eigen(M)$vectors          # W=matrix of right eigenvectors 
  w=abs(W[,1])					      # w=stable age distribution, unscaled
  sad = w/sum(w)              # w=stable age distribution (sad), scaled
  if(sum(HV)>0){
    sad0 = sad
    # Iteratively re-calculate sad, acounting for harvrest mortality 
    for (y in 1:20){
      nt0 = Nt*sad0
      nt0 = M%*%nt0 - HV
      sad0 = nt0/sum(nt0)
    }
    sad = sad0
    nt = Nt*sad
    nt1 = M%*%nt - HV
    Nt1 = sum(nt1)
    lam = Nt1/Nt
  }else{
    nt = Nt*sad
    nt1 = M%*%nt 
    Nt1 = sum(nt1)
  }
  # Calculate expected pup counts, by breeding Area
  pups = numeric(length = 3)
  pups[1] = sum((0.5*nt*Fi)*PA[1]*sqrt(SiceA[1]))*PupDetect
  pups[2] = sum((0.5*nt*Fi)*PA[2]*sqrt(SiceA[2]))*PupDetect
  pups[3] = sum((0.5*nt*Fi)*PA[3]*sqrt(SiceA[3]))*PupDetect
  result <- list(lambda=lam,sad=sad,mat=M,pups=pups,nt1=nt1,Nt1=Nt1,Fecund=Fi)
  return(result)
}