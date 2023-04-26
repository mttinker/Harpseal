Do_sims <- function(Pars,Nyrs,Nsec,r){

  # Initialize variables for year 1
  # Otters (Ott array = sections by year by sex):
  Ott = array(0,dim-c(Nsec,Nyrs,2))
  Ott[,1,1] = Pars$Otters$Ot.inits[,1,r]
  Ott[,1,2] = Pars$Otters$Ot.inits[,2,r] 
  # Urchins:
  #
  # Markov:
  #
  # Loop through years
  for(y in 2:Nyrs){
    # Otter dynamics:
    Nvals = Otter_pop_dynamics(Ott[,y-1,1],Ott[,y-1,2],
                               Pars$Otters$Ot.vars$gamma,
                               Pars$Otters$Ot.vars$epsY[,y,r],
                               Pars$Otters$Ot.vars$sigY[[r]],
                               Pars$Otters$Ot.vars$K[[r]],
                               Pars$Otters$Ot.vars$delta[[r]],
                               Pars$Otters$Ot.vars$logemmF[[r]],
                               Pars$Otters$Ot.vars$logemmM[[r]],
                               Pars$Otters$Ot.vars$dspF[[r]],
                               Pars$Otters$Ot.vars$dspM[[r]])
    Ott[,y,1] = Nvals[,1]
    Ott[,y,2] = Nvals[,2]
    #
    # Urchin dynamics:
    #
    # Markov dynamiocs:
    #
    #
  }
  
  Rslts = list(Ott=Ott)
    
  return(Rslts)  
}
