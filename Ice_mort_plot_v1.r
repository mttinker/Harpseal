Ice_mort_plot <- function(psipri1a,psipri1b,psipri2a,psipri2b) {
  # Function to Plot effects of ice cover on pup mortality
  require(ggplot2)
  gamm0 = -7
  Ice_cover = seq(-1,1,by=0.02)
  gammIce = 8 * inv.logit(psipri1a - psipri2a*Ice_cover) # psipri1a * (1/exp(Ice_cover))^psipri2a - 1
  # S_ice = inv.logit(-psipri1a + 7*(Ice_cover+1)^psipri2a)
  S_ice = exp(-1 * exp(gamm0 + gammIce))
  S_ice_r = matrix(0,nrow = 1000,ncol = length(Ice_cover))
  for(r in 1:1000){
    psi1r = rnorm(1,psipri1a,psipri1b); psi2r = rnorm(1,psipri2a,psipri2b)  
    gammIce = 8 * inv.logit(psi1r - psi2r * Ice_cover) #  psi1r * (1/exp(Ice_cover))^psi2r - 1 ;
    S_ice_r[r,1:length(Ice_cover)] = exp(-1 * exp(gamm0 + gammIce)) ;
    # S_ice_r[r,1:length(Ice_cover)] = inv.logit(-psi1r +7*(Ice_cover+1)^psi2r)
  }
  dfplt = data.frame(Ice_cover=Ice_cover,
                     Mean_S_ice = S_ice,
                     Lo = apply(S_ice_r,2,quantile,0.025),
                     Hi = apply(S_ice_r,2,quantile,0.975))
  plt = ggplot(dfplt,aes(x=Ice_cover,y=Mean_S_ice)) +
    geom_ribbon(aes(ymin=Lo,ymax=Hi),alpha=0.3) +
    geom_line() + geom_vline(xintercept = 0) +
    labs(x="Ice Index (% Cover Anomaly)",y="Relative pup survival") +
    ggtitle("Effect of ice cover on pup survival") +
    theme_classic()
  # print(plt)
  return(plt)
}

