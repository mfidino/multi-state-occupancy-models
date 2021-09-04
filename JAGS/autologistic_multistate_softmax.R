model{
  #
  #THE LATENT STATE MODEL
  #
  for(site in 1:nsite){
    # LATENT STATE LINEAR PREDICTORS.
    #   Note: I am assuming here that while the model is dynamic,
    #   The probabilities do not vary by year (because the
    #   linear predictors do not vary through time). If you want
    #   to add this then all of the linear predictors need
    #   to be indexed by site AND year (e.g., copy all the psi stuff
    #   and add it into the year loop as well, indexing psi by year
    #   (e.g., psi[site,year,1:3,1:3]). Directly below you'd add 
    #   a fourth dimension to psi as well, but index it as 1 for
    #   first year (i.e., psi[site,1,1:3,1:3]).
    # 
    # latent state probabilities given state == 1
    psi[site,1,1] <- 1
    psi[site,1,2] <- exp( inprod(beta2, x[site,] ) )
    psi[site,1,3] <- exp( inprod(beta3, u[site,] ) )
    # latent state probabilities given state == 2
    psi[site,2,1] <- 1
    psi[site,2,2] <- exp( inprod(beta2, x[site,] ) + theta2g2 ) 
    psi[site,2,3] <- exp( inprod(beta3, u[site,] ) + theta3g2 ) 
    # latent state probabilities given state == 3
    psi[site,3,1] <- 1
    psi[site,3,2] <- exp( inprod(beta2, x[site,] ) + theta2g3 )
    psi[site,3,3] <- exp( inprod(beta3, u[site,] ) + theta3g3 ) 
    # Estimate latent state year == 1. Setting to first
    #  row of psi because we have no prior knowledge
    #  on state of site before we started sampling.
    z[site,1] ~ dcat(
      psi[site,1,] / sum(psi[site,1,])
    )
    # Do the remaining years. We grab the correct row
    #  of the transition matrix based on the state in the 
    #  previous time step.
    for(year in 2:nyear){
      z[site,year] ~ dcat(
        psi[site,z[year-1],] / sum(psi[site,z[year-1],])
      )
    }
  }
  #
  # THE DATA MODEL
  #
  for(site in 1:nsite){
    for(survey in 1:nsurvey){
      # Set up the logit linear predictors
      #   Note: I am assuming here that while the model is dynamic,
      #   The probabilities do not vary by year (because the
      #   linear predictors do not vary through time). If you want
      #   to add this then all of the linear predictors need
      #   to be indexed by site AND year.
      # Fill in detection probability matrix
      # First row: TS = 1
      eta[site,survey,1,1] <- 1 # -------------------------------------- OS = 1
      eta[site,survey,1,2] <- 0 # -------------------------------------- OS = 2
      eta[site,survey,1,3] <- 0 # -------------------------------------- OS = 3
      # Second row: TS = 2
      eta[site,survey,2,1] <- 1 # -------------------------------------- OS = 1
      eta[site,survey,2,2] <- exp( inprod( rho2g2, k[site,survey,] ) ) # OS = 2
      eta[site,survey,2,3] <- 0 # -------------------------------------- OS = 3
      # Third row: TS = 3
      eta[site,survey,3,1] <- 1 # -------------------------------------- OS = 1
      eta[site,survey,3,2] <- exp( inprod( rho2g3, k[site,survey,] ) ) # OS = 2
      eta[site,survey,3,3] <- exp( inprod( rho3g3, k[site,survey,] ) ) # OS = 3
      for(yr in 1:nyear){
        # Index the appropriate row of eta based on the current latent state.
        # Again, we are assuming there is no variation among years or sampling.
        y[site,survey,yr] ~ dcat(
          eta[site,survey,z[yr],] / sum(eta[site,survey,z[yr],])
        )
      }
    }
  }
  #
  # Priors
  #
  # Pr(latent state 2)
  for(b2 in 1:nbeta2){
    beta2[b2] ~ dlogis(0,1)
  }
  # Pr(latent state 3)
  for(b3 in 1:nbeta2){
    beta3[b3] ~ dlogis(0,1)
  }
  # Autologistic terms
  theta2g2 ~ dlogis(0,1)
  theta3g2 ~ dlogis(0,1)
  theta2g3 ~ dlogis(0,1)
  theta3g3 ~ dlogis(0,1)
  # Pr(OS = 2 | TS = 2 or 3)
  for(ts2 in 1:nts2){
    rho2g3[ts2] ~ dlogis(0,1)
  }
  # Pr(OS = 3 | TS = 3) 
  for(ts3 in 1:nts3){
    rho3g3[ts3] ~ dlogis(0,1)
  } 
}