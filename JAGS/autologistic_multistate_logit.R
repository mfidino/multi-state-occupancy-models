model{
  #
  # THE LATENT STATE MODEL
  #
  for(site in 1:nsite){
    # Set up the logit linear predictors
    #   Note: I am assuming here that while the model is dynamic,
    #   the probabilities do not vary by year (because the logit
    #   linear predictors do not vary through time).
    # LATENT STATE LINEAR PREDICTORS.
    # Probability of either state 2 or 3
    #  given state at year-1 == 1 (or first year).
    logit(omega[site]) <- inprod(
      beta23, x[site,]
    )
    # Conditional probability of state 3
    logit(delta[site]) <- inprod(
      beta3, u[site,]
    )
    # Probability of either state 2 or 3
    #  given state at year-1 != 1
    logit(omega_cond[site]) <- inprod(
      beta23, x[site,]
    ) + theta23
    # Conditional probability of state 3 
    #  given state at year-1 = 3.
    logit(delta_cond[site]) <- inprod(
      beta3, u[site,]
    ) + theta3
    # Fill in the transition matrix
    # latent state probabilities given state == 1
    psi[site,1,1] <- 1 - omega[site]
    psi[site,1,2] <- omega[site] * (1-delta[site])
    psi[site,1,3] <- omega[site] * delta[site]
    # latent state probabilities given state == 2
    psi[site,2,1] <- 1 - omega_cond[site]
    psi[site,2,2] <- omega_cond[site] * (1-delta[site])
    psi[site,2,3] <- omega_cond[site] * delta[site]
    # latent state probabilities given state == 3
    psi[site,3,1] <- 1 - omega_cond[site]
    psi[site,3,2] <- omega_cond[site] * (1-delta_cond[site])
    psi[site,3,3] <- omega_cond[site] * delta_cond[site]
    # Estimate latent state year == 1. Setting to first
    #  row of psi because we have no prior knowledge
    #  on state of site before we started sampling.
    z[site,1] ~ dcat(
      psi[site,1,]
    )
    # Do the remaining years. We grab the correct row
    #  of the transition matrix based on the state in the 
    #  previous time step.
    for(year in 2:nyear){
      z[site,year] ~ dcat(
        psi[site,z[site,year-1],]
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
      #   The probabilities do not vary by year (because the logit
      #   linear predictors do not vary through time).
      # Probability of detecting either state 2 or 3
      logit(eta23[site,survey]) <- inprod(
        rho23, k[site,survey,]
      )
      # Probability of detecting state 3 given 2
      logit(eta3[site,survey]) <- inprod(
        rho3, k[site,survey,]
      )
      # Fill in detection probability matrix
      # First row: TS = 1
      eta[site,survey,1,1] <- 1 # -------------------------------------- OS = 1
      eta[site,survey,1,2] <- 0 # -------------------------------------- OS = 2
      eta[site,survey,1,3] <- 0 # -------------------------------------- OS = 3
      # Second row: TS = 2
      eta[site,survey,2,1] <- 1 - eta23[site,survey] # ----------------- OS = 1
      eta[site,survey,2,2] <- eta23[site,survey] # --------------------- OS = 2
      eta[site,survey,2,3] <- 0 # -------------------------------------- OS = 3
      # Third row: TS = 3
      eta[site,survey,3,1] <- 1 - eta23[site,survey] # ----------------- OS = 1
      eta[site,survey,3,2] <- eta23[site,survey]*(1-eta3[site,survey]) # OS = 2
      eta[site,survey,3,3] <- eta23[site,survey]*eta3[site,survey] # --- OS = 3
      for(yr in 1:nyear){
        # Index the appropriate row of eta based on the current latent state.
        # Again, we are assuming there is no variation among years or sampling.
        y[site,survey,yr] ~ dcat(
          eta[site,survey,z[site,yr],]
        )
      }
    }
  }
  #
  # Priors
  #
  # Pr(latent state 2 or 3)
  for(b2 in 1:nbeta23){
    beta23[b2] ~ dlogis(0,1)
  }
  # Pr(latent state 3)
  for(b3 in 1:nbeta3){
    beta3[b3] ~ dlogis(0,1)
  }
  # Autologistic terms
  theta23 ~ dlogis(0,1)
  theta3 ~ dlogis(0,1)
  # Pr(OS = 2 or 3 | TS = 2 or 3)
  for(ts2 in 1:nts2){
    rho23[ts2] ~ dlogis(0,1)
  }
  # Pr(OS = 3 | TS = 3) 
  for(ts3 in 1:nts3){
    rho3[ts3] ~ dlogis(0,1)
  } 
}