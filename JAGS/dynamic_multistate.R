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
    # For season = 1
    logit(psi23[site]) <- inprod(
      a23, x[site,]
    )
    logit(psi3[site]) <- inprod(
      a3, x[site,]
    )
    # For the remaining seasons
    
    # colonization of state 2 or 3
    logit(gamma[site]) <- inprod(
      b23, x[site,]
    )
    # Conditional colonization probability of state 3
    logit(lambda[site]) <- inprod(
      b3, u[site,]
    )
    # Persistence of state 2 or 3
    logit(phi23[site]) <- inprod(
      d23, x[site,]
    )
    # conditional persistence probability of state 3
    logit(phi3[site]) <- inprod(
      d3, u[site,]
    ) 
    # Fill out occupancy vector for first season
    s1_psi[site,1] <- 1 - psi23[site]
    s1_psi[site,2] <- psi23[site] * (1 - psi3[site])
    s1_psi[site,3] <- psi23[site] * psi3[site]
    # Fill in the transition matrix
    # latent state probabilities given state == 1
    psi[site,1,1] <- 1 - gamma[site]
    psi[site,1,2] <- gamma[site] * (1-lambda[site])
    psi[site,1,3] <- gamma[site] * lambda[site]
    # latent state probabilities given state == 2
    psi[site,2,1] <- 1 - phi23[site]
    psi[site,2,2] <- phi23[site] * (1-lambda[site])
    psi[site,2,3] <- phi23[site] * lambda[site]
    # latent state probabilities given state == 3
    psi[site,3,1] <- 1 - phi23[site]
    psi[site,3,2] <- phi23[site] * (1-phi3[site])
    psi[site,3,3] <- phi23[site] * phi3[site]
    # Estimate latent state year == 1. Setting to first
    #  row of psi because we have no prior knowledge
    #  on state of site before we started sampling.
    z[site,1] ~ dcat(
      s1_psi[site,]
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
  # initial occupancy
  for(ia23 in 1:na23){
    a23[ia23] ~ dlogis(0,1)
    
  }
  for(ia3 in 1:na3){
    a3[ia3] ~ dlogis(0,1)
  }
  # colonization
  for(ib23 in 1:nb23){
    b23[ib23] ~ dlogis(0,1)
  }
  # Pr(latent state 3)
  for(ib3 in 1:nb3){
    b3[ib3] ~ dlogis(0,1)
  }
  # persistence
  for(id23 in 1:nd23){
    d23[id23] ~ dlogis(0,1)
  }
  for(id3 in 1:nd3){
    d3[id3] ~ dlogis(0,1)
  }
  # Pr(OS = 2 or 3 | TS = 2 or 3)
  for(irho23 in 1:nrho23){
    rho23[irho23] ~ dlogis(0,1)
  }
  # Pr(OS = 3 | TS = 3) 
  for(irho3 in 1:nrho3){
    rho3[irho3] ~ dlogis(0,1)
  } 
}