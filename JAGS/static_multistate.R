model{
  for(site in 1:nsite){
    # latent state model linear predictors
    psi[site,1] <- 1
    psi[site,2] <- exp( inprod( beta2, x[site,] ) )
    psi[site,3] <- exp( inprod( beta3, u[site,] ) )
    # latent state model, dividing by sum of psi
    #  to complete softmax function
    z[site] ~ dcat(
      psi[site,] / sum(psi[site,])
    )
    # data model linear predictors. I'm assuming we used
    #  the same design matrix for all of the detection
    #  probabilities. 
    for(survey in 1:nsurvey){
    # TS = True state, OS = observed state
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
    # data model, we use the latent state z[site] to
    #   index the appropriate row of eta.
      y[site,survey] ~ dcat(
        eta[site, survey,z[site],] / sum(eta[site,survey,z[site],])
      )
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
  # Pr(OS = 2 | TS = 2)
  for(ts2 in 1:nts2){
    rho2g2[ts2] ~ dlogis(0,1)
  }
  # Pr(OS = 2 | TS = 3) & Pr(OS = 3 | TS = 3)
  for(ts3 in 1:nts3){
    rho2g3[ts3] ~ dlogis(0,1)
    rho3g3[ts3] ~ dlogis(0,1)
  } 
}