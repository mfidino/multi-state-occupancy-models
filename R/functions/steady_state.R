# Function to calculate the expected occupancy from:
#  "autologistic_logit_example.R"
# A function to calculate the expected occupancy
#  of the autologistic logit model:
#  Can be found at: "./R/functions/steady_state.R"
# Returns:
#  Named list.
#  The first four elements are nrow(pred_covs[[1]]) by 3 matrices
#    that contain the 0.025,0.5, and 0.975 quantiles.
#  --------------------------------------------------------------
#  omega = Pr(site is in state 2 or 3)
#  delta = Pr(site is in state 3 | site is in state 2 or 3)
#  omega_cond = Pr(site is in state 2 or 3 | site was in 
#                  state 2 or 3 at t-1)
#  delta_cond = Pr(site is in state 3 | site is in state 2 or 3
#                  at t and t-1)
#  --------------------------------------------------------------
#  The next element is a 3 by nrow(covs[[1]]) by 3 array.
#  Organized by 0.025,0.5, and 0.975 quantiles, number of
#  sites in the prediction design matrix, and number of states.
#  --------------------------------------------------------------
#  steady_states = Expected occupancy of each state. 
#  --------------------------------------------------------------
#  The final element is a 3 by nrow(covs[[1]]) by 3 by 3 array.
#  It is estimates for the 3 by 3 transition matrix. The first
#  dimension is the 0.025,0.5,and 0.975 quantiles.
#  --------------------------------------------------------------
#  psi = the transition probability matrix
#  --------------------------------------------------------------
ex_occ_logit <- function(mcmc,covs){
  omega <- plogis(
    grab(mcmc,"beta23") %*% t(covs$beta23)
  )
  delta <- plogis(
    grab(mcmc, "beta3") %*% t(covs$beta3)
  )
  omega_cond <- plogis(
    grab(mcmc,"beta23") %*% t(covs$beta23) + grab(mcmc,"theta23")
  )
  delta_cond <- plogis(
    grab(mcmc,"beta3") %*% t(covs$beta3) + grab(mcmc, "theta3")
  )
  # fill up transition matrix
  psi <-  array(NA, dim = c(nrow(mcmc), nrow(covs$beta23), 3,3))
  psi[,,1,1] <- 1 - omega
  psi[,,1,2] <- omega * (1 - delta)
  psi[,,1,3] <- omega * delta
  psi[,,2,1] <- 1 - omega_cond
  psi[,,2,2] <- omega_cond * (1 - delta)
  psi[,,2,3] <- omega_cond * delta
  psi[,,3,1] <- 1 - omega_cond
  psi[,,3,2] <- omega_cond * (1 - delta_cond)
  psi[,,3,3] <- omega_cond * delta_cond
  
  # Get quantiles of psi
  locs <- expand.grid(1:3, 1:3)
  psi_quantiles <- array(NA, dim= c(3, nrow(covs[[1]]), 3,3))
  for(i in 1:9){
    tmp <- psi[,,locs[i,1], locs[i,2]]
    tmp <-apply(
      tmp, 2, quantile, probs = c(0.025,0.5,0.975)
    )
    psi_quantiles[,,locs[i,1], locs[i,2]] <- tmp
  }
  
  # Create a matrix to store the steady states
  steady_states <- array(
    NA,
    dim = c(nrow(mcmc), nrow(covs[[1]]), 3)
  )
  
  # loop through each site and mcmc step
  pb <- txtProgressBar(max = nrow(mcmc))
  for(step in 1:nrow(mcmc)){
    setTxtProgressBar(pb,step)
    for(site in 1:nrow(covs[[1]])){
      mark <- new(
        "markovchain",
        states = as.character(1:3),
        byrow = TRUE,
        psi[step,site,,]
      )
      steady_states[step,site,] <- markovchain::steadyStates(mark)
    }
  }
    
  to_return <- list(
    omega = t(
      apply(
        omega, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    delta = t(
      apply(
        delta, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    omega_cond = t(
      apply(
        omega_cond, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    delta_cond = t(
      apply(
        delta_cond, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    steady_states = apply(
      steady_states, c(2,3), quantile, probs = c(0.025,0.5,0.975)
    ),
    psi = psi_quantiles
  )
  return(to_return)
}




# Function to calculate the expected occupancy from:
#  "autologistic_logit_example.R"
# A function to calculate the expected occupancy
#  of the autologistic logit model:
#  Can be found at: "./R/functions/steady_state.R"
# Returns:
#  Named list.
#  --------------------------------------------------------------
#  The first element is a 3 by nrow(covs[[1]]) by 3 array.
#  Organized by 0.025,0.5, and 0.975 quantiles, number of
#  sites in the prediction design matrix, and number of states.
#  --------------------------------------------------------------
#  steady_states = Expected occupancy of each state. 
#  --------------------------------------------------------------
#  The final element is a 3 by nrow(covs[[1]]) by 3 by 3 array.
#  It is estimates for the 3 by 3 transition matrix. The first
#  dimension is the 0.025,0.5,and 0.975 quantiles.
#  --------------------------------------------------------------
#  psi = the transition probability matrix
#  --------------------------------------------------------------
ex_occ_softmax <- function(mcmc,covs){

  # fill up transition matrix
  psi <-  array(NA, dim = c(nrow(mcmc), nrow(covs[[1]]), 3,3))
  psi[,,1,1] <- 1
  psi[,,1,2] <- exp(grab(mcmc,"beta2") %*% t(covs$beta2))
  psi[,,1,3] <- exp(grab(mcmc,"beta3") %*% t(covs$beta3))
  psi[,,2,1] <- 1
  psi[,,2,2] <- exp(
    grab(mcmc,"beta2") %*% t(covs$beta2) + grab(mcmc, "theta2g2")
  )
  psi[,,2,3] <- exp(
    grab(mcmc,"beta3") %*% t(covs$beta3) + grab(mcmc, "theta3g2")
  )
  psi[,,3,1] <- 1
  psi[,,3,2] <- exp(
    grab(mcmc,"beta2") %*% t(covs$beta2) + grab(mcmc, "theta2g3")
  )
  psi[,,3,3] <- exp(
    grab(mcmc,"beta3") %*% t(covs$beta3) + grab(mcmc, "theta3g3")
  )
  # Convert to probability
  for(i in 1:nrow(mcmc)){
    for(j in 1:nrow(covs[[1]])){
      psi[i,j,1,] <- psi[i,j,1,] / sum(psi[i,j,1,])
      psi[i,j,2,] <- psi[i,j,2,] / sum(psi[i,j,2,])
      psi[i,j,3,] <- psi[i,j,3,] / sum(psi[i,j,3,])
    }
  }
  
  # Get quantiles of psi
  locs <- cbind(rep(1:3, each = 3), rep(1:3, 3))
  psi_quantiles <- array(NA, dim= c(3, nrow(covs[[1]]), 3,3))
  for(i in 1:9){
    tmp <- psi[,,locs[i,1], locs[i,2]]
    tmp <-apply(
      tmp, 2, quantile, probs = c(0.025,0.5,0.975)
    )
    psi_quantiles[,,locs[i,1], locs[i,2]] <- tmp
  }
  
  # Create a matrix to store the steady states
  steady_states <- array(
    NA,
    dim = c(nrow(mcmc), nrow(covs[[1]]), 3)
  )
  
  # loop through each site and mcmc step
  pb <- txtProgressBar(max = nrow(mcmc))
  for(step in 1:nrow(mcmc)){
    setTxtProgressBar(pb,step)
    for(site in 1:nrow(covs[[1]])){
      mark <- new(
        "markovchain",
        states = as.character(1:3),
        byrow = TRUE,
        psi[step,site,,]
      )
      steady_states[step,site,] <- markovchain::steadyStates(mark)
    }
  }
  
  to_return <- list(
    steady_states = apply(
      steady_states, c(2,3), quantile, probs = c(0.025,0.5,0.975)
    ),
    psi = psi_quantiles
  )
  return(to_return)
}

ex_occ_dynamic <- function(mcmc,covs){
  
  psi23 <- plogis(covs$a23 %*% params$a23)
  psi3 <- plogis(covs$a3 %*% params$a3)
  gamma <- plogis(covs$b23 %*% params$b23)
  lambda <- plogis(covs$b3 %*% params$b3)
  phi23 <- plogis(covs$d23 %*% params$d23)
  phi3 <- plogis(covs$d3 %*% params$d3)
  
  # for first season
  s1_psi <- matrix(NA, nrow = nrow(covs$a23), ncol = 3)
  s1_psi[,1] <- (1 - psi23)
  s1_psi[,2] <- psi23 * (1 - psi3)
  s1_psi[,3] <- psi23 * psi3
  
  
  # for rest of the seasons
  psi <- array(NA, dim = c(nrow(covs$a23), 3, 3))
  
  
  # previous state == 1
  psi[,1,1] <- (1 - gamma)
  psi[,1,2] <- gamma * (1 - lambda)
  psi[,1,3] <- gamma * lambda
  # previous state == 2
  psi[,2,1] <- (1 - phi23)
  psi[,2,2] <- phi23 * (1 - lambda)
  psi[,2,3] <- phi23 * lambda
  # previous state == 3
  psi[,3,1] <- (1 - phi23)
  psi[,3,2] <- phi23 * (1 - phi3)
  psi[,3,3] <- phi23 * phi3
  
  
  gamma <- plogis(
    grab(mcmc,"b23") %*% t(covs$a23)
  )
  lambda <- plogis(
    grab(mcmc,"b3") %*% t(covs$a3)
  )
  phi23 <- plogis(
    grab(mcmc,"d23") %*% t(covs$a23)
  )
  phi3 <- plogis(
    grab(mcmc, "d3") %*% t(covs$a3)
  )
  
  # fill up transition matrix
  psi <-  array(NA, dim = c(nrow(mcmc), nrow(covs$a23), 3,3))
  # previous state == 1
  psi[,,1,1] <- (1 - gamma)
  psi[,,1,2] <- gamma * (1 - lambda)
  psi[,,1,3] <- gamma * lambda
  # previous state == 2
  psi[,,2,1] <- (1 - phi23)
  psi[,,2,2] <- phi23 * (1 - lambda)
  psi[,,2,3] <- phi23 * lambda
  # previous state == 3
  psi[,,3,1] <- (1 - phi23)
  psi[,,3,2] <- phi23 * (1 - phi3)
  psi[,,3,3] <- phi23 * phi3
  
  # Get quantiles of psi
  locs <- expand.grid(1:3, 1:3)
  psi_quantiles <- array(NA, dim= c(3, nrow(covs[[1]]), 3,3))
  for(i in 1:9){
    tmp <- psi[,,locs[i,1], locs[i,2]]
    tmp <-apply(
      tmp, 2, quantile, probs = c(0.025,0.5,0.975)
    )
    psi_quantiles[,,locs[i,1], locs[i,2]] <- tmp
  }
  
  # Create a matrix to store the steady states
  steady_states <- array(
    NA,
    dim = c(nrow(mcmc), nrow(covs[[1]]), 3)
  )
  
  # loop through each site and mcmc step
  pb <- txtProgressBar(max = nrow(mcmc))
  for(step in 1:nrow(mcmc)){
    setTxtProgressBar(pb,step)
    for(site in 1:nrow(covs[[1]])){
      mark <- new(
        "markovchain",
        states = as.character(1:3),
        byrow = TRUE,
        psi[step,site,,]
      )
      steady_states[step,site,] <- markovchain::steadyStates(mark)
    }
  }
  
  to_return <- list(
    gamma = t(
      apply(
        gamma, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    lambda = t(
      apply(
        lambda, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    phi23 = t(
      apply(
        phi23, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    phi3 = t(
      apply(
        phi3, 2, quantile, probs = c(0.025,0.5,0.975)
      )
    ),
    steady_states = apply(
      steady_states, c(2,3), quantile, probs = c(0.025,0.5,0.975)
    ),
    psi = psi_quantiles
  )
  return(to_return)
}

