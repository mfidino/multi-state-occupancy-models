simulate_static <- function(params, covs){
  if(!is.list(params)){
    stop("params must be a list")
  }
  if(!is.list(covs)){
    stop("covs must be a list")
  }
  # latent state model
  psi <- matrix(0, ncol = 3, nrow = nrow(covs$beta2))
  psi[,2] <- covs$beta2 %*% params$beta2
  psi[,3] <- covs$beta3 %*% params$beta3
  # convert to probability
  psi <- t(
    apply(
      psi,
      1,
      softmax
    )
  )
  # simulate latent state
  z <- apply(
    psi,
    1,
    function(x) sample(1:3, 1, prob = x)
  )
  # data model
  eta <- array(0, dim = c(nrow(covs$beta2), dim(covs$rho2g2)[2],3,3))
  # TS = 1
  eta[,,1,1] <- 1  # OS = 1
  eta[,,1,2] <- 0  # OS = 2
  eta[,,1,3] <- 0  # OS = 3
  # TS = 2
  eta[,,2,1] <- 1  # OS = 1
  for(i in 1:dim(eta)[2]){
  eta[,i,2,2] <- exp(covs$rho2g2[,i,] %*% params$rho2g2) # OS = 2
  }
  eta[,,2,3] <- 0 # OS = 3
  # TS = 3
  eta[,,3,1] <- 1 # OS = 1
  for(i in 1:dim(eta)[2]){
    eta[,i,3,2] <- exp(covs$rho2g3[,i,] %*% params$rho2g3) # OS = 2
    eta[,i,3,3] <- exp(covs$rho3g3[,i,] %*% params$rho3g3) # OS = 2
  }
  # convert to probability
    for(i in 1:dim(eta)[2]){
      eta[,i,2,] <- eta[,i,2,] / rowSums(eta[,i,2,])
      eta[,i,3,] <- eta[,i,3,] / rowSums(eta[,i,3,])
    }

# sample y.
y <- matrix(
  NA,
  ncol = dim(eta)[2],
  nrow = dim(eta)[1]
)
for(site in 1:nrow(y)){
  for(survey in 1:ncol(y)){
    y[site,survey] <- sample(1:3, 1, prob = eta[site,survey,z[site],])
  }
}

return(y)
}

grab <- function(x,y){
  x[,grep(y, colnames(x))]
}

simulate_autologistic <- function(params, covs, link){
  if(link == "logit"){
    response <- .simulate_autologistic_logit(params = params, covs = covs)
  }
  if(link == "softmax"){
    response <- .simulate_autologistic_softmax(params = params, covs = covs)
  }
  return(response)
}


.simulate_autologistic_logit <- function(params, covs){
  if(!is.list(params)){
    stop("params must be a list")
  }
  if(!is.list(covs)){
    stop("covs must be a list")
  }
  # latent state model
  omega <- plogis(covs$beta23 %*% params$beta23)
  delta <- plogis(covs$beta3 %*% params$beta3)
  omega_cond <- plogis(covs$beta23 %*% params$beta23 + params$theta23)
  delta_cond <- plogis(covs$beta3 %*% params$beta3 + params$theta3)
  
  psi <- array(NA, dim = c(nrow(covs$beta23), 3, 3))
  
  # previous state == 1
  psi[,1,1] <- (1 - omega)
  psi[,1,2] <- omega * (1 - delta)
  psi[,1,3] <- omega * delta
  # previous state == 2
  psi[,2,1] <- (1 - omega_cond)
  psi[,2,2] <- omega_cond * (1 - delta)
  psi[,2,3] <- omega_cond * delta
  # previous state == 3
  psi[,3,1] <- (1 - omega_cond)
  psi[,3,2] <- omega_cond * (1 - delta_cond)
  psi[,3,3] <- omega_cond * delta_cond
  # convert to probability
  # simulate latent state
  z <- matrix(NA, ncol = params$nyear, nrow = nrow(covs$beta23))
  z[,1] <- apply(
    psi[,1,],
    1,
    function(x) sample(1:3, 1, prob = x)
  )
  for(year in 2:params$nyear){
    # grab the correct probabilities based on previous state
    for(site in 1:nrow(covs$beta23)){
      z[site,year] <- sample(
        1:3, 1, prob = psi[site,z[site,year-1],]
      )
    }
  }
  # data model
  eta <- array(0, dim = c(nrow(covs$beta2), dim(covs$rho23)[2],3,3))
  # TS = 1
  eta[,,1,1] <- 1  # OS = 1
  eta[,,1,2] <- 0  # OS = 2
  eta[,,1,3] <- 0  # OS = 3
  # TS = 2
  for(i in 1:dim(eta)[2]){
    eta[,i,2,1] <- 1 - plogis(covs$rho23[,i,] %*% params$rho23)
    eta[,i,2,2] <- plogis(covs$rho23[,i,] %*% params$rho23)
  }
  eta[,,2,3] <- 0 # OS = 3
  # TS = 3
  eta[,,3,1] <- 1 # OS = 1
  for(i in 1:dim(eta)[2]){
    eta[,i,3,1] <- 1 - plogis(covs$rho23[,i,] %*% params$rho23)
    eta[,i,3,2] <- plogis(covs$rho23[,i,] %*% params$rho23) * 
                  (1 - plogis(covs$rho3[,i,] %*% params$rho3))
    eta[,i,3,3] <- plogis(covs$rho23[,i,] %*% params$rho23) * 
                   plogis(covs$rho3[,i,] %*% params$rho3)
  }
  # sample y.
  y <- array(
    NA,
    dim = c(dim(eta)[1], dim(eta)[2], params$nyear)
  )
  for(site in 1:dim(y)[1]){
    for(survey in 1:dim(y)[2]){
      for(year in 1:dim(y)[3]){
        y[site,survey,year] <- sample(
          1:3,
          1, prob = eta[site,survey,z[site,year],])
      }
    }
  }
  return(y)
}

.simulate_autologistic_softmax <- function(params, covs){
  if(!is.list(params)){
    stop("params must be a list")
  }
  if(!is.list(covs)){
    stop("covs must be a list")
  }
  # latent state model
  psi <- array(NA, dim = c(nrow(covs$beta2), 3, 3))
  
  # previous state == 1
  psi[,1,1] <- 1
  psi[,1,2] <- exp(covs$beta2 %*% params$beta2)
  psi[,1,3] <- exp(covs$beta3 %*% params$beta3)
  # previous state == 2
  psi[,2,1] <- 1
  psi[,2,2] <- exp(covs$beta2 %*% params$beta2 + params$theta2g2)
  psi[,2,3] <- exp(covs$beta3 %*% params$beta3 + params$theta3g2)
  # previous state == 3
  psi[,3,1] <- 1
  psi[,3,2] <- exp(covs$beta2 %*% params$beta2 + params$theta2g3)
  psi[,3,3] <- exp(covs$beta3 %*% params$beta3 + params$theta3g3)
  # convert to probability
  
  # simulate latent state
  z <- matrix(NA, ncol = params$nyear, nrow = nrow(covs$beta2))
  z[,1] <- apply(
    psi[,1,],
    1,
    function(x) sample(1:3, 1, prob = x / sum(x))
  )
  for(year in 2:params$nyear){
    # grab the correct probabilities based on previous state
    for(site in 1:nrow(covs$beta2)){
      z[site,year] <- sample(
        1:3,
        1,
        prob = psi[site,z[site,year-1],] / sum(psi[site,z[site,year-1],])
      )
    }
  }
  # data model
  eta <- array(0, dim = c(nrow(covs$beta2), dim(covs$rho2g2)[2],3,3))
  # TS = 1
  eta[,,1,1] <- 1  # OS = 1
  eta[,,1,2] <- 0  # OS = 2
  eta[,,1,3] <- 0  # OS = 3
  # TS = 2
  for(i in 1:dim(eta)[2]){
    eta[,i,2,1] <- 1 
    eta[,i,2,2] <- exp(covs$rho2g2[,i,] %*% params$rho2g2)
  }
  eta[,,2,3] <- 0 # OS = 3
  # TS = 3
  eta[,,3,1] <- 1
  for(i in 1:dim(eta)[2]){
    eta[,i,3,2] <- exp(covs$rho2g3[,i,] %*% params$rho2g3)
    eta[,i,3,3] <- exp(covs$rho3g3[,i,] %*% params$rho3g3)
  }
  # sample y.
  y <- array(
    NA,
    dim = c(dim(eta)[1], dim(eta)[2], params$nyear)
  )
  for(site in 1:dim(y)[1]){
    for(survey in 1:dim(y)[2]){
      for(year in 1:dim(y)[3]){
        y[site,survey,year] <- sample(
          1:3,
          1, 
          prob = eta[site,survey,z[site,year],] /
            sum(eta[site,survey,z[site,year],])
        )
      }
    }
  }
  return(y)
}

transitions <- function(y){
    tmp <- apply(y, c(1,3), max)
  ans <- matrix(NA, ncol = ncol(tmp)-1, nrow = nrow(tmp))
  for(i in 1:ncol(ans)){
    ans[,i] <- paste0(tmp[,i],"_to_", tmp[,i+1])
  }
  return(table(ans))
}

simulate_dynamic <- function(params, covs){
  if(!is.list(params)){
    stop("params must be a list")
  }
  if(!is.list(covs)){
    stop("covs must be a list")
  }
  # latent state model
  psi23 <- plogis(covs$psi23 %*% params$psi23)
  psi3 <- plogis(covs$psi3 %*% params$psi3)
  gamma <- plogis(covs$gamma %*% params$gamma)
  lambda <- plogis(covs$lambda %*% params$lambda)
  phi23 <- plogis(covs$phi23 %*% params$phi23)
  phi3 <- plogis(covs$phi3 %*% params$phi3)
  
  # for first season
  s1_psi <- matrix(NA, nrow = nrow(covs$psi23), ncol = 3)
  s1_psi[,1] <- (1 - psi23)
  s1_psi[,2] <- psi_23 * (1 - psi3)
  s1_psi[,3] <- psi_23 * psi3
  
  
  # for rest of the seasons
  psi <- array(NA, dim = c(nrow(covs$psi), 3, 3))
  
  
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
  # convert to probability
  # simulate latent state
  z <- matrix(NA, ncol = params$nyear, nrow = nrow(covs$psi23))
  z[,1] <- apply(
    s1_psi,
    1,
    function(x) sample(1:3, 1, prob = x)
  )
  for(year in 2:params$nyear){
    # grab the correct probabilities based on previous state
    for(site in 1:nrow(covs$psi23)){
      z[site,year] <- sample(
        1:3, 1, prob = psi[site,z[site,year-1],]
      )
    }
  }
  # data model
  eta <- array(0, dim = c(nrow(covs$psi23), dim(covs$rho23)[2],3,3))
  # TS = 1
  eta[,,1,1] <- 1  # OS = 1
  eta[,,1,2] <- 0  # OS = 2
  eta[,,1,3] <- 0  # OS = 3
  # TS = 2
  for(i in 1:dim(eta)[2]){
    eta[,i,2,1] <- 1 - plogis(covs$rho23[,i,] %*% params$rho23)
    eta[,i,2,2] <- plogis(covs$rho23[,i,] %*% params$rho23)
  }
  eta[,,2,3] <- 0 # OS = 3
  # TS = 3
  eta[,,3,1] <- 1 # OS = 1
  for(i in 1:dim(eta)[2]){
    eta[,i,3,1] <- 1 - plogis(covs$rho23[,i,] %*% params$rho23)
    eta[,i,3,2] <- plogis(covs$rho23[,i,] %*% params$rho23) * 
      (1 - plogis(covs$rho3[,i,] %*% params$rho3))
    eta[,i,3,3] <- plogis(covs$rho23[,i,] %*% params$rho23) * 
      plogis(covs$rho3[,i,] %*% params$rho3)
  }
  # sample y.
  y <- array(
    NA,
    dim = c(dim(eta)[1], dim(eta)[2], params$nyear)
  )
  for(site in 1:dim(y)[1]){
    for(survey in 1:dim(y)[2]){
      for(year in 1:dim(y)[3]){
        y[site,survey,year] <- sample(
          1:3,
          1, prob = eta[site,survey,z[site,year],])
      }
    }
  }
  return(y)
}
