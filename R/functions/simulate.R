


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
  eta[,,3,2] <- 0 # OS = 3
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
