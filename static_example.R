# A simulated example of the static multistate occupancy model

library(runjags)
library(mcmcplots)
library(scales)


# Load all the functions we need to simulate the data
sapply(
  list.files("./R/functions", full.names = TRUE),
  source, verbose = FALSE
)

# General bookkeeping
nsite <- 100
nsurvey <- 8
set.seed(134)

# make a named list for the parameters in the model.
#  For our owl example, we are assuming that the
#  areas where they breed and do not breed vary
#  in opposite directions along some environmental
#  gradient (e.g., urban intensity, forest cover, etc.)
my_params <- list(
  beta2 =  c(0.5, -1),    # latent state = 2
  beta3 =  c(-1,1),       # latent state = 3
  rho2g2 = c(-0.5, 0.75), # observed state = 2 | true state = 2
  rho2g3 = c(0, -1),      # observed state = 2 | true state = 3
  rho3g3 = c(0.5,0.5)     # observed state = 3 | true state = 3
)

# make design matrice
#  for latent state = 2
x <- matrix(1, ncol = 2, nrow = nsite)
x[,2] <- rnorm(nsite)

# for latent state = 3, assuming same covariate
u <- x

# for detection data model. Assuming the same
#  design matrix for all parameters.
k <- array(1, dim = c(nsite, nsurvey,2))

# Some sort of covariate that varies across
#  surveys. Note that if you also have covariates
#  that do not vary across surveys you just repeat
#  the values across the second dimension of k.
k[,,2] <- rnorm(nsite * nsurvey)

# combine them into a named list
my_covs <- list(
  beta2 = x,
  beta3 = x,
  rho2g2 = k,
  rho2g3 = k,
  rho3g3 = k
)

# simulate the observed data. Check out
#  this function in ./R/functions/simulate.R
#  if you are interested in how I did this.
y <- simulate_static(my_params, my_covs)

# make the data list for modeling
data_list <- list(
  y = y,
  x = x,
  u = u,
  k = k,
  nsite = nsite,
  nsurvey = nsurvey,
  nbeta2 = length(my_params$beta2),
  nbeta3 = length(my_params$beta3),
  nts2 = length(my_params$rho2g2),
  nts3 = length(my_params$rho2g3)
)

# fit the model
m1 <- run.jags(
  model = "./JAGS/static_multistate.R",
  monitor = c("beta2", "beta3", "rho2g2", "rho2g3", "rho3g3"),
  data = data_list,
  n.chains = 4,
  inits = init_static,
  adapt = 1000,
  burnin = 20000,
  sample = 20000,
  modules = "glm",
  method = "parallel"
)

# summarise model
m1sum <- summary(
  m1
)

saveRDS(m1, "static_mcmc.RDS")

# plot out the model coeffcients, compare
#  to true values
mcmcplots::caterplot(
  m1, reorder = FALSE
)
# and overlay the true simualted values
points(x = unlist(my_params), y = rev(1:10), pch = 19)
legend("topleft",
       c("Estimate", "Truth"),
       pch = 19,
       col = c(mcmcplots::mcmcplotsPalette(1), "black")
)

# Plot out some of the results.
#  make a design matrix with the covariate
#  you want to predict with. 
for_pred <- cbind(1, seq(-2,2, 0.05))

# get mcmc 
my_mcmc <- do.call("rbind", m1$mcmc)

# calculate latent state linear predictors/
#  grab is a utility function I put together
#  it's in ./R/functions/simulate.R. It just
#  pulls out the columns of a given linear
#  predictor with regex.
tmp1 <- for_pred %*% t(grab(my_mcmc, "beta2"))
tmp2 <- for_pred %*% t(grab(my_mcmc, "beta3"))

psi_preds <- list(
  state1 = matrix(NA, ncol = 3, nrow = nrow(for_pred)),
  state2 = matrix(NA, ncol = 3, nrow = nrow(for_pred)),
  state3 = matrix(NA, ncol = 3, nrow = nrow(for_pred)),
  marginal_occupancy = matrix(NA, ncol = 3, nrow = nrow(for_pred)),
  cond_breeding = matrix(NA, ncol = 3, nrow = nrow(for_pred))
)
# could write in a way to do this faster,
#  but this is easier to read.
pb <- txtProgressBar(max = nrow(for_pred))
for(i in 1:nrow(for_pred)){
  setTxtProgressBar(pb, i)
  tmp_pred <- t(
    apply(
      cbind(0, tmp1[i,], tmp2[i,]),
      1,
      softmax
    )
  )
  # marginal occupancy
  marg_occ <- tmp_pred[,2] + tmp_pred[,3]
  cond_occ <- tmp_pred[,3] / marg_occ
  # calculate quantiles of the 3 states
  tmp_pred <- apply(
    tmp_pred,
    2,
    quantile, 
    probs = c(0.025,0.5,0.975)
  )
  psi_preds$state1[i,] <- tmp_pred[,1]
  psi_preds$state2[i,] <- tmp_pred[,2]
  psi_preds$state3[i,] <- tmp_pred[,3]
  psi_preds$marginal_occupancy[i,] <- quantile(
    marg_occ,
    probs = c(0.025,0.5,0.975)
  )
  psi_preds$cond_breeding[i,] <- quantile(
    cond_occ,
    probs = c(0.025,0.5,0.975)
  )
}

# plot it out

plot(
  1~1, 
  type = "n", 
  xlim = c(-2,2),
  ylim = c(0,1),
  xlab = "Environmental covariate",
  ylab = "Occupancy probability",
  bty = "l",
  las = 1,
  cex.lab = 1.5
)
# 95% CI for state 2
polygon(
  x = c(for_pred[,2], rev(for_pred[,2])),
  y = c(psi_preds$state2[,1], rev(psi_preds$state2[,3])),
  col = scales::alpha("purple", 0.5),
  border = NA
)
# 95% CI for state 3
polygon(
  x = c(for_pred[,2], rev(for_pred[,2])),
  y = c(psi_preds$state3[,1], rev(psi_preds$state3[,3])),
  col = scales::alpha("gray50", 0.5),
  border = NA
)
# add median prediction
lines(
  x = for_pred[,2],
  y = psi_preds$state2[,2],
  lwd = 3,
  col = "purple"
)

lines(
  x = for_pred[,2],
  y = psi_preds$state3[,2],
  lwd = 3,
  lty = 3,
  col = "black"
)

legend(
  "topright",
  legend = c(
    "E(Owls present, no breeding)",
    "95% CI Owls present, no breeding",
    "E(Owls present, breeding)",
    "95% CI Owls present, breeding"
  ),
  fill = c(
    NA,
    scales::alpha("purple", 0.5),
    NA,
    scales::alpha("gray50", 0.5)
    ),
  border = c(NA, "black", NA, "black"),
  lty = c(1, NA, 3, NA),
  lwd = 3,
  col = c("purple", NA, "black", NA),
  bty = "n",
  cex = 1.2,
  seg.len = 1.5,
  x.intersp = c(2.2,1,2.2,1)
)

# plot out conditional probability of breeding | presence.

plot(
  1~1, 
  type = "n", 
  xlim = c(-2,2),
  ylim = c(0,1),
  xlab = "Environmental covariate",
  ylab = "Probability of breeding | owls present",
  bty = "l",
  las = 1,
  cex.lab = 1.5
)
# 95% CI for state 2
polygon(
  x = c(for_pred[,2], rev(for_pred[,2])),
  y = c(psi_preds$cond_breeding[,1], rev(psi_preds$cond_breeding[,3])),
  col = scales::alpha("purple", 0.5),
  border = NA
)

# add median prediction
lines(
  x = for_pred[,2],
  y = psi_preds$cond_breeding[,2],
  lwd = 3,
  col = "purple"
)



