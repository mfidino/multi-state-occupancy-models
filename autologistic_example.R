# Two examples of the autologistic multistate occupancy model

library(runjags)
library(mcmcplots)
library(scales)

# Load all the functions we need to simulate the data
sapply(
  list.files("./R/functions", full.names = TRUE),
  source, verbose = FALSE
)

# -------------------------------------------
# Autologistic model that uses the logit link
# -------------------------------------------

# General bookkeeping
nsite <- 75
nsurvey <- 8
nyear <- 4
set.seed(22025)

my_params <- list(
  beta23 = c(0.5,0.5),
  beta3 = c(0,1),
  theta23 = 0.75,
  theta3 = 0.5,
  rho23 = c(0.5,1),
  rho3 = c(-0.5,0.5),
  nyear = nyear
)

# For latent state being either 2 or 3
x <- cbind(1, rnorm(nsite))
# For conditional probability site is in state 3
u <- x

# for detection data model. Assuming the same
#  design matrix for all parameters.
k <- array(1, dim = c(nsite, nsurvey,2))

# Some sort of covariate that varies across
#  surveys. Note that if you also have covariates
#  that do not vary across surveys you just repeat
#  the values across the second dimension of k.
k[,,2] <- rnorm(nsite * nsurvey)

my_covs <- list(
 beta23 = x,
 beta3 = u,
 rho23 = k,
 rho3 = k
)

y <- simulate_autologistic(my_params, my_covs, "logit")

# make the data list for modeling
data_list <- list(
  y = y,
  x = x,
  u = u,
  k = k,
  nsite = nsite,
  nsurvey = nsurvey,
  nyear = nyear,
  nbeta23 = length(my_params$beta23),
  nbeta3 = length(my_params$beta3),
  nts2 = length(my_params$rho23),
  nts3 = length(my_params$rho3)
)

# fit the model
m1 <- run.jags(
  model = "./JAGS/autologistic_multistate_logit.R",
  monitor = c("beta23", "beta3", "rho23", "rho3", "theta23", "theta3"),
  data = data_list,
  n.chains = 4,
  inits = init_autologistic_logit,
  adapt = 1000,
  burnin = 20000,
  sample = 20000,
  modules = "glm",
  method = "parallel"
)

msum <- summary(m1)
round(msum,2)[,1:3]
