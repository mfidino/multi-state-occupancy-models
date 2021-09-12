# The autologistic multistate occupancy model using softmax

library(runjags)
library(mcmcplots)
library(markovchain)

# Note: I parameterized this model complex than the aulogistic model
#       that used the logit link. I only did this to show how
#       there are different assumptions you can make with this model.
#       This does, however, make this model much more data hungry.
#       And it's not so much overall sample size, it's the number of
#       transitions among different states!

# Load all the functions we need to simulate the data
sapply(
  list.files("./R/functions", full.names = TRUE),
  source, verbose = FALSE
)

# -------------------------------------------
# Autologistic model that uses the logit link
# -------------------------------------------

# General bookkeeping
nsite <- 100
nsurvey <- 8
nyear <- 6
set.seed(188)

# OS = observed state
# TS = true state

my_params <- list(
  # latent state 2 parameters
  beta2 = c(-0.5,0.5),
  # latent state 3 parameters
  beta3 = c(-1,2),
  # detection, OS=2|TS=2
  rho2g2 = c(0.5,1),
  # detection, OS=2|TS=3
  rho2g3 = c(0.75,0),
  # detection, OS=3|TS=3
  rho3g3 = c(1,-0.5),
  # autologistic, 2 given 2 at t-1
  theta2g2 = 0.75,
  # autologistic, 3 given 2 at t-1
  theta3g2 = 2,
  # autologistic, 2 given 3 at t-1
  theta2g3 = -0.5,
  # autologistic, 3 given 3 at t-1
  theta3g3 = 0.5,
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

# Make a list for covariates
my_covs <- list(
  beta2 = x,
  beta3 = u,
  rho2g2 = k,
  rho2g3 = k,
  rho3g3 = k
)

# simulate data
y <- simulate_autologistic(
  params = my_params,
  covs = my_covs,
  link = "softmax"
)

transitions(y)

# make the data list for modeling
data_list <- list(
  # observed data
  y = y,
  # design matrices
  x = x,
  u = u,
  k = k,
  # for looping
  nsite = nsite,
  nsurvey = nsurvey,
  nyear = nyear,
  # number of parameters for different
  #  linear predictors
  nbeta2 = length(my_params$beta2),
  nbeta3 = length(my_params$beta3),
  nts2 = length(my_params$rho2g2),
  nts3 = length(my_params$rho3g3)
)

# fit the model
m1 <- run.jags(
  model = "./JAGS/autologistic_multistate_softmax.R",
  monitor = c("beta2", "beta3", "rho2g2", "rho2g3",
              "rho3g3", "theta2g2", "theta3g2",
              "theta2g3", "theta3g3"),
  data = data_list,
  n.chains = 4,
  inits = init_autologistic_softmax,
  adapt = 1000,
  burnin = 20000,
  sample = 10000,
  modules = "glm",
  method = "parallel"
)

# summarise model
msum <- summary(m1)

# save output for later
saveRDS(m1, "autologistic_softmax_mcmc.rds")


med_ests <- apply(
  do.call("rbind", m1$mcmc),
  2,
  median
)
med_ests <- med_ests[order(names(med_ests))]
pnames <- names(med_ests)

# Compare mcmc output to truth
jpeg("./figures/alm_softmax_mcmc.jpeg", quality = 100,
     width = 600, height = 480)
par(mar = c(5,6,2,2))

mcmcplots::caterplot(
  m1, 
  parms = pnames,
  reorder = FALSE,
  style = "plain",
  bty = "l",
  lwd = c(5,8),
  cex.labels = 1.25,
  cex.axis = 1.25
)

points(
  x = med_ests,
  y = rev(1:length(med_ests)),
  bg = mcmcplots::mcmcplotsPalette(1),
  pch = 21,
  cex = 2.25
)
# and overlay the true simualted values
points(
  x = c(
    unlist(
      my_params[c("beta2", "beta3",
                  "rho2g2", "rho2g3",
                  "rho3g3", "theta2g2",
                  "theta2g3", "theta3g2",
                  "theta3g3")]
    )
  ),
  y = rev(1:14), pch = 21,
  cex = 2.25,
  bg = "gray40"
)
legend("topright",
       c("Estimate", "Truth"),
       pch = 21,
       pt.bg = c(mcmcplots::mcmcplotsPalette(1), "gray40"),
       cex = 1.3,
       pt.cex = 2.25,
       bty = "n"
)
mtext("Coefficients", 1,line = 3, cex = 2)

dev.off()

# Calculate and plot expected occupancy.

# The environmental gradient used in the analysis
pred_covs <- list(
  beta2 = cbind(1, seq(-2.75, 2.75,0.05)),
  beta3 = cbind(1, seq(-2.75, 2.75,0.05))
)

# 5000 random samples from the mcmc
#  sub-sampling because we need to calculate
#  the steady state for each mcmc step and 
#  "site" (the number of rows in the prediction
#   design matrix).
my_mcmc <- do.call("rbind", m1$mcmc)
my_mcmc <- my_mcmc[sample(1:nrow(my_mcmc), 5000),]


# A function to calculate the expected occupancy
#  of the autologistic logit model:
#  Can be found at: "./R/functions/steady_state.R"
# Returns:
#  Named list.
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
to_plot <- ex_occ_softmax(
  mcmc = my_mcmc,
  covs = pred_covs
)

# Colors for plotting our the three states
my_cols <- c("#1b9e77", "#d95f02","#7570b3")
# the environmental gradient
x <- pred_covs$beta2[,2]

jpeg(
  "./figures/alm_softmax_states.jpeg",
  width = 750,
  height = 550,
  quality = 100,
  pointsize = 17
)
plot(
  1~1,
  type = "n", bty = "l", xlim = c(-2.75,2.75),
  ylim = c(0,1), xlab = "Environmental gradient",
  ylab = "Pr(Occupancy)", las = 1,
  cex.lab = 1.2,
  cex.axis = 1.2,
  xaxs = "i",
  yaxs = "i"
)

# State 1 95% CI. This is a wrapper function I
#  wrote for polygon(). Can be found in
#  "./R/functions/plotting_utilities.R
ribbon(
  x = x,
  y = t(to_plot$steady_states[c(1,3),,1]),
  col = my_cols[1],
  alpha = 0.5
)
# state 2 95% CI
ribbon(
  x = x,
  y = t(to_plot$steady_states[c(1,3),,2]),
  col = my_cols[2],
  alpha = 0.5
)
# State 3 95% CI
ribbon(
  x = x,
  y = t(to_plot$steady_states[c(1,3),,3]),
  col = my_cols[3],
  alpha = 0.5
)
# Median estimate state 1
lines(
  x = x,
  y = to_plot$steady_states[2,,1],
  col = my_cols[1], lwd = 3
)
# Median estimate state 2
lines(
  x = x,
  y = to_plot$steady_states[2,,2],
  col = my_cols[2], lwd = 3, lty = 2
)
# Median estimate state 3
lines(
  x = x,
  y = to_plot$steady_states[2,,3],
  col = my_cols[3], lwd = 3, lty = 3
)

# A little trick to have colored boxes
# in the legend
# Do a big wide line for 95% CI,
#  And write text in white
par(fg = "white",lend = 2)
legend(
  "right",
  c("Absent", "Present, not breeding", "Present & breeding"),
  lwd = 18,
  col = c(
    alpha(my_cols[1], 0.5),
    alpha(my_cols[2], 0.5),
    alpha(my_cols[3], 0.5)
  ),
  lty = c(1,1,1),
  cex = 1,
  seg.len = 2,
  bty = "n"
)
# Now do the median lines, and actually put the
#  text in.
par(fg = "black",lend = 0)
legend(
  "right", c("Absent", "Present, not breeding", "Present & breeding"),
  lwd = 3,
  col = c(
    my_cols[1],
    my_cols[2],
    my_cols[3]
  ),
  lty = c(1,2,3),
  cex = 1,
  seg.len = 2,
  bty = "n"
)
dev.off()

# Plot out transition matrix

jpeg("./figures/alm_softmax_transitions.jpeg",
     width = 8, height = 8, units = "in", quality = 100, res = 300
)
nr <- 2
hm <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
m <- matrix(
  c(rep(0,( nr * 4)),
    0, rep(hm[1:3], each = nr),0,
    0, rep(hm[1:3], each = nr),0,
    0, rep(hm[4:6], each = nr),0,
    0, rep(hm[4:6], each = nr),0,
    0, rep(hm[7:9], each = nr),0,
    0, rep(hm[7:9], each = nr),0,
    rep(0,( nr * 4))),
  ncol = (nr * 4) ,
  nrow = 8,
  byrow = TRUE
)
layout(m)

par(mar = c(2.5,1,0.5,2.5), xpd = NA)

# Get the locations of each element in the transition matrix
locs <- cbind(rep(1:3, each = 3), rep(1:3,3))

# plot out each of them
for(i in 1:9){
  plot(1~1, ylim = c(0,1), xlim = c(-2.75,2.75),
       xlab = "", ylab = "",
       bty = "l", type = "l", xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n")
  axis(1, seq(-2.75,2.75, length.out = 5), labels = FALSE,  tck = -0.025)
  
  
  if(i %in% c(1:3)){
    mtext(
      sprintf("%0.2f", seq(0, 1, 0.25)), 2, at = seq(0,1,0.25), 
      cex = 1, line = 1 ,las = 1)
  }
  if(i %in% c(3,6,9)){
    mtext(c(-2.75, 0, 2.75), 1, at = c(-2.75, 0, 2.75),
          line = 1.2)
  }
  axis(2, seq(0,1,0.25), labels = FALSE,  tck = -0.025)
  axis(2, seq(0,1,0.25/2), labels = FALSE,  tck = -0.025/2)
  
  ribbon(
    x = x,
    y = t(to_plot$psi[c(1,3),,locs[i,1],locs[i,2]]),
    col = my_cols[3],
    alpha = 0.5
  )
  
  lines(
    to_plot$psi[2,,locs[i,1], locs[i,2]] ~ x,
    lwd = 3, col = my_cols[3]
  )
  if(i == 1){
    par(xpd = NA)
    ty <- 1.25
    tc <- 1.6
    text(0, y = ty, "'No use' to...", cex = tc)
    mtext("Probability of transition", 2, line = 4.75, at = -0.7 ,cex = 1.6)
  }
  if(i == 4){
    par(xpd = NA)
    text(0, y = ty, "'Present, no breeding' to...", cex = tc)
  }
  if(i == 7){
    par(xpd = NA)
    text(0, y = ty, "'Present & breeding' to...", cex = tc)
    text(4, 0.5, "'No use'", srt = 270, cex = tc)
  }
  if(i == 8){
    text(4, 0.5, "'Present, no breeding'", srt = 270, cex = tc)
  }
  if(i == 9){
    text(4, 0.5, "'Present & breeding'", srt = 270, cex = tc)
    mtext("Environmental gradient", 1, line = 4.5, at = -6.75 ,cex = 1.6)
  }
  
}
dev.off()








