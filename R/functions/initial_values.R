init_static <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = rep(3, data_list$nsite),
      beta2 = rnorm(data_list$nbeta2),
      beta3 = rnorm(data_list$nbeta3),
      rho2g2 = rnorm(data_list$nts2),
      rho2g3 = rnorm(data_list$nts3),
      rho3g3 = rnorm(data_list$nts3),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"
      ),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

init_autologistic_logit <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(3, data_list$nsite, data_list$nyear),
      beta23 = rnorm(data_list$nbeta2),
      beta3 = rnorm(data_list$nbeta3),
      rho23 = rnorm(data_list$nts2),
      rho3 = rnorm(data_list$nts3),
      theta23 = rnorm(1),
      theta3 = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"
      ),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

init_autologistic_softmax <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(3, data_list$nsite, data_list$nyear),
      beta2 = rnorm(data_list$nbeta2),
      beta3 = rnorm(data_list$nbeta3),
      rho2g2 = rnorm(data_list$nts2),
      rho2g3 = rnorm(data_list$nts3),
      rho3g3 = rnorm(data_list$nts3),
      theta2g2 = rnorm(1),
      theta3g2 = rnorm(1),
      theta2g3 = rnorm(1),
      theta3g3 = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"
      ),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

init_dynamic <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(3, data_list$nsite, data_list$nyear),
      a23 = rnorm(data_list$na23),
      a3 = rnorm(data_list$na3),
      b23 = rnorm(data_list$nb23),
      b3 = rnorm(data_list$nb3),
      d23 = rnorm(data_list$nd23),
      d3 = rnorm(data_list$nd3),
      rho23 = rnorm(data_list$nrho23),
      rho3 = rnorm(data_list$nrho3),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"
      ),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}
