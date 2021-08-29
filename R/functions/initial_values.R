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
