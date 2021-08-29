# The softmax function. Convert any number of linear
#  predictors to a probability!
#
# arguments
#
# x = The log-linear predictors. One of these must be the
#       the baseline category (which equals zero.)  
#
softmax <- function(x){
  if(sum(x == 0) == 0){
    stop("No baseline category. One value in x must equal zero.")
  }
  exp(x) / sum(exp(x))
}
