ribbon <- function(x, y, density=NULL, angle=45, border=NA,
                   col=NA, lty= par("lty"),...,fillOddEven=FALSE,
                   alpha = NULL
){
  # error checks
  if(!any(is.na(col)) & length(col)>1){
    warning("Two values input to col. Only first element used.")
    col <- col[1]
  }
  # check if y is a matrix
  if(is.matrix(y)|is.data.frame(y)){
    y <- c(y[,1], rev(y[,2]))
    # check if x is half the length of y
    if(length(y)/length(x) == 2){
      x <- c(x, rev(x))
    }
  }
  # evaluate color and alpha channel
  if( is.na(col) ){
    my_col <- NA
  } else { # otherwise go through color process
    # 
    if(
      length(grep("^#", col)) == 1 & # if start with hash
      nchar(col)>7 # & alpha channel is present
    ){
      if(!is.null(alpha)){
        warning("col already has alpha channel, ignoring alpha argument.")
        my_col <- col
      } else {
        my_col <- col
      }
    } else {
      # get rgb
      my_rgbs <- col2rgb(col)
      # set color
      my_col <- rgb(
        my_rgbs[1],my_rgbs[2],my_rgbs[3],max = 255,alpha = 255 * alpha
      )
    } 
  }
  polygon(
    x = x, y = y, density = density, angle = angle,
    border = border, col = my_col, lty = lty,
    fillOddEven = fillOddEven, ...
  )
}

alpha <- function(col, alpha){
    # get rgb
    my_rgbs <- col2rgb(col)
    # set color
    my_col <- rgb(
      my_rgbs[1],my_rgbs[2],my_rgbs[3],max = 255,alpha = 255 * alpha
    )
    return(my_col)
}
