
#---------------------------------------------------
#
# Binary Segmentation (\theta = E(x), P = Gauss)
#
#---------------------------------------------------


BinSeg <- function(ss, ee, xx, lambda)
{
  if (ee - ss < 1) return(NULL)
  
  nn <- ee - ss + 1
  
  SS <- sapply(ss:(ee-1), function(bb) abs(sqrt((ee-bb)/(nn*(bb-ss+1))) * sum(xx[ss:bb]) - sqrt((bb-ss+1)/(nn*(ee-bb))) * sum(xx[(bb+1):ee])))
  
  bb_ind <- which.max(SS)
  
  bb <- (ss:(ee-1))[bb_ind]
  
  if (abs(SS[bb_ind]) > lambda)
  {
    SS <- (SS-min(SS))/(max(SS)-min(SS)) * (max(xx)-min(xx)) + min(xx)
    
    lines(SS ~ c(ss:(ee-1)), col = bb)
    
    abline(v = bb, col = bb, lty = 2)  
    
    return(c(
      BinSeg(ss,bb,xx,lambda),
      bb,
      BinSeg(bb+1,ee,xx,lambda)))
    
  } else {
    
    return(NULL)

  }
}

doBinSeg <- function(xx)
{
  nn <- length(xx)
  
  lambda <- sqrt(3*mad(xx)*log(nn))
  
  return(BinSeg(1, nn, xx, lambda))
}

#----------------------
#
# A small example
#
#----------------------

## One change point
##

theta <- c(rep(0,50), rep(3,50))

xx <- theta + rnorm(length(theta))

xx |> plot(col = "grey", lwd = 2)

theta |> lines(lwd = 2, lty = 2)

doBinSeg(xx)



## Two change points
##

theta <- c(rep(0,30), rep(3,30), rep(0,30))

xx <- theta + rnorm(length(theta))

xx |> plot(col = "grey", lwd = 2)

theta |> lines(lwd = 2, lty = 2)

doBinSeg(xx)

