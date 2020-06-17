#' muhat
#' @description muhat for phi = 0
#' @param y a vector of observations
#' @param kis a list of k(i) and k
#' @param omega a float
#' @return a vector of the estimator of mu for phi=0

muhat <- function(y, kis, omega)
{
  n <- length(y)
  res <- rep(0, n)
  ki <- kis$ki
  k <- kis$k

  for(i in 1:n)
  {
    s1 <- 0
    for(j in 1:i)
    {
      s1 <- s1 + ki[j]*y[j]
    }
    if(i < n)
    {
      s2 <- 0
      for(j in (i + 1):n)
      {
        s2 <- s2 + ki[n - j + 1]*y[j]
      }
    }

    res[i] <- (ki[n - i + 1]*s1 + ki[i]*s2)/(omega*k*ki[1])
  }

  return(res)
}

#' muhat1
#' @description muhat for phi = 0
#' @param y a vector of observations
#' @param kis a list of k(i) and k
#' @param omega a float
#' @return a vector of the estimator of mu for phi=0

muhat1 <- function(y, u, rp)
{
  n <- length(y)
  res <- rep(0, n)

  facteur <- 1/(u*(1 - rp^(-2*n)))
  for(i in 1:n)
  {
    s1 <- 0
    for(j in 1:i)
    {
      s1 <- s1 + (rp^(j - 1/2) + rp^(1/2 - j))*y[j]
    }
    if(i < n)
    {
      s2 <- 0
      for(j in (i + 1):n)
      {
        s2 <- s2 + (rp^(1/2 - j) + rp^(j - 2*n - 1/2))*y[j]
      }
    }

    res[i] <- facteur*((rp^(1/2 - i) + rp^(i - 2*n - 1/2))*s1 + (rp^(i - 1/2) + rp^(1/2 - i))*s2)
  }

  return(res)
}



