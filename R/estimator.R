
#' def
#' @description Return the values of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param sdEta standard deviation in Random Walk
#' @param sdNu standard deviation in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @return a list of rPlus, rMinus, varphiPlus, varphiMinus, omega

def <- function(sdEta = 0.4, sdNu = 0.5, phi = 0.6)
{
  omega <- (sdNu/sdEta)^2
  u <- sqrt((1 + phi)^2 + 4*omega)
  rPlus <- (u + 1- phi)/(u - 1 + phi)
  rMinus <- 1/rPlus
  varphiPlus <- (2*(phi - 1)*(u + 1 + phi))/((u - phi + 1)*(u + phi - 1))
  varphiMinus <- (2*(phi - 1)*(- u + 1 + phi))/((u - phi + 1)*(u + phi - 1))

  return(list(rp = rPlus, rm = rMinus, phip = varphiPlus, phim = varphiMinus, om = omega))
}

#' ki
#' @description Return the value of ki
#' @param val a list of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param n an integer
#' @return a list of the values of ki and the value of k

ki <- function(val, n)
{
  phip <- val$phip
  phim <- val$phim
  rp <- val$rp
  rm <- val$rm

  res <- rep(0, n + 1)
  for(i in 1:(n + 1))
  {
    res[i] <- phip*rp^(i-1) - phim*rm^(i-1)
  }

  k <- phip^2*rp^(n - 1) - phim^2*rm^(n - 1)


  return(list(ki = res, k = k))
}


#' mu_hat
#' @description Estimator of mu
#' @param val a list of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param n an integer
#' @param y a vector of observations
#' @return a vector of the estimator of mu

mu_hat <- function(val, n=500,y = datag)
{
  muhat <- rep(0,n)

  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  k0 <- ki(val,0)
  omega <- val$om

  for(i in 1:n)
  {
    s1 <- 0
    if(i > 1)
    {
      for(j in 1:(i - 1))
      {
        knj <- ki(val, n - j)
        s1 <- s1 + y[j]*knj
      }
    }
    s2 <- 0
    if(i < n)
    {
      for(j in (i + 1):(n))
      {
        kj1 <- ki(val, j - 1)
        s2 <- s2 + y[j]*kj1
      }
    }
    ki1 <-  ki(val, i - 1)
    kni <- ki(val, n - i)
    estimateuri <- (ki1*s1 + ki1*kni*y[i] + kni*s2)/(omega*k0*k)
    muhat[i] <- estimateuri

  }
  return(muhat)
}



