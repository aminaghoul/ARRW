library(ARRW)

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
#' @param i an integer
#' @return the value of ki

ki <- function(val, i)
{
  phip <- esti$phip
  phim <- esti$phim
  rp <- esti$rp
  rm <- esti$rm

  res <- phip*rp^(i) - phim*rm^(i)
  return(res)
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

  k <- esti$phip^2*esti$rp^(n - 1) - esti$phim^2*esti$rm^(n - 1)
  k0 <- ki(esti,0)
  omega <- esti$om

  for(i in 1:n)
    {
    s1 <- 0
    if(i >1)
    {
      for(j in 1:(i-1))
      {
        knj <- ki(esti, n - j)
        s1 <- s1 + y[j]*knj
      }
    }
    s2 <- 0
    if(i<n)
    {
      for(j in (i+1):(n))
      {
        kj1 <- ki(esti, n-j)
        s2 <- s2 + y[j]*kj1
      }
    }
    ki1 <-  ki(esti, i-1)
    kni <- ki(esti, n-i)
    estimateuri <- (ki1*s1+ki1*kni*y[i]+ki1*s2)/(omega*k0*k)
    muhat[i] <- estimateuri

  }
  return(muhat)
}



