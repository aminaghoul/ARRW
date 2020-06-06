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


#' var
#' @description Return the values of the variance of the estimator of the signal
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

var <- function(sigma,n)
{
  val <- def(sdEta = sigma, sdNu = 0.5, phi = 0)
  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  u <- sqrt(1 + 4*val$om)

  res <- rep(0,n)
  for(i in 1:n)
  {
    terme1 <- cosh(1)*(4*i*sinh(n - 2*i + 1)*cosh(n) + sinh(2*n))
    print(terme1)
    terme2 <- 2*n*sinh(1)*(1 + cosh(2*i - 1))
    terme3 <- 2*sinh(n)*cosh(2*i - n)
    facteur <- ((4*sigma^2)/(sinh(1)*(u^4)*(k*ki(val,0))^2))
    print(sinh(1)*(u^4))
    res[i] <- facteur*(terme1 + terme2 + terme3)
  }

  return(res)

}
res <- var(6,10)
plot(res)
