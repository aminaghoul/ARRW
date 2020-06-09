library(tidyverse)

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

#' var
#' @description Return the values of the variance of the estimator of the signal
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

var <- function(val, sigma, n)
{
  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  k0 <- val$phip - val$phim
  u <- sqrt(1 + 4*val$om)

  res <- rep(0,n)
  for(i in 1:n)
  {
    terme1 <- 2*sinh(1)*(n + n*cosh(2*i - 1)+2*i*sinh(n - 2*i + 1)*sinh(n))
    terme2 <- 2*sinh(n)*cosh(n - 2*i) + cosh(1)*sinh(2*n)
    facteur <- (4*sigma^2)/(sinh(1)*(u^4)*(k*k0)^2)

    res[i] <- facteur*(terme1 + terme2)
  }
  res <- data.frame(1:n,res)
  colnames(res) = c("i", "variance")
  return(res)

}
sigma <- 10.0
val <- def(sdEta = sigma, sdNu = 0.5, phi = 0.0)
val
val$om
res <- var(val, sigma,100)
ggplot(res) + ggtitle("Variance de l'estimateur du modèle \n marche aléatoire plus bruit") + geom_point(aes(x = i, y = variance))

