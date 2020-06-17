#' var
#' @description Return the values of the variance of the estimator of the signal
#' @param val a list of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

vari <- function(val, sigma, n)
{
  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  k0 <- val$phip - val$phim
  u <- sqrt(1 + 4*val$om)
  R <- log(val$rp)

  res <- rep(0,n)
  for(i in 1:n)
  {
    terme1 <- 2*sinh(R)*(n + n*cosh(R*(2*i - 1))+2*i*sinh(R*(n - 2*i + 1))*sinh(R*n))
    terme2 <- 2*sinh(R*n)*cosh(R*(n - 2*i)) + cosh(R)*sinh(2*R*n)
    facteur <- (2*sigma^2)/(sinh(R)*(u^4)*(k*k0)^2)

    res[i] <- facteur*(terme1 + terme2)
  }
  res <- data.frame(1:n,res)
  colnames(res) = c("i", "variance")
  return(res)

}


#' var1
#' @description Return the values of the variance of the estimator of the signal
#' @param val a list of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

var1 <- function(val, sigma, n)
{
  rp <-val$rp
  omega <- val$om
  u <- sqrt(1 + 4*val$om)

  res <- rep(0,n)

  terme1 <- (2*omega + 1)*(1+rp^(-2*n))
  terme2 <- 2*u*(2*n*rp^(-2*n))/(1-rp^(-2*n))
  facteur <- (omega^2*sigma^2)/((1-rp^(-2*n))*u^3)
  for(i in 1:n)
  {
    terme3 <- 2*u*n*(rp^(2*i - 1-2*n) + rp^(1 - 2*i-2*n))/(1-rp^(-2*n))
    terme4 <- 2*u*i*(rp^(- 2*i + 1) - rp^(-2*n + 2*i - 1))
    terme5 <- 2*omega*(rp^(-2*i)+rp^(2*i-2*n))

    res[i] <- facteur*(terme1 + terme2+ terme3 + terme4 + terme5)
  }

  return(res)
}

