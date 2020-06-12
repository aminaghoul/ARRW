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
    facteur <- (2*sigma^2)/(sinh(1)*(u^4)*(k*k0)^2)

    res[i] <- facteur*(terme1 + terme2)
  }
  res <- data.frame(1:n,res)
  colnames(res) = c("i", "variance")
  return(res)

}


#' var1
#' @description Return the values of the variance of the estimator of the signal
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

var1 <- function(val, sigma, n)
{
  rp <-val$rp
  omega <- val$om
  u <- sqrt(1 + 4*val$om)

  res <- rep(0,n)
  for(i in 1:n)
  {
    terme1 <- 4*(n + n*(rp^(2*i - 1) + rp^(1 - 2*i))/2)/((rp-rp^(-1))^(3)*(rp^n - rp^(-n))^2)
    terme2 <- 2*(i*(rp^(n - 2*i + 1) - rp^(-n + 2*i - 1)))/((rp-rp^(-1))^(3)*(rp^n - rp^(-n)))
    terme3 <- (rp^(2*n) - rp^(-2*n))*(rp+rp^(-1))/((rp-rp^(-1))^(4)*(rp^n - rp^(-n))^2)
    terme4 <- 2*(rp^(n - 2*i) + rp^(2*i - n))/((rp-rp^(-1))^(4)*(rp^n - rp^(-n)))
    facteur <- (sigma^2)/(omega^2)

    res[i] <- facteur*(terme1 + terme2 + terme3 + terme4)
  }
  res <- data.frame(1:n,res)
  colnames(res) = c("i", "variance")
  return(res)

}

