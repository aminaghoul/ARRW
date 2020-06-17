#' var_log
#' @description Return the values of the variance of the estimator of the signal
#' @param val a list of rPlus, rMinus, varphiPlus, varphiMinus, omega
#' @param sigma variance of epsilon
#' @param n number of observations
#' @return res a vector of the value of the variance

var_log <- function(val, sigma, n)
{
  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  k0 <- val$phip - val$phim
  u <- val$u
  omega <- val$om
  R <- log(val$rp)

  res <- rep(0,n)
  facteur <- (2*sigma^2)/(4*sinh(R)*(omega^4)*(k*k0)^2)
  for(i in 1:n)
  {
    terme1 <- 2*sinh(R)*(n + n*cosh(R*(2*i - 1)) + 2*i*sinh(R*(n - 2*i + 1))*sinh(R*n))
    terme2 <- 2*sinh(R*n)*cosh(R*(n - 2*i)) + cosh(R)*sinh(2*R*n)

    res[i] <- facteur*(terme1 + terme2)
  }

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
  u <- val$u

  res <- rep(0,n)

  terme1 <- (2*omega + 1)*(1+rp^(-2*n))
  terme2 <- 2*u*(2*n*rp^(-2*n))/(1-rp^(-2*n))
  facteur <- (sigma^2)/((1-rp^(-2*n))*u^3)
  for(i in 1:n)
  {
    terme3 <- 2*u*n*(rp^(2*i - 1-2*n) + rp^(1 - 2*i-2*n))/(1-rp^(-2*n))
    terme4 <- 2*u*i*(rp^(- 2*i + 1) - rp^(-2*n + 2*i - 1))
    terme5 <- 2*omega*(rp^(-2*i)+rp^(2*i-2*n))

    res[i] <- facteur*(terme1 + terme2+ terme3 + terme4 + terme5)
  }

  return(res)
}

#' var_emp
#' @description empirical variance of generated observations
#' @param n an integer of the number of observations
#' @param nb_rep an integer the number of repetitions of the simulation
#' @param sigma a float the standard deviation of y
#' @return a vector of the empirical variance
var_emp <- function(n, nb_rep, sigma,sdEta)
{
  val <- def(sdEta, sdNu = sigma, phi = 0)
  u <- val$u
  rp <- val$rp

  vari <- matrix(0, nrow = nb_rep, ncol = n)
  mu <- rnorm(n, 0, sdEta)

  for(i in 1:nb_rep)
  {
    epsilon <- rnorm(n, 0, sigma)
    y <- epsilon + mu
    estim <- muhat1(y, u, rp)
    vari[i,] <- estim
  }
  varempi <- rep(0, n)
  for(j in 1:n)
  {
    varempi[j] <- sd(vari[,j]) ^2
  }
  return(varempi)
}

