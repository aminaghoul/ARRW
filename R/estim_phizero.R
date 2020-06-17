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


#' var_emp
#' @description empirical variance of generated observations
#' @param n an integer of the number of observations
#' @param nb_rep an integer the number of repetitions of the simulation
#' @param sigma a float the standard deviation of y
#' @return a vector of the empirical variance
var_emp <- function(n, nb_rep, sigma)
{
  val <- def(sdEta = 0.8^2, sdNu = sigma^2, phi = 0)
  kisi <- ki(val, n)
  omega <- val$om

  vari <- matrix(0, nrow = nb_rep, ncol = n)

  for(i in 1:nb_rep)
  {
    #generate y
    Y <- dataRWAR(n = n, poisParam = .01, meanGap = 15, phi = 0, sdEta = 0.8^2, sdNu = sigma^2)
    datai= Y$y
    # calculate the estimator
    estim <- muhat(y = datai, kis = kisi, omega)
    vari[i,] <- estim
  }
  varempi <- rep(0, n)
  for(j in 1:n)
  {
    varempi[j] <- sd(vari[,j]) ^2
  }
  return(varempi)
}

