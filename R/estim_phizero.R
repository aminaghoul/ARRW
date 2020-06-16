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

#' cost
#' @description Return the values of the cost function
#' @param y a vector of observations
#' @param estim a vector of the estimator
#' @return a vector of the value of the cost function for y

cost <- function(y, estim)
{
  n <- length(y)
  cout <- rep(0, n)

  cout[1] <-  (y[1] - estim[1])^2

  for(t in 2:n)
  {
    cout[t] <- (y[t] - estim[t])^2 + (estim[t] - estim[t - 1])^2
  }

  return((cout))
}

##################################################################################################################
# A FINIR

#' OP
#' @description Return the values of the changepoints
#' @param cost cost function of the model
#' @param beta penalty constant
#' @param y observations data
#' @return a list of cp


OP <- function(cost = cout, beta = 0.5, y = data, muhat = estim)
{
  n <- length(y)
  f <- rep(- beta,n + 1)
  cp <- rep(NA,n + 1)

  for (taustar in 1:n)
  {
    tab <- rep(0, taustar)
    for(tau in 2:(taustar -1))
    {
      tab[tau] <- f[tau] + cost(y[(tau):taustar], estim[(tau):taustar]) + beta
    }
    print(tab)

  }
  f[taustar] <- min(tab)
  tab[taustar] <- f[taustar] + cost(y[1:1], estim[1:1]) + beta
  tauprim = which.min(tab)
  cp[taustar] = tauprim

  return(cp)
}



