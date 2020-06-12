
#' muhat
#' @description muhat for phi = 0
#' @param sdEta standard deviation in Random Walk
#' @param sdNu standard deviation in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @return a list of rPlus, rMinus, varphiPlus, varphiMinus, omega

muhat <- function(y, kis, omega)
{
  n <- length(y)

  res <- rep(0, n)

  ki <- kis$ki
  k <- kis$k

  for(i in 1:n)
  {
    s1 <- 0
    if(i > 1)
    {
      for(j in 1:(i - 1))
      {
        s1 <- s1 + ki[j]*y[j]
      }
    }
    if(i < n)
    {
      s2 <- 0
      for(j in (i + 1):n)
      {
        s2 <- s2 + ki[n - j + 1]*y[j]
      }
    }
    print(ki[i])

    res[i] <- (1/(omega*k*ki[1]))*(ki[n - i + 1]*s1 + ki[n - i + 1]*ki[i]*y[i] + ki[i]*s2)

  }

  return(res)
}



#' cost
#' @description Return the values of the cost function
#' @param ki float value of k(i)
#' @param y observations data
#' @return a float the value of the cost function for y

cost <- function(y, estim)
{
  n <- length(y)
  cout <- rep(0, n)

  cout[1] <-  (y[1] - estim[1])^2

  for(t in 2:n)
  {
    cout[t] <- (y[t] - estim[t])^2 + (estim[t] - estim[t - 1])^2
  }

  return(sum(cout))
}

estim
y

cout <- cost(data[1:2], estim)

cout

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
OP()


