
#' muhat
#' @description muhat for phi = 0
#' @param y standard deviation in Random Walk
#' @param kis standard deviation in AR(1)
#' @param omega AR(1) autocorrelation parameter
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

  return(cout)
}

#' f
#' @description Return the values of the function F(t) + cost + beta
#' @param t
#' @param tstar penalty constant
#' @param F
#' @param y observation data
#' @param cost cost function
#' @param beta penalty constant
#' @return a list of res

f <- function(t,tstar,F,cost,y, beta)
{
  res = rep(0,tstar)
  for(i in 1:(tstar-1))
  {
    res[i] = F(t)+cost(ki, y[i+1:tstar])+beta
  }

  return(res)
}



#' OP
#' @description Return the values of the changepoints
#' @param cost cost function of the model
#' @param beta penalty constant
#' @param y observations data
#' @return a list of cp


OP <- function(cost, beta = 0.5, y)
{
  n = length(y)
  F = rep(- beta,n + 1)
  cp = rep(NA,n + 1)

  for (taustar in 2:(n + 1)){
    res = (f(t,tstar,F,cost,y, beta))
    F(taustar) = min(res)
    tauprim = which.min(res)
    cp[taustar] = tauprim

  }

  return(cp)
}


#' cost
#' @description Return the values of the cost function
#' @param ki float value of k(i)
#' @param y observations data
#' @return a float the value of the cost function for y

cost <- function(y)
{
  n = length(y)

  val <- def(sdEta = 0.3, sdNu = 0.5, phi = 0)
  k <- val$phip^2*val$rp^(n - 1) - val$phim^2*val$rm^(n - 1)
  omega <- val$om
  u <- sqrt(1 + 4*omega)
  print(k)

  # somme de (mut - mut-1)^2
  somme1 <- 0

  # somme de (yt - mut)^2
  somme2 <- 0
  for (t in 1:n)
  {

    sum1 <- 0
    sum2 <- 0
    sum3 <- 0
    for(j in 1:(t - 1))
    {
      sum1 <- sum1 + (y[j] * ki(val, j - 1))^2
      sum2 <- sum2 + y[j] * ki(val, j - 1)

      for(i in 1:(j - 1))
      {
        sum3 <- sum3 + y[i] * y[j] * ki(val, n - j) * ki(val, n - i)
      }
    }

    sum4 <- 0
    sum5 <- 0
    sum6 <- 0

    for(j in t:n)
    {
      sum4 <- sum4 + (y[j] * ki(val, n - j))^2
      sum5 <- sum5 + y[j] * ki(val, n - j)

      for(i in t:(j - 1))
      {
        sum6 <- sum6 + y[i] * y[j] * ki(val, n - j) * ki(val, n - i)
      }
    }

    terme1 <- (cosh(1) - 1) * (cosh(2 * (n - t + 1)) - 1) * (sum1 + 2 * sum3)
    terme2 <- (cosh(2*(t - 1)) - 1) * (2*sum6 + sum4)
    terme3 <- 2*(cosh(n - 2*t + 2) - cosh(n))*(sum2 + sum5)
    r <- (4/k*ki(val, 0)*u^2)*(terme1 + terme2 + terme3)
    somme1 <- somme1 + r

    sum7 <- sum1 +  (y[t] * ki(val, t - 1))^2
    sum8 <- sum3 + sum2
    sum9 <- sum2 + y[t] * ki(val, t - 1)
    sum10 <- sum5 - y[t] * ki(val, n - t)

    terme4 <- (cosh(n - t + 1/2))^2 * (sum7 + 2 * sum8)
    terme5 <- cosh(n - t + 1/2) * cosh(t - 1/2) * (sum9 * sum10)
    terme6 <- (cosh(t - 1/2))^2*(sum4 + 2*sum6)
    terme7 <- omega*k*ki(val, 0)* y[t] * (2*ki(val, n - t)*sum9 + 2*ki(val, t - 1)*sum10 - omega*k*ki(val, 0)* y[t])

    r1 <- (1/(k*ki(val,0)))*((4/u^2)*terme4 + (8/u^2)*terme5 + (4/u^2)*terme6 + (k*ki(val, 0)*y[t])*terme7)
    somme2 <- somme2 + r1

  }
  res <- somme1 + somme2 +

  return(res)
}

data = c(12,45,12,88,78,0.2,13,14,68,59,48)
res <- cost(y=data)

res
#' f
#' @description Return the values of the function F(t) + cost + beta
#' @param t
#' @param tstar penalty constant
#' @param F
#' @param y observation data
#' @param cost cost function
#' @param beta penalty constant
#' @return a list of res

f <- function(t,tstar,F,cost,y, beta)
{
  res = rep(0,tstar)
  for(i in 1:(tstar-1))
  {
    res[i] = F(t)+cost(ki, y[i+1:tstar])+beta
  }

  return(res)
}



#' OP
#' @description Return the values of the changepoints
#' @param cost cost function of the model
#' @param beta penalty constant
#' @param y observations data
#' @return a list of cp


OP <- function(cost, beta = 0.5, y)
{
  n = length(y)
  F = rep(- beta,n + 1)
  cp = rep(NA,n + 1)

  for (taustar in 2:(n + 1)){
    res = (f(t,tstar,F,cost,y, beta))
    F(taustar) = min(res)
    tauprim = which.min(res)
    cp[taustar] = tauprim

  }

  return(cp)
}
