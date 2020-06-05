
#' cost
#' @description Return the values of the cost function
#' @param ki float value of k(i)
#' @param y observations data
#' @return a float the value of the cost function for y

cost <- function(ki, y)
{

  n = length(y)
  # somme de (mut - mut-1)^2
  somme1 = 0

  # somme de (yt - mut)^2
  somme2 = 0
  for (t in 2:n)
  {

    sum1 = 0
    sum2 = 0
    sum3 = 0
    for(j in 1:(t - 1))
    {
      sum1 = sum1 + (y[j] * k(j - 1))^2
      sum2 = sum2 + y[j] * k(j - 1)

      for(i in 1:(j - 1))
      {
        sum3 = sum3 + y[i] * y[j] * k(n - j) * k(n - i)
      }
    }

    sum4 = 0
    sum5 = 0
    sum6 = 0
    for(j in t:n)
    {
      sum4 = sum4 + (y[j] * k(j - 1))^2
      sum5 = sum5 + y[j] * k(j - 1)

      for(i in t:(j - 1))
      {
        sum6 = sum6 + y[i] * y[j] * k(n - j) * k(n - i)
      }
    }

    terme1 = (cosh(1) - 1) * (cosh(2 * (n - t + 1)) - 1) * (sum1 + 2 * sum3)
    terme2 = (cosh(2(t - 1)) - 1) * (2*sum6 + sum4)
    terme3 = 2*(cosh(n - 2*t + 2) - cosh(n))*(sum2 + sum5)
    r = (4/k*k(0)*u^2)*(terme1 + terme2 + terme3)
    somme1 = somme1 + r


    # TODO À FINIR
    terme4 = (cosh(n - t + 1/2))^2 * (sum4 + 2 * sum6)
    terme5 = cosh(n - t + 1/2) * cosh(t - 1/2) * (sum5 * sum2)
    terme6 = (cosh(t - 1/2))^2
    terme7 = omega*k*k(0)* y[t] * ()
    r1 = terme4 + terme5 + terme6 + terme7
    somme2 = somme2 + r1

  }
  terme8 = #y1 - mu1
  res = somme1 + somme2 + terme8

  return(res)
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
