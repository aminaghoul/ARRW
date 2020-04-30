
estimator <- function(sdEta = 0.4, sdNu = 0.5, phi = 0.6)
{
  omega <- (sdNu/sdEta)^2
  u <- sqrt((1+phi)^2 + 4*omega)
  rPlus <- (u + 1- phi)/(u - 1 + phi)
  rMinus <- 1/rPlus
  varphiPlus <- (2*(phi - 1)*(u + 1 + phi))/((u - phi + 1)*(u + phi - 1))
  varphiMinus <- (2*(phi - 1)*(-u + 1 + phi))/((u - phi + 1)*(u + phi - 1))

  return(list(rp = rPlus, rm = rMinus, phip = varphiPlus, phim = varphiMinus))
}
e <- estimator(1,2,0.5)

