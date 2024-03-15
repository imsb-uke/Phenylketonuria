library(stats)

x = rnorm(100)
y = rnorm(100)

gaussian2D <- function(x, y, params) {
  mux = params[1]
  muy = params[2]
  sigmax = params[3]
  sigmay = params[4]
  rho = params[5] # Correlation coefficient
  
  z = 1/(2*pi*sigmax*sigmay*sqrt(1-rho^2)) * exp(-1/(2*(1-rho^2)) * (((x-mux)^2/sigmax^2) + ((y-muy)^2/sigmay^2) - (2*rho*(x-mux)*(y-muy)/(sigmax*sigmay))))
  return(z)
}

nll = function(params){
  negLogLk = -sum(log(gaussian2D(x, y, params)))
  return(negLogLk)
}

err = function(params){
  z_hat = gaussian2D(x, y, params)
  mse = mean((z - z_hat)^2)
  return(mse)
}

params_initial = c(mean(x), mean(y), sd(x), sd(y), 0)

fit = optim(params_initial, nll, method = "L-BFGS-B", lower = c(-Inf, -Inf, 0, 0, -1), upper = c(Inf, Inf, Inf, Inf, 1))

print(fit$par)

plot(x, y, xlim = c(-3, 3), ylim = c(-3, 3))



x = 0:1000/100
y = exp(- (x-5)^2)
x = log(x)
plot(x, y)


