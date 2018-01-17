SimRegDat <- function(n = 100, p = 200, coef, data.type = "indep", miss.type="MCAR", rate = 0.1)
{
  mu_x <- rep(0,p)
  sigma_x <- 2
  sigma_eps <- 1
  beta <- coef
  beta_0 <- rep(1,n)
  if(data.type=="indep"){
    cat("Independent covariates generated.\n")
    x <- rmvnorm(n,mu_x,diag(sigma_x^2,p))
    eps <- rnorm(n,0,sigma_eps)
    y <- beta_0 +x%*%beta+eps
    sec <- sample(1:length(x),size=rate*length(x))
    x[sec] <- NA
    result <- list(x=x,y=y,coef= beta)
  }else if(data.type=="dep"){
    cat("Dependent covariates generated.\n")
    C <- diag(1, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (abs(j - i) == 1) 
          C[i, j] = 0.5 else if (abs(j - i) == 2) 
            C[i, j] = 0.25
      }
    }
    A <- diag(0, p)
    sigma <- solve(C)
    sigma2 <- sigma
    for (i in 1:p) {
      for (j in 1:p) {
        if (sigma[i, j] <= sigma[j, i]) {
          sigma2[i, j] <- sigma[j, i]
        } else {
          sigma2[i, j] <- sigma[i, j]
        }
        if (C[i, j] != 0 && i != j) 
          A[i, j] <- 1
      }
    }
    x <- rmvnorm(n, mu_x, sigma2)
    eps <- rnorm(n,0,sigma_eps)
    y <- beta_0 +x%*%beta+eps
    if(data.type=="dep" && miss.type=="MAR"){
      beta <- c(1,beta)
      sigma_eps <- 1
      mu_hat <- matrix(1:(n*p),ncol=p)
      for(h in 1:n)
      {
        for(k in 1:p)
        {
          x_sub <- cbind(x[,k],x[,A[,k]==1]) 
          mu_x <- rep(0,ncol(x_sub))
          theta <- C[c(k,which(A[,k]==1)),c(k,which(A[,k]==1))]
          mu <- (beta[k+1]*(y[h]-sum(x[h,]*beta[-1])+x[h,k]*beta[k+1]-beta[1])+sigma_eps^2*mu_x[1]*theta[1,1]-
                   sum(((x_sub[h,]-mu_x)*theta[1,])[-1])*sigma_eps^2)/(beta[k+1]^2+sigma_eps^2*theta[1,1])
          mu_hat[h,k] <- abs(mu)
        }
        
      }
      sec <- sample(matrix(1:(n*p),ncol=p),size=rate*length(x),prob=mu_hat)
      x[sec] <- NA
    }else if(data.type=="indep" && miss.type=="MAR"){
      beta <- c(1,beta)
      sigma_eps <- 1
      sigma_xk <- 2
      mu_hat <- matrix(1:(n*p),ncol=p)
      for(h in 1:n)
      {
        for(k in 1:p)
        {
          mu <- (beta[k + 1] * sigma_xk^2 * (y[h] - sum(x[h, ] * beta[-1]) + x[h, k] * beta[k + 1] - beta[1]))/(beta[k + 1]^2 * sigma_xk^2 + sigma_eps^2)
          mu_hat[h,k] <- abs(mu)
        }
        
      }
      sec <- sample(matrix(1:(n*p),ncol=p),size=rate*length(x),prob=mu_hat)
      x[sec] <- NA
    }else{
      sec <- sample(1:length(x),size=rate*length(x))
      x[sec] <- NA
    }
    beta <- coef
    result <- list(x=x,y=y,coef= beta)
  } else{
    stop("Only 'dep' or 'indep' types are provided !")
  }
  return(result)
}