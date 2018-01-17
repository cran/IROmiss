RegICRO <- function(x, y, coef, type = "indep", alpha1 = 0.1, alpha2 = 0.05, iteration = 30, warm = 20)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  leng <- iteration - warm
  miss <- which(is.na(x), arr.ind = TRUE)
  x_new <- x
  for (i in 1:length(miss[, 1]))
  {
    x_new[miss[i, 1], miss[i, 2]] <- median(na.omit(x[, miss[i, 2]]))
  }
  x_result <- x_new
  if (type == "dep")
  {
    beta_hat <- NULL
    for (iter in 1:iteration)
    {
      GraRes <- equSAR(x_result, ALPHA1 = alpha1, ALPHA2 = alpha2)
      A2 <- GraRes$Adj
      cvfit <- cv.ncvreg(x_result, y)
      fit <- cvfit$fit
      beta <- as.numeric(fit$beta[, cvfit$min])
      beta_hat <- cbind(beta_hat, beta)
      eps_hat <- y - rep(beta[1], n) - x_result %*% as.matrix(beta[-1])
      sigma_eps <- sqrt(sum(eps_hat^2)/(n - length(which(beta != 0))))
      for (i in 1:length(miss[, 1]))
      {
        h <- miss[i, 1]
        k <- miss[i, 2]
        x_sub <- cbind(x_result[, k], x_result[, A2[, k] == 1])
        mu_x <- apply(x_sub, 2, mean)
        theta <- solve(cov(x_sub))
        mu <- (beta[k + 1] * (y[h] - sum(x_result[h, ] * beta[-1]) + x_result[h, k] * beta[k + 1] - beta[1]) + sigma_eps^2 * mu_x[1] * theta[1, 1] - sum(((x_sub[h, ] - mu_x) * theta[1, ])[-1]) * sigma_eps^2)/(beta[k + 1]^2 + sigma_eps^2 * theta[1, 1])
        sigma <- sqrt(sigma_eps^2/(beta[k + 1]^2 + sigma_eps^2 * theta[1, 1]))
        x_result[miss[i, 1], miss[i, 2]] <- rnorm(1, mu, sigma)
      }
    }
  }
  else if (type == "indep")
  {
    beta_hat <- NULL
    for (iter in 1:iteration)
    {
      cvfit <- cv.ncvreg(x_result, y)
      fit <- cvfit$fit
      beta <- as.numeric(fit$beta[, cvfit$min])
      beta_hat <- cbind(beta_hat, beta)
      eps_hat <- y - rep(beta[1], n) - x_result %*% as.matrix(beta[-1])
      sigma_eps <- sqrt(sum(eps_hat^2)/(n - length(which(beta != 0))))
      dep <- which(beta != 0) - 1
      for (i in 1:length(miss[, 1]))
      {
        h <- miss[i, 1]
        k <- miss[i, 2]
        if (k %in% dep)
        {
          mu_xk <- mean(x_result[, k])
          sigma_xk <- sd(x_result[, k])
          mu <- (beta[k + 1] * sigma_xk^2 * (y[h] - sum(x_result[h, ] * beta[-1]) + x_result[h, k] * beta[k + 1] - beta[1]) + sigma_eps^2 * mu_xk)/(beta[k + 1]^2 * sigma_xk^2 + sigma_eps^2)
          sigma <- sqrt(sigma_eps^2 * sigma_xk^2/(beta[k + 1]^2 * sigma_xk^2 + sigma_eps^2))
          x_result[miss[i, 1], miss[i, 2]] <- rnorm(1, mu, sigma)
        }
      }
    }
  }
  else
  {
    stop("Only 'dep' or 'indep' types are provided !")
  }
  ave <- NULL
  for (i in 1:(p + 1))
  {
    if (length(which(beta_hat[i, ] == 0)) >= leng/2) 
      beta_hat[i, ] <- 0
      ave <- cbind(ave, mean(as.numeric(beta_hat[i, (iteration - leng):iteration])))
  }
  true_beta <- c(1,which(coef!=0)+1)
  diff <- function(x)
  {
    res <- matrix(1:3, nrow = 1)
    fsr <- length(setdiff(which(x != 0), true_beta))/length(which(x != 0))
    nsr <- length(setdiff(true_beta, which(x != 0)))/length(true_beta)
    bias <- sum((c(1, coef) - x)^2)
    res[1] <- bias
    res[2] <- fsr
    res[3] <- nsr
    colnames(res) <- c("bias", "fsr", "nsr")
    return(res)
  }
  Var <- t(as.matrix(ave[which(ave!=0)][-1]))
  colnames(Var) <- which(ave[-1]!=0)
  result <- list()
  result$Var <- Var
  result$table <- diff(ave)
  return(result)
}