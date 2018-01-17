EyeICRO <- function(x, y, rate = 0.05, alpha1 = 0.1, alpha2 = 0.1, iteration = 30, warm = 20)
{
  leng <- iteration - warm
  n <- dim(x)[1]
  p <- dim(x)[2]
  sec <- sample(1:length(x), size = rate * length(x))
  x[sec] <- NA
  miss <- which(is.na(x), arr.ind = TRUE)
  x_new <- x
  for (i in 1:length(miss[, 1]))
  {
    x_new[miss[i, 1], miss[i, 2]] <- median(na.omit(x[, miss[i, 2]]))
  }
  x_result <- x_new
  beta_hat <- NULL
  for (iter in 1:iteration)
  {
    GraRes <- equSAR(huge.npn(x_result), ALPHA1 = alpha1, ALPHA2 = alpha2)
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
      mu <- (beta[k + 1] * (y[h] - sum(x_result[h, ] * beta[-1]) + x_result[h, k] * beta[k + 1] - 
                              beta[1]) + sigma_eps^2 * mu_x[1] * theta[1, 1] - sum(((x_sub[h, ] - mu_x) * theta[1, 
                                                                                                                ])[-1]) * sigma_eps^2)/(beta[k + 1]^2 + sigma_eps^2 * theta[1, 1])
      sigma <- sqrt(sigma_eps^2/(beta[k + 1]^2 + sigma_eps^2 * theta[1, 1]))
      x_result[miss[i, 1], miss[i, 2]] <- rnorm(1, mu, sigma)
    }
  }
  ave_hat <- as.matrix(beta_hat[, (iteration - leng + 1):iteration])
  ave_order <- matrix(1:(2 * p), ncol = 2)
  for (i in 1:p)
  {
    ave_order[i, ] <- c(i, length(which(ave_hat[i + 1, ] != 0)))
  }
  ave_rank <- ave_order[order(ave_order[, 2], decreasing = TRUE), 1]
  topVar <- ave_rank
  return(topVar)
}