YeastIRO <- function(data, alpha1 = 0.05, alpha2 = 0.01, alpha3 = 0.01, iteration = 30, warm = 20)
{
  n <- dim(data)[1]
  p <- dim(data)[2]
  name <- colnames(data)
  leng <- iteration - warm
  miss <- which(is.na(data), arr.ind = TRUE)
  data_result <- data
  for (i in 1:length(miss[, 1]))
  {
    data_result[miss[i, 1], miss[i, 2]] <- median(na.omit(data[, miss[i, 2]]))
  }
  U <- NULL
  for (j in 1:iteration)
  {
    GraRes <- equSAR(data_result, ALPHA1 = alpha1, ALPHA2 = alpha2)
    U1 <- GraRes$score
    A2 <- GraRes$Adj
    thresh <- n/log(n)
    censor <- which(apply(A2, 1, sum) >= thresh)
    for (k in censor)
    {
      score_sub <- U1[U1[, 1] == k | U1[, 2] == k, ]
      score_ord <- score_sub[order(-score_sub[, 3]), ]
      ind <- as.numeric(apply(score_ord[, 1:2], 1, function(x) ifelse(x[1] == k, x[2], x[1])))
      A2[, k] <- 0
      A2[k, ] <- 0
      A2[ind[1:thresh], k] <- 1
      A2[k, ind[1:thresh]] <- 1
    }
    U <- cbind(U, U1[, 3])
    for (i in 1:length(miss[, 1]))
    {
      dep <- which(A2[miss[i, 2], ] == 1)
      if (length(dep) > 0)
      {
        combine <- data.frame(data_result[, miss[i, 2]], data_result[, dep])
        cov <- cov(combine)
        mu1 <- mean(data_result[, miss[i, 2]])
        mu2 <- apply(data_result[, dep, drop = FALSE], 2, mean)
        mu <- mu1 + cov[-1, 1] %*% solve(cov[-1, -1]) %*% t(data_result[miss[i, 1], dep, drop = FALSE] - mu2)
        sigma <- cov[1, 1] - cov[-1, 1] %*% solve(cov[-1, -1]) %*% cov[1, -1]
        data_result[miss[i, 1], miss[i, 2]] <- rnorm(1, mu, sqrt(sigma))
      }
      else
      {
        next
      }
    }
  }
  U <- cbind(U1[, 1:2], U)
  z <- -apply(U[, -(1:(iteration - leng + 2))], 1, sum)/leng
  q <- pnorm(-abs(z), log.p = TRUE)
  q <- q + log(2)
  s <- qnorm(q, log.p = TRUE)
  s <- (-1) * s
  U <- cbind(U[, 1:2], s)
  psiscore <- U
  U <- U[order(U[, 3]), 1:3]
  N <- p * (p - 1)/2
  ratio <- ceiling(N/1e+05)
  m <- floor(N/ratio)
  m0 <- N - m * ratio
  s <- sample.int(ratio, m, replace = TRUE)
  for (i in 1:length(s)) s[i] <- s[i] + (i - 1) * ratio
  if (m0 > 0)
  {
    s0 <- sample.int(m0, 1) + length(s) * ratio
    s <- c(s, s0)
  }
  Us <- U[s, ]
  qqqscore <- pcorselR(Us, ALPHA2 = alpha3)
  U <- psiscore
  U <- U[U[, 3] >= qqqscore, ]
  A <- matrix(rep(0, p * p), ncol = p)
  for (i in 1:nrow(U))
  {
    k1 <- U[i, 1]
    k2 <- U[i, 2]
    A[k1, k2] <- 1
    A[k2, k1] <- 1
  }
  colnames(A) <- name
  rownames(A) <- name
  return(A)
}