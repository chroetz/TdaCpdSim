cusum_mat_ts <- function(y) {
  y <- as.matrix(y)
  n <- nrow(y)
  s <- apply(y, 2, \(x) cumsum(x[-n] - mean(x)))
  s_norm <- rowSums(abs(s)^2)
  k <- 1:(n-1)
  corr <- k/n*(n-k)
  s_corr <- s_norm / corr
  s_corr
}

cusum_mat <- function(y) {
  which.max(cusum_mat_ts(y))
}

