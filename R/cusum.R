cusum_ts <- function(y) {
  n <- length(y)
  s <- cumsum(y[-n] - mean(y))
  s_norm2 <- s^2
  k <- 1:(n-1)
  corr <- k/n*(n-k)
  s_corr <- s_norm2 / corr
  s_corr
}

cusum <- function(y) {
  which.max(cusum_ts(y))
}
