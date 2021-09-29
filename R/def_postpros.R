cusum_mat_ts <- function(y) {
  n <- nrow(y)
  s <- apply(y, 2, \(x) cumsum(x[-n] - mean(x)))
  s_norm <- rowSums(s^2)
  k <- 1:(n-1)
  corr <- k/n*(n-k)
  s_corr <- s_norm / corr
  s_corr
}

# signature: double(l, p) -> double(k)
postpros <- list(
  id = \(X) as.vector(X),
  cusum = cusum_mat_ts
)

register_postpro <- function(name, fun) {
  postpros[[name]] <- fun
}

get_postpro_names <- function() {
  names(postpros)
}

get_postpros <- function() {
  postpros
}

