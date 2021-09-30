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

#' @export
register_postpro <- function(name, fun) {
  postpros[[name]] <- fun
}

#' @export
get_postpro_names <- function() {
  names(postpros)
}

#' @export
get_postpro <- function(name) {
  postpros[[name]]
}

