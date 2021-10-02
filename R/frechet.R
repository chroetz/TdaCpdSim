restricted_fm <- function(x, dist, ...) {
  D <- dist_mat(x, dist, ...)
  i <- which.min(rowMeans(D^2))
  x[[i]]
}

#' Fréchet Change Point Detection Stats
#'
#' @param dist_mat_sq nxn matrix of squared distances
fcp_stats <- function(dist_mat_sq) {
  D <- dist_mat_sq

  stopifnot(length(D) > 0)
  stopifnot(is.matrix(D))
  stopifnot(nrow(D) == ncol(D))

  n <- nrow(D)

  D_cumsum_left <- apply(D, 2, \(x) c(0, cumsum(x)))
  mu_i_left <- apply(D_cumsum_left, 1, which.min)
  D_cumsum_right <- apply(D, 2, \(x) rev(c(0, cumsum(rev(x)))))
  mu_i_right <- apply(D_cumsum_right, 1, which.min)

  sep <- 0:n
  u <- sep / n
  sel <- -c(1, n+1)
  corr <- u*(1-u) # correction factor

  # statistic 1
  mu_dist <- D[cbind(mu_i_left[sel], mu_i_right[sel])] # * corr[sel]

  # statistic 2
  v_trans_left <- D_cumsum_left[cbind(sep+1, mu_i_right)]/sep
  v_cis_left <- D_cumsum_left[cbind(sep+1, mu_i_left)]/sep
  v_trans_right <- D_cumsum_right[cbind(sep+1, mu_i_left)]/rev(sep)
  v_cis_right <- D_cumsum_right[cbind(sep+1, mu_i_right)]/rev(sep)
  v_mean_dist <- corr*(v_trans_left - v_cis_left +
                         v_trans_right - v_cis_right)^2

  # statistic 3
  v_var_dist <- corr*(v_cis_left - v_cis_right)^2

  cbind(mu_dist, v_mean_dist[sel], v_var_dist[sel])
}

#' A Fréchet variance statistic without need for calculating a Fréchet mean.
#'
#' @param dist_mat_sq nxn matrix of squared distances
inco_var <- function(dist_mat_sq) {
  D <- dist_mat_sq
  n <- nrow(D)
  u <- (0:n) / n
  prod_var <- sapply(0:n, \(k) {
    i <- seq_len(k)
    j <- setdiff(1:n, i)
    sum(D[i,i]) / k^2 - sum(D[j,j])/(n-k)^2})^2 * u*(1-u)
  as.matrix(prod_var[-c(1, n)])
}
