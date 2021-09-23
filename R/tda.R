# extract persistence diagram from point cloud (matrix Y)
# uses Rips filtration or alpha complex filtration
pd <- function(Y, max_dim=1, filtration=c("Rips", "Alpha")) {
  filtration <- match.arg(filtration)
  if (filtration == "Rips")
    pd <- TDA::ripsDiag(Y,
                        maxdimension=max_dim,
                        maxscale=sqrt(ncol(Y))*(max(Y)-min(Y))/2,
                        library="Dionysus")$diagram
  else if(filtration=="Alpha")
    pd <- TDA::alphaComplexDiag(Y,
                                maxdimension=max_dim)$diagram
  pd
}

dist_mat_wasserstein <- function(pds, power, dim) {
  if (power == Inf) {
    D <- dist_mat(pds, TDA::bottleneck, dimension=dim)
  } else {
    D <- dist_mat(pds, TDA::wasserstein, p=power, dimension=dim)
  }
  D
}

restricted_fm <- function(x, dist, ...) {
  D <- dist_mat(x, dist, ...)
  i <- which.min(rowMeans(D^2))
  x[[i]]
}

# Fréchet Change Point Detection Stats
#' @param dist_mat nxn matrix of distances (squared distances)
fcp_stats <- function(dist_mat) {
  D <- dist_mat

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
  mu_dist <- D[cbind(mu_i_left[sel], mu_i_right[sel])] * corr[sel]

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

inco_var <- function(dist_mat) {
  D <- dist_mat
  n <- nrow(D)
  u <- (0:n) / n
  # TODO: is factor u*(1-u) legit?
  prod_var <- sapply(0:n, \(k) {
    i <- seq_len(k)
    j <- setdiff(1:n, i)
    sum(D[i,i]) / k^2 - sum(D[j,j])/(n-k)^2}) * u*(1-u)
  prod_var[-c(1, n)]
}

eval_betti <- function(pd, epsilon, dimension) {
  sel <- pd[, "dimension"] == dimension
  gt_b <- outer(epsilon, pd[sel,"Birth"], `>=`)
  gt_d <- outer(epsilon, pd[sel,"Death"], `>=`)
  rowSums(gt_b) - rowSums(gt_d)
}

projected_betty <- function(pds, epsilon, num_pc, dim) {
  beta0 <- sapply(pds, eval_betti, epsilon=epsilon, dimension=dim)
  pca <- prcomp(beta0, center=FALSE, scale.=FALSE)
  pca$rotation[,1:num_pc]
}

make_projected_betty <- function(epsilon, num_pc, dim) {
  force(epsilon)
  force(num_pc)
  force(dim)
  function(pds) projected_betty(pds, epsilon, num_pc=num_pc, dim=dim)
}

pd_values <- function(pds, dim) {
  lst <- lapply(pds, \(x) sort(x[x[,"dimension"]==dim,"Death"]-
                                 x[x[,"dimension"]==dim,"Birth"],
                               decreasing = TRUE))
  nr <- max(sapply(lst, NROW))
  X <- sapply(lst, \(x) x[1:nr])
  X[is.na(X)] <- 0
  t(X)
}


