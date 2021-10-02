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


eval_betti <- function(pd, epsilon, dimension) {
  sel <- pd[, "dimension"] == dimension
  gt_b <- outer(epsilon, pd[sel,"Birth"], `>=`)
  gt_d <- outer(epsilon, pd[sel,"Death"], `>=`)
  rowSums(gt_b) - rowSums(gt_d)
}


projected_betty <- function(pds, epsilon, num_pc, dim) {
  beta0 <- sapply(pds, eval_betti, epsilon=epsilon, dimension=dim)
  pca <- stats::prcomp(beta0, center=FALSE, scale.=FALSE)
  pca$rotation[,1:num_pc]
}


pd_values <- function(pds, dim) {
  lst <- lapply(pds, \(x) sort(x[x[,"dimension"]==dim,"Death"]-
                                 x[x[,"dimension"]==dim,"Birth"],
                               decreasing = TRUE))
  nr <- max(sapply(lst, NROW))
  if (nr == 0) return(matrix(0, nrow=length(pds), ncol=1))
  X <- sapply(lst, \(x) x[1:nr])
  X <- matrix(X, nrow=nr)
  X[is.na(X)] <- 0
  t(X)
}


