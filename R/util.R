# apply function fun to each consecutive subsequence in x of length l
rollapply <- function(x, l, fun, ...) {
  n <- length(x)
  out <- list()
  for (i in 0:(n-l))
    out[[i+1]] <- fun(x[i+1:l], ...)
  simplify2array(out)
}

dist_mat <- function(x, dist, ...) {
  n <- length(x)
  D <- matrix(0, nrow=n, ncol=n)
  for (i in 2:n) for (j in 1:(i-1)) {
    D[i,j] <- D[j,i] <- dist(x[[i]], x[[j]], ...)
  }
  D
}

which_max_pad <- function(x, pad=0) {
  n <- length(x)
  y <- x[(1+pad):(n-pad)]
  mx <- max(y)
  i <- mean((1:length(y))[y == mx])
  pad + i
}
