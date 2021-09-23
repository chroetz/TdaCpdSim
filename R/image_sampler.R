normalized_coor <- function(x) {
  t(t(arrayInd(1:length(x), dim(x))-1)/(dim(x)-1)-0.5)*2
}

load_img <- function(file) {
  x <- png::readPNG(file)
  d <- dim(x)
  if (length(d) == 3) {
    x <- apply(x^2, c(1,2), \(x) sqrt(sum(x)))
  }
  stopifnot(is.matrix(x))
  x <- x / max(x)
  x <- t(x)
  x <- x[,ncol(x):1]
  x
}

sample_img <- function(n, x, sd=1/sqrt(length(x)), trans=NULL, rot=NULL) {
  X <- normalized_coor(x)
  idx <- (1:length(x))[x > 0]
  i <- sample(idx, size=n, prob=x[x>0], replace=TRUE)
  v <- X[i,,drop=FALSE]
  r <- v + rnorm(length(v), sd=sd)
  if (!is.null(trans)) {
    r <- r - rep(trans, each=nrow(r))
  }
  if (!is.null(rot)) {
    r <- r %*% rot
  }
  r
}

make_img_sampler <- function(file) {
  x <- load_img(file)

  Y <- sample_img(1e7, x)
  trans <- colMeans(Y)
  en <- eigen(cov(Y))
  rot <- en$vectors %*% diag(en$values^(-0.5)) %*% t(en$vectors)

  function(n) sample_img(n, x, trans=trans, rot=rot)
}

create_img_sampler <- function(img_names) {
  lst <- lapply(img_names, \(x) make_img_sampler(paste0("img/",x,".png")))
  nms <- names(img_names)
  if (is.null(nms)) nms <- img_names
  names(lst) <- nms
  lst
}

