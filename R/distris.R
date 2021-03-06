normalized_coor <- function(x) {
  t(t(arrayInd(1:length(x), dim(x))-1)/(dim(x)-1)-0.5)*2
}

image_folder_path <- "img/"

#' @export
get_image_folder_path <- function() {
  image_folder_path
}

#' @export
set_image_folder_path <- function(p) {
  force(p)
  image_folder_path <- p
  assignInMyNamespace("image_folder_path", image_folder_path)
  invisible(NULL)
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

sample_img <- function(n, img_mat, sd=1/sqrt(length(img_mat)), trans=NULL, rot=NULL) {
  pix_coor <- normalized_coor(img_mat)
  idx <- (1:length(img_mat))[img_mat > 0]
  i <- sample(idx, size=n, prob=img_mat[img_mat>0], replace=TRUE)
  v <- pix_coor[i,,drop=FALSE]
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
  img_mat <- load_img(file)

  # calculate transformations for zero mean and identity covariance
  Y <- sample_img(1e7, img_mat)
  trans <- colMeans(Y)
  en <- eigen(cov(Y))
  rot <- en$vectors %*% diag(en$values^(-0.5)) %*% t(en$vectors)
  rm(Y, en, file)

  function(n) sample_img(n, img_mat, trans=trans, rot=rot)
}

#' @export
create_distris <- function(names, path=NULL, type="image") {
  if (is.null(path)) path <- image_folder_path
  switch(
    type,
    image = {
      lst <- lapply(names, \(x) make_img_sampler(paste0(path,x,".png")))
      nms <- names(names)
      if (is.null(nms)) nms <- names
      names(lst) <- nms
      lst
    },
    stop("unknown type of distri ", type))
}

