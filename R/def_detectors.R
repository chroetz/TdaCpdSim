# signature: double(k), integer(1) -> double(1)
detectors <- list(
  argmax = \(y, n) which.max(y) + (n-length(y))/2
)

#' @export
register_detector <- function(name, fun) {
  force(fun)
  detectors[[name]] <- fun
  assignInMyNamespace("detectors", detectors)
  invisible(NULL)
}

#' @export
get_detector_names <- function() {
  names(detectors)
}

#' @export
get_detector <- function(name) {
  detectors[[name]]
}
