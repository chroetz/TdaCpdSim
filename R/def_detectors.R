# signature: double(k), integer(1) -> double(1)
detectors <- list(
  argmax = \(y, n) which.max(y) + (n-length(y))/2
)

register_cpd <- function(name, fun) {
  detectors[[name]] <- fun
}

get_detector_names <- function() {
  names(detectors)
}

get_detectors <- function() {
  detectors
}
