# signature: double(k) -> double(1)
cpds <- list(
  argmax = which.max,
  cusum = cusum)

register_cpd <- function(name, fun) {
  cpds[[name]] <- fun
}

get_cpd_names <- function() {
  names(cpds)
}

get_cpds <- function() {
  cpds
}
