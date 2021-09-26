# signature: double(l, p) -> double(k)
postpros <- list(
  id = \(X) X,
  norm = \(X) sqrt(rowSums(X^2)),
  sel1_norm = \(X) abs(X[,1]),
  sel123_norm = \(X) sqrt(rowSums(X[,1:3]^2)),
  sel1_pad20 = \(X) X[21:(nrow(X)-20),1],
  sel2_pad20 = \(X) X[21:(nrow(X)-20),2],
  sel3_pad20 = \(X) X[21:(nrow(X)-20),3]
)

register_postpro <- function(name, fun) {
  postpros[[name]] <- fun
}

get_postpro_names <- function() {
  names(postpros)
}

get_postpros <- function() {
  postpros
}

