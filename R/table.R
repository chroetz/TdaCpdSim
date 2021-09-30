#' @export
pretty_table <- function(x) {
  kableExtra::kable_styling(
    kableExtra::kbl(x),
    c("striped", "hover", "condensed"),
    full_width = F)
}

bkgnd_cols <- colorRampPalette(c(
  rgb(0.5,1,0.5),
  rgb(1,1,0.5),
  rgb(1,0.5,0.5)))(100)

rmae_color <- function(rmae) {
  i <- 1 + round(rmae*100)
  i[i > 100] <- 100
  bkgnd_cols[i]
}

#' @export
pretty_color_table <- function(tb, color_cols, color_tb) {
  k <- pretty_table(tb)
  for (i in seq_along(color_cols)) {
    k <- kableExtra::column_spec(k, color_cols[i],
                                 background = rmae_color(color_tb[[i]]))
  }
  k
}
