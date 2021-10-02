# signature: list(n), list(n), list(n) -> double(l, q)
prepros <- list(
  FR_RM_D0 = \(point_ts, pd_ts, dist_mats) fcp_stats(dist_mats[["0-Inf"]]^2),
  FR_RM_D1 = \(point_ts, pd_ts, dist_mats) fcp_stats(dist_mats[["1-Inf"]]^2),
  FR_VV_D0 = \(point_ts, pd_ts, dist_mats) inco_var(dist_mats[["0-Inf"]]^2),
  FR_VV_D1 = \(point_ts, pd_ts, dist_mats) inco_var(dist_mats[["1-Inf"]]^2),
  IYG_D0 = \(point_ts, pd_ts, dist_mats)
    projected_betty(pd_ts, epsilon=(0:49)*0.01, num_pc=3, dim=0),
  IYG_D1 = \(point_ts, pd_ts, dist_mats)
    projected_betty(pd_ts, epsilon=(0:49)*0.01, num_pc=3, dim=1),
  LT_D0 = \(point_ts, pd_ts, dist_mats) pd_values(pd_ts, dim=0),
  LT_D1 = \(point_ts, pd_ts, dist_mats) pd_values(pd_ts, dim=1)
)

#' @export
register_prepro <- function(name, fun) {
  force(fun)
  prepros[[name]] <- fun
  assignInMyNamespace("prepros", prepros)
  invisible(NULL)
}

#' @export
get_prepro_names <- function() {
  names(prepros)
}

#' @export
get_prepro <- function(name) {
  prepros[[name]]
}
