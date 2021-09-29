# MAE of a uniform distribution on a centered support of length
# support_len given that the original length of the time series is len
# and the true change point is located at true.
# The function is vectorized.
mae_guess <- function(len, true, support_len) {
  pad <- (len-support_len)/2
  mae <- double(length(len))
  for (i in seq_along(len)) {
    mae[i] <- mean(abs(pad[i] + 1:support_len[i] - true[i]))
  }
  mae
}

#' Evaluate results of simulation by calculating the relative mean absolute error.
#'
#' @param results The results of a simulation.
#' @param samplers The table of samplers used in the simulation.
#' @param estimators The table of estimators used in the simulation.
#' @return A table containing the relative mean absolute error (rmae) for each run.
#' @export
rmae_eval <- function(results, samplers, estimators) {
  results %>%
    left_join(samplers %>% mutate(true = n1+0.5) %>% select(s_name, n, true),
              by="s_name") %>%
    mutate(ae = abs(estimation-true)) %>%
    group_by(s_name, e_name, len, true, n) %>%
    summarize(mae = mean(ae), .groups="drop") ->
    mea

  mea %>%
    mutate(ref = mae_guess(n, true, len)) %>%
    mutate(rmae = mae / ref) ->
    rmae

  rmae %>%
    select(s_name, e_name, rmae)
}
