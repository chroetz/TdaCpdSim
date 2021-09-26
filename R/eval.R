mae_guess <- function(total, true, support_len) {
  pad <- (total-support_len)/2
  mae <- double(length(total))
  for (i in seq_along(total)) {
    mae[i] <- mean(abs(pad[i] + 1:support_len[i] - true[i]))
  }
  mae
}

rmae_eval <- function(results, samplers, estimators) {
  results %>%
    left_join(samplers %>% select(s_name, n1, n2), by="s_name") %>%
    mutate(true = n1+0.5, total=n1+n2) %>%
    mutate(estim_centered = (total-len)/2 + estimation) %>%
    mutate(ae = abs(estim_centered-true)) %>%
    group_by(s_name, e_name, len, true, total) %>%
    summarize(mae = mean(ae), .groups="drop") ->
    mea

  mea %>%
    mutate(ref = mae_guess(total, true, len)) %>%
    mutate(rmae = mae / ref) ->
    rmae

  rmae %>%
    select(s_name, e_name, rmae)
}
