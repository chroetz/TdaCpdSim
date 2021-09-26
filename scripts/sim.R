library(tidyverse)

# s_name must be unique
samplers <- tribble(
  ~s_name, ~distri1, ~distri2,
  "o-.", "circ", "dot",
  "o-n", "circ", "horseshoe") %>%
  #mutate(m=20, n1=220, n2=80, rate=0.05,
  mutate(m=2, n1=21, n2=21, rate=0.05,
         kind = "birth_death",
         filt = "Rips",
         reps=2)

# d_name must be unique
required_dists <- tribble(
  ~d_name, ~power, ~dim,
  "0-Inf", Inf, 0,
  "1-Inf", Inf, 1)

# e_name must be unique
estimators <- tribble(
  ~e_name, ~prepro, ~postpro, ~cpd,
  "IYG_PC1_D0", "IYG_D0", "sel1_norm", "cusum",
  "IYG_PC3_D0", "IYG_D0", "sel123_norm", "cusum",
  "IYG_PC1_D1", "IYG_D1", "sel1_norm", "cusum",
  "IYG_PC3_D1", "IYG_D1", "sel123_norm", "cusum",
  "LT_D0", "LT_D0", "norm", "cusum",
  "LT_D1", "LT_D1", "norm", "cusum",
  "FR_M_D0", "FR_RM_D0", "sel1_pad20", "argmax",
  "FR_MV_D0", "FR_RM_D0", "sel2_pad20", "argmax",
  "FR_V_D0", "FR_RM_D0", "sel3_pad20", "argmax",
  "FR_VV_D0", "FR_VV_D0", "sel1_pad20", "argmax",
  "FR_M_D0", "FR_RM_D1", "sel1_pad20", "argmax",
  "FR_MV_D0", "FR_RM_D1", "sel2_pad20", "argmax",
  "FR_V_D0", "FR_RM_D1", "sel3_pad20", "argmax",
  "FR_VV_D0", "FR_VV_D1", "sel1_pad20", "argmax")

set.seed(0)
#results <- simulate(samplers, estimators, required_dists)
results <- simulate_parallel(samplers, estimators, required_dists, 1)

write_rds(
  list(results=results, samplers=samplers, estimators=estimators),
  paste0("sim_", format(Sys.time(), "%Y%m%d-%H%M%S")))

rmae <- rmae_eval(results, samplers, estimators)

