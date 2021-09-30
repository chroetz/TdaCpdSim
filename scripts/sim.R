library(tidyverse)
library(TdaCpdSim)

# s_name must be unique
samplers <- tribble(
  ~s_name, ~distri1, ~distri2,
  "o-.", "circ", "dot",
  "o-n", "circ", "horseshoe") %>%
  mutate(m=20, n=300, n1=200, rate=0.05,
  #mutate(m=2, n=42, n1=21, rate=0.05,
         kind = "birth_death",
         filt = "Rips",
         reps=3)

# d_name must be unique
required_dists <- tribble(
  ~d_name, ~power, ~dim,
  "0-Inf", Inf, 0,
  "1-Inf", Inf, 1)

# e_name must be unique
estimators <- tribble(
  ~e_name, ~prepro, ~select, ~trim, ~postpro, ~detector,
  "IYG_PC1_D0", "IYG_D0", "1",   0, "cusum", "argmax",
  "IYG_PC3_D0", "IYG_D0", "1:3", 0, "cusum", "argmax",
  "IYG_PC1_D1", "IYG_D1", "1",   0, "cusum", "argmax",
  "IYG_PC3_D1", "IYG_D1", "1:3", 0, "cusum", "argmax",
  "LT_D0",      "LT_D0",  "",    0, "cusum", "argmax",
  "LT_D1",      "LT_D1",  "",    0, "cusum", "argmax",
  "FR_M_D0",    "FR_RM_D0", "1", 20, "id", "argmax",
  "FR_MV_D0",   "FR_RM_D0", "2", 20, "id", "argmax",
  "FR_V_D0",    "FR_RM_D0", "3", 20, "id", "argmax",
  "FR_VV_D0",   "FR_VV_D0", "",  20, "id", "argmax",
  "FR_M_D1",    "FR_RM_D1", "1", 20, "id", "argmax",
  "FR_MV_D1",   "FR_RM_D1", "2", 20, "id", "argmax",
  "FR_V_D1",    "FR_RM_D1", "3", 20, "id", "argmax",
  "FR_VV_D1",   "FR_VV_D1", "",  20, "id", "argmax")

set.seed(1)
results <- simulate_ts(samplers, estimators, required_dists, 4)

write_rds(
  list(results=results,
       samplers=samplers,
       estimators=estimators,
       required_dists=required_dists),
  paste0("sim_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".RDS"))

# lst <- read_rds("sim_20210928-142547")
# list2env(lst, rlang::global_env())

# rmae <- rmae_eval(results, samplers, estimators)

