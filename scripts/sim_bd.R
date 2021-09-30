library(tidyverse)
library(TdaCpdSim)

estimators <- read_csv("estimators.csv")
required_dists <- read_csv("required_dists.csv")

read_csv("samplers_base.csv") %>%
  mutate(
    n = 300, n1 = 200,
    filt = "Rips",
    kind = "birth_death",
    m = 20, rate = 0.05,
    reps = 4) ->
  samplers

seed <- 1
set.seed(seed)
results <- simulate_ts(samplers, estimators, required_dists, 4)

write_rds(
  list(results=results,
       samplers=samplers,
       estimators=estimators,
       required_dists=required_dists,
       seed=seed),
  paste0("sim_bd_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".RDS"))
