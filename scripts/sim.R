library(tidyverse)

# TODO:
# * no id, only name
# * move defs to package, add register stuff function
# * eval
# * save simulation to file
# * parallel

# s_name must be unique
samplers <- tribble(
  ~s_name, ~distri1, ~distri2,
  "o-.", "circ", "dot",
  "o-n", "circ", "horseshoe") %>%
  mutate(m=20, n1=220, n2=80, rate=0.05,
  #mutate(m=20, n1=40, n2=40, rate=0.05,
         kind = "birth_death",
         filt="Rips", reps=10)

distri_names <- unique(c(samplers$distri1, samplers$distri2))
distributions <- create_img_sampler(distri_names)

# d_name must be unique
required_dists <- tribble(
  ~d_name, ~power, ~dim,
  "0-Inf", Inf, 0,
  "1-Inf", Inf, 1)

# signature: list(n), list(n), list(n) -> double(l, p)
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

# signature: double(k) -> double(1)
cpds <- list(
  argmax = which.max,
  cusum = cusum)

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

def <- list(
  distributions=distributions,
  prepros=prepros,
  postpros=postpros,
  cpds=cpds)

set.seed(0)
results <- simulate(samplers, estimators, required_dists, def)

rmae <- rmae_eval(results, samplers, estimators)
