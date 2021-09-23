library(tidyverse)

pre_pros <- list(
  fcp2 = \(dist_mat) fcp_stats(dist_mat^2),
  inco_var2 = \(dist_mat) inco_var(dist_mat^2),
  islambekov_dim0 = make_projected_betty(epsilon=(0:49)*0.01, num_pc=3, dim=0),
  islambekov_dim1 = make_projected_betty(epsilon=(0:49)*0.01, num_pc=3, dim=1),
  pd_values0 = \(pds) pd_values(pds, dim=0),
  pd_values1 = \(pds) pd_values(pds, dim=1))
post_pro_funs <- list(
  which_max_pad_20 = \(x) which_max_pad(x, pad=20),
  cusum_mat = cusum_mat)
post_pros <- list(
  fcp20 = tibble(
    name = c("mean_diff", "v_mean_diff", "v_var_diff"),
    column_sel = list(1, 2, 3),
    fun = "which_max_pad_20",
    pad = 20),
  whichmax20 = tibble(
    name = NA,
    column_sel = list(NULL),
    fun = "which_max_pad_20",
    pad = 20),
  cusumPC3 = tibble(
    name = c("PC1", "PC1&2&3"),
    column_sel = list(1, 1:3),
    fun = "cusum_mat",
    pad = 0),
  cusum = tibble(
    name = NA,
    column_sel = list(NULL),
    fun = "cusum_mat",
    pad = 0))
post_pros <- lapply(post_pros, \(x) rowid_to_column(x, "post_id"))

tibble(
  pre_pro = c("islambekov_dim0", "islambekov_dim1",
              "pd_values0", "pd_values1"),
  arg_dist = FALSE,
  arg_dist_power = NA,
  arg_dist_dim = NA,
  arg_pds = TRUE,
  post_pro = c("cusumPC3", "cusumPC3", "cusum", "cusum")
) -> direct_processing
tibble(pre_pro = c("fcp2", "inco_var2")) %>%
  crossing(expand.grid(arg_dist_power = Inf, arg_dist_dim = c(0, 1))) %>%
  mutate(arg_dist = TRUE, arg_pds=FALSE) %>%
  mutate(post_pro = ifelse(pre_pro == "fcp2", "fcp20", "whichmax20")) ->
  processing_frechet
bind_rows(direct_processing, processing_frechet) %>%
  rowid_to_column("proc_id") ->
  pc_opt

tribble(
  ~distri1, ~distri2,
  "circ", "dot",
  "circ", "horseshoe") %>%
  rowid_to_column("s_id") %>%
  mutate(m=20, n1=220, n2=80, rate=0.05,
         kind = "birth_death",
         filt="Rips", reps=2) ->
  samplers
distri_names <- unique(c(samplers$distri1,samplers$distri2))
distributions <- create_img_sampler(distri_names)

def <- list(
  distri=distributions,
  pre=pre_pros,
  post_fun=post_pro_funs,
  post=post_pros)

results <- simulate(samplers, pc_opt, def)
