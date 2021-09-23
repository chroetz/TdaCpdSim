# sample n1 points via distri1(m1) and n2 points via distri2(m12)
sample_iid <- function(distri1, m1, n1, distri2, m2, n2) {
  c(replicate(n1, distri1(m1), simplify=FALSE),
    replicate(n2, distri2(m2), simplify=FALSE))
}

sample_birth_death <- function(distri1, distri2, n1, n2, m, rate) {
  n <- n1 + n2
  locations <- distri1(m)
  lifetimes <- rexp(m, rate=rate)
  # If CP should be between n1 and n1+1, then threshold for changing form distri1 to distri2
  # for new births should be at n1+0.5 - 0.5*(1/rate). Then on average there will be
  # 50% points form distri1 and 50% points from distri2 at n1+0.5.
  threshold <- n1+0.5 - 0.5*(1/rate)
  res <- list()
  res[[1]] <- locations
  for (i in 2:n) {
    lifetimes <- lifetimes-1
    sel <- lifetimes<0
    cnt <- sum(sel)
    locations[sel,] <- if (i <= threshold) distri1(cnt) else distri2(cnt)
    lifetimes[sel] <- rexp(cnt, rate=rate)
    res[[i]] <- locations
  }
  res
}

sample_point_ts <- function(sampler, def) {
  switch(
    sampler$kind,
    birth_death = sample_birth_death(
      def$distri[[sampler$distri1]],
      def$distri[[sampler$distri2]],
      sampler$n1, sampler$n2,
      sampler$m,
      sampler$rate),
    iid = sample_iid(
      def$distri[[sampler$distri1]],
      sampler$m1,
      sampler$n1,
      def$distri[[sampler$distri2]],
      sampler$m2,
      sampler$n2),
    stop("unknown kind of sampler ", sampler$kind))
}


call_proc <- function(proc, pre_pros, dist_mats, pd_ts) {
  proc_fun <- pre_pros[[proc$pre_pro]]
  proc_fun_args <- list()
  if (proc$arg_dist) {
    proc_fun_args$dist_mat <- semi_join(
      dist_mats,
      proc,
      by=setdiff(names(dist_mats), "matrix"))$matrix[[1]]
  }
  if (proc$arg_pds)
    proc_fun_args$pds <- pd_ts
  do.call(proc_fun, proc_fun_args)
}

#' @export
simulate <- function(samplers, estimators, def) {
  sim_pt <- proc.time()
  cat("sim: start\n")

  estimators %>%
    filter(arg_dist) %>%
    select(starts_with("arg_dist_")) %>%
    distinct() ->
    required_dist

  lapply(seq_len(nrow(samplers)), \(i) {
    samp_pt <- proc.time()
    cat("  sampler:", i, "/", nrow(samplers), "start\n")
    s <- (samplers %>% filter(s_id == i))[1,]
    lapply(seq_len(s$reps), \(j) {
      rep_pt <- proc.time()
      cat("    rep:", j, "/", s$reps)
      point_ts <- sample_point_ts(s, def)
      pd_ts <- map(point_ts, pd, filtration=s$filt)

      required_dist %>%
        rowwise() %>%
        mutate(matrix = list(dist_mat_wasserstein(pd_ts, arg_dist_power, arg_dist_dim))) ->
        dist_mats

      # TODO:
      # from estimators, extract preprocessors, postprocessors, cpds


      preps <- list()
      values <- list()
      for (l in seq_len(nrow(estimators))) {
        prep <- call_proc(estimators[l, ], def$pre, dist_mats, pd_ts)
        preps[[l]] <- prep
        post <- def$post[[estimators$post_pro[[l]]]]
        v <- double(nrow(post))
        for (k in seq_len(nrow(post))) {
          column_sel <- post$column_sel[[k]]
          x <- if (length(column_sel) > 0) prep[,column_sel] else prep
          v[k] <- def$post_fun[[post$fun[[k]]]](x)
        }
        values[[l]] <- tibble(estim=v, post_id=seq_len(nrow(post)))
      }

      cat(", duration:", (proc.time()-rep_pt)[3], "s\n")
      tb <- tibble(proc_id=estimators$proc_id, rep_nr=j, values=values)
      tb$preps <- preps
      tb$len <- map_int(preps, NROW)
      tb
    }) -> r
    r %>%
      bind_rows() %>%
      mutate(opt_id = i, .before=1) ->
      res
    cat("  opt:", i, "/", nrow(samplers), "end, duration:", (proc.time()-samp_pt)[3], "s\n")
    res
  }) %>%
    bind_rows() ->
    results
  cat("sim: end, duration:", (proc.time()-sim_pt)[3], "s\n")
  results
}
