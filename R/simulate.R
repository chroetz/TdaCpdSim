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
simulate <- function(samplers, estimators, required_dists, def) {
  sim_pt <- proc.time()
  cat("sim: start\n")

  prepros <- unique(estimators$prepro)

  lapply(seq_len(nrow(samplers)), \(i) {
    s <- (samplers %>% filter(s_id == i))[1,]
    samp_pt <- proc.time()
    cat("  sampler:", i, "/", nrow(samplers), "(",s$s_name,") start\n")
    lapply(seq_len(s$reps), \(j) {
      rep_pt <- proc.time()
      cat("    rep:", j, "/", s$reps)
      point_ts <- sample_point_ts(s, def)
      pd_ts <- map(point_ts, pd, filtration=s$filt)

      required_dists %>%
        rowwise() %>%
        mutate(matrix = list(dist_mat_wasserstein(pd_ts, power, dim))) ->
        dist_mats_tb
      dist_mats <- dist_mats_tb$matrix
      names(dist_mats) <- dist_mats_tb$d_name

      prepro_tss <- lapply(
        prepros,
        \(prep) def$prepros[[prep]](point_ts, pd_ts, dist_mats))
      names(prepro_tss) <- prepros

      postpro_tss <- list()
      estimations <- double()
      for (l in seq_len(nrow(estimators))) {
        e <- (estimators %>% filter(e_id == l))[1,]
        if (e$e_name == "LT_D1") browser() # something produces NA
        postpro_tss[[l]] <- def$postpros[[e$postpro]](prepro_tss[[e$prepro]])
        estimations[[l]] <- def$cpds[[e$cpd]](prepro_tss[[l]])
      }

      cat(", duration:", (proc.time()-rep_pt)[3], "s\n")
      tb <- tibble(e_name=estimators$e_name, rep_nr=j, estimation=estimations)
      tb$postpro_ts <- postpro_tss
      tb$len <- map_int(postpro_tss, length)
      tb
    }) -> r
    r %>%
      bind_rows() %>%
      mutate(s_id = i, .before=1) ->
      res
    cat("  sampler:", i, "/", nrow(samplers), "end, duration:", (proc.time()-samp_pt)[3], "s\n")
    res
  }) %>%
    bind_rows() ->
    results
  cat("sim: end, duration:", (proc.time()-sim_pt)[3], "s\n")
  results
}
