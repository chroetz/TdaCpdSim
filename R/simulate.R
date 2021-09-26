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

sample_point_ts <- function(sampler, distris) {
  switch(
    sampler$kind,
    birth_death = sample_birth_death(
      distris[[sampler$distri1]],
      distris[[sampler$distri2]],
      sampler$n1, sampler$n2,
      sampler$m,
      sampler$rate),
    iid = sample_iid(
      distris[[sampler$distri1]],
      sampler$m1,
      sampler$n1,
      distris[[sampler$distri2]],
      sampler$m2,
      sampler$n2),
    stop("unknown kind of sampler ", sampler$kind))
}

#' @export
simulate <- function(samplers, estimators, required_dists) {
  sim_pt <- proc.time()

  cat("prepare simulation\n")

  distri_names <- unique(c(samplers$distri1, samplers$distri2))
  distributions <- create_img_sampler(distri_names)
  unique_prepro_names <- unique(estimators$prepro)

  cat("sim: start\n")

  lapply(seq_along(samplers$s_name), \(i) {
    sn <- samplers$s_name[i]
    s <- (samplers %>% filter(s_name == sn))[1,]
    samp_pt <- proc.time()
    cat("  sampler:", i, "/", nrow(samplers), "(",sn,") start\n")
    lapply(seq_len(s$reps), \(j) {
      rep_pt <- proc.time()
      cat("    rep:", j, "/", s$reps)
      point_ts <- sample_point_ts(s, distributions)
      pd_ts <- map(point_ts, pd, filtration=s$filt)

      required_dists %>%
        rowwise() %>%
        mutate(matrix = list(dist_mat_wasserstein(pd_ts, power, dim))) ->
        dist_mats_tb
      dist_mats <- dist_mats_tb$matrix
      names(dist_mats) <- dist_mats_tb$d_name

      prepro_tss <- lapply(
        unique_prepro_names,
        \(pn) prepros[[pn]](point_ts, pd_ts, dist_mats))
      names(prepro_tss) <- unique_prepro_names

      postpro_tss <- list()
      estimations <- double()
      for (l in seq_along(estimators$e_name)) {
        en <- estimators$e_name[l]
        e <- (estimators %>% filter(e_name == en))[1,]
        postpro_tss[[l]] <- postpros[[e$postpro]](prepro_tss[[e$prepro]])
        estimations[[l]] <- cpds[[e$cpd]](postpro_tss[[l]])
      }

      cat(", duration:", (proc.time()-rep_pt)[3], "s\n")
      tb <- tibble(e_name=estimators$e_name, rep_nr=j, estimation=estimations)
      tb$postpro_ts <- postpro_tss
      tb$len <- map_int(postpro_tss, length)
      tb
    }) -> r
    r %>%
      bind_rows() %>%
      mutate(s_name = sn, .before=1) ->
      res
    cat("  sampler:", i, "/", nrow(samplers), "end, duration:", (proc.time()-samp_pt)[3], "s\n")
    res
  }) %>%
    bind_rows() ->
    results
  cat("sim: end, duration:", (proc.time()-sim_pt)[3], "s\n")
  results
}

simulate_parallel <- function(samplers, estimators, required_dists, n_cores=parallel::detectCores()-1) {
  sim_pt <- proc.time()

  cat("prepare simulation\n")

  distri_names <- unique(c(samplers$distri1, samplers$distri2))
  distributions <- create_img_sampler(distri_names)
  unique_prepro_names <- unique(estimators$prepro)

  clstr <- parallel::makeCluster(n_cores, type = "PSOCK")
  parallel::clusterExport(
    clstr,
    c("samplers", "estimators", "required_dists", "distributions"),
    envir=rlang::current_env())
  parallel::clusterExport(
    clstr,
    rlang::env_names(rlang::env_parent()),
    envir=rlang::env_parent())

  cat("sim: start\n")

  lapply(seq_along(samplers$s_name), \(i) {
    sn <- samplers$s_name[i]
    s <- (samplers %>% filter(s_name == sn))[1,]
    samp_pt <- proc.time()
    cat("  sampler:", i, "/", nrow(samplers), "(",sn,") start\n")
    parallel::parLapply(clstr, seq_len(s$reps), s=s, fun=\(j, s) {
      point_ts <- sample_point_ts(s, distributions)
      pd_ts <- purrr::map(point_ts, pd, filtration=s$filt)

      dist_mats_tb <- dplyr::mutate(
        dplyr::rowwise(required_dists),
        matrix = list(dist_mat_wasserstein(pd_ts, power, dim)))
      dist_mats <- dist_mats_tb$matrix
      names(dist_mats) <- dist_mats_tb$d_name

      prepro_tss <- lapply(
        unique_prepro_names,
        \(pn) prepros[[pn]](point_ts, pd_ts, dist_mats))
      names(prepro_tss) <- unique_prepro_names

      postpro_tss <- list()
      estimations <- double()
      for (l in seq_along(estimators$e_name)) {
        en <- estimators$e_name[l]
        e <- (dplyr::filter(estimators, e_name == en))[1,]
        postpro_tss[[l]] <- postpros[[e$postpro]](prepro_tss[[e$prepro]])
        estimations[[l]] <- cpds[[e$cpd]](postpro_tss[[l]])
      }

      tb <- tibble::tibble(e_name=estimators$e_name, rep_nr=j, estimation=estimations)
      tb$postpro_ts <- postpro_tss
      tb$len <- purrr::map_int(postpro_tss, length)
      tb
    }) -> r
    r %>%
      bind_rows() %>%
      mutate(s_name = sn, .before=1) ->
      res
    cat("  sampler:", i, "/", nrow(samplers), "end, duration:", (proc.time()-samp_pt)[3], "s\n")
    res
  }) %>%
    bind_rows() ->
    results
  parallel::stopCluster(clstr)
  cat("sim: end, duration:", (proc.time()-sim_pt)[3], "s\n")
  results
}


