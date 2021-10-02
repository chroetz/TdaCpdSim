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
      sampler$n1, sampler$n-sampler$n1,
      sampler$m,
      sampler$rate),
    iid = sample_iid(
      distris[[sampler$distri1]],
      sampler$m1,
      sampler$n1,
      distris[[sampler$distri2]],
      sampler$m2,
      sampler$n-sampler$n1),
    stop("unknown kind of sampler ", sampler$kind))
}

run_experiment <- function(j, s, seed, verbose,
                           distributions, required_prepros, estimators) {
  ppdiff <- setdiff(required_prepros, names(prepros))
  if (length(ppdiff) > 0) stop("Unknow Pre-Processors ", paste(ppdiff,collapse=","))

  if (verbose) {
    rep_pt <- proc.time()
    cat("    rep:", j, "/", s$reps)
  }
  set.seed(seed)
  point_ts <- sample_point_ts(s, distributions)
  pd_ts <- purrr::map(point_ts, pd, filtration=s$filt)

  if (nrow(required_dists) > 0) {
    dist_mats_tb <- dplyr::mutate(
      dplyr::rowwise(required_dists),
      matrix = list(dist_mat_wasserstein(pd_ts, power, dim)))
    dist_mats <- dist_mats_tb$matrix
    names(dist_mats) <- dist_mats_tb$d_name
  } else {
    dist_mats_tb <- list()
  }

  prepro_tss <- lapply(
    required_prepros,
    \(pn) prepros[[pn]](point_ts, pd_ts, dist_mats))
  names(prepro_tss) <- required_prepros

  postpro_tss <- list()
  estimations <- double()
  for (l in seq_along(estimators$e_name)) {
    en <- estimators$e_name[l]
    e <- (estimators %>% filter(e_name == en))[1,]
    ts <- prepro_tss[[e$prepro]]
    if (isTRUE(nchar(e$select) > 0))
      ts <- ts[, eval(str2lang(e$select)), drop=FALSE]
    ts <- ts[(1+e$trim):(nrow(ts)-e$trim), , drop=FALSE]
    postpro_tss[[l]] <- postpros[[e$postpro]](ts)
    estimations[[l]] <- detectors[[e$detector]](postpro_tss[[l]], s$n)
  }

  tb <- tibble::tibble(e_name=estimators$e_name, rep_nr=j, estimation=estimations)
  tb$postpro_ts <- postpro_tss
  tb$len <- purrr::map_int(postpro_tss, length)

  if (verbose) cat(", duration:", (proc.time()-rep_pt)[3], "s\n")

  tb
}

#' Run a Simulation of a Time Series
#'
#' @param samplers A tibble. Each row describes one distribution for a time
#'   series of sets of points.
#' @param estimators A tibble. Each row describes an estimator of the change
#'   point in the time series.
#' @param required_dists A tibble. Each row describes a distance of persistence
#'   diagrams.
#' @param n_cores An integer. The number of cores to use. If greater 1, than R
#'   package parallel is used to parallelize repetitions of experiments.
#' @return A tibble. Each row describes the results of applying one estimator to
#'   a time series created by one of the samplers.
#' @export
simulate_ts <- function(samplers, estimators, required_dists, n_cores=1) {
  sim_pt <- proc.time()

  cat("prepare simulation\n")

  distri_names <- unique(c(samplers$distri1, samplers$distri2))
  distributions <- create_distris(distri_names)
  required_prepros <- unique(estimators$prepro)

  if (n_cores > 1) {
    clstr <- parallel::makeCluster(n_cores, type = "PSOCK")
    parallel::clusterExport(
      clstr,
      c("samplers", "estimators", "required_dists"),
      envir=rlang::current_env())
    parallel::clusterExport(
      clstr,
      rlang::env_names(rlang::env_parent()),
      envir=rlang::env_parent())
  }

  cat("sim: start\n")

  # creating seeds to make parallel simulation reproducible
  seeds <- sample.int(2^31-1, size=nrow(samplers))

  lapply(seq_along(samplers$s_name), \(i) {
    set.seed(seeds[i])
    sn <- samplers$s_name[i]
    s <- (samplers %>% filter(s_name == sn))[1,]

    # creating seeds to make parallel simulation reproducible
    seedss <- sample.int(2^31-1, size=s$reps)

    samp_pt <- proc.time()
    cat("  sampler:", i, "/", nrow(samplers), "start,", sn, "\n")
    if (n_cores > 1) {
      r <- parallel::parLapply(
        clstr,
        seq_len(s$reps),
        fun=\(j) run_experiment(j, s, seedss[j], verbose=FALSE,
                                distributions=distributions,
                                required_prepros=required_prepros,
                                estimators=estimators))
    } else {
      r <- lapply(
        seq_len(s$reps),
        FUN=\(j) run_experiment(j, s, seedss[j], verbose=TRUE,
                                distributions=distributions,
                                required_prepros=required_prepros,
                                estimators=estimators))
    }
    r %>%
      bind_rows() %>%
      mutate(s_name = sn, .before=1) ->
      res
    cat("  sampler:", i, "/", nrow(samplers), "end, duration:", (proc.time()-samp_pt)[3], "s\n")
    res
  }) %>%
    bind_rows() ->
    results
  if (n_cores > 1) parallel::stopCluster(clstr)
  cat("sim: end, duration:", (proc.time()-sim_pt)[3], "s\n")
  results
}


