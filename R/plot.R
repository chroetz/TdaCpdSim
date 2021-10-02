plot_ts_summary <- function(X, n, p=0.5) {
  qs <- apply(X, 1, \(x) quantile(x, probs=c((1-p)/2, (1+p)/2), na.rm=TRUE))
  means <- rowMeans(X)
  k <- length(means)
  x <- (n-k)/2+(1:k)
  plot(NA, xlim=c(1,n), ylim=range(qs), frame.plot=FALSE, xlab="t", ylab="post-processed")
  grid()
  lines(x, means, col=2, lwd=2)
  lines(x, qs[1,], col=2)
  lines(x, qs[2,], col=2)
}

plot_ts_all <- function(X, n, colors) {
  k <- nrow(X)
  x <- (n-k)/2+(1:k)
  plot(NA, xlim=c(1,n), ylim=range(X), frame.plot=FALSE, xlab="t", ylab="post-processed")
  grid()
  for (i in 1:ncol(X)) lines(x, X[,i], col=colors[i])
}

plot_reps_summary <- function(postpro_ts, estimation, n, true, title="", p=0.5) {
  S <- sapply(postpro_ts, \(x) as.matrix(x))
  S <- S[!apply(S, 1, \(x) any(is.na(x))),]
  plot_ts_summary(S, n=n, p=p)
  abline(v=mean(estimation), col=3, lwd=2)
  abline(v=quantile(estimation, probs=(1-p)/2), col=3)
  abline(v=quantile(estimation, probs=(1+p)/2), col=3)
  abline(v=true, col=4)
  pad <- (n-nrow(S))/2
  abline(v=pad, lty="dashed")
  abline(v=n-pad, lty="dashed")
  abline(v=1)
  abline(v=n)
  title(title)
}


plot_reps_all <- function(postpro_ts, estimation, n, true, title="") {
  S <- sapply(postpro_ts, \(x) as.matrix(x))
  S <- S[!apply(S, 1, \(x) any(is.na(x))),]
  colors <- rainbow(ncol(S), alpha=0.3)
  plot_ts_all(S, n=n, colors=colors)
  # for (i in 1:ncol(S)) abline(v=estimation[i], col=colors[i])
  # abline(v=true, col=4)
  pad <- (n-nrow(S))/2
  abline(v=pad, lty="dashed")
  abline(v=n-pad, lty="dashed")
  abline(v=1)
  abline(v=n)
  title(title)
}

#' @export
plot_result_summary <- function(results, samplers, sn, en, p=0.5) {
  r <- results %>% filter(s_name==sn, e_name==en)
  s <- samplers %>% filter(s_name==sn)
  plot_reps_summary(
    postpro_ts = r$postpro_ts,
    estimation = r$estimation,
    n = s$n,
    true = s$n1 + 0.5,
    title = paste0("Sampler ", sn, "    Estimator ", en),
    p = p)
}

#' @export
plot_result_all <- function(results, samplers, sn, en) {
  r <- results %>% filter(s_name==sn, e_name==en)
  s <- samplers %>% filter(s_name==sn)
  plot_reps_all(
    postpro_ts = r$postpro_ts,
    estimation = r$estimation,
    n = s$n,
    true = s$n1 + 0.5,
    title = paste0("Sampler ", sn, "    Estimator ", en))
}
