plot_reps <- function(X, p=0.5) {
  qs <- apply(X, 1, \(x) quantile(x, probs=c((1-p)/2, (1+p)/2), na.rm=TRUE))
  means <- rowMeans(X)
  plot(means, ylim=range(qs), type="l", col=2, lwd=2)
  lines(qs[1,], col=2)
  lines(qs[2,], col=2)
}

plot_cusum_reps <- function(Z, title="", p=0.5) {
  S <- sapply(Z, cusum_ts)
  I <- apply(S, 2, \(s) which.max(s))
  plot_reps(S, p=p)
  abline(v=mean(I), col=3, lwd=2)
  abline(v=quantile(I, probs=(1-p)/2), col=3)
  abline(v=quantile(I, probs=(1+p)/2), col=3)
  title(title)
}

plot_max_reps <- function(Z, title="", p=0.5) {
  S <- sapply(Z, \(x) as.matrix(x))
  S <- S[!apply(S,1,\(x) any(is.na(x))),]
  I <- apply(S, 2, \(s) which.max(s))
  plot_reps(S, p=p)
  abline(v=mean(I), col=3, lwd=2)
  abline(v=quantile(I, probs=(1-p)/2), col=3)
  abline(v=quantile(I, probs=(1+p)/2), col=3)
  title(title)
}

plot_cusum_result <- function(results, sn, en) {
  plot_cusum_reps(
    results %>% filter(s_name==sn, e_name==en) %>% pluck("postpro_ts"),
    title = paste0("sampler ", sn, ", estimator ", en))
}

plot_max_result <- function(results, sn, en) {
  plot_max_reps(
    results %>% filter(s_name==sn, e_name==en) %>% pluck("postpro_ts"),
    title = paste0("sampler ", sn, ", estimator ", en))
}
