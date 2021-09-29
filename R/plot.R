# TODO: true length n as argument, visualize

plot_reps <- function(X, n, p=0.5) {
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

plot_max_reps <- function(Z, n, true, title="", p=0.5) {
  S <- sapply(Z, \(x) as.matrix(x))
  S <- S[!apply(S, 1, \(x) any(is.na(x))),]
  I <- apply(S, 2, \(s) which.max(s))
  plot_reps(S, n=n, p=p)
  abline(v=mean(I), col=3, lwd=2)
  abline(v=quantile(I, probs=(1-p)/2), col=3)
  abline(v=quantile(I, probs=(1+p)/2), col=3)
  abline(v=true, col=4)
  pad <- (n-nrow(S))/2
  abline(v=pad, lty="dashed")
  abline(v=n-pad, lty="dashed")
  abline(v=1)
  abline(v=n)
  title(title)
}

#' @export
plot_result <- function(results, samplers, sn, en, p=0.5) {
  plot_max_reps(
    results %>% filter(s_name==sn, e_name==en) %>% pluck("postpro_ts"),
    n = samplers %>% filter(s_name==sn) %>% pluck("n"),
    true = samplers %>% filter(s_name==sn) %>% pluck("n1") + 0.5,
    title = paste0("Sampler ", sn, "    Estimator ", en),
    p = p)
}
