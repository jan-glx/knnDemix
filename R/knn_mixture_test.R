#' mixture.test
#'
#' @description Tests against the true mixing coefficient of \eqn{F_0} in the mixture \eqn{F_m} to be smaller than \code{alpha}
#' using a knn statistic calculated from one sample from F0 and one sample from the mixture \eqn{F_m}.
#' @details \eqn{F_m=\alpha F_0 + (1-\alpha) F_1} thus \eqn{f_m=\alpha f_0 + (1-\alpha) f_1},
#' \eqn{f_1 >=0} thus \eqn{f_m >= \alpha f_0} ... TODO: finish explaination
#' @import stats
#' @param X0 numeric matrix with rows beeing data points from a sample from \eqn{F_0}. Vector inputs are coerced to nx1 matrices.
#' @param Xm numeric matrix with rows beeing data points from a sample from \eqn{F_m}. Vector inputs are coerced to nx1 matrices.
#' @param alpha numeric value of the fraction of \eqn{F_0} in \eqn{F_m} that is tested against
#' @param alternative a character string specifying the alternative hypothesis - only 'lower' is allowed and only lower makes sense
#' @param conf.level confidence level of the interval.
#' @param k integer value indicating the which nearest neighbour is used for the testing density ratios.
#' @param sub_sample numeric value specifing factor subsamplamping -- only n*sub_sample_exp of all n pvalues will be used in order to reduce their depency
#' @param seed_ seed to be used for subsampling. Set to \code{NULL} if want reduce variance by averaging several subsamplings.
#' @param calc.CI logical value indicating if confidence intervals for alpha schould be computed.
#' Can be turend of to save computation time in case it is not needed.
#' @examples
#' n0 <- 100 # 100 sample points from \eqn{F_0}
#' nm <- 100 # 100 sample points from \eqn{F_0}
#' n1 <- 50  # of which 50 are from \eqn{F_1}
#' p <- 2    # two dimensions
#' X0 <- matrix(rnorm(n0*p), ncol=p) # \eqn{F_0} is standard normal, \eqn{F_1} has mean c(10,10):
#' Xm <- matrix(c(rnorm((nm-n1)*p), rnorm(n1*p)+10), ncol=p, byrow = TRUE)
#'
#' mixture.test(X0, Xm, alpha=0.5, calc.CI=FALSE)$p.value #should not be significant
#' mixture.test(X0, Xm, alpha=1, calc.CI=FALSE)$p.value #should be significant
#' mixture.test(X0, Xm) # computes confidence intervals
#' @export
mixture.test <- function(X0, Xm, alpha = 1, alternative = "lower", conf.level = 0.95, k = 1, sub_sample = 1/4, seed_ = 1, calc.CI=TRUE) {
  alternative <- match.arg(alternative)
  if (!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) || alpha < 0 || alpha > 1))
    stop("'alpha' must be a single number between 0 and 1")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!missing(sub_sample) && (length(sub_sample) != 1 || !is.finite(sub_sample) || sub_sample < 0 || sub_sample > 1))
    stop("'sub_sample_exp' must be a single number between 0 and 1")
  dname <- paste(deparse(substitute(X0)), "and", deparse(substitute(Xm)))
  X0 <- as.matrix(X0)
  Xm <- as.matrix(Xm)
  n0 <- nrow(X0)
  nm <- nrow(Xm)
  p <- ncol(X0)
  if (ncol(Xm)!=p)
    stop("X0 and Xm must have the same dimesions")

  n <- ceiling(n0*sub_sample)
  if (!exists(".Random.seed")) set.seed(NULL)
  old_seed <- .Random.seed
  set.seed(seed_)
  ss <- sample(seq_len(n0), n) # subsampling bc. p.vals are not independent
  .Random.seed <- old_seed

  dnn0 <- FNN::knn.dist(X0,k)[ss,k]
  dnnm <- FNN::knnx.dist(Xm,X0[ss,],k)[,k]

  r <- (k/(dnnm^p*nm))/(k/(dnn0^p*(n0-1)))
  p.vals <- extraDistr::pbetapr(alpha/r,k,k,lower.tail = FALSE) # == 1-1/(r/alpha+1) for k=1
  S <- -2*sum(log(p.vals))

  pval <- pchisq(S, lower.tail = FALSE, df = 2*n)

  S_upper <- qchisq(conf.level, df=2*n)
  cint <-  if (!calc.CI) c(0,1) else c(0,  pmax(optimize(function(alpha) {
    p.vals <- extraDistr::pbetapr(alpha/r,k,k,lower.tail = FALSE)
    S <- -2*sum(log(p.vals))
    abs(S-S_upper)
  },c(0,1))$minimum,0) ) #there should be a smarter way...
  attr(cint, "conf.level") <- conf.level

  S_med <- qchisq(0.5, df=2*n)
  estimate <- if (!calc.CI) 0 else min(1, pmax(optimize(function(alpha) {
    p.vals <- extraDistr::pbetapr(alpha/r,k,k,lower.tail = FALSE)
    S <- -2*sum(log(p.vals))
    abs(S-S_med)
  },c(0,1))$minimum,0) )
  names(estimate) <- "alpha_0"

  rval <- list(statistic = S, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = alpha,
               alternative = alternative, method = "knn mixture test", data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
