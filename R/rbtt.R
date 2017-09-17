library(parallel)
library(data.table)

get.std.err <- function(x, y, var.equal = FALSE)
{
  n1 <- length(x)
  n2 <- length(y)

  if(var.equal)
  {
    s_p <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    return(s_p * sqrt(1/n1 + 1/n2))
  }
  else
  {
    s_d = sqrt(var(x)/n1 + var(y)/n2)
    return(s_d)
  }
}

get.t.stat <- function(x, y, var.equal = FALSE)
{
  if (var(x) == 0 && var(y) == 0)
  {
    stop("Data are essentially constant")
  }

  std.err <- get.std.err(x, y, var.equal)
  return((mean(x) - mean(y)) / std.err)
}

get.ttest.output <- function(x, y, method, t, conf.int, p.val, conf.level, x.name = NULL, y.name = NULL)
{
  if (is.null(x.name) || is.null(y.name))
  {
    dname <- paste(deparse(substitute(x)),"and",
                   deparse(substitute(y)))
  }
  else
  {
    dname <- paste(x.name, "and", y.name)
  }

  mx <- mean(x)
  my <- mean(y)
  mu <- 0
  names(mu) <- "difference in means"
  estimate <- c(mx, my)
  names(estimate) <- c("mean of x","mean of y")
  if (method == 1) method.str = "Two-sided robust bootstrapped t-test assuming equal variance"
  else method.str = "One-sided robust bootstrapped t-test not assuming equal variance"

  if (method == 1)
  {
    alternative <- "two.sided"
  }
  else if (method == 2)
  {
    if (t > 0)
    {
      alternative <- "greater"
    }
    else
    {
      alternative <- "less"
    }
  }

  names(t) <- "t"
  attr(conf.int, "conf.level") <- conf.level

  rval <- list(statistic = t, p.value = p.val,
               conf.int = conf.int, estimate = estimate, null.value = mu,
               alternative = alternative,
               method = method.str, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}

rbtt.1 <- function(x, y, n.boot, n.cores = 1, conf.level)
{
  arguments <- as.list(match.call())

  x.star <- x / mean(x)
  y.star <- y / mean(y)

  boot.t.vals.list <- mclapply(seq(1:n.boot), function(i)
  {
    boot.x <- sample(x.star, size = length(x.star), replace = TRUE)
    boot.y <- sample(y.star, size = length(y.star), replace = TRUE)

    t.boot <- get.t.stat(boot.x, boot.y, var.equal = TRUE)
    return(t.boot)
  }, mc.cores = n.cores)

  boot.t.vals <- unlist(boot.t.vals.list)

  t <- get.t.stat(x, y, var.equal = TRUE)
  p.val <- sum(abs(boot.t.vals) > abs(t)) / n.boot
  conf.int <- quantile(boot.t.vals, probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2))
  conf.int <- mean(x) - mean(y) + conf.int * get.std.err(x, y, var.equal = TRUE)

  return(get.ttest.output(x = x, y = y, method = 1, t = t, conf.int = conf.int,
                          p.val = p.val, conf.level = conf.level))
}

rbtt.2 <- function(x, y, n.boot, n.cores = 1, conf.level)
{
  arguments <- as.list(match.call())

  x.star <- x / mean(x)
  y.star <- y / mean(y)

  x.y.star <- c(x.star, y.star)

  boot.t.vals.list <- mclapply(seq(1:n.boot), function(i)
  {
    boot.x <- sample(x.y.star, size = length(x.star), replace = TRUE)
    boot.y <- sample(x.y.star, size = length(y.star), replace = TRUE)

    t.boot <- get.t.stat(boot.x, boot.y, var.equal = FALSE)
    return(t.boot)
  }, mc.cores = n.cores)

  boot.t.vals <- unlist(boot.t.vals.list)

  t <- get.t.stat(x, y, var.equal = FALSE)

  if (t > 0)
  {
    p.val <- 2 * (sum(boot.t.vals > t) / n.boot)
    q <- quantile(boot.t.vals, probs = (1 - (1 - conf.level) / 2))
    conf.int <- c(-q, q)
  }
  else
  {
    p.val <- 2 * (sum(boot.t.vals < t) / n.boot)
    q <- quantile(boot.t.vals, probs = c((1 - conf.level) / 2))
    conf.int <- c(q, -q)
  }

  conf.int <- mean(x) - mean(y) + conf.int * get.std.err(x, y, var.equal = FALSE)
  return(get.ttest.output(x = x, y = y, method = 2, t = t, conf.int = conf.int,
                          p.val = p.val, conf.level = conf.level))
}

rbtt.combined <- function(x, y, n.boot, n.cores = 1, conf.level = 0.95, x.name, y.name)
{
  arguments <- as.list(match.call())

  x.star <- x / mean(x)
  y.star <- y / mean(y)

  x.y.star <- c(x.star, y.star)

  boot.t.vals.list <- mclapply(seq(1:n.boot), function(i)
  {
    boot.x.1 <- sample(x.star, size = length(x.star), replace = TRUE)
    boot.y.1 <- sample(y.star, size = length(y.star), replace = TRUE)
    boot.x.2 <- sample(x.y.star, size = length(x.star), replace = TRUE)
    boot.y.2 <- sample(x.y.star, size = length(y.star), replace = TRUE)

    t.boot.1 <- get.t.stat(boot.x.1, boot.y.1, var.equal = TRUE)
    t.boot.2 <- get.t.stat(boot.x.2, boot.y.2, var.equal = FALSE)

    return(c(t.boot.1, t.boot.2))
  }, mc.cores = n.cores)

  boot.t.vals <- do.call(rbind, boot.t.vals.list)

  t.1 <- get.t.stat(x, y, var.equal = TRUE)
  t.2 <- get.t.stat(x, y, var.equal = FALSE)

  p.val.1 <- sum(abs(boot.t.vals[,1]) > abs(t.1)) / n.boot
  conf.int.1 <- quantile(boot.t.vals, probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2))
  conf.int.1 <- mean(x) - mean(y) + conf.int.1 * get.std.err(x, y, var.equal = TRUE)

  if (t.2 > 0)
  {
    p.val.2 <- 2 * (sum(boot.t.vals[,2] > t.2) / n.boot)
    q.2 <- quantile(boot.t.vals, probs = (1 - (1 - conf.level) / 2))
    conf.int.2 <- c(-q.2, q.2)
  }
  else
  {
    p.val.2 <- 2 * (sum(boot.t.vals[,2] < t.2) / n.boot)
    q.2 <- quantile(boot.t.vals, probs = c((1 - conf.level) / 2))
    conf.int.2 <- c(q.2, -q.2)
  }

  conf.int.2 <- mean(x) - mean(y) + conf.int.2 * get.std.err(x, y, var.equal = FALSE)

  output.1 <- get.ttest.output(x = x, y = y, method = 1, t = t.1, conf.int = conf.int.1, p.val = p.val.1,
                               conf.level = conf.level, x.name = x.name, y.name = y.name)
  output.2 <- get.ttest.output(x = x, y = y, method = 2, t = t.2, conf.int = conf.int.2, p.val = p.val.2,
                               conf.level = conf.level, x.name = x.name, y.name = y.name)

  return(list(output.1, output.2))
}

#' Perform robust bootstrapped t-tests
#'
#' Perform robust bootstrapped two-sample t-tests that aim to better control type-I error rates
#' when comparing means of non-negative distributions with excess zero observations.
#'
#' @importFrom parallel mclapply
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom data.table setnames
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param n.boot number of bootstrap resamples to perform
#' @param n.cores number of cores to use for parallelization. Defaults to 1. If using Windows, set n.cores = 1.
#' @param method Which robust bootstrapped t-test to perform. Set `method=1’ for a two-sample t-test under the equal variance assumption, 'method = 2' for a two-sample t-test without the equal variance assumption, and 'method = "both"' to perform both methods simultaneously.
#' @param conf.level Desired confidence level for computing confidence intervals: a number between 0 and 1.
#' @return A list (or two lists in the case of method = "combined") containing the following components:\cr
#' \item{statistic}{the value of the t-statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a bootstrap-based confidence interval for the difference in means.}
#' \item{estimate}{the estimated difference in means.}
#' \item{null.value}{the hypothesized value of the mean difference, zero.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string describing the type of two-sample bootstrapped t-test used}
#' \item{data.name}{a character string giving the names of the data}
#' @export
#' @examples
#' x=rbinom(50,1,0.5)*rlnorm(50,0,1)
#' y=rbinom(150,1,0.3)*rlnorm(150,2,1)
#'
#' rbtt(x, y, n.boot=999)
#'
#' # Perform bootstrap resamples on 2 cores
#' rbtt(x, y, n.boot=999, n.cores=2)
#'
#' # Use methods 1 or 2 individually
#' rbtt(x, y, n.boot = 999, method = 1)
#' rbtt(x, y, n.boot = 999, method = 2)
#'
#' # Use a confidence level of 0.99
#' rbtt(x, y, n.boot = 999, conf.level = 0.99)
rbtt <- function(x, y, n.boot, n.cores = 1, method = "combined", conf.level = 0.95)
{
  if (Sys.info()['sysname'] == "Windows" && n.cores > 1)
  {
    warning("Multi-core processing is not supported on Windows. Setting n.cores to 1.")
    n.cores <- 1
  }

  if (method == "combined")
  {
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    output <- rbtt.combined(x, y, n.boot = n.boot, n.cores = n.cores,
                            conf.level = conf.level, x.name = x.name, y.name = y.name)
  }
  else if (method == 1 || method == 2)
  {
    if (method == 1)
    {
      output <- rbtt.1(x, y, n.boot = n.boot, n.cores = n.cores, conf.level = conf.level)
    }
    else
    {
      output <- rbtt.2(x, y, n.boot = n.boot, n.cores = n.cores, conf.level = conf.level)
    }

    dname <- paste(deparse(substitute(x)),"and",
                   deparse(substitute(y)))

    output$data.name <- dname
  }
  else
  {
    stop("Invalid method specification.\n
         Use \"method = 1\" or \"method = 2\" or \"method = \'combined\'\"")
  }

  return(output)
}

