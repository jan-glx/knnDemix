
context("mixture.test")
if (!exists(".Random.seed")) set.seed(NULL)
old_seed <- .Random.seed
set.seed(1)

n0 <- 500 # 100 sample points from \eqn{F_0}
nm <- 500 # 100 sample points from \eqn{F_0}
n1 <- 250  # of which 50 are from \eqn{F_1}
p <- 2    # two dimensions
X0 <- matrix(rnorm(n0*p), ncol=p) # \eqn{F_0} is standard normal, \eqn{F_1} has mean c(10,10):
Xm <- matrix(c(rnorm((nm-n1)*p), rnorm(n1*p)+10), ncol=p, byrow = TRUE)

test_that("mixture.test runs smooth and gives reasonable results", {
  expect_gt(
    mixture.test(X0, Xm, alpha=0.5, calc.CI=FALSE)$p.value, #should not be significant
    0.05)
  expect_lte(
    mixture.test(X0, Xm, alpha=1, calc.CI=FALSE)$p.value, #should be significant
    0.05)
  expect_lt(
    mixture.test(X0, Xm)$conf.int[2],
    1) #should exclude 1
})

n0 <- 5000 # 100 sample points from \eqn{F_0}
nm <- 500 # 100 sample points from \eqn{F_0}
n1 <- 250  # of which 50 are from \eqn{F_1}
p <- 2    # two dimensions
X0 <- matrix(rnorm(n0*p), ncol=p) # \eqn{F_0} is standard normal, \eqn{F_1} has mean c(10,10):
Xm <- matrix(c(rnorm((nm-n1)*p), rnorm(n1*p)+10), ncol=p, byrow = TRUE)

test_that("mixture.test for scewed sample numbers", {
  expect_gt(
    mixture.test(X0, Xm, alpha=0.5, calc.CI=FALSE)$p.value, #should not be significant
    0.05)
  expect_lte(
    mixture.test(X0, Xm, alpha=1, calc.CI=FALSE)$p.value, #should be significant
    0.05)
  expect_lt(
    mixture.test(X0, Xm)$conf.int[2],
    1) #should exclude 1
})


n0 <- 5000 # 100 sample points from \eqn{F_0}
nm <- 500 # 100 sample points from \eqn{F_0}
n1 <- 250  # of which 50 are from \eqn{F_1}
X0 <- rnorm(n0) # \eqn{F_0} is standard normal, \eqn{F_1} has mean 10:
Xm <- c(rnorm((nm-n1)), rnorm(n1)+10)

test_that("mixture.test with vectors", {
  expect_gt(
    mixture.test(X0, Xm, alpha=0.5, calc.CI=FALSE)$p.value, #should not be significant
    0.05)
  expect_lte(
    mixture.test(X0, Xm, alpha=1, calc.CI=FALSE)$p.value, #should be significant
    0.05)
  expect_lt(
    mixture.test(X0, Xm)$conf.int[2],
    1) #should exclude 1
})

.Random.seed <- old_seed
