##############################
## Test functions for ARS #
#############################
## Last edit: Dec 15 by Shuyao #

library(testthat)
library(ars)

# Main function

test_that("Check if the user's input includes a valid function", {

  # when g is not a function.
  expect_error(ars(1,100))
  expect_error(ars(x^3, 100, -1, 1))

})

test_that("Check if the user's input of bound is valid or not", {

  # In order to see the difference, I put all normal function.
  # n is not numeric.
  expect_error(ars(dnorm, a, -1, 10))

  # Not both values are numeric.
  expect_error(ars(dnorm, 100, a, 4))

  # Lower and upper are equal.
  expect_error(ars(dnorm, 100, -1, -1))

  # Function with good bounds.
  expect_length(ars(dnorm, 100, -1, 10), 100)

  # Do not specify bounds.
  expect_warning(ars(dnorm, 100, upper=1))

  # Lower are greater than Upper.
  expect_warning(ars(dnorm, 100))

})

test_that("Check if the function is concave or not", {

  # Function that is concave
  log.norm <- function(x) log(dnorm(x))
  expect_true(isConcave(log.norm, 0, 10))

  # Function that is not concave
  log.cauchy <- function(x) log(dcauchy(x))
  expect_false(isConcave(log.cauchy, 0, 10))
})

test_that("check for log-concavity of final results", {
  # Function is not concave.
  expect_error(ars(function(x) exp(x), 100, -10, 10))
  expect_error(ars(dcauchy, 100, -1, 10))

  # Function is concave.
  expect_length(ars(dnorm, 100, -1, 10), 100)
})

test_that("check if initializing the abscissae function successfully", {
  # Result should be in form of numeric.
  expect_is(init_abscissae(dnorm, -1, 10), 'numeric')

  # Result should be length of 3.
  expect_length(init_abscissae(dnorm, -1, 10), 3)

})

test_that("check if utility function z_fun return the coorrect result.", {

  dx <- 1e-8
  g<- dnorm
  h <- function(x) log(g(x))
  x <- init_abscissae(h, -1, 10)
  hx <- h(x)
  hpx <- (h(x + dx) - hx)/dx

  # Result should be in form of numeric.
  expect_is(z_fun(hx, hpx, x, -1, 10), 'numeric')

  # zfindj function return the output between 1 and 3.
  expect_gt(zfindj(5,z_fun(hx, hpx, x, -1, 10)), 0)

})

test_that("check if u_fun return the correct result.", {

  dx <- 1e-8
  g<- dnorm
  h <- function(x) log(g(x))
  x <- init_abscissae(h, -1, 10)
  hx <- h(x)
  hpx <- (h(x + dx) - hx)/dx
  z <- z_fun(hx, hpx, x, -1, 10)

  # Result should be in form of numeric.
  expect_is(u_fun(5, hx, hpx, x, z), 'numeric')

  # Result should return a single value.
  expect_length(u_fun(1, hx, hpx, x, z),1)

})


test_that("check if sinvcdf_fun return the correct result.", {

  dx <- 1e-8
  g<- dnorm
  h <- function(x) log(g(x))
  x <- init_abscissae(h, -1, 10)
  hx <- h(x)
  hpx <- (h(x + dx) - hx)/dx
  z <- z_fun(hx, hpx, x, -1, 10)
  u <- u_fun(5, hx, hpx, x, z)
  T_u <- function(x_val) {u_fun(x_val, hx, hpx, x, z)}

  # Result should be in form of numeric.
  expect_is(sinvcdf_fun(0.5, hx, hpx, x, z, T_u), 'numeric')

  # Result should return a single value.
  expect_length(sinvcdf_fun(0.5, hx, hpx, x, z, T_u), 1)

  # p should be in range of the [0,1].
  expect_error(sinvcdf_fun(2, hx, hpx, x, z, T_u))

})

test_that("check if l_fun return the correct result.", {

  dx <- 1e-8
  g<- dnorm
  h <- function(x) log(g(x))
  x <- init_abscissae(h, -1, 10)
  hx <- h(x)
  l_fun(5, hx, x)

  # Result should be in form of numeric.
  expect_is(l_fun(5, hx, x), 'numeric')

  # Result should return a single value.
  expect_length(l_fun(5, hx, x), 1)

})

test_that("check if the distribution we got finally is close to the original one", {
  expect_equal(ks.test(ars(dnorm, 100), rnorm(100))$p.value > 0.05, TRUE)
})


