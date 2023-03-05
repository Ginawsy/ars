#' @importFrom assertthat assert_that
#' @importFrom numDeriv grad

library(testthat)
if(!exists('z_fun', mode='function')) source('./R/utility.R')

#' Draw a sample from s_k(x) using inverse CDF method and perform rejection test.
#' Returns a list:
#' accept: {TRUE, FALSE} TRUE means x* accepted.
#' x_star: value of the sample drawn (x*).
#' hEval: {TRUE, FALSE} TRUE means x* should be included in T_k+1.

ars_sample <- function(h, Tk) {
  lower <- Tk$z[1]
  upper <- Tk$z[length(Tk$z)]
  p <- runif(1)
  x_star <- Tk$sinvcdf(p)
  if ((x_star < lower) | (x_star > upper)){
    warning('Sampled x* is not inside bounds.')
  }
  w <- runif(1)
  if(w <= exp(Tk$l(x_star) - Tk$u(x_star))) {
    return(list(accept = TRUE, x_star = x_star, hEval = FALSE))
  } else if (w <= exp(h(x_star) - Tk$u(x_star))) {
    return(list(accept = TRUE, x_star = x_star, hEval = TRUE))
  }
  return(list(accept = FALSE, x_star = x_star, hEval = TRUE))
}


#' Helper function to initialize the abscissae of k = 3.

init_abscissae <- function(h, lower, upper) {
  dx <- 1e-8
  if ( (h(lower + dx) - h(lower))/dx == 0) {
    stop("Derivative at lower bound is equal to 0.")
  }
  if ( (h(upper + dx) - h(upper))/dx == 0) {
    stop("Derivative at upper bound is equal to 0.")
  }

  x <- rep(0, 3)

  x[1] <- lower
  x[3] <- upper

  x[2] <- runif(1, min = lower/3*2 + upper/3, max = upper/3*2 + lower/3)
  if ((h(x[2] + dx) - h(x[2]))/dx == 0) {
    x[2] <- runif(1, min = lower/8*7 + upper/8, max = upper/8*7 + lower/8)
  }
  if ((h(x[2] + dx) - h(x[2]))/dx == 0) {
    stop("Derivative is equal to 0 for a large portion of domain.")
  }

  return(x)
}


#' Update T_k to T_k+1 by inserting x* into the right place.
#' Update h(x), h'(x) at the same time.

update_Tk <- function(Tk, h, x_star) {
  k <- length(Tk$x)
  x_new <- rep(0, k + 1)
  hx_new <- rep(0, k + 1)
  hpx_new <- rep(0, k + 1)

  if (x_star <= Tk$x[1]) {

    # If x* is smaller than the smallest value of x_j, insert at beginning.

    x_new <- c(x_star, Tk$x)
    hx_new <- c(h(x_star), Tk$hx)
    hpx_new <- c(grad(h, x_star), Tk$hpx)
  } else if (x_star >= Tk$x[k]) {

    # If x* is larger than the largest value of x_j, append at end.

    x_new <- c(Tk$x, x_star)
    hx_new <- c(Tk$hx, h(x_star))
    hpx_new <- c(Tk$hpx, grad(h, x_star))
  } else {

    # else insert in the middle.

    pos <- 0
    for (i in 1:length(Tk$x)) {
      if (Tk$x[i] > x_star) { pos <- i; break }
    }
    x_new <- c(Tk$x[1:pos - 1], x_star, Tk$x[pos:k])
    hx_new <- c(Tk$hx[1:pos - 1], h(x_star), Tk$hx[pos:k])
    hpx_new <- c(Tk$hpx[1:pos - 1], grad(h, x_star), Tk$hpx[pos:k])
  }

  # Calculate the new z
  z_new <- z_fun(hx_new, hpx_new, x_new, Tk$z[1], Tk$z[k])

  Tk$x <- x_new
  Tk$hx <- hx_new
  Tk$hpx <- hpx_new
  Tk$z <- z_new

  return(Tk)
}


#' The main ars function callable for users.
#' @param g  a probability density that is log concave of class \code{function} to be drawn from
#' @param n an positive integer denoting the sample size
#' @param lower numeric denoting the lower bound of g to draw from
#' @param upper numeric denoting the upper bound of g to draw from
#' @export
ars <- function(g, n, lower = -Inf, upper = Inf) {

  # First, we need to check the input.
  # Check for the bounds:
  # If lower and upper are valid.

  assert_that((is.numeric(lower) == TRUE) & (is.numeric(upper) == TRUE) & (is.numeric(n)))

  # Find valid domain if user doesn't input one.

  if(lower == -Inf | upper == Inf) {
    warning('Setting bounds by default.')
    domain <- findDomain(g)
    lower <- max(domain[1], lower)
    upper <- min(domain[2], upper)
  }

  # Reverse upper and lower values if upper is smaller than lower.

  if (lower > upper){
    warning('Automatically reverse the order of bounds.')
    report <- list(lower, upper)
    report[1] <- upper
    report[2] <- lower
    upper <- report[2]
    lower <- report[1]
  }

  # Check if upper and lower are equal.

  assert_that (lower!=upper)

  # Check that g is valid probability density.

  assert_that(is.function(g))
  assert_that(all(g(seq(from = lower, to = upper, by = n)) >= 0))

  h <- function(x) log(g(x))

  # Check concavity of h here.

  assert_that(isConcave(FUN = h, x1 = lower, x2 = upper))

  x_samples <- rep(NA, n)

  # Initialize T_k.

  Tk <- list()
  Tk$x <- init_abscissae(h, lower, upper)
  if (any(!is.finite(Tk$x))) {
    stop("initialize the abscissae function returns invalid values of x (infinite / NA).")
  }

  # Store the function values of h(x_j) and h'(x_j) to avoid repeated evaluation.

  Tk$hx <- h(Tk$x)
  dx <- 1e-8

  #Tk$hpx <- grad(h, Tk$x)

  Tk$hpx <- (h(Tk$x + dx) - Tk$hx)/dx

  # Check if upper and lower bounds are not actually bounding the density.

  if (any(!is.finite(Tk$hx))) {
    stop("h(x) has invalid values (infinite / NA).")
  }
  if (any(!is.finite(Tk$hpx))) {
    stop("h'(x) has invalid values (infinite / NA).")
  }

  # Calculate z_j's using utility function and store.
  # Here j = 1,...,k+1 instead of 0,...,k since R vectors are 1-indexed.

  Tk$z <- z_fun(Tk$hx, Tk$hpx, Tk$x, lower, upper)

  # Store u_k, inverse CDF of s_k, l_k as functions to avoid redundant calculation.

  Tk$u <- function(x_val) { u_fun(x_val, Tk$hx, Tk$hpx, Tk$x, Tk$z) }
  Tk$sinvcdf <- function(p) { sinvcdf_fun(p, Tk$hx, Tk$hpx, Tk$x, Tk$z, Tk$u) }
  Tk$l <- function(x_val) { l_fun(x_val, Tk$hx, Tk$x) }

  # Keep drawing samples until reach n points.

  curr <- 1
  while (curr <= n) {
    sample_result <- ars_sample(h, Tk)
    if (sample_result$accept) {

      # If x* is accepted during sampling, add to output vector.

      x_samples[curr] <- sample_result$x_star
      curr <- curr + 1
    }

    # If h(x*) is evaluated, update T_k to T_k+1.

    if (sample_result$hEval) Tk <- update_Tk(Tk, h, sample_result$x_star)
  }
  return(x_samples)
}
