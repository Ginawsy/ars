#' @importFrom assertthat assert_that

#' Calculate the intersection points z_j's of the tangents at x_j's.
#' x: vector of k abscissae.
#' hx: vector of function evals at x.
#' hpx: vector of derivative evals at x.
#' lower and upper: the lower bound and upper bound of domain D.

z_fun <- function(hx, hpx, x, lower, upper) {

  k <- length(x)

  z <- rep(0, k+1)
  z[1] <- lower
  z[k+1] <- upper
  assert_that(length(hx) == k)
  assert_that(length(hpx) == k)

  # The original formula for z_j's are for j in 0,...,k.
  # Since R vectors are 1-indexed instead of 0-indexed,
  # the range of j is shifted to 1,...,k+1.
  # So instead of filling in z_1,...,z_k-1, we fill in z_2,...,z_k.
  #for (j in 1:k-1) {
  #  z[j+1] <- (hx[j+1] - hx[j] - x[j+1] * hpx[j+1] + x[j] * hpx[j])/
  #            (hpx[j] - hpx[j+1])
  #}

  z[2:k] <- (hx[2:k] - hx[1:k-1] - x[2:k] * hpx[2:k] + x[1:k-1] * hpx[1:k-1]) / (hpx[1:k-1] - hpx[2:k])
  return(z)
}

#' Helper function that returns j such that z[j] <= x_val < z[j+1].
#' i.e. z_j-1 <= x_val < z_j (since index of z is shifted).

zfindj <- function(x_val, z) {

  # Verify that x_val is within D.
  # assert_that(x_val >= z[1] & x_val <= z[length(z)])

  for (j in 1:(length(z)-1)) {
    if(z[j+1] > x_val) return(j)
  }
  return(length(z) - 1)
}

#' Return u_k(x) given x (denoted x_val here, similar below).

u_fun <- function(x_val, hx, hpx, x, z){
  # assert_that(z[1]<= x_val && z[length(z)]>= x_val)
  k <- length(x)

  j <- zfindj(x_val, z)

  return(hx[j] + (x_val-x[j]) * hpx[j])
}


#' Helper function that calculate the integral of exp u_k(x) over D.

calculateC <- function(u, z, hpx) {
  k <- length(hpx)
  cdf <- 0
  for (j in 1:k) {
    cdf <- cdf + (exp(u(z[j+1])) - exp(u(z[j]))) / hpx[j]
  }
  return(cdf)
}


#' Helper function that find j such that the CDF up to z[j] <= p < CDF up to z[j+1].
#' Returns a vector of: [1] j [2] CDF up to z[j] (z[1] to z[j]).

pfindj <- function(p, C, hpx, z, u) {
  cdf_val <- 0
  for (j in 1:(length(z)-1)) {
    cdf_curr <- (exp(u(z[j + 1])) - exp(u(z[j]))) / (C * hpx[j])
    if(cdf_val + cdf_curr > p){
      return(c(j, cdf_val))
    }
    else{
      cdf_val <- cdf_val + cdf_curr
    }
  }
  return(c(length(z), 1))
}

#' Calculate the inverse CDF of s_k(x), giving the probability p, p in [0,1].
#' More explanation in helper document.

sinvcdf_fun <- function(p, hx, hpx, x, z, u){
  assert_that(p >= 0 && p<=1)
  k <- length(x)
  C <- calculateC(u, z, hpx)

  cdf_result <- pfindj(p, C, hpx, z, u)
  j <- cdf_result[1]
  cdf_val <- cdf_result[2]

  ux <- log( exp(u(z[j])) + (p - cdf_val) * C * hpx[j] )
  x <- (ux - hx[j] + x[j] * hpx[j]) / hpx[j]
  return(x)
}

#' Helper function that returns j such that x_j <= x_val < x_j+1.
#' Returns 0 if x_val < x_1 or x_val > x_k.

xfindj <- function(x_val, x) {
  if (x_val < x[1] || x_val >= x[length(x)]) return(0)
  for (j in 1:length(x)-1) {
    if(x[j+1] > x_val) return(j)
  }
  return(length(x))
}

#' Return l_k(x) given x (denoted x_val)

l_fun <- function(x_val, hx, x){

  j <- xfindj(x_val, x)
  if (j == 0) return(-Inf)

  return(((x[j + 1]-x_val) * hx[j] + (x_val-x[j]) * hx[j + 1]) / (x[j + 1]-x[j]))
}

#' Source: https://stackoverflow.com/questions/58735415/check-convexity-of-any-function-in-r.
#' x1 is the lower bound, x2 is the upper bound, must be finite.

isConcave <- function(FUN, x1, x2) { # Don't need lambda as an argument.
  lambda <- seq(0, 1, 0.01)
  test <- FUN(lambda * x1 + (1-lambda) * x2 ) - lambda * FUN(x1) - (1-lambda) * FUN(x2) >= 0

  return(all(test))
}

#' find a domain if unspecified by the user, assuming the function
#' is a valid probability density function

findDomain <- function(FUN){
  domain <- seq(-500,500,1)

  vals <- FUN(domain)

  valid <- which(vals > 0)

  bounds <- c(domain[valid[1]], domain[valid[length(valid)]])

  return(bounds)
}
