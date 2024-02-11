## ==============================================================================
# library_simulation.R

# Key functions for generating simulated data


# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca


# ============================================================================== ==============================================================================
print("Loading library_simulation.R")

# Continous approximation of the sign function

sigmoid <- function(x) 1 / (1 + exp(-x))
dsigmoid <- function(x) sigmoid(x) * (1 - sigmoid(x))
d2sigmoid <- function(x) dsigmoid(x) * (1 - 2 * sigmoid(x))
approx_sign <- function(t, lam) sigmoid(lam * t)
dapprox_sign <- function(t, lam) lam * dsigmoid(lam * t)
ddapprox_sign <- function(t, lam) lam^2 * d2sigmoid(lam * t)
##
# ==============================================================================
# Continous approximation to the Heavisine function
# See Donoho and Johnstone (1994)
# ==============================================================================
heavisine <- function(t, lam) {
    term1 <- 4 * sin(4 * pi * t)
    term2 <- -approx_sign(t - 0.3, lam)
    term3 <- -approx_sign(0.72 - t, lam)
    return(term1 + term2 + term3)
}


dheavisine <- function(t, lam) {
    dterm1 <- 4 * 4 * pi * cos(4 * pi * t)
    dterm2 <- -dapprox_sign(t - 0.3, lam)
    dterm3 <- dapprox_sign(0.72 - t, lam)
    return(dterm1 + dterm2 + dterm3)
}

ddheavisine <- function(t, lam) {
    ddterm1 <- -4 * pi * 4 * 4 * pi * sin(4 * pi * t)
    ddterm2 <- -ddapprox_sign(t - 0.3, lam)
    ddterm3 <- -ddapprox_sign(0.72 - t, lam)
    return(ddterm1 + ddterm2 + ddterm3)
}

# ==============================================================================
# Continous approximation to the Blocks
# See Donoho and Johnstone (1994)
# ==============================================================================
blocks <- function(t, lambda = 500,
                   tj = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.82),
                   hj = c(4, -5, 4, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)) {
    f <- rep(0, length(t))
    for (i in 1:length(tj)) {
        f <- f + hj[i] * ((1 + approx_sign(t - tj[i], lambda)) / 2)
    }
    return(f)
}

# 1st deriative
dblocks <- function(t, lambda = 500,
                    tj = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.82),
                    hj = c(4, -5, 4, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)) {
    f <- rep(0, length(t))
    for (i in 1:length(tj)) {
        f <- f + hj[i] * ((1 + dapprox_sign(t - tj[i], lambda)) / 2)
    }
    return(f)
}

# 2nd deriative
ddblocks <- function(t, lambda = 500,
                     tj = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.82),
                     hj = c(4, -5, 4, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)) {
    f <- rep(0, length(t))
    for (i in 1:length(tj)) {
        f <- f + hj[i] * ((1 + ddapprox_sign(t - tj[i], lambda)) / 2)
    }
    return(f)
}


# ==============================================================================
# Doppler function
# See Donoho and Johnstone (1994)
# ==============================================================================
# Doppler

doppler <- function(t, lambda = 0.05) {
    return((t * (1 - t))^(0.5) * sin((2 * pi * (1 + lambda)) / (t + lambda)))
}

ddoppler <- function(t, lambda = 0.05) {
    output <- (0.5 * (t * (1 - t))^(-0.5) * (1 - 2 * t) * sin((2 * pi * (1 + lambda)) / (t + lambda)) +
        (t * (1 - t))^(0.5) * cos((2 * pi * (1 + lambda)) / (t + lambda)) * (-2 * pi * (1 + lambda)) / (t + lambda)^2)
    return(output)
}

dddoppler <- function(t, lambda = 0.05) {
    lam1 <- lambda + 1
    twopi <- 2 * pi
    output <- ((-((1 - 2 * t)^2) / (4 * ((1 - t) * t)^(3 / 2)) - ((1 - t) * t)^(-1 / 2)) * sin((2 * pi * (lambda + 1)) / (lambda + t))) -
        ((twopi * (lam1) * (1 - 2 * t) * cos((twopi * lam1) / (lambda + t))) / (((1 - t) * t)^(1 / 2) * (lambda + t)^2)) +
        (((1 - t) * t)^(1 / 2) * (((4 * pi * lam1 * cos((twopi * lam1) / (lambda + t))) / (lambda + t)^3) - ((4 * pi^2 * lam1^2 * (sin((2 * pi * lam1) / (lambda + t)))) / ((lambda + t)^4))))
    return(output)
}
