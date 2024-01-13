# ==============================================================================
# Functions for generating simulated data
# Used for simulations performed in
# Adapative P-Splines for challenging filtering problems in biomechanics
# Pohl et al. (2024)
# ==============================================================================

print("Loading library_simulation.R")

# Continuous approximation of the sgn function
approx_sign <- function(t, lam) {
    # Approximate the sign function via atan

    return(2 / pi * atan(lam * t))
}

dapprox_sign <- function(t, lam) {
    # Approximate the sign function via atan
    return((2 * lam) / (pi * (1 + (lam * t)^2)))
}

ddapprox_sign <- function(t, lam) {
    # Approximate the sign function via atan

    return(-4 * lam^2 * t / (pi * (1 + (lam * t)^2)^2))
}

# Continuous approximation of the absolute value function
approx_abs <- function(t, lam) {
    return(t^2 / (sqrt(t^2 + lam^2)))
}

dapprox_abs <- function(t, lam) {
    return((2 * t * (t^2 + lam^2) - t^3) / ((t^2 + lam^2)^(3 / 2)))
}

ddapprox_abs <- function(t, lam) {
    return(
        (2 * (t^2 + lam^2)^(3 / 2) - 3 * t^2 * (t^2 + lam^2)^(1 / 2)) /
            ((t^2 + lam^2)^3)
    )
}

# ==============================================================================
# Continous approximation to the Heavisine function
# See Donoho and Johnstone (1994)
# ==============================================================================
heavisine <- function(t, lam) {
    return(
        4 * sin(4 * pi * t) -
            approx_sign(t - 0.3, lam) -
            approx_sign(0.72 - t, lam)
    )
}

# 1st deriative
dheavisine <- function(t, lam) {
    return(
        16 * pi * cos(4 * pi * t) -
            dapprox_sign(t - 0.3, lam) +
            dapprox_sign(0.72 - t, lam)
    )
}

# 2nd deriative
ddheavisine <- function(t, lam) {
    return(
        -64 * pi^2 * sin(4 * pi * t) -
            ddapprox_sign(t - 0.3, lam) -
            ddapprox_sign(0.72 - t, lam)
    )
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
# Continous approximation to the Bumps function
# See Donoho and Johnstone (1994)
# ==============================================================================
# TODO:  Note currently buggy ... - need to fix


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
