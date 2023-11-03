if (!require("truncnorm")) install.packages("truncnorm")

library("truncnorm")


source("./src/library.R")
seed <- 123
# Set RNG seed
set.seed(seed)
# ============================= Load Data ======================================

print("    Prepping data...")

# Process data
y <- as.numeric(read.delim2("data/dowling/Dowling_Disp.txt")[, 1])
ddy <- as.numeric(read.delim2("data/dowling/Dowling_Accn.txt")[1:length(y), 1])
fs <- 512
t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)

# standardize
y_mu <- mean(y)
y_sd <- sd(y)
y_std <- (y - y_mu) / y_sd

print("    Prepping spline basis...")

# Main Basis
deg <- 5
pDeg <- 2
nK <- 100 # ceiling(100*max(t))  #floor(length(t) / 5) we defined 50 basis functions per second
B <- basis(t, nK, deg = deg)
Bm <- B$matrix
K <- dim(Bm)[2]
ddB <- dbasis(B, 2)
ddBm <- ddB$matrix
P <- penalty(dim(B$matrix)[2], pDeg)
MM <- basis_mm(P)

# Basis for tau.
nK_sm <- floor(nK / 5)
bDeg_sm <- 3
C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)
Cm <- C$matrix
K_C <- dim(Cm)[2]


# ===================== 1) P-Spline  ===========================
# Determine scale xi for half t prior on tau^2 via prior predictive simulation
# across a range of xi values.

d_xi <- seq(1e-5, 1e-2, length.out = 50)
qderv <- rep(NA, length(d_xi))

for (l in 1:length(d_xi)) {
    print(paste0("l = ", l, " xi = ", d_xi[l]))

    xi <- d_xi[l]
    max_derv <- rep(NA, 20000)

    for (i in 1:length(max_derv)) {
        # Generate tau2 ~ t_3+(0, xi^2)
        x <- rchisq(1, df = 3)
        v <- (3 * xi^2) / x # v ~ invchi2(nu, xi^2)
        tau2 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt(v)) # tau2 ~ t_3+(0, xi^2)

        tau <- sqrt(tau2)

        # theta | . is 2nd order RW
        theta <- rep(0, K)
        for (k in (pDeg + 1):K) {
            theta[k] <- 2 * theta[k - 1] - theta[k - 2] + tau * rnorm(1)
        }

        # Compute max of 2nd deriative
        max_derv[i] <- max(abs(ddBm %*% theta))
    }

    # Compute 90th quantile of max(|y''(t)|)
    qderv[l] <- quantile(max_derv, 0.90, na.rm = TRUE)
    print(paste0(".   q = ", qderv[l]))
}

# Visualize relationship between xi and 90% quantile of max(|y''(t)|)
plot(d_xi, qderv,
    type = "l",
    ylab = "",
    xlab = TeX(r"($\tilde{\xi}$)"),
    ylim = c(0, 2000)
)
abline(h = 1000, col = "red", lty = 2)
title(ylab = TeX(r"(90% quantile: $max(|y''(t)|)$)"), line = 2, cex.lab = 1)

# Determine which xi meets the 90% quantile threshold
xi <- d_xi[which.min(abs(qderv - 1000))]

# Suggests that xi = 0.002 meets criteria for half-t prior on tau^2

# ===================== 2) APS_ind ===========================
# Determine scale xi for half t prior on tau^2 via prior predictive simulation
# across a grid of prior scale values d_xi


d_xis <- seq(1e-3, 2e-2, length.out = 50) # grid of d_xis
qderv <- matrix(NA, length(d_xis))

for (l in 1:length(d_xis)) {
    print(paste0("l = ", l, " xi = ", d_xis[l]))

    dxi <- d_xis[l]
    max_derv <- rep(NA, 20000)

    nu <- 3 # 1/runif(1)  # Can also simulate nu from its inverse uniform prior

    for (i in 1:length(max_derv)) {
        # Generate xi ~ t_3+(0, dxi^2)
        x <- rchisq(1, df = 3)
        v_xi <- (3 * dxi^2) / x # tau2_k ~ invchi2(nu, xi^2)
        xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt(v_xi))

        # tau2_k ~ invchi2(nu, xi^2)
        x <- rchisq(K - pDeg, df = 3)
        tau2_k <- (nu * xi^2) / x # tau2_k ~ invchi2(nu, xi^2)

        # theta | . is 2nd order RW
        theta <- rep(0, K)
        for (k in (pDeg + 1):K) {
            theta[k] <- 2 * theta[k - 1] - theta[k - 2] + sqrt(tau2_k[k - pDeg]) * rnorm(1)
        }

        # Compute max of 2nd deriative
        max_derv[i] <- max(abs(ddBm %*% theta))
    }

    # Compute 90th quantile of max(|y''(t)|)
    qderv[l] <- quantile(max_derv, 0.90)
    print(paste0(".   q = ", qderv[l]))
}

# visualize relationship between scale of prior and 90% quantile of max(|y''(t)|)
plot(d_xis, qderv,
    type = "l",
    ylab = TeX(r"(90% quantile of $max(|f''(t)|) \, \, (rad\cdot s^{-2})$)"),
    xlab = TeX(r"(Scale of prior for $\xi$)"),
    ylim = c(0, 2000)
)
abline(h = 1000, col = "red", lty = 2)

# Determine which prior scale meets the 90% quantile threshold
d_xis[min(which(qderv > 1000))]


# ===================== 3) APS_ar ===========================
# Determine scale xi for half t prior on tau^2 via prior predictive simulation
# Note prior for alpha_0 is fixed so that P-Spline is within space of models

# Define grid of prior scales
delta_xi0 <- seq(1, 1.5, length.out = 50)
delta_a0 <- 0.002 # Set to ensure that standard P-Spline is within space of models

hparms <- expand.grid(xi0 = delta_xi0, a0 = delta_a0)
qderv <- rep(NA, nrow(hparms))

for (l in 1:nrow(hparms)) {
    print(paste0("l = ", l, "/", nrow(hparms), ": delta_xi0 = ", hparms[l, "xi0"], " - delta_a0 = ", hparms[l, "a0"]))

    delta_ao <- hparms[l, "a0"]
    d_xi0 <- hparms[l, "xi0"]
    max_derv <- rep(NA, 20000)
    for (i in 1:length(max_derv)) {
        # generate \alpha_0 from half-t prior with delta_a0 = 0.002
        x <- rchisq(1, df = 3)
        a0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * delta_ao^2) / x))

        # generate xi from half t priro with scale d_xi
        x <- rchisq(1, df = 3)
        xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * d_xi0^2) / x))

        # generate phi from uniform prior
        phi <- runif(1)

        # ltau2 | phi, xi, a0 ~ AR(1)
        ltau2_k <- rep(NA, K - pDeg)
        ltau2_k[1] <- log(a0^2) + xi * rnorm(1)
        for (k in 2:(K - pDeg)) {
            ltau2_k[k] <- log(a0^2) + phi * (ltau2_k[k - 1] - log(a0^2)) + xi * rnorm(1)
        }

        # theta | . is 2nd order RW
        theta <- rep(0, K)
        for (k in (pDeg + 1):K) {
            theta[k] <- 2 * theta[k - 1] - theta[k - 2] + sqrt(exp(ltau2_k[k - pDeg])) * rnorm(1)
        }

        # Compute max of 2nd deriative
        max_derv[i] <- max(abs(ddBm %*% theta))
    }

    # generate 90th quantile of deriative
    qderv[l] <- quantile(max_derv, 0.90, na.rm = T)
    print(paste0(".   q = ", qderv[l]))
}

# visualize relationship between scale of prior and 90% quantile of max(|y''(t)|)
plot(hparms[, "xi0"], qderv,
    type = "l",
    ylab = TeX(r"(90% quantile of $max(|f''(t)|) \, \, (rad\cdot s^{-2})$)"),
    xlab = TeX(r"(Scale of prior for $\xi$)"),
    ylim = c(0, 2000)
)
abline(h = 1000, col = "red", lty = 2)

# Determine which prior scale meets the 90% quantile threshold
hparms[min(which(qderv > 1000)), "xi0"]

# Suggests 1.25 sufficient

# ===================== 3) APS_spl ===========================
# Determine scale xi for half t prior on tau^2 via prior predictive simulation
# Note prior for alpha_0 is fixed so that P-Spline is within space of models

# Define grid of prior scales
delta_xi0 <- seq(1, 1.5, length.out = 50)
delta_a0 <- 0.002
hparms <- expand.grid(xi0 = delta_xi0, a0 = delta_a0)

qderv <- rep(NA, nrow(hparms))
for (l in 1:nrow(hparms)) {
    print(paste0("l = ", l, " delta_xi0 = ", hparms[l, "xi0"]))

    delta_xi0 <- hparms[l, "xi0"]
    delta_a0 <- hparms[l, "a0"]
    max_derv <- rep(NA, 20000)
    for (i in 1:length(max_derv)) {
        # generate \alpha_0 from half-t prior with delta_a0 = 0.002
        x <- rchisq(1, df = 3)
        a0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * delta_a0^2) / x)) # tau2_k ~ invchi2(nu, xi^2)

        # generate xi from half t priro with scale d_xi
        x <- rchisq(1, df = 3)
        xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * delta_xi0^2) / x))

        # gamma ~ N(log(a0^2), xi^2)
        # ltau2 | phi, xi, a0 ~ Cubic Spline given basis matrix C
        gamma <- log(a0^2) + xi * rnorm(K_C, 1)
        ltau2_k <- Cm %*% gamma

        # theta | . is 2nd order RW
        theta <- rep(0, K)
        for (k in (pDeg + 1):K) {
            theta[k] <- 2 * theta[k - 1] - theta[k - 2] + sqrt(exp(ltau2_k[k - pDeg])) * rnorm(1)
        }
        # Compute max of 2nd deriative
        max_derv[i] <- max(abs(ddBm %*% theta))
    }
    # generate 90th quantile of deriative
    qderv[l] <- quantile(max_derv, 0.95, na.rm = T)
    print(paste0(".   q = ", qderv[l]))
}

# visualize relationship between scale of prior and 90% quantile of max(|y''(t)|)
plot(hparms[, "xi0"], qderv,
    type = "l",
    ylab = TeX(r"(90% quantile of $max(|f''(t)|) \, \, (rad\cdot s^{-2})$)"),
    xlab = TeX(r"(Scale of prior for $\xi$)"),
    ylim = c(0, 2000)
)
abline(h = 1000, col = "red", lty = 2)

# Determine which prior scale meets the 90% quantile threshold
hparms[min(which(qderv > 1000)), "xi0"]

# Suggests 1.15 sufficient
