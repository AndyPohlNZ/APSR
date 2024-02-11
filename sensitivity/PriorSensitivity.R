# ==============================================================================
# PriorSensitivty.R

# Performs sensitivity studies for parameters of the prior distribution for APS
# models.


# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================
setwd(".")

library("rjags")
library("bayesplot")
library("coda")
library("snow")

source("./src/library.R")
source("./src/library_filter.R")

# ============================= Preliminaries ======================================
# Get command line arguments
seed <- 123
K <- 100

# Set RNG seed
set.seed(seed)
# ============================= Load Data ======================================

print("Running AdPSpline analysis on Dowling data...")
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

# Set MCMC parameters
nchains <- 4
nadapt <- 5000
nsampling <- 10000
thin <- 1

# ====================== Initialize spline basis ==============================================
print("    Prepping spline basis...")

# Main Basis
deg <- 5
pDeg <- 2
nK <- K # ceiling(100*max(t))  #floor(length(t) / 5) we defined 50 basis functions per second
B <- basis(t, nK, deg = deg)
ddB <- dbasis(B, 2)
P <- penalty(dim(B$matrix)[2], pDeg)
MM <- basis_mm(P)

# Basis for tau.
nK_sm <- floor(nK / 5)
bDeg_sm <- 3
C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)

# ============================= P-Spline ================================
print("    Computing Bayes Models:")
print("        1: PSpline")
xi_ins <- 0.002 * c(5, 10, 100)
for (i in 1:length(xi_ins)) {
    xi_in <- xi_ins[i]
    print(paste0("        2: PSpline - xi = ", xi_in))


    coda.samples.wrapper <- function(j) {
        init_fn <- function() {
            mle <- mle_pspl(B, MM, y_std)
            theta <- as.numeric(mle$theta_hat)
            list(
                theta = theta,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = j
            )
        }
        mod <- jags.model("./mdl/jags/PSpline.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                ord = B$deg,
                xi = xi_in
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("xi_in", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
    out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
    for (i in 1:length(out)) {
        out[[i]] <- out[[i]][[1]]
    }
    class(out) <- "mcmc.list"
    stopCluster(cl)
    t2 <- proc.time()

    # compute y and ddy
    nchains <- length(out)
    nsamps <- dim(out[[1]])[1]
    theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
    for (i in 1:nchains) {
        theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
    }
    y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    for (i in 1:nrow(theta_samps)) {
        y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
        ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
    }

    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "tau"]
    }

    saveRDS(
        list(
            samples = out,
            y = y_samps,
            ddy = ddy_samps,
            tau = tau_samps,
            time = t2 - t1
        ),
        file = paste0("./PriorSensitivity/mdl/PSpline_xi-", xi_in, "_seed-", seed, ".RDS")
    )
}
# ========================= Prior sensitivity ==================================

xis <- 0.002 * c(5, 10, 100)
ddy_hat <- y_hat <- matrix(data = NA, ncol = length(t), nrow = length(xis))
y_rmse <- ddy_rmse <- rep(NA, nrow(y_hat))

for (i in 1:length(xis)) {
    mdl <- readRDS(paste0("./PriorSensitivity/mdl/PSpline_xi-", xis[i], "_seed-", 123, ".RDS"))
    y_hat[i, ] <- apply(mdl$y, 2, median)

    ddy_hat[i, ] <- apply(mdl$ddy, 2, median)

    y_rmse[i] <- rmse(y, y_hat[i, ])
    ddy_rmse[i] <- rmse(ddy, ddy_hat[i, ])
}

plot(t, ddy, type = "l")
lines(t, est[1, ], col = adjustcolor(pal[1], alpha = 0.5))
lines(t, est[2, ], col = adjustcolor(pal[2], alpha = 0.5))
lines(t, est[3, ], col = adjustcolor(pal[3], alpha = 0.5))
lines(t, est[4, ], col = adjustcolor(pal[4], alpha = 0.5))

# Compute plots

rmse(ddy, est[1, ])
rmse(ddy, est[2, ])
rmse(ddy, est[3, ])
rmse(ddy, est[4, ])

# ============================== APS_ind ======================================
xi_ins <- 1.65 * c(5, 10, 100)
for (i in 1:length(xi_ins)) {
    xi_in <- xi_ins[i]

    print(paste0("        2: APS_ind xi = ", xi_in))
    coda.samples.wrapper <- function(j) {
        init_fn <- function() {
            mle <- mle_pspl(B, MM, y_std)
            theta <- as.numeric(mle$theta_hat)
            list(
                theta = theta,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = j
            )
        }
        mod <- jags.model("./mdl/jags/APS_ind.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                d_a0 = 0.002,
                d_xi0 = xi_in
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("xi_in", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
    out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
    for (i in 1:length(out)) {
        out[[i]] <- out[[i]][[1]]
    }
    class(out) <- "mcmc.list"
    stopCluster(cl)
    t2 <- proc.time()

    # compute y and ddy
    nchains <- length(out)
    nsamps <- dim(out[[1]])[1]
    theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
    for (i in 1:nchains) {
        theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
    }
    y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    for (i in 1:nrow(theta_samps)) {
        y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
        ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
    }

    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2] - pDeg)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(dim(B$matrix)[2] - pDeg), "]")]
    }
    #
    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0("./PriorSensitivity/mdl/APS_ind_xi-", xi_in, "_seed-", seed, ".RDS"))
}



xis <- 1.65 * c(5, 10, 100)
ddy_hat <- y_hat <- matrix(data = NA, ncol = length(t), nrow = length(xis))
y_rmse <- ddy_rmse <- rep(NA, nrow(y_hat))

for (i in 1:length(xis)) {
    mdl <- readRDS(paste0("./PriorSensitivity/mdl/APS_ind_xi-", xis[i], "_seed-", 123, ".RDS"))
    y_hat[i, ] <- apply(mdl$y, 2, median)

    ddy_hat[i, ] <- apply(mdl$ddy, 2, median)

    y_rmse[i] <- rmse(y, y_hat[i, ])
    ddy_rmse[i] <- rmse(ddy, ddy_hat[i, ])
}

plot(t, ddy, type = "l")
lines(t, ddy_hat[1, ], col = adjustcolor(pal[1], alpha = 0.5))
lines(t, ddy_hat[2, ], col = adjustcolor(pal[2], alpha = 0.5))
lines(t, ddy_hat[3, ], col = adjustcolor(pal[3], alpha = 0.5))
lines(t, ddy_hat[4, ], col = adjustcolor(pal[4], alpha = 0.5))

# Compute plots
rmse(y, y_hat[1, ])
rmse(y, y_hat[2, ])
rmse(y, y_hat[3, ])
rmse(y, y_hat[4, ])


rmse(ddy, ddy_hat[1, ])
rmse(ddy, ddy_hat[2, ])
rmse(ddy, ddy_hat[3, ])
rmse(ddy, ddy_hat[4, ])

# ============================== APS_ar ======================================

xi_ins <- 1.25 * c(5, 10, 100)
for (i in 1:length(xi_ins)) {
    xi_in <- xi_ins[i]
    print(paste0("        2: APS_ar - xi = ", xi_in))
    coda.samples.wrapper <- function(j) {
        init_fn <- function() {
            mle <- mle_pspl(B, MM, y_std)
            theta <- as.numeric(mle$theta_hat)
            list(
                theta = theta,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = j
            )
        }
        mod <- jags.model("./mdl/jags/APS_ar.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                d_a0 = 0.002,
                d_xi0 = xi_in
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "sigma", "xi", "phi", "alpha"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("xi_in", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
    out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
    for (i in 1:length(out)) {
        out[[i]] <- out[[i]][[1]]
    }
    class(out) <- "mcmc.list"
    stopCluster(cl)
    t2 <- proc.time()

    # compute y and ddy
    nchains <- length(out)
    nsamps <- dim(out[[1]])[1]
    theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
    for (i in 1:nchains) {
        theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
    }
    y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    for (i in 1:nrow(theta_samps)) {
        y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
        ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
    }

    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2] - pDeg)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(dim(B$matrix)[2] - pDeg), "]")]
    }

    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0("./PriorSensitivity/mdl/APS_ar_xi-", xi_in, "_seed-", seed, ".RDS"))
}

xis <- 1.25 * c(5, 10, 100)
ddy_hat <- y_hat <- matrix(data = NA, ncol = length(t), nrow = length(xis))
y_rmse <- ddy_rmse <- rep(NA, nrow(y_hat))

for (i in 1:length(xis)) {
    mdl <- readRDS(paste0("./PriorSensitivity/mdl/APS_ar_xi-", xis[i], "_seed-", 123, ".RDS"))
    y_hat[i, ] <- apply(mdl$y, 2, median)

    ddy_hat[i, ] <- apply(mdl$ddy, 2, median)

    y_rmse[i] <- rmse(y, y_hat[i, ])
    ddy_rmse[i] <- rmse(ddy, ddy_hat[i, ])
}

plot(t, ddy, type = "l")
lines(t, ddy_hat[1, ], col = adjustcolor(pal[1], alpha = 0.5))
lines(t, ddy_hat[2, ], col = adjustcolor(pal[2], alpha = 0.5))
lines(t, ddy_hat[3, ], col = adjustcolor(pal[3], alpha = 0.5))
lines(t, ddy_hat[4, ], col = adjustcolor(pal[4], alpha = 0.5))

# Compute plots


rmse(ddy, est[1, ])
rmse(ddy, est[2, ])
rmse(ddy, est[3, ])
rmse(ddy, est[4, ])

# ============================== APS_spl ======================================
xi_ins <- 1.15 * c(1, 5, 10, 100)
for (i in 1:length(xi_ins)) {
    xi_in <- xi_ins[i]
    print(paste0("        2: APS_spl - xi = ", xi_in))
    coda.samples.wrapper <- function(j) {
        init_fn <- function() {
            mle <- mle_pspl(B, MM, y_std)
            theta <- as.numeric(mle$theta_hat)

            list(
                theta = theta,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = j
            )
        }
        mod <- jags.model("./mdl/jags/APS_spl.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                Q_sm = dim(C$matrix)[2],
                B = B$matrix,
                C = C$matrix,
                y = y_std,
                d_a0 = 0.002,
                d_xi0 = xi_in
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "sigma", "xi", "phi", "alpha"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("xi_in", "MM_likelihood", "mle_pspl", "t", "B", "C", "MM", "y_std", "nadapt", "nsampling", "thin"))
    out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
    for (i in 1:length(out)) {
        out[[i]] <- out[[i]][[1]]
    }
    class(out) <- "mcmc.list"
    stopCluster(cl)
    t2 <- proc.time()

    # compute y and ddy
    nchains <- length(out)
    nsamps <- dim(out[[1]])[1]
    theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
    for (i in 1:nchains) {
        theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
    }
    y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
    for (i in 1:nrow(theta_samps)) {
        y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
        ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
    }

    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2] - pDeg)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(dim(B$matrix)[2] - pDeg), "]")]
    }


    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0("./PriorSensitivity/mdl/APS_spl_xi-", xi_in, "_seed-", seed, ".RDS"))
}


xis <- 1.15 * c(1, 5, 10, 100)
ddy_hat <- y_hat <- matrix(data = NA, ncol = length(t), nrow = length(xis))
y_rmse <- ddy_rmse <- rep(NA, nrow(y_hat))

for (i in 1:length(xis)) {
    mdl <- readRDS(paste0("./PriorSensitivity/mdl/APS_spl_xi-", xis[i], "_seed-", 123, ".RDS"))
    y_hat[i, ] <- apply(mdl$y, 2, median)

    ddy_hat[i, ] <- apply(mdl$ddy, 2, median)

    y_rmse[i] <- rmse(y, y_hat[i, ])
    ddy_rmse[i] <- rmse(ddy, ddy_hat[i, ])
}

plot(t, ddy, type = "l")
lines(t, ddy_hat[1, ], col = adjustcolor(pal[1], alpha = 0.5))
lines(t, ddy_hat[2, ], col = adjustcolor(pal[2], alpha = 0.5))
lines(t, ddy_hat[3, ], col = adjustcolor(pal[3], alpha = 0.5))
lines(t, ddy_hat[4, ], col = adjustcolor(pal[4], alpha = 0.5))

# Compute plots
rmse(ddy, est[1, ])
rmse(ddy, est[2, ])
rmse(ddy, est[3, ])
rmse(ddy, est[4, ])
