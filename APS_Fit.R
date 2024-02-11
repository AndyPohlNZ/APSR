# ==============================================================================
# APS_Fit.R
#
# Steps through applying the APS model to the dowling dataset.
#
# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================

if (!require("rjags")) install.packages("rjags")
if (!require("coda")) install.packages("coda")
if (!require("snow")) install.packages("snow")

library("rjags")
library("coda")
library("snow")

source("./src/library.R")
source("./src/library_filter.R")

# ============================= Preliminaries ======================================


seed <- 123
model_type <- "All" # One of "FIlter", "PSpline", "APS_ind", "APS_spl", "APS_ar" or "All"
save_folder <- "./results/"

# Set RNG seed
set.seed(seed)

# Set MCMC parameters
nchains <- 4 # Number of independent chains
nadapt <- 5000 # Size of adaption phase
nsampling <- 10000 # Number of samples per chain
thin <- 1 # Thinning rate per chain.

# ============================= Load Data ======================================

print("Running AdPSpline analysis on Dowling data...")
print("    Prepping data...")

# Process data
y <- as.numeric(read.delim2("data/Dowling_Disp.txt")[, 1])
ddy <- as.numeric(read.delim2("data/Dowling_Accn.txt")[1:length(y), 1])

# Set time vector
fs <- 512 # Sampling frequency of Dowling Dataset
t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)

# standardize
y_mu <- mean(y)
y_sd <- sd(y)
y_std <- (y - y_mu) / y_sd

# ====================== Initialize spline basis ==============================================
print("    Prepping spline basis...")

# Main Basis
deg <- 5
pDeg <- 2
nK <- 100
B <- basis(t, nK, deg = deg)
ddB <- dbasis(B, 2)
P <- penalty(dim(B$matrix)[2], pDeg)
MM <- basis_mm(P)

# Basis for tau.
nK_sm <- floor(nK / 5)
bDeg_sm <- 3
C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)

# ============================= Compute filters ================================
if (model_type %in% c("Filter", "All")) {
    print("    Computing filters...")
    t1 <- proc.time()
    fc_acf <- cutoff_acf(y_std, fs)
    y_acf <- filtbutter(y_std, butter(2, fc_acf$fc / (fs / 2), type = "low")) * y_sd + y_mu
    ddy_acf <- finite_diff(y_acf, t, 2)
    t2 <- proc.time()
    acf_time <- t2 - t1

    t1 <- proc.time()
    fc_resid <- cutoff_residual(y_std, fs)
    y_resid <- filtbutter(y_std, butter(2, fc_resid$fc / (fs / 2), type = "low")) * y_sd + y_mu
    ddy_resid <- finite_diff(y_resid, t, 2)
    t2 <- proc.time()
    resid_time <- t2 - t1

    t1 <- proc.time()
    fc_opt <- dowling_cutoff_optimal(y, ddy, t, fs)
    y_opt <- filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")) * y_sd + y_mu
    ddy_opt <- finite_diff(y_opt, t, 2)
    t2 <- proc.time()
    opt_time <- t2 - t1

    t1 <- proc.time()
    fc_pow <- cutoff_power(y_std, 0.9995, fs)
    y_pow <- filtbutter(y_std, butter(2, fc_pow$fc / (fs / 2), type = "low")) * y_sd + y_mu
    ddy_pow <- finite_diff(y_pow, t, 2)
    t2 <- proc.time()
    pow_time <- t2 - t1

    t1 <- proc.time()
    fc <- 15
    fc_fix <- list(
        fc = fc,
        resid = sqrt(sum((y_std - filtbutter(y_std, butter(2, fc / (fs / 2), type = "low")))^2))
    )
    y_fix <- filtbutter(y_std, butter(2, fc_fix$fc / (fs / 2), type = "low")) * y_sd + y_mu
    ddy_fix <- finite_diff(y_fix, t, 2)
    t2 <- proc.time()
    ten_time <- t2 - t1

    saveRDS(list(
        fc_acf = fc_acf,
        fc_resid = fc_resid,
        fc_opt = fc_opt,
        fc_pow = fc_pow,
        fc_fix = fc_fix,
        y_acf = y_acf,
        y_resid = y_resid,
        y_opt = y_opt,
        y_pow = y_pow,
        y_fix = y_fix,
        ddy_acf = ddy_acf,
        ddy_resid = ddy_resid,
        ddy_opt = ddy_opt,
        ddy_pow = ddy_pow,
        ddy_fix = ddy_fix,
        time = list(acf = acf_time, resid = resid_time, opt = opt_time, pow = pow_time, ten = ten_time)
    ), file = paste0(save_folder, "FilterOutput.RDS"))
}
# ============================= P-Spline ================================
if (model_type %in% c("PSpline", "All")) {
    print("    Computing Bayes Models:")
    print("        1: PSpline")
    
    hyperparms <- get_prior_hyperparms_opt(ddB, 1000, "PSpline", pDeg = pDeg)

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
        mod <- jags.model("./mdl/PSpline.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                ord = B$deg,
                xi = hyperparms$xi
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
    clusterExport(cl, list("hyperparms", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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

    # Plot estimates
    par(mfrow = c(1, 1))
    plot(t, ddy, type = "l", main = "PSpline")
    for (i in sample(1:nrow(ddy_samps), 100)) {
        lines(t, ddy_samps[i, ], col = adjustcolor(pal[1], alpha = 0.05))
    }
    lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[1], lwd = 3)

    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0(save_folder, "PSpline.RDS"))
}
# ====================== APS_ind =============================================
if (model_type %in% c("APS_ind", "All")) {
    print("        2: APS_ind")
        hyperparms <- get_prior_hyperparms_opt(ddB, 1000, "APS_ind", pDeg = pDeg)
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
        mod <- jags.model("./mdl/APS_ind.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                d_a0 = hyperparms$a0,
                d_xi0 = hyperparms$xi
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
    clusterExport(cl, list("hyperparms", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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


    # Plot estimates
    par(mfrow = c(1, 1))
    plot(t, ddy, type = "l", main = "APS_ind")
    for (i in sample(1:nrow(ddy_samps), 100)) {
        lines(t, ddy_samps[i, ], col = adjustcolor(pal[2], alpha = 0.05))
    }
    lines(t, apply(ddy_samps, 2, mean), type = "l", col = pal[2], lwd = 3)

    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0(paste0(save_folder, "APS_ind.RDS")))
}

# ====================== AdPSpline_ar =========================================

if (model_type %in% c("APS_ar", "All")) {
    print("        3: APS_ar")
    
    hyperparms <- get_prior_hyperparms_opt(ddB, 1000, "APS_ar", pDeg = pDeg)

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
        mod <- jags.model("./mdl/APS_ar.jags",
            data = list(
                T = length(B$x),
                Q = dim(B$matrix)[2],
                B = B$matrix,
                y = y_std,
                ord = B$deg,
                d_a0 = hyperparms$a0,
                d_xi0 = hyperparms$xi
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "alpha", "xi", "phi", "sigma"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("hyperparms", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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
    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
    }

    # Plot estimates
    par(mfrow = c(1, 1))
    plot(t, ddy, type = "l", main = "AdPSpline_ar")
    for (i in sample(1:nrow(ddy_samps), 100)) {
        lines(t, ddy_samps[i, ], col = adjustcolor(pal[3], alpha = 0.05))
    }
    lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[3], lwd = 3)
    text(0.1, -300, paste0("Min: ", round(min(apply(ddy_samps, 2, mean)), 2)))


    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0(save_folder, "APS_ar.RDS"))
}
# ====================== AdPSpline_spl =========================================

if (model_type %in% c("APS_spl", "All")) {
    
    hyperparms <- get_prior_hyperparms_opt(ddB, 1000, "APS_spl", pDeg = pDeg)
    print("        4: APS_spl")
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
        mod <- jags.model("./mdl/APS_spl.jags",
            data = list(
                T = length(t),
                Q = dim(B$matrix)[2],
                Q_sm = dim(C$matrix)[2],
                B = B$matrix,
                C = C$matrix,
                y = y_std,
                ord = B$deg,
                d_a0 = hyperparms$a0,
                d_xi0 = hyperparms$xi
            ),
            inits = init_fn,
            n.chains = 1,
            n.adapt = nadapt
        )
        coda.samples(mod, c("theta", "tau", "alpha", "xi", "sigma"),
            n.iter = nsampling,
            thin = thin
        )
    }

    t1 <- proc.time()
    cl <- makeCluster(nchains, "SOCK")
    clusterEvalQ(cl, library(rjags))
    clusterExport(cl, list("hyperparms", "MM_likelihood", "mle_pspl", "t", "B", "MM", "C", "y_std", "nadapt", "nsampling", "thin"))
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
    tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
    for (i in 1:nchains) {
        tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
    }

    # Plot estimates
    par(mfrow = c(1, 1))
    plot(t, ddy, type = "l", main = "APS_spl", ylim = c(-450, 150))
    for (i in sample(1:nrow(ddy_samps), 500)) {
        lines(t, ddy_samps[i, ], col = adjustcolor(pal[3], alpha = 0.05))
    }
    lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[3], lwd = 3)
    lines(t, ddy)

    saveRDS(list(
        samples = out,
        y = y_samps,
        ddy = ddy_samps,
        tau = tau_samps,
        time = t2 - t1
    ), file = paste0(save_folder, "APS_spl.RDS"))
}


# ====================== End =========================================
print("ANALYSIS COMPLETE!")
