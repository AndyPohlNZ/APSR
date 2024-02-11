# ==============================================================================
# APS_Fit_JSON.R

# Performs a simulation study of APS given a synthetic dataset or fits APS to
# the dowling dataset. Outputs results as a JSON file.
#
# Datasets condidered:
#   - Heaisine
#   - Blocks
#   - Impact simulation
#   - Dowling dataset

# Models considered:
#   - Filter (ACF, Optimal, Power, Fixed)
#   - P-Spline
#   - APS_ind
#   - APS_spl
#   - APS_ar

# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================

if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("rjags")) install.packages("rjags")
if (!require("coda")) install.packages("coda")
if (!require("snow")) install.packages("snow")
if (!require("argparse")) install.packages("argparse")
if (!require("truncnorm")) install.packages("truncnorm")
if (!require("posterior")) install.packages("posterior")
if (!require("RJSONIO")) install.packages("RJSONIO")
library("argparse")

# ============================= Set Argparser ===================================
parser <- ArgumentParser(
    description = "Run APS analysis on a synthetic dataset."
)
parser$add_argument("-s", "--seed",
    type = "integer",
    default = 1, help = "Set seed for simulation"
)

parser$add_argument("-wd", "--working_directory",
    type = "character",
    default = "./simulations/",
    help = "Set working directory"
)

parser$add_argument("-m", "--model",
    type = "character",
    default = "APS_spl",
    help = "Select model to apply")

parser$add_argument("-a", "--append",
    type = "logical",
    default = FALSE,
    help = "Append to existing results"
)

parser$add_argument("-f", "--file",
    type = "character",
    default = "./results/",
    help = "Set file to save results"
)

parser$add_argument("-mf", "--model_folder",
    type = "character",
    default = "./mdl/",
    help = "Set folder to save results"
)

parser$add_argument("-ds", "--dataset",
    type = "character",
    default = "heavisine",
    help = "Set dataset to use. One of 'heavisine', 'blocks', 'impact' or dowling"
)

parser$add_argument("-sig", "--sigma",
    type = "numeric",
    default = 0.1,
    help = "Set noise level"
)

parser$add_argument("-fs", "--sampling_freq",
    type = "integer",
    default = 512,
    help = "Set sampling frequency"
)

parser$add_argument("-bDeg", "--basis_degree",
    type = "integer",
    default = 5,
    help = "Set degree of spline basis"
)

parser$add_argument("-nk", "--basis_dim",
    type = "integer",
    default = 100,
    help = "Set dimension of Spline basis"
)

parser$add_argument("-pDeg", "--penalty_degree",
    type = "integer",
    default = 2,
    help = "Set penalty degree"
)

parser$add_argument("-nc", "--nchains",
    type = "integer",
    default = 4,
    help = "Set number of MCMC chains"
)

parser$add_argument("-nw", "--nwarmup",
    type = "integer",
    default = 5000,
    help = "Set number of warmup iterations per chain"
)
parser$add_argument("-ns", "--nsamples",
    type = "integer",
    default = 10000,
    help = "Set number of samples per chain"
)

parser$add_argument("-thin", "--thin",
    type = "integer",
    default = 1,
    help = "Set thinning rate"
)

parser$add_argument("-crit", "--criterion",
    type = "numeric",
    default = 100,
    help = "Set criterion for prior 90th quantile of absolute value of second derivative"
)

parser$add_argument("-alp", "--alpha",
    type = "numeric",
    default = 0.1,
    help = "Set alpha level for credible intervals"
)

# ============================= Set parameters ======================================
# get and check arguments
args <- parser$parse_args()

if (args$seed < 1) {
    stop("Seed must be a positive integer.")
} else {
    seed <- args$seed
}
save_file <- args$file
to_append <- args$append
mdl_folder <- args$model_folder

if (args$model %in% c("Filt_acf", "Filt_opt", "Filt_pow", "Filt_fix", "PSpline", "APS_ind", "APS_spl", "APS_ar")) {
    model_type <- args$model
} else {
    stop("Model not supported")
}

if (args$dataset %in% c("heavisine", "blocks", "impact", "dowling")) {
    dataset <- args$dataset
} else {
    stop("Dataset must be one of 'heavisine', 'blocks', 'impact' or 'dowling'")
}

if (args$sigma < 0) {
    stop("Sigma must be a positive number.")
} else {
    sigma <- args$sigma
}

if (args$sampling_freq < 1) {
    stop("Sampling frequency must be a positive integer.")
} else {
    fs <- args$sampling_freq
}

if (args$basis_degree < 3) {
    stop("Basis degree must be a positive integer > 3.")
} else {
    deg <- args$basis_degree
}

if (args$basis_dim < 1) {
    stop("Basis dimension must be a postiive integer.")
} else {
    nK <- args$basis_dim
}

if (!(args$penalty_degree %in% c(2, 3, 4))) {
    stop("Penalty degree must be eithr 2, 3 or 4.")
} else {
    pDeg <- args$penalty_degree
}

if (args$nchains < 1) {
    stop("Number of chains must be a positive integer.")
} else {
    nchains <- args$nchains
}

if (args$nwarmup < 1) {
    stop("Number of warmup iterations must be a positive integer.")
} else {
    nadapt <- args$nwarmup
}

if (args$nsamples < 1) {
    stop("Number of samples must be a positive integer.")
} else {
    nsampling <- args$nsamples
}

if (args$thin < 1) {
    stop("Thinning rate must be a positive integer.")
} else {
    thin <- args$thin
}

if (args$criterion < 0) {
    stop("Criterion must be a positive number.")
} else {
    crit <- args$criterion
}

if (args$alpha < 0 || args$alpha > 1) {
    stop("Alpha must be a number between 0 and 1.")
} else {
    alpha <- args$alpha
}

# ========== Preliminaries ======================================
# Set RNG seed
print(paste0("Setting seed to ", seed))
set.seed(seed)

print(paste0("Setting Working directory to ", args$working_directory))
setwd(args$working_directory)

library("rjags")
library("coda")
library("snow")
library("RColorBrewer")
library("truncnorm")
library("posterior")
library("RJSONIO")

source("./src/library.R")
source("./src/library_analysis.R")
source("./src/library_filter.R")
source("./src/library_simulation.R")

# ============================= Prep Data ======================================
print(paste0("Running AdPSpline analysis on ", dataset, " simulation..."))
# Set time vector
t <- seq(0, 1, length.out = fs)
# prep data
switch(dataset,
    heavisine = {
        print("    Dataset: Heavisine")
        y <- heavisine(t, 200)
        dy <- dheavisine(t, 200)
        ddy <- ddheavisine(t, 200)
    },
    blocks = {
        print("    Dataset: Blocks")
        y <- blocks(t, 100)
        dy <- dblocks(t, 100)
        ddy <- ddblocks(t, 100)
    },
    impact = {
        print("    Dataset: Impact Sim")
        obj <- readRDS("./data/SimulatedImpact.rds")
        y <- obj$y
        dy <- obj$dy
        ddy <- obj$ddy
        t <- obj$time
        fs <- 1 / (t[2] - t[1])
    },
    dowling = {
        print("    Dataset: Dowling")
        y <- as.numeric(read.delim2("./data/Dowling_Disp.txt")[, 1])
        dy <- rep(NA, length(y))
        ddy <- as.numeric(read.delim2("./data/Dowling_Accn.txt")[1:length(y), 1])
        fs <- 512
        t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)
    }
)

if (dataset == "dowling") {
    y_obs <- y
} else {
    # Add noise
    y_obs <- y + rnorm(length(t), 0, sigma)
}
# standardize
y_mu <- mean(y_obs)
y_sd <- sd(y_obs)
y_std <- (y_obs - y_mu) / y_sd
t_std <- t - mean(t)
snr <- var(y) / sigma^2

results_function <- function(model, output) {
    print(paste0("            Generating output for model: ", model))
    if (model %in% c("PSpline", "APS_ind", "APS_ar", "APS_spl")) {
        y_hat <- apply(output$y_hat, 2, median)
        y_lwr <- apply(output$y_hat, 2, quantile, probs = alpha / 2)
        y_upr <- apply(output$y_hat, 2, quantile, probs = 1 - alpha / 2)
        dy_hat <- apply(output$dy_hat, 2, median)
        dy_lwr <- apply(output$dy_hat, 2, quantile, probs = alpha / 2)
        dy_upr <- apply(output$dy_hat, 2, quantile, probs = 1 - alpha / 2)
        ddy_hat <- apply(output$ddy_hat, 2, median)
        ddy_lwr <- apply(output$ddy_hat, 2, quantile, probs = alpha / 2)
        ddy_upr <- apply(output$ddy_hat, 2, quantile, probs = 1 - alpha / 2)

        sigma_hat <- rep(median(output$sigma), length(t))
        sigma_lwr <- rep(quantile(output$sigma, probs = alpha / 2), length(t))
        sigma_upr <- rep(quantile(output$sigma, probs = 1 - alpha / 2), length(t))

        ess_y <- effectiveSize(output$y_hat)
        ess_dy <- effectiveSize(output$dy_hat)
        ess_ddy <- effectiveSize(output$ddy_hat)

        ess_sigma <- effectiveSize(output$sigma)
        rhat_sigma <- rhat(matrix(output$sigma, nrow = nsampling, ncol = nchains))

        tmp <- array(output$y_hat, dim = c(nsampling, nchains, length(t)))
        rhat_y <- apply(tmp, 3, rhat)
        tmp <- array(output$dy_hat, dim = c(nsampling, nchains, length(t)))
        rhat_dy <- apply(tmp, 3, rhat)
        tmp <- array(output$ddy_hat, dim = c(nsampling, nchains, length(t)))
        rhat_ddy <- apply(tmp, 3, rhat)

        t_tau <- tau_hat <- tau_lwr <- tau_upr <- rep(NA, length(t))
        ess_tau <- rhat_tau <- rep(NA, length(t))

        if (model == "PSpline") {
            tau_hat <- mean(output$tau_hat)
            tau_lwr <- quantile(output$tau_hat, probs = alpha / 2)
            tau_upr <- quantile(output$tau_hat, probs = 1 - alpha / 2)
            ess_tau <- effectiveSize(output$tau_hat)
            rhat_tau <- rhat(matrix(output$tau_hat, nrow = nsampling, ncol = nchains))
        } else {
            tau_hat <- apply(output$tau_hat, 2, median)
            tau_lwr <- apply(output$tau_hat, 2, quantile, probs = alpha / 2)
            tau_upr <- apply(output$tau_hat, 2, quantile, probs = 1 - alpha / 2)
            ess_tau <- effectiveSize(output$tau_hat)
            tmp <- array(output$tau_hat, dim = c(nsampling, nchains, dim(output$tau_hat)[2]))
            rhat_tau <- apply(tmp, 3, rhat)
        }
    } else {
        y_hat <- output$y_hat
        dy_hat <- output$dy_hat
        ddy_hat <- output$ddy_hat
    }


    if (model %in% c("PSpline", "APS_ind", "APS_ar", "APS_spl")) {
        output_ls <- list(
            model = model,
            dataset = dataset,
            sigma = sigma,
            seed = seed,
            basis_parm = list(nk = nK, deg = deg, pdeg = pDeg),
            true_vals = list(t = t, y = y, y_obs = y_obs, dy = dy, ddy = ddy),
            estimates = list(
                y_hat = y_hat, dy_hat = dy_hat, ddy_hat = ddy_hat,
                tau_hat = tau_hat, sigma_hat = sigma_hat
            ),
            CI = list(
                y_lwr = as.vector(y_lwr), y_upr = as.vector(y_upr),
                dy_lwr = as.vector(dy_lwr), dy_upr = as.vector(dy_upr),
                ddy_lwr = as.vector(ddy_lwr), ddy_upr = as.vector(ddy_upr),
                tau_lwr = as.vector(tau_lwr), tau_upr = as.vector(tau_upr),
                sigma_lwr = as.vector(sigma_lwr), sigma_upr = as.vector(sigma_upr)
            ),
            ess = list(
                ess_y = as.vector(ess_y), ess_dy = as.vector(ess_dy), ess_ddy = as.vector(ess_ddy),
                ess_tau = as.vector(ess_tau), ess_sigma = as.vector(ess_sigma)
            ),
            rhat = list(
                rhat_y = as.vector(rhat_y), rhat_dy = as.vector(rhat_dy), rhat_ddy = as.vector(rhat_ddy),
                rhat_tau = as.vector(rhat_tau), rhat_sigma = as.vector(rhat_sigma)
            ),
            run_time = as.numeric(output$run_time["elapsed"]),
            call = args
        )
    } else {
        output_ls <- list(
            model = model,
            dataset = dataset,
            sigma = sigma,
            seed = seed,
            basis_parm = NA,
            true_vals = list(t = t, y = y, y_obs = y_obs, dy = dy, ddy = ddy),
            estimates = list(y_hat = as.vector(y_hat), dy_hat = as.vector(dy_hat), ddy_hat = as.vector(ddy_hat)),
            CI = NA,
            ess = NA,
            rhat = NA,
            run_time = as.numeric(output$run_time["elapsed"]),
            call = args
        )
    }

    json_output <- toJSON(output_ls)
    write(json_output, file = save_file, append = to_append)
}
# ====================== Initialize spline basis ==============================
print("    Prepping spline basis...")
# Main Basis
B <- basis(t, nK, deg = deg)
dB <- dbasis(B, 1)
ddB <- dbasis(B, 2)
P <- penalty(dim(B$matrix)[2], pDeg)
MM <- basis_mm(P)
nK_sm <- floor(nK / 5)
bDeg_sm <- 3
C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)

# ============================= Compute filters ================================
switch(model_type,
    Filt_acf = {
        print("    Computing filters...")
        t1 <- proc.time()
        fc_acf <- cutoff_acf(y_std, fs)
        y_acf <- filtbutter(y_std, butter(1, fc_acf$fc / (fs / 2), type = "low")) * y_sd + y_mu
        dy_acf <- finite_diff(y_acf, t, 1)
        ddy_acf <- finite_diff(y_acf, t, 2)
        t2 <- proc.time()
        acf_time <- t2 - t1
        output <- list(
            y_hat = y_acf,
            dy_hat = dy_acf,
            ddy_hat = ddy_acf,
            run_time = acf_time,
            fc = fc_acf
        )
        results_function("Filt_acf", output)
    },
    Filt_opt = {
        if (dataset == "dowling") {
            t1 <- proc.time()
            fc_opt <- cutoff_optimal(y_obs, ddy, 2, t, fs)
            y_opt <- filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")) * y_sd + y_mu
            dy_opt <- finite_diff(filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")), t, 1) * y_sd
            ddy_opt <- finite_diff(filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")), t, 2) * y_sd
            t2 <- proc.time()
            opt_time_y <- t2 - t1
            opt_time_dy <- opt_time_y
            opt_time_ddy <- opt_time_y
        } else {
            t1 <- proc.time()
            fc_opt <- cutoff_optimal(y_obs, y, 0, t, fs)
            y_opt <- filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")) * y_sd + y_mu
            t2 <- proc.time()
            opt_time_y <- t2 - t1

            t1 <- proc.time()
            fc_opt <- cutoff_optimal(y_obs, dy, 1, t, fs)
            dy_opt <- finite_diff(filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")), t, 1) * y_sd
            t2 <- proc.time()
            opt_time_dy <- t2 - t1

            t1 <- proc.time()
            fc_opt <- cutoff_optimal(y_obs, ddy, 2, t, fs)
            ddy_opt <- finite_diff(filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")), t, 2) * y_sd
            t2 <- proc.time()
            opt_time_ddy <- t2 - t1
        }

        output <- list(
            y_hat = y_opt,
            dy_hat = dy_opt,
            ddy_hat = ddy_opt,
            run_time = c(opt_time_y, opt_time_dy, opt_time_ddy),
            fc = fc_opt
        )
        results_function("Filt_opt", output)
    },
    Filt_pow = {
        t1 <- proc.time()
        fc_pow <- cutoff_power(y_std, 0.9995, fs)
        y_pow <- filtbutter(y_std, butter(2, fc_pow$fc / (fs / 2), type = "low")) * y_sd + y_mu
        dy_pow <- finite_diff(y_pow, t, 1)
        ddy_pow <- finite_diff(y_pow, t, 2)
        t2 <- proc.time()
        pow_time <- t2 - t1
        output <- list(
            y_hat = y_pow,
            dy_hat = dy_pow,
            ddy_hat = ddy_pow,
            run_time = pow_time,
            fc = fc_pow
        )
        results_function("Filt_pow", output)
    },
    Filt_fix = {
        t1 <- proc.time()
        fc <- 15
        fc_fix <- list(
            fc = fc,
            resid = sqrt(sum((y_std - filtbutter(y_std, butter(2, fc / (fs / 2), type = "low")))^2))
        )
        y_fix <- filtbutter(y_std, butter(2, fc_fix$fc / (fs / 2), type = "low")) * y_sd + y_mu
        dy_fix <- finite_diff(y_fix, t, 1)
        ddy_fix <- finite_diff(y_fix, t, 2)
        t2 <- proc.time()
        ten_time <- t2 - t1
        output <- list(
            y_hat = y_fix,
            dy_hat = dy_fix,
            ddy_hat = ddy_fix,
            run_time = ten_time,
            fc = fc_fix
        )
        results_function("Filt_fix", output)
    },
    # ============================= P-Spline ================================
    PSpline = {
        print("         PSpline")

        hyperparms <- get_prior_hyperparms_opt(ddB, crit, "PSpline", pDeg = pDeg)

        print("            Fitting PSpline")
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
            mod <- jags.model(paste0(mdl_folder, "/PSpline.jags"),
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
        clusterExport(cl, list("hyperparms", "mdl_folder", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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
        dy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        for (i in 1:nrow(theta_samps)) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
        }

        tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
        for (i in 1:nchains) {
            tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "tau"]
        }

        sigma_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
        for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "sigma"]
        }

        output <- list(
            y_hat = y_samps,
            dy_hat = dy_samps,
            ddy_hat = ddy_samps,
            tau_hat = tau_samps,
            sigma_hat = sigma_samps,
            run_time = t2 - t1
        )
        results_function("PSpline", output)
    },
    APS_ind = {
        # ====================== APS_ind =============================================
        print("        APS_ind")

        hyperparms <- get_prior_hyperparms_opt(ddB, crit, "APS_ind", pDeg = pDeg)
        print("             Fitting APS_ind")
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
            mod <- jags.model(paste0(mdl_folder, "/APS_ind.jags"),
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
        clusterExport(cl, list("hyperparms", "mdl_folder", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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
        dy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        for (i in 1:nrow(theta_samps)) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
        }
        tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2] - pDeg)
        for (i in 1:nchains) {
            tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(dim(B$matrix)[2] - pDeg), "]")]
        }

        sigma_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
        for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "sigma"]
        }

        output <- list(
            y_hat = y_samps,
            dy_hat = dy_samps,
            ddy_hat = ddy_samps,
            tau_hat = tau_samps,
            sigma_hat = sigma_samps,
            run_time = t2 - t1
        )


        results_function("APS_ind", output)
    },
    APS_ar = {
        # ====================== AdPSpline_ar =========================================

        print("        APS_ar")

        hyperparms <- get_prior_hyperparms_opt(ddB, crit, "APS_ar", pDeg = pDeg)
        print("            Fitting APS_ar")
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
            mod <- jags.model(paste0(mdl_folder, "APS_ar.jags"),
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
        clusterExport(cl, list("hyperparms", "mdl_folder", "MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
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
        dy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        for (i in 1:nrow(theta_samps)) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
        }
        tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
        for (i in 1:nchains) {
            tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
        }

        sigma_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
        for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "sigma"]
        }

        output <- list(
            y_hat = y_samps,
            dy_hat = dy_samps,
            ddy_hat = ddy_samps,
            tau_hat = tau_samps,
            sigma_hat = sigma_samps,
            run_time = t2 - t1
        )

        results_function("APS_ar", output)
    },
    APS_spl = {
        # ====================== AdPSpline_spl =========================================

        print("        APS_spl")

        hyperparms <- get_prior_hyperparms_opt(ddB, crit, "APS_spl", pDeg = pDeg)
        print("            Fitting APS_spl")
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
            mod <- jags.model(paste0(mdl_folder, "/APS_spl.jags"),
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
        clusterExport(cl, list("hyperparms", "mdl_folder", "MM_likelihood", "mle_pspl", "t", "B", "MM", "C", "y_std", "nadapt", "nsampling", "thin"))
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
        dy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
        for (i in 1:nrow(theta_samps)) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
        }
        tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
        for (i in 1:nchains) {
            tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
        }

        sigma_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
        for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "sigma"]
        }

        output <- list(
            y_hat = y_samps,
            dy_hat = dy_samps,
            ddy_hat = ddy_samps,
            tau_hat = tau_samps,
            sigma_hat = sigma_samps,
            run_time = t2 - t1
        )

        results_function("APS_spl", output)
    }
)

# ==========================================================================
# End
print("Analysis Complete")
# ==========================================================================







#
#
# if (!require("rjags")) install.packages("rjags")
# if (!require("coda")) install.packages("coda")
# if (!require("snow")) install.packages("snow")
#
# library("rjags")
# library("coda")
# library("snow")
#
# source("./src/library.R")
# source("./src/library_filter.R")
#
# # ============================= Preliminaries ======================================
#
#
# seed <- 123
# model_type <- "All" # One of "FIlter", "PSpline", "APS_ind", "APS_spl", "APS_ar" or "All"
# save_folder <- "./results/"
#
# # Set RNG seed
# set.seed(seed)
#
# # Set MCMC parameters
# nchains <- 4 # Number of independent chains
# nadapt <- 10000 # Size of adaption phase
# nsampling <- 50000 # Number of samples per chain
# thin <- 1 # Thinning rate per chain.
#
# # ============================= Load Data ======================================
#
# print("Running AdPSpline analysis on Dowling data...")
# print("    Prepping data...")
#
# # Process data
# y <- as.numeric(read.delim2("data/Dowling_Disp.txt")[, 1])
# ddy <- as.numeric(read.delim2("data/Dowling_Accn.txt")[1:length(y), 1])
#
# # Set time vector
# fs <- 512 # Sampling frequency of Dowling Dataset
# t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)
#
# # standardize
# y_mu <- mean(y)
# y_sd <- sd(y)
# y_std <- (y - y_mu) / y_sd
#
# # ====================== Initialize spline basis ==============================================
# print("    Prepping spline basis...")
#
# # Main Basis
# deg <- 5
# pDeg <- 2
# nK <- 100
# B <- basis(t, nK, deg = deg)
# ddB <- dbasis(B, 2)
# P <- penalty(dim(B$matrix)[2], pDeg)
# MM <- basis_mm(P)
#
# # Basis for tau.
# nK_sm <- floor(nK / 5)
# bDeg_sm <- 3
# C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)
#
# # ============================= Compute filters ================================
# if (model_type %in% c("Filter", "All")) {
#     print("    Computing filters...")
#     t1 <- proc.time()
#     fc_acf <- cutoff_acf(y_std, fs)
#     y_acf <- filtbutter(y_std, butter(2, fc_acf$fc / (fs / 2), type = "low")) * y_sd + y_mu
#     ddy_acf <- finite_diff(y_acf, t, 2)
#     t2 <- proc.time()
#     acf_time <- t2 - t1
#
#     t1 <- proc.time()
#     fc_resid <- cutoff_residual(y_std, fs)
#     y_resid <- filtbutter(y_std, butter(2, fc_resid$fc / (fs / 2), type = "low")) * y_sd + y_mu
#     ddy_resid <- finite_diff(y_resid, t, 2)
#     t2 <- proc.time()
#     resid_time <- t2 - t1
#
#     t1 <- proc.time()
#     fc_opt <- dowling_cutoff_optimal(y, ddy, t, fs)
#     y_opt <- filtbutter(y_std, butter(2, fc_opt$fc / (fs / 2), type = "low")) * y_sd + y_mu
#     ddy_opt <- finite_diff(y_opt, t, 2)
#     t2 <- proc.time()
#     opt_time <- t2 - t1
#
#     t1 <- proc.time()
#     fc_pow <- cutoff_power(y_std, 0.9995, fs)
#     y_pow <- filtbutter(y_std, butter(2, fc_pow$fc / (fs / 2), type = "low")) * y_sd + y_mu
#     ddy_pow <- finite_diff(y_pow, t, 2)
#     t2 <- proc.time()
#     pow_time <- t2 - t1
#
#     t1 <- proc.time()
#     fc <- 15
#     fc_fix <- list(
#         fc = fc,
#         resid = sqrt(sum((y_std - filtbutter(y_std, butter(2, fc / (fs / 2), type = "low")))^2))
#     )
#     y_fix <- filtbutter(y_std, butter(2, fc_fix$fc / (fs / 2), type = "low")) * y_sd + y_mu
#     ddy_fix <- finite_diff(y_fix, t, 2)
#     t2 <- proc.time()
#     ten_time <- t2 - t1
#
#     saveRDS(list(
#         fc_acf = fc_acf,
#         fc_resid = fc_resid,
#         fc_opt = fc_opt,
#         fc_pow = fc_pow,
#         fc_fix = fc_fix,
#         y_acf = y_acf,
#         y_resid = y_resid,
#         y_opt = y_opt,
#         y_pow = y_pow,
#         y_fix = y_fix,
#         ddy_acf = ddy_acf,
#         ddy_resid = ddy_resid,
#         ddy_opt = ddy_opt,
#         ddy_pow = ddy_pow,
#         ddy_fix = ddy_fix,
#         time = list(acf = acf_time, resid = resid_time, opt = opt_time, pow = pow_time, ten = ten_time)
#     ), file = paste0(save_folder, "FilterOutput.RDS"))
# }
# # ============================= P-Spline ================================
# if (model_type %in% c("PSpline", "All")) {
#     print("    Computing Bayes Models:")
#     print("        1: PSpline")
#     coda.samples.wrapper <- function(j) {
#         init_fn <- function() {
#             mle <- mle_pspl(B, MM, y_std)
#             theta <- as.numeric(mle$theta_hat)
#             list(
#                 theta = theta,
#                 .RNG.name = "base::Wichmann-Hill",
#                 .RNG.seed = j
#             )
#         }
#         mod <- jags.model("./mdl/PSpline.jags",
#             data = list(
#                 T = length(t),
#                 Q = dim(B$matrix)[2],
#                 B = B$matrix,
#                 y = y_std,
#                 ord = B$deg,
#                 xi = 0.002
#             ),
#             inits = init_fn,
#             n.chains = 1,
#             n.adapt = nadapt
#         )
#         coda.samples(mod, c("theta", "tau", "sigma"),
#             n.iter = nsampling,
#             thin = thin
#         )
#     }
#
#     t1 <- proc.time()
#     cl <- makeCluster(nchains, "SOCK")
#     clusterEvalQ(cl, library(rjags))
#     clusterExport(cl, list("MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
#     out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
#     for (i in 1:length(out)) {
#         out[[i]] <- out[[i]][[1]]
#     }
#     class(out) <- "mcmc.list"
#     stopCluster(cl)
#     t2 <- proc.time()
#
#     # compute y and ddy
#     nchains <- length(out)
#     nsamps <- dim(out[[1]])[1]
#     theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
#     for (i in 1:nchains) {
#         theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
#     }
#     y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     for (i in 1:nrow(theta_samps)) {
#         y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
#         ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
#     }
#
#     tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = 1)
#     for (i in 1:nchains) {
#         tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, "tau"]
#     }
#
#     # Plot estimates
#     par(mfrow = c(1, 1))
#     plot(t, ddy, type = "l", main = "PSpline")
#     for (i in sample(1:nrow(ddy_samps), 100)) {
#         lines(t, ddy_samps[i, ], col = adjustcolor(pal[1], alpha = 0.05))
#     }
#     lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[1], lwd = 3)
#
#     saveRDS(list(
#         samples = out,
#         y = y_samps,
#         ddy = ddy_samps,
#         tau = tau_samps,
#         time = t2 - t1
#     ), file = paste0(save_folder, "PSpline.RDS"))
# }
# # ====================== APS_ind =============================================
# if (model_type %in% c("APS_ind", "All")) {
#     print("        2: APS_ind")
#     coda.samples.wrapper <- function(j) {
#         init_fn <- function() {
#             mle <- mle_pspl(B, MM, y_std)
#             theta <- as.numeric(mle$theta_hat)
#             list(
#                 theta = theta,
#                 .RNG.name = "base::Wichmann-Hill",
#                 .RNG.seed = j
#             )
#         }
#         mod <- jags.model("./mdl/APS_ind.jags",
#             data = list(
#                 T = length(t),
#                 Q = dim(B$matrix)[2],
#                 B = B$matrix,
#                 y = y_std,
#                 d_a0 = 0.002,
#                 d_xi0 = 1.65
#             ),
#             inits = init_fn,
#             n.chains = 1,
#             n.adapt = nadapt
#         )
#         coda.samples(mod, c("theta", "tau", "sigma"),
#             n.iter = nsampling,
#             thin = thin
#         )
#     }
#
#     t1 <- proc.time()
#     cl <- makeCluster(nchains, "SOCK")
#     clusterEvalQ(cl, library(rjags))
#     clusterExport(cl, list("MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
#     out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
#     for (i in 1:length(out)) {
#         out[[i]] <- out[[i]][[1]]
#     }
#     class(out) <- "mcmc.list"
#     stopCluster(cl)
#     t2 <- proc.time()
#
#     # compute y and ddy
#     nchains <- length(out)
#     nsamps <- dim(out[[1]])[1]
#     theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
#     for (i in 1:nchains) {
#         theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
#     }
#     y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     for (i in 1:nrow(theta_samps)) {
#         y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
#         ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
#     }
#     tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2] - pDeg)
#     for (i in 1:nchains) {
#         tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(dim(B$matrix)[2] - pDeg), "]")]
#     }
#
#
#     # Plot estimates
#     par(mfrow = c(1, 1))
#     plot(t, ddy, type = "l", main = "APS_ind")
#     for (i in sample(1:nrow(ddy_samps), 100)) {
#         lines(t, ddy_samps[i, ], col = adjustcolor(pal[2], alpha = 0.05))
#     }
#     lines(t, apply(ddy_samps, 2, mean), type = "l", col = pal[2], lwd = 3)
#
#     saveRDS(list(
#         samples = out,
#         y = y_samps,
#         ddy = ddy_samps,
#         tau = tau_samps,
#         time = t2 - t1
#     ), file = paste0(paste0(save_folder, "APS_ind.RDS")))
# }
#
# # ====================== AdPSpline_ar =========================================
#
# if (model_type %in% c("APS_ar", "All")) {
#     print("        3: APS_ar")
#     coda.samples.wrapper <- function(j) {
#         init_fn <- function() {
#             mle <- mle_pspl(B, MM, y_std)
#             theta <- as.numeric(mle$theta_hat)
#             list(
#                 theta = theta,
#                 .RNG.name = "base::Wichmann-Hill",
#                 .RNG.seed = j
#             )
#         }
#         mod <- jags.model("./mdl/APS_ar.jags",
#             data = list(
#                 T = length(B$x),
#                 Q = dim(B$matrix)[2],
#                 B = B$matrix,
#                 y = y_std,
#                 ord = B$deg,
#                 d_a0 = 0.002,
#                 d_xi0 = 1.25
#             ),
#             inits = init_fn,
#             n.chains = 1,
#             n.adapt = nadapt
#         )
#         coda.samples(mod, c("theta", "tau", "alpha", "xi", "phi", "sigma"),
#             n.iter = nsampling,
#             thin = thin
#         )
#     }
#
#     t1 <- proc.time()
#     cl <- makeCluster(nchains, "SOCK")
#     clusterEvalQ(cl, library(rjags))
#     clusterExport(cl, list("MM_likelihood", "mle_pspl", "t", "B", "MM", "y_std", "nadapt", "nsampling", "thin"))
#     out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
#     for (i in 1:length(out)) {
#         out[[i]] <- out[[i]][[1]]
#     }
#     class(out) <- "mcmc.list"
#     stopCluster(cl)
#     t2 <- proc.time()
#
#     # compute y and ddy
#     nchains <- length(out)
#     nsamps <- dim(out[[1]])[1]
#     theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
#     for (i in 1:nchains) {
#         theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
#     }
#     y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     for (i in 1:nrow(theta_samps)) {
#         y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
#         ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
#     }
#     tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
#     for (i in 1:nchains) {
#         tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
#     }
#
#     # Plot estimates
#     par(mfrow = c(1, 1))
#     plot(t, ddy, type = "l", main = "AdPSpline_ar")
#     for (i in sample(1:nrow(ddy_samps), 100)) {
#         lines(t, ddy_samps[i, ], col = adjustcolor(pal[3], alpha = 0.05))
#     }
#     lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[3], lwd = 3)
#     text(0.1, -300, paste0("Min: ", round(min(apply(ddy_samps, 2, mean)), 2)))
#
#
#     saveRDS(list(
#         samples = out,
#         y = y_samps,
#         ddy = ddy_samps,
#         tau = tau_samps,
#         time = t2 - t1
#     ), file = paste0(save_folder, "APS_ar.RDS"))
# }
# # ====================== AdPSpline_spl =========================================
#
# if (model_type %in% c("APS_spl", "All")) {
#     print("        4: APS_spl")
#     coda.samples.wrapper <- function(j) {
#         init_fn <- function() {
#             mle <- mle_pspl(B, MM, y_std)
#             theta <- as.numeric(mle$theta_hat)
#             list(
#                 theta = theta,
#                 .RNG.name = "base::Wichmann-Hill",
#                 .RNG.seed = j
#             )
#         }
#         mod <- jags.model("./mdl/APS_spl.jags",
#             data = list(
#                 T = length(t),
#                 Q = dim(B$matrix)[2],
#                 Q_sm = dim(C$matrix)[2],
#                 B = B$matrix,
#                 C = C$matrix,
#                 y = y_std,
#                 ord = B$deg,
#                 d_a0 = 0.002,
#                 d_xi0 = 1.15
#             ),
#             inits = init_fn,
#             n.chains = 1,
#             n.adapt = nadapt
#         )
#         coda.samples(mod, c("theta", "tau", "alpha", "xi", "sigma"),
#             n.iter = nsampling,
#             thin = thin
#         )
#     }
#
#     t1 <- proc.time()
#     cl <- makeCluster(nchains, "SOCK")
#     clusterEvalQ(cl, library(rjags))
#     clusterExport(cl, list("MM_likelihood", "mle_pspl", "t", "B", "MM", "C", "y_std", "nadapt", "nsampling", "thin"))
#     out <- clusterApply(cl, 1:nchains, coda.samples.wrapper)
#     for (i in 1:length(out)) {
#         out[[i]] <- out[[i]][[1]]
#     }
#     class(out) <- "mcmc.list"
#     stopCluster(cl)
#     t2 <- proc.time()
#
#     # compute y and ddy
#     nchains <- length(out)
#     nsamps <- dim(out[[1]])[1]
#     theta_samps <- matrix(NA, nrow = nsamps * nchains, ncol = dim(B$matrix)[2])
#     for (i in 1:nchains) {
#         theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
#     }
#     y_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     ddy_samps <- matrix(NA, nrow = nrow(theta_samps), ncol = length(t))
#     for (i in 1:nrow(theta_samps)) {
#         y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
#         ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
#     }
#     tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = ncol(B$matrix) - 2)
#     for (i in 1:nchains) {
#         tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("tau[", 1:(ncol(B$matrix) - 2), "]")]
#     }
#
#     # Plot estimates
#     par(mfrow = c(1, 1))
#     plot(t, ddy, type = "l", main = "APS_spl", ylim = c(-450, 150))
#     for (i in sample(1:nrow(ddy_samps), 500)) {
#         lines(t, ddy_samps[i, ], col = adjustcolor(pal[3], alpha = 0.05))
#     }
#     lines(t, apply(ddy_samps, 2, median), type = "l", col = pal[3], lwd = 3)
#     lines(t, ddy)
#
#     saveRDS(list(
#         samples = out,
#         y = y_samps,
#         ddy = ddy_samps,
#         tau = tau_samps,
#         time = t2 - t1
#     ), file = paste0(save_folder, "APS_spl.RDS"))
# }
#
#
# # ====================== End =========================================
# print("ANALYSIS COMPLETE!")
