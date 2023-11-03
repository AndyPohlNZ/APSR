# Run analysis of Impact dataset

if (!require("coda")) install.packages("coda")
if (!require("bayesplot")) install.packages("bayesplot")
if (!require("latex2exp")) install.packages("latex2exp")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("xtable")) install.packages("xtable")


library("bayesplot")
library("coda")
library("latex2exp")
library("RColorBrewer")
library("xtable")

source("./src/library.R")
source("./src/library_analysis.R")

pal <- brewer.pal(5, "Set1") # Set color pallette for plots


# ============================= Load Data ======================================

print("Running AdPSpline analysis on Dowling data...")
print("    Prepping data...")

# Process data
y <- as.numeric(read.delim2("data/Dowling_Disp.txt")[, 1])
ddy <- as.numeric(read.delim2("data/Dowling_Accn.txt")[1:length(y), 1])
fs <- 512
t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)

# Main Basis
nK <- 100
deg <- 5
pDeg <- 2
B <- basis(t, nK, deg = deg)
ddB <- dbasis(B, 2)
P <- penalty(dim(B$matrix)[2], pDeg)

# Basis for tau.
nK_sm <- floor(nK / 5)
bDeg_sm <- 3
C <- basis(seq(0, 1, length.out = dim(B$matrix)[2] - P$deg), nK_sm, deg = bDeg_sm)

# ============================= Load Fitted models ======================================

Filter_fit <- readRDS("./results/FilterOutput.RDS")
PSpline_fit <- readRDS("./results/PSpline.RDS")
APS_ind_fit <- readRDS("./results/APS_ind.RDS")
APS_spl_fit <- readRDS("./results/APS_spl.RDS")
APS_ar_fit <- readRDS("./results/APS_ar.RDS")

# ============================= Compute Bayes estimates ======================================
PSpline_est <- bayes_est(PSpline_fit)
APS_ind_est <- bayes_est(APS_ind_fit)
APS_spl_est <- bayes_est(APS_spl_fit)
APS_ar_est <- bayes_est(APS_ar_fit)

# Compute diagonstics for MCMC chains
PSpline_mcmcdiag <- ess_est(PSpline_fit, B, ddB)
APS_ind_mcmcdiag <- ess_est(APS_ind_fit, B, ddB)
APS_spl_mcmcdiag <- ess_est(APS_spl_fit, B, ddB)
APS_ar_mcmcdiag <- ess_est(APS_ar_fit, B, ddB)


# Specify number of models
nmodels <- 8

# ============================= Compute RMSE ======================================
rmse_results <- data.frame(type = rep("", nmodels), name = rep("", nmodels), y_rmse = as.numeric(rep(NA, nmodels)), ddy_rmse = as.numeric(rep(NA, nmodels)))
rmse_results[1, ] <- c("Filter", "15", rmse(y, Filter_fit$y_fix), rmse(ddy, Filter_fit$ddy_fix))

rmse_results[2, ] <- c("Filter", "ACF", rmse(y, Filter_fit$y_acf), rmse(ddy, Filter_fit$ddy_acf))
rmse_results[3, ] <- c("Filter", "Resid", rmse(y, Filter_fit$y_resid), rmse(ddy, Filter_fit$ddy_resid))
rmse_results[4, ] <- c("Filter", "Optimal", rmse(y, Filter_fit$y_opt), rmse(ddy, Filter_fit$ddy_opt))

rmse_results[5, ] <- c("PSpline", "Regular", rmse(y, PSpline_est$y_hat), rmse(ddy, PSpline_est$ddy_hat))
rmse_results[6, ] <- c("APS_ind", "Regular", rmse(y, APS_ind_est$y_hat), rmse(ddy, APS_ind_est$ddy_hat))

rmse_results[7, ] <- c("APS_spl", "Regular", rmse(y, APS_spl_est$y_hat), rmse(ddy, APS_spl_est$ddy_hat))

rmse_results[8, ] <- c("APS_ar", "Regular", rmse(y, APS_ar_est$y_hat), rmse(ddy, APS_ar_est$ddy_hat))

rmse_results[, 3] <- round(as.numeric(rmse_results$y_rmse), 4)
rmse_results[, 4] <- round(as.numeric(rmse_results$ddy_rmse), 3)

print(xtable(rmse_results, digits = 4), include.rownames = FALSE)

# =============================  Effective Sample Size  ======================================
ess_results <- data.frame(
    type = rep("", nmodels), name = rep("", nmodels),
    y_ess = rep(NA, nmodels), y_ess_sd = rep(NA, nmodels),
    ddy_ess = rep(NA, nmodels), ddy_ess_sd = rep(NA, nmodels),
    tau_ess = rep(NA, nmodels), tau_ess_sd = rep(NA, nmodels),
    sigma_ess = rep(NA, nmodels)
)

ess_results[1, ] <- c("Filter", "15", rep(NA, 7))
ess_results[2, ] <- c("Filter", "ACF", rep(NA, 7))
ess_results[3, ] <- c("Filter", "Resid", rep(NA, 7))
ess_results[4, ] <- c("Filter", "Optimal", rep(NA, 7))

ess_results[5, ] <- c(
    "PSpline", "Regular", mean(PSpline_mcmcdiag$ess_y), sd(PSpline_mcmcdiag$ess_y),
    mean(PSpline_mcmcdiag$ess_ddy), sd(PSpline_mcmcdiag$ess_ddy),
    PSpline_mcmcdiag$ess_tau, NA,
    PSpline_mcmcdiag$ess_sigma
)

ess_results[6, ] <- c(
    "APS_ind", "Regular", mean(APS_ind_mcmcdiag$ess_y), sd(APS_ind_mcmcdiag$ess_y),
    mean(APS_ind_mcmcdiag$ess_ddy), sd(APS_ind_mcmcdiag$ess_ddy),
    APS_ind_mcmcdiag$ess_tau, NA,
    APS_ind_mcmcdiag$ess_sigma
)

ess_results[7, ] <- c(
    "APS_spl", "Regular", mean(APS_spl_mcmcdiag$ess_y), sd(APS_spl_mcmcdiag$ess_y),
    mean(APS_spl_mcmcdiag$ess_ddy), sd(APS_spl_mcmcdiag$ess_ddy),
    mean(APS_spl_mcmcdiag$ess_tau), sd(APS_spl_mcmcdiag$ess_tau),
    APS_spl_mcmcdiag$ess_sigma
)


ess_results[8, ] <- c(
    "APS_ar", "Regular", mean(APS_ar_mcmcdiag$ess_y), sd(APS_ar_mcmcdiag$ess_y),
    mean(APS_ar_mcmcdiag$ess_ddy), sd(APS_ar_mcmcdiag$ess_ddy),
    mean(APS_ar_mcmcdiag$ess_tau), sd(APS_ar_mcmcdiag$ess_tau),
    APS_ar_mcmcdiag$ess_sigma
)

for (i in 3:ncol(ess_results)) {
    ess_results[, i] <- as.numeric(ess_results[, i])
}
print(xtable(ess_results[, c(1, 2, 3, 5, 6, 7, 8)], digits = 0), include.rownames = FALSE)

# ============================= Rhat  ======================================

rhat_results <- data.frame(
    type = rep("", nmodels), name = rep("", nmodels),
    y_rhat = rep(NA, nmodels), y_rhat_sd = rep(NA, nmodels),
    ddy_rhat = rep(NA, nmodels), ddy_rhat_sd = rep(NA, nmodels),
    tau_rhat = rep(NA, nmodels), tau_rhat_sd = rep(NA, nmodels),
    sigma_rhat = rep(NA, nmodels)
)

rhat_results[1, ] <- c("Filter", "15", rep(NA, 7))
rhat_results[2, ] <- c("Filter", "ACF", rep(NA, 7))
rhat_results[3, ] <- c("Filter", "Resid", rep(NA, 7))
rhat_results[4, ] <- c("Filter", "Optimal", rep(NA, 7))

rhat_results[5, ] <- c(
    "PSpline", "Regular", mean(PSpline_mcmcdiag$rhat_y$psrf[, 1]), sd(PSpline_mcmcdiag$rhat_y$psrf[, 1]),
    mean(PSpline_mcmcdiag$rhat_ddy$psrf[, 1]), sd(PSpline_mcmcdiag$rhat_ddy$psrf[, 1]),
    mean(PSpline_mcmcdiag$rhat_tau$psrf[, 1]), sd(PSpline_mcmcdiag$rhat_tau$psrf[, 1]),
    mean(PSpline_mcmcdiag$rhat_sigma$psrf[, 1])
)

rhat_results[6, ] <- c(
    "APS_ind", "Regular", mean(APS_ind_mcmcdiag$rhat_y$psrf[, 1]), sd(APS_ind_mcmcdiag$rhat_y$psrf[, 1]),
    mean(APS_ind_mcmcdiag$rhat_ddy$psrf[, 1]), sd(APS_ind_mcmcdiag$rhat_ddy$psrf[, 1]),
    mean(APS_ind_mcmcdiag$rhat_tau$psrf[, 1]), sd(APS_ind_mcmcdiag$rhat_tau$psrf[, 1]),
    mean(APS_ind_mcmcdiag$rhat_sigma$psrf[, 1])
)

rhat_results[7, ] <- c(
    "APS_spl", "Regular", mean(APS_spl_mcmcdiag$rhat_y$psrf[, 1]), sd(APS_spl_mcmcdiag$rhat_y$psrf[, 1]),
    mean(APS_spl_mcmcdiag$rhat_ddy$psrf[, 1]), sd(APS_spl_mcmcdiag$rhat_ddy$psrf[, 1]),
    mean(APS_spl_mcmcdiag$rhat_tau$psrf[, 1]), sd(APS_spl_mcmcdiag$rhat_tau$psrf[, 1]),
    mean(APS_spl_mcmcdiag$rhat_sigma$psrf[, 1])
)



rhat_results[8, ] <- c(
    "APS_ar", "Regular", mean(APS_ar_mcmcdiag$rhat_y$psrf[, 1]), sd(APS_ar_mcmcdiag$rhat_y$psrf[, 1]),
    mean(APS_ar_mcmcdiag$rhat_ddy$psrf[, 1]), sd(APS_ar_mcmcdiag$rhat_ddy$psrf[, 1]),
    mean(APS_ar_mcmcdiag$rhat_tau$psrf[, 1]), sd(APS_ar_mcmcdiag$rhat_tau$psrf[, 1]),
    mean(APS_ar_mcmcdiag$rhat_sigma$psrf[, 1])
)



# ===================== Computational time  ====================================
comp_time <- data.frame(
    type = rep("", nmodels), name = rep("", nmodels),
    time = rep(NA, nmodels), y_ess_s = rep(NA, nmodels),
    ddy_ess_s = rep(NA, nmodels), tau_ess_s = rep(NA, nmodels),
    sigma_ess_s = rep(NA, nmodels)
)

comp_time[1, ] <- c("Filter", "ACF", Filter_fit$time$acf[3], rep(NA, 4))
comp_time[2, ] <- c("Filter", "Resid", Filter_fit$time$resid[3], rep(NA, 4))
comp_time[3, ] <- c("Filter", "Optimal", Filter_fit$time$opt[3], rep(NA, 4))
comp_time[4, ] <- c("PSpline", "Regular", PSpline_fit$time[3], ess_results[4, 3] / PSpline_fit$time[3], ess_results[4, 5] / PSpline_fit$time[3], ess_results[4, 7] / PSpline_fit$time[3], ess_results[4, 9] / PSpline_fit$time[3])
comp_time[5, ] <- c("APS_ind", "Regular", APS_ind_fit$time[3], ess_results[5, 3] / APS_ind_fit$time[3], ess_results[5, 5] / APS_ind_fit$time[3], ess_results[5, 7] / APS_ind_fit$time[3], ess_results[5, 9] / APS_ind_fit$time[3])
comp_time[6, ] <- c("APS_spl", "Regular", APS_spl_fit$time[3], ess_results[6, 3] / APS_spl_fit$time[3], ess_results[6, 5] / APS_spl_fit$time[3], ess_results[6, 7] / APS_spl_fit$time[3], ess_results[6, 9] / APS_spl_fit$time[3])
comp_time[7, ] <- c("APS_ar", "Regular", APS_ar_fit$time[3], ess_results[7, 3] / APS_ar_fit$time[3], ess_results[7, 5] / APS_ar_fit$time[3], ess_results[7, 7] / APS_ar_fit$time[3], ess_results[7, 9] / APS_ar_fit$time[3])

for (i in 3:7) {
    comp_time[, i] <- as.numeric(comp_time[, i])
}
print(xtable(comp_time, digits = 2), include.rownames = FALSE)

# ============================= Save processed results ======================================
saveRDS(
    list(
        rmse_results = rmse_results,
        ess_results = ess_results,
        rhat_results = rhat_results,
        comp_time = comp_time
    ),
    file = "./results/APS_results.rds"
)


# ============================= Plot Figure 1 (Fit) ======================================

# Plot fit
pdf(file = "results/figures/Impact_ddy.pdf", width = 11, height = 8.5)
par(mfrow = c(1, 1))
plot(t, ddy,
    type = "l", col = "black",
    ylim = c(-450, 150),
    xlab = TeX(r'(Time $(s)$)'),
    ylab = TeX(r'(Angular Acceleration $(rad\cdot s^{-2})$)')
)
lines(t, Filter_fit$ddy_acf, col = pal[1])
lines(t, PSpline_est$ddy_hat, col = pal[2])
lines(t, APS_ind_est$ddy_hat, col = pal[3])
lines(t, APS_spl_est$ddy_hat, col = pal[4])
lines(t, APS_ar_est$ddy_hat, col = pal[5])
legend(
    x = 0.95, y = -300,
    legend = c("Criterion", TeX(r'(Filter$_{acf}$)'), "P-Spline", TeX(r'(APS$_{ind}$)'), TeX(r'(APS$_{spl}$)'), TeX(r'(APS$_{ar})')),
    col = c("black", pal[1:5]),
    lty = 1,
    bty = "n"
)

# insert A (peak)
v <- c(0.2, 0.3, 0.2, 0.55)
par(fig = v, new = TRUE, mar = c(0, 0, 0, 0))

idx <- which(t > 0.365 & t < 0.41)
ttmp <- t[idx]
# par(mfrow=c(1,1))
plot(ttmp, ddy[idx],
    type = "l", col = "black",
    ylim = c(-450, 20),
    xlab = "", ylab = ""
    # bty='n'#, xaxt='n', yaxt='n'
)
lines(ttmp, Filter_fit$ddy_acf[idx], col = pal[1])
lines(ttmp, PSpline_est$ddy_hat[idx], col = pal[2])
lines(ttmp, APS_ind_est$ddy_hat[idx], col = pal[3])
lines(ttmp, APS_spl_est$ddy_hat[idx], col = pal[4])
lines(ttmp, APS_ar_est$ddy_hat[idx], col = pal[5])

# insert B (stationary)
v <- c(0.45, 0.65, 0.2, 0.55)
par(fig = v, new = TRUE, mar = c(0, 0, 0, 0))
idx <- which(t > 0.42 & t < 0.65)
ttmp <- t[idx]
# par(mfrow=c(1,1))
plot(ttmp, ddy[idx],
    type = "l", col = "black",
    ylim = c(-50, 50),
    xlab = "", ylab = ""
    # bty='n'#, xaxt='n', yaxt='n'
)
lines(ttmp, Filter_fit$ddy_acf[idx], col = pal[1])
lines(ttmp, PSpline_est$ddy_hat[idx], col = pal[2])
lines(ttmp, APS_ind_est$ddy_hat[idx], col = pal[3])
lines(ttmp, APS_spl_est$ddy_hat[idx], col = pal[4])
lines(ttmp, APS_ar_est$ddy_hat[idx], col = pal[5])
dev.off()
# ========================= Plot Figure 2 Tau ===============================

# Plot tau
pdf(file = "results/figures/Impact_tau.pdf", width = 11, height = 8.5)
t_tau <- seq(from = min(t), to = max(t), length.out = length(APS_spl_est$tau_hat))
plot(t_tau, rep(PSpline_est$tau_hat, length(t_tau)),
    col = pal[2], type = "l",
    xlab = TeX(r'(Time $(s)$)'), ylab = TeX(r'($\tau$)'),
    ylim = c(0, 0.12)
)

lines(t_tau, APS_ind_est$tau_hat, col = pal[3])
lines(t_tau, APS_spl_est$tau_hat, col = pal[4])
lines(t_tau, APS_ar_est$tau_hat, col = pal[5])
legend(
    "topright",
    legend = c("P-Spline", TeX(r'(APS$_{ind}$)'), TeX(r'(APS$_{spl}$)'), TeX(r'(APS$_{ar})')),
    col = pal[2:5],
    lty = 1,
    bty = "n"
)
dev.off()

# ============================= Plot Figure 3 90% CI ======================================

nsims <- 100
simidx <- sample(1:nrow(PSpline_fit$ddy), nsims, replace = FALSE)
peak_idx <- 196
stat_idx <- which(abs(t - 0.6) == min(abs(t - 0.6)))

pdf(paste0("./results/figures/Impact_CI.pdf"), width = 11, height = 8.5)
par(mfcol = c(2, 2), mar = c(4, 4.1, 1.1, 0.5), mgp = c(2, 1, 0))
plot(t, ddy,
    type = "l", col = "black",
    ylim = c(-450, 200),
    xlab = "",
    ylab = TeX(r'(Angular Acceleration $(rad\cdot s^{-2})$)'),
    main = "P-Spline"
)

for (i in simidx) {
    lines(t, PSpline_fit$ddy[i, ], col = adjustcolor(pal[2], alpha = 0.025))
}
lines(t, PSpline_est$ddy_hat, col = pal[2])
lines(t, PSpline_est$ddy_ci[, 1], col = pal[2], lty = 3)
lines(t, PSpline_est$ddy_ci[, 2], col = pal[2], lty = 3)
text(x = 0.8, y = -300, TeX(paste0("Avg CI width: ", round(mean(PSpline_est$ddy_ci[, 2] - PSpline_est$ddy_ci[, 1]), 2), r'($rads^{-2}$)')), col = pal[2])

plot(t, ddy,
    type = "l", col = "black",
    ylim = c(-450, 200),
    xlab = TeX(r'(Time $(s)$)'),
    ylab = TeX(r'(Angular Acceleration $(rad\cdot s^{-2})$)'),
    main = TeX(r'($APS_{ind}$)', bold = TRUE)
)

for (i in simidx) {
    lines(t, APS_ind_fit$ddy[i, ], col = adjustcolor(pal[3], alpha = 0.025))
}
lines(t, APS_ind_est$ddy_hat, col = pal[3])
lines(t, APS_ind_est$ddy_ci[, 1], col = pal[3], lty = 3)
lines(t, APS_ind_est$ddy_ci[, 2], col = pal[3], lty = 3)
text(x = 0.8, y = -300, TeX(paste0("Avg CI width: ", round(mean(APS_ind_est$ddy_ci[, 2] - APS_ind_est$ddy_ci[, 1]), 2), r'($rads^{-2}$)')), col = pal[3])


plot(t, ddy,
    type = "l", col = "black",
    ylim = c(-450, 200),
    ylab = "",
    main = TeX(r'($APS_{spl}$)', bold = TRUE)
)

for (i in simidx) {
    lines(t, APS_spl_fit$ddy[i, ], col = adjustcolor(pal[4], alpha = 0.025))
}
lines(t, APS_spl_est$ddy_hat, col = pal[4])
lines(t, APS_spl_est$ddy_ci[, 1], col = pal[4], lty = 3)
lines(t, APS_spl_est$ddy_ci[, 2], col = pal[4], lty = 3)
text(x = 0.8, y = -300, TeX(paste0("Avg CI width: ", round(mean(APS_spl_est$ddy_ci[, 2] - APS_spl_est$ddy_ci[, 1]), 2), r'($rads^{-2}$)')), col = pal[4])


plot(t, ddy,
    type = "l", col = "black",
    ylim = c(-450, 200),
    xlab = TeX(r'(Time $(s)$)'),
    ylab = "",
    main = TeX(r'($APS_{ar}$)', bold = TRUE)
)

for (i in simidx) {
    lines(t, APS_ar_fit$ddy[i, ], col = adjustcolor(pal[5], alpha = 0.025))
}
lines(t, APS_ar_est$ddy_hat, col = pal[5])
lines(t, APS_ar_est$ddy_ci[, 1], col = pal[5], lty = 3)
lines(t, APS_ar_est$ddy_ci[, 2], col = pal[5], lty = 3)
text(x = 0.8, y = -300, TeX(paste0("Avg CI width: ", round(mean(APS_ar_est$ddy_ci[, 2] - APS_ar_est$ddy_ci[, 1]), 2), r'($rads^{-2}$)')), col = pal[5])
dev.off()

# ========================= Compute COverage ===============================

coverage <- rep(NA, 4)
coverage[1] <- mean(ddy < PSpline_est$ddy_ci[, 2] & ddy > PSpline_est$ddy_ci[, 1])
coverage[2] <- mean(ddy < APS_ind_est$ddy_ci[, 2] & ddy > APS_ind_est$ddy_ci[, 1])

coverage[3] <- mean(ddy < APS_spl_est$ddy_ci[, 2] & ddy > APS_spl_est$ddy_ci[, 1])
coverage[4] <- mean(ddy < APS_ar_est$ddy_ci[, 2] & ddy > APS_ar_est$ddy_ci[, 1])

# ======================= Plot Figure 4: Generated Quatntites =======================

# Generate estimates of peak acceleration
nsamps <- dim(APS_spl_fit$ddy)[1]
peak_accn_PSpline <- rep(NA, nsamps)
peak_accn_APS_ind <- rep(NA, nsamps)
peak_accn_APSspl <- rep(NA, nsamps)
peak_accn_APSar <- rep(NA, nsamps)

for (i in 1:nsamps) {
    peak_accn_PSpline[i] <- max(abs(PSpline_fit$ddy[i, ]))
    peak_accn_APS_ind[i] <- max(abs(APS_ind_fit$ddy[i, ]))

    peak_accn_APSspl[i] <- max(abs(APS_spl_fit$ddy[i, ]))
    peak_accn_APSar[i] <- max(abs(APS_ar_fit$ddy[i, ]))
}



# Generate estiamtes of stationary duration
stationary_duration_criterion <- 0
tidx <- which(t > 0.4 & t < 0.8)
l <- 0
for (j in 1:length(tidx)) {
    if (abs(ddy[tidx[j]]) < 5) {
        l <- l + 1
    } else {
        l <- 0
    }
    if (l > stationary_duration_criterion) {
        stationary_duration_criterion <- l
    }
}
stationary_duration_criterion <- stationary_duration_criterion / fs

stationary_duration_filter_acf <- 0
l <- 0
for (j in 1:length(tidx)) {
    if (abs(Filter_fit$ddy_acf[tidx[j]]) < 5) {
        l <- l + 1
    } else {
        l <- 0
    }
    if (l > stationary_duration_filter_acf) {
        stationary_duration_filter_acf <- l
    }
}
stationary_duration_filter_acf <- stationary_duration_filter_acf / fs

stationary_duration_PSpline <- rep(0, nsamps)
stationary_duration_APS_ind <- rep(0, nsamps)
stationary_duration_APS_spl <- rep(0, nsamps)
stationary_duration_APS_ar <- rep(0, nsamps)

for (i in 1:nsamps) {
    # Pspline
    l <- 0
    for (j in 1:length(tidx)) {
        if (abs(PSpline_fit$ddy[i, tidx[j]]) < 5) {
            l <- l + 1
        } else {
            l <- 0
        }
        if (l > stationary_duration_PSpline[i]) {
            stationary_duration_PSpline[i] <- l
        }
    }

    # Pspline heavy
    l <- 0
    for (j in 1:length(tidx)) {
        if (abs(APS_ind_fit$ddy[i, tidx[j]]) < 5) {
            l <- l + 1
        } else {
            l <- 0
        }
        if (l > stationary_duration_APS_ind[i]) {
            stationary_duration_APS_ind[i] <- l
        }
    }


    # APS_spl_fit
    l <- 0
    for (j in 1:length(tidx)) {
        if (abs(APS_spl_fit$ddy[i, tidx[j]]) < 5) {
            l <- l + 1
        } else {
            l <- 0
        }
        if (l > stationary_duration_APS_spl[i]) {
            stationary_duration_APS_spl[i] <- l
        }
    }

    # APS_ar_fit
    l <- 0
    for (j in 1:length(tidx)) {
        if (abs(APS_ar_fit$ddy[i, tidx[j]]) < 5) {
            l <- l + 1
        } else {
            l <- 0
        }
        if (l > stationary_duration_APS_ar[i]) {
            stationary_duration_APS_ar[i] <- l
        }
    }
}

stationary_duration_PSpline <- stationary_duration_PSpline / fs
stationary_duration_APS_ind <- stationary_duration_APS_ind / fs
stationary_duration_APS_spl <- stationary_duration_APS_spl / fs
stationary_duration_APS_ar <- stationary_duration_APS_ar / fs

# Make Plot
pdf(paste0("./results/figures/GenQuantities.pdf"), width = 11, height = 6.5)
par(mfrow = c(1, 2))
plot(density(peak_accn_APS_ind)$x, density(peak_accn_APS_ind)$y,
    col = pal[3], type = "l",
    ylim = c(-0.0005, 0.02), xlim = c(150, 500), xlab = TeX(r'(Impact Peak ($rad \cdot s^{-2}$))'), ylab = "Posterior Density"
)

h <- density(peak_accn_APS_ind)$y[which.min(abs(density(peak_accn_APS_ind)$x - median(peak_accn_APS_ind)))]
segments(x0 = median(peak_accn_APS_ind), y0 = 0, x1 = median(peak_accn_APS_ind), y1 = h, col = pal[3])
segments(x0 = quantile(peak_accn_APS_ind, 0.05), y0 = 0, x1 = quantile(peak_accn_APS_ind, 0.95), y1 = 0, col = pal[3], lty = 1, lwd = 2)

lines(density(peak_accn_PSpline)$x, density(peak_accn_PSpline)$y,
    col = pal[2]
)
h <- density(peak_accn_PSpline)$y[which.min(abs(density(peak_accn_PSpline)$x - median(peak_accn_PSpline)))]

segments(x0 = median(peak_accn_PSpline), y0 = 0.0002, x1 = median(peak_accn_PSpline), y1 = h, col = pal[2])
segments(x0 = quantile(peak_accn_PSpline, 0.05), y0 = 0.0002, x1 = quantile(peak_accn_PSpline, 0.95), y1 = 0.0002, col = pal[2], lty = 1, lwd = 2)


lines(density(peak_accn_APSspl)$x, density(peak_accn_APSspl)$y,
    col = pal[4]
)
h <- density(peak_accn_APSspl)$y[which.min(abs(density(peak_accn_APSspl)$x - median(peak_accn_APSspl)))]

segments(x0 = median(peak_accn_APSspl), y0 = -0.0002, x1 = median(peak_accn_APSspl), y1 = h, col = pal[4])
segments(x0 = quantile(peak_accn_APSspl, 0.05), y0 = -0.0002, x1 = quantile(peak_accn_APSspl, 0.95), y1 = -0.0002, col = pal[4], lty = 1, lwd = 2)

lines(density(peak_accn_APSar)$x, density(peak_accn_APSar)$y, col = pal[5])
h <- density(peak_accn_APSar)$y[which.min(abs(density(peak_accn_APSar)$x - median(peak_accn_APSar)))]

segments(x0 = median(peak_accn_APSar), y0 = -0.0004, x1 = median(peak_accn_APSar), y1 = h, col = pal[5])
segments(x0 = quantile(peak_accn_APSar, 0.05), y0 = -0.0004, x1 = quantile(peak_accn_APSar, 0.95), y1 = -0.0004, col = pal[5], lty = 1, lwd = 2)
abline(v = (max(abs(Filter_fit$ddy_acf))), col = pal[1], lty = 2)
abline(v = max(abs(ddy)), col = "black", lwd = 2)



plot(density(stationary_duration_PSpline, adjust = 3)$x, density(stationary_duration_PSpline, adjust = 3)$y,
    col = pal[2], type = "l",
    xlab = TeX(r'(Stationary Duration(s))'), ylab = "Posterior Density",
    ylim = c(-5, 100), xlim = c(0, 0.3)
)
h <- density(stationary_duration_PSpline, adjust = 3)$y[which.min(abs(density(stationary_duration_PSpline, adjust = 3)$x - median(stationary_duration_PSpline)))]
segments(x0 = median(stationary_duration_PSpline), y0 = -1.5, x1 = median(stationary_duration_PSpline), y1 = h, col = pal[2])
segments(x0 = quantile(stationary_duration_PSpline, 0.05), y0 = -1.5, x1 = quantile(stationary_duration_PSpline, 0.95), y1 = -1.5, col = pal[2], lty = 1, lwd = 2)

lines(density(stationary_duration_APS_ind, adjust = 3)$x, density(stationary_duration_APS_ind, adjust = 3)$y,
    col = pal[3]
)
h <- density(stationary_duration_APS_ind, adjust = 3)$y[which.min(abs(density(stationary_duration_APS_ind, adjust = 3)$x - median(stationary_duration_APS_ind)))]
segments(x0 = median(stationary_duration_APS_ind), y0 = -1, x1 = median(stationary_duration_APS_ind), y1 = h, col = pal[3])
segments(x0 = quantile(stationary_duration_APS_ind, 0.05), y0 = -1, x1 = quantile(stationary_duration_APS_ind, 0.95), y1 = -1, col = pal[3], lty = 1, lwd = 2)

h <- density(stationary_duration_APS_spl)$y[which.min(abs(density(stationary_duration_APS_spl)$x - median(stationary_duration_APS_spl)))]
lines(density(stationary_duration_APS_spl)$x, density(stationary_duration_APS_spl)$y,
    col = pal[4]
)
segments(x0 = median(stationary_duration_APS_spl), y0 = -0.5, x1 = median(stationary_duration_APS_spl), y1 = h, col = pal[4])
segments(x0 = quantile(stationary_duration_APS_spl, 0.05), y0 = -0.5, x1 = quantile(stationary_duration_APS_spl, 0.95), y1 = -0.5, col = pal[4], lty = 1, lwd = 2)


lines(density(stationary_duration_APS_ar)$x, density(stationary_duration_APS_ar)$y,
    col = pal[5]
)
h <- density(stationary_duration_APS_ar)$y[which.min(abs(density(stationary_duration_APS_ar)$x - median(stationary_duration_APS_ar)))]
segments(x0 = median(stationary_duration_APS_ar), y0 = 0, x1 = median(stationary_duration_APS_ar), y1 = h, col = pal[5])
segments(x0 = quantile(stationary_duration_APS_ar, 0.05), y0 = 0, x1 = quantile(stationary_duration_APS_ar, 0.95), y1 = 0, col = pal[5], lty = 1, lwd = 2)

abline(v = stationary_duration_criterion, col = "black", lwd = 2)
abline(v = stationary_duration_filter_acf, col = pal[1], lty = 2)

legend("topright", c("Criterion", TeX(r'(Filter$_{acf}$)'), TeX(r"(P-Spline)"), TeX(r"(APS$_{ind}$)"), TeX(r"(APS$_{spl}$)"), TeX(r"(APS$_{ar}$)")),
    col = c("black", adjustcolor(pal[1:5], 1.0)),
    lty = c(1, 2, rep(1, 4)),
    bty = "n"
)
dev.off()
