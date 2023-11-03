if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("latex2exp")) install.packages("latex2exp")

library("RColorBrewer")
library("latex2exp")

setwd("./sensitivity/")
source("../src/library.R")
source("../src/library_analysis.R")


# ============================= Generate Data ======================================
set.seed(123)
print("...Read data...")

# Process data
y <- as.numeric(read.delim2("../data/Dowling_Disp.txt")[, 1])
ddy <- as.numeric(read.delim2("../data/Dowling_Accn.txt")[1:length(y), 1])
fs <- 512
t <- seq(from = 0, to = (length(y) - 1) * (1 / fs), by = 1 / fs)

# standardize
y_mu <- mean(y)
y_sd <- sd(y)
y_std <- (y - y_mu) / y_sd


# ============================= Load Models ======================================

pspl25 <- bayes_est(readRDS("./mdl/PSpline_K25_seed123.RDS"))
pspl50 <- bayes_est(readRDS("./mdl/PSpline_K50_seed123.RDS"))
pspl100 <- bayes_est(readRDS("./mdl/PSpline_K100_seed123.RDS"))
pspl200 <- bayes_est(readRDS("./mdl/PSpline_K200_seed123.RDS"))
pspl400 <- bayes_est(readRDS("./mdl/PSpline_K400_seed123.RDS"))
pspl <- list(pspl25, pspl50, pspl100, pspl200, pspl400)
rm(pspl25, pspl50, pspl100, pspl200, pspl400)

aps_ind25 <- bayes_est(readRDS("./mdl/APS_ind_K25_seed123.RDS"))
aps_ind50 <- bayes_est(readRDS("./mdl/APS_ind_K50_seed123.RDS"))
aps_ind100 <- bayes_est(readRDS("./mdl/APS_ind_K100_seed123.RDS"))
aps_ind200 <- bayes_est(readRDS("./mdl/APS_ind_K200_seed123.RDS"))
aps_ind400 <- bayes_est(readRDS("./mdl/APS_ind_K400_seed123.RDS"))
aps_ind <- list(aps_ind25, aps_ind50, aps_ind100, aps_ind200, aps_ind400)
rm(aps_ind25, aps_ind50, aps_ind100, aps_ind200, aps_ind400)

aps_spl25 <- bayes_est(readRDS("./mdl/APS_spl_K25_seed123.RDS"))
aps_spl50 <- bayes_est(readRDS("./mdl/APS_spl_K50_seed123.RDS"))
aps_spl100 <- bayes_est(readRDS("./mdl/APS_spl_K100_seed123.RDS"))
aps_spl200 <- bayes_est(readRDS("./mdl/APS_spl_K200_seed123.RDS"))
aps_spl400 <- bayes_est(readRDS("./mdl/APS_spl_K400_seed123.RDS"))
aps_spl <- list(aps_spl25, aps_spl50, aps_spl100, aps_spl200, aps_spl400)
rm(aps_spl25, aps_spl50, aps_spl100, aps_spl200, aps_spl400)

aps_ar25 <- bayes_est(readRDS("./mdl/APS_ar_K25_seed123.RDS"))
aps_ar50 <- bayes_est(readRDS("./mdl/APS_ar_K50_seed123.RDS"))
aps_ar100 <- bayes_est(readRDS("./mdl/APS_ar_K100_seed123.RDS"))
aps_ar200 <- bayes_est(readRDS("./mdl/APS_ar_K200_seed123.RDS"))
aps_ar400 <- bayes_est(readRDS("./mdl/APS_ar_K400_seed123.RDS"))
aps_ar <- list(aps_ar25, aps_ar50, aps_ar100, aps_ar200, aps_ar400)
rm(aps_ar25, aps_ar50, aps_ar100, aps_ar200, aps_ar400)


# ============================= Plot 1 Residuals ======================================

# Generate Residuals
pdf(paste0("./figures/residual_BasisSize.pdf"), width = 11, height = 8.5)
Bcheck <- basis(t, 50, deg = 3)
par(mfrow = c(5, 4), mar = c(4.1, 4.8, 1.1, 0.5), mgp = c(2, 1, 0))
k_labs <- c(
    substitute(paste(bold("K=25"))), substitute(paste(bold("K=50"))),
    substitute(paste(bold("K=100"))), substitute(paste(bold("K=200"))),
    substitute(paste(bold("K=400")))
)
mdl_labs <- c(
    substitute(paste(bold("P-Spline"))), TeX(r"(APS$_{ind}$)", bold = T),
    TeX(r"(APS$_{spl}$)", bold = T), TeX(r"(APS$_{ar}$)", bold = T)
)
for (i in 1:length(pspl)) {
    if (i == 1) {
        plot(t, y - pspl[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[2], alpha = 0.5),
            ylab = k_labs[i], xlab = NA,
            cex.lab = 1.2,
            main = mdl_labs[1],
            # xlab = "Time (s)", ylab = "Residual (rad)",
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - pspl[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_ind[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[3], alpha = 0.5),
            main = mdl_labs[2],
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_ind[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_spl[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[4], alpha = 0.5),
            main = mdl_labs[3],
            cex.lab = 1.2,
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_spl[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_ar[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[5], alpha = 0.5),
            main = mdl_labs[4],
            cex.lab = 1.2,
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_ar[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)
    } else {
        plot(t, y - pspl[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[2], alpha = 0.5),
            ylab = k_labs[i],
            cex.lab = 1.2,
            xlab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - pspl[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_ind[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[3], alpha = 0.5),
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_ind[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_spl[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[4], alpha = 0.5),
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_spl[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)

        plot(t, y - aps_ar[[i]]$y_hat,
            pch = 16, col = adjustcolor(pal[5], alpha = 0.5),
            xlab = NA, ylab = NA,
            ylim = c(-0.04, 0.04)
        )
        pspl_coef <- lm((y - aps_ar[[i]]$y_hat) ~ Bcheck$matrix - 1)$coef
        lines(t, Bcheck$matrix %*% pspl_coef)
    }
}
mtext(text = "Time (s)", side = 1, line = -1.5, outer = TRUE)
mtext(text = "Residual (rad)", side = 2, line = -1.5, outer = TRUE)
dev.off()


# ======================== Create RMSE figure ========================

# RMSE
rmse_y <- matrix(NA, nrow = 5, ncol = 4)
rmse_ddy <- matrix(NA, nrow = 5, ncol = 4)

edf <- edf2 <- matrix(NA, nrow = 5, ncol = 4)



basis_size <- c(25, 50, 100, 200, 400)
for (i in 1:length(pspl)) {
    B <- basis(t, basis_size[i], 5)$matrix
    Bp <- t(B)
    Q <- dim(B)[2]
    P <- penalty(Q, 2)
    D <- P$D
    Dp <- t(D)

    # Pspline
    rmse_y[i, 1] <- rmse(pspl[[i]]$y_hat, y)
    rmse_ddy[i, 1] <- rmse(pspl[[i]]$ddy_hat, ddy)
    S <- Dp %*% diag(rep(pspl[[i]]$tau_hat^2, Q - 2)) %*% D
    A <- B %*% solve(Bp %*% B + S) %*% Bp
    edf[i, 1] <- sum(diag(A))
    edf2[i, 1] <- 2 * sum(diag(A)) - sum(diag(A %*% A))
    # PSpline Heavy
    rmse_y[i, 2] <- rmse(aps_ind[[i]]$y_hat, y)
    rmse_ddy[i, 2] <- rmse(aps_ind[[i]]$ddy_hat, ddy)
    S <- Dp %*% diag(rep(pspl[[i]]$tau_hat^2, Q - 2)) %*% D
    A <- B %*% solve(Bp %*% B + S) %*% Bp
    edf[i, 2] <- sum(diag(A))
    edf2[i, 2] <- 2 * sum(diag(A)) - sum(diag(A %*% A))

    # APS_spl
    rmse_y[i, 3] <- rmse(aps_spl[[i]]$y_hat, y)
    rmse_ddy[i, 3] <- rmse(aps_spl[[i]]$ddy_hat, ddy)
    S <- Dp %*% diag(aps_spl[[i]]$tau_hat^2) %*% D
    A <- B %*% solve(Bp %*% B + S) %*% Bp
    edf[i, 3] <- sum(diag(A))
    edf2[i, 3] <- 2 * sum(diag(A)) - sum(diag(A %*% A))


    rmse_y[i, 4] <- rmse(aps_ar[[i]]$y_hat, y)
    rmse_ddy[i, 4] <- rmse(aps_ar[[i]]$ddy_hat, ddy)

    S <- Dp %*% diag(aps_ar[[i]]$tau_hat^2) %*% D
    A <- B %*% solve(Bp %*% B + S) %*% Bp
    edf[i, 4] <- sum(diag(A))
    edf2[i, 4] <- 2 * sum(diag(A)) - sum(diag(A %*% A))
}

pdf(paste0("./figures/RMSE_BasisSize.pdf"), width = 11, height = 4.5)
par(mfrow = c(1, 2))
plot(basis_size, rmse_y[, 1],
    type = "o", col = pal[2], lty = 2,
    xlab = "Number of knot points", ylab = "RMSE Angle (rad)",
    ylim = range(rmse_y), xaxt = "n"
)
axis(side = 1, at = c(25, 50, 100, 200, 400))
lines(basis_size, rmse_y[, 2], type = "o", col = pal[3], lty = 2)
lines(basis_size, rmse_y[, 3], type = "o", col = pal[4], lty = 2)
lines(basis_size, rmse_y[, 4], type = "o", col = pal[5], lty = 2)

plot(basis_size, rmse_ddy[, 1],
    type = "o", col = pal[2], lty = 2,
    xlab = "Number of knot points",
    ylab = TeX(r'(RMSE Accn $(rad\cdot s^{-2})$)'),
    ylim = range(rmse_ddy), xaxt = "n"
)
axis(side = 1, at = c(25, 50, 100, 200, 400))
lines(basis_size, rmse_ddy[, 2], type = "o", col = pal[3], lty = 2)
lines(basis_size, rmse_ddy[, 3], type = "o", col = pal[4], lty = 2)
lines(basis_size, rmse_ddy[, 4], type = "o", col = pal[5], lty = 2)
legend("topleft",
    legend = c("P-Spline", TeX(r"(APS$_{ind}$)"), TeX(r"(APS$_{spl}$)"), TeX(r"(APS$_{ar}$)")),
    col = pal[2:5],
    pch = c(1, 1, 1, 1),
    lty = 2,
    bty = "n"
)
dev.off()
