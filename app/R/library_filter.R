# Filter functions for AdPSpline Analysis.
library("signal")
library("fftw")

logit <- function(x) {
    return(log(x / (1 - x)))
}
inv_logit <- function(x) {
    return(1 / (1 + exp(-x)))
}

bounded_transform <- function(x, a, b) {
    return(logit((x - a) / (b - a)))
}
inv_bounded_transform <- function(x, a, b) {
    return(a + (b - a) * inv_logit(x))
}

filtbutter <- function(x, b, padtype = "circular", pad_length = NA) {
    if (padtype == "circular") {
        firsthalf <- -rev(x[1:floor(length(x) / 2)])
        secondhalf <- -rev(x[(floor(length(x) / 2) + 1):length(x)])

        ypad <- c(firsthalf + (x[1] - firsthalf[length(firsthalf)]), x, secondhalf + (x[length(x)] - secondhalf[1]))

        fit <- lm(ypad ~ c(1:length(ypad)))

        ypad <- ypad - fit$fitted.values

        yhat <- filtfilt(b, ypad)
        yhat <- yhat[(length(firsthalf) + 1):(length(firsthalf) + length(x))]
        yhat <- yhat + fit$fitted.values[(length(firsthalf) + 1):(length(firsthalf) + length(x))]
    } else if (padtype == "extend") {
        ypad <- c(rep(x[1], pad_length), x, rep(x[length(x)], pad_length))
        yhat <- filtfilt(b, ypad)
        yhat <- yhat[(pad_length + 1):(pad_length + length(x))]
    }
    return(yhat)
}

cutoff_acf <- function(x, fs) {
    cost <- function(lfc, x, fs) {
        fc <- inv_bounded_transform(lfc, 0, fs / 2)
        xhat <- filtbutter(x, butter(2, fc / (fs / 2), type = "low"))
        return(sum(acf(x - xhat, plot = FALSE, lag.max = length(x))$acf^2))
    }
    f <- function(lfc) {
        return(cost(lfc, x, fs))
    }
    res <- optim(bounded_transform(10, 0, fs / 2), f, method = "BFGS")
    fc <- inv_bounded_transform(res$par, 0, fs / 2)
    xf <- filtbutter(x, butter(2, fc / (fs / 2), type = "low"))

    return(list(
        fc = fc,
        resid = sqrt(sum((x - xf)^2))
    ))
}

cutoff_residual <- function(x, fs, tune = 0.7) {
    fcs <- seq(0, fs / 2, length.out = 10000)
    resid <- rep(0, length(fcs))
    for (i in 1:length(fcs)) {
        xf <- filtbutter(x, butter(2, fcs[i] / (fs / 2), type = "low"))
        resid[i] <- sqrt(sum((x - xf)^2))
    }

    n_fcs <- length(fcs) - floor(tune * length(fcs))
    R2 <- rep(0, n_fcs)
    fit_hat <- list()
    for (i in 1:n_fcs) {
        rsamp <- resid[i:length(resid)]
        fit_hat[[i]] <- lm(rsamp ~ fcs[i:length(fcs)])
        y_hat <- predict(fit_hat[[i]])
        R2[i] <- 1 - sum((rsamp - y_hat)^2) / sum((rsamp - mean(rsamp))^2)
    }
    fit <- fit_hat[[which(R2 == max(R2))]]

    cutoff <- fcs[which(log((resid - fit$coefficients[1])^2) == min(log((resid - fit$coefficients[1])^2)))]
    xf <- filtbutter(x, butter(2, cutoff / (fs / 2), type = "low"))
    resid <- sqrt(sum((x - xf)^2))
    return(list(
        fc = cutoff,
        resid = resid
    ))
}


central_finite_diff <- function(x, t, pad = TRUE) {
    if (!is.vector(x)) {
        stop("x must be 1d")
    }
    if (pad) {
        x_pad <- c(x[1], x, x[length(x)])
        t_pad <- c(t[1], t, t[length(t)])
    } else {
        x_pad <- x
        t_pad <- t
    }
    x <- rep(0, length(x))
    for (i in 2:(length(x_pad) - 1)) {
        x[i - 1] <- (x_pad[i + 1] - x_pad[i - 1]) / (t_pad[i + 1] - t_pad[i - 1])
    }
    return(x)
}


finite_diff <- function(x, t, ord, pad = TRUE) {
    xtmp <- x
    if (ord == 0) {
        return(xtmp)
    } else {
        for (i in 1:ord) {
            xtmp <- central_finite_diff(xtmp, t, pad = pad)
        }
        return(xtmp)
    }
}


# Optimal
opt_cost <- function(fcstar, y_kin, y_frc, mom, fmg, fs, fc, to, criteria = "mom") {
    cut_freq <- inv_bounded_transform((fcstar), 0, fs / 2)
    y_kin_tmp <- matrix(NA, nrow = nrow(y_kin), ncol = ncol(y_kin))
    for (i in 1:ncol(y_kin)) {
        y_kin_tmp[, i] <- filtbutter(y_kin[, i], butter(2, cut_freq[i] / (fs / 2), type = "low"), padtype = "circular")
    }
    y_frc_tmp <- matrix(NA, nrow = nrow(y_frc), ncol = ncol(y_frc))
    for (i in 1:ncol(y_frc)) {
        y_frc_tmp[, i] <- filtbutter(y_frc[, i], butter(2, cut_freq[i + ncol(y_frc)] / (fs / 2), type = "low"), padtype = "circular")
    }
    mom_tmp <- compute_mom(y_kin_tmp, y_frc_tmp)
    if (criteria == "mom") {
        cst <- sum((cbind(mom_tmp$Mhip[fc:to], mom_tmp$Mknee[fc:to], mom_tmp$Mankle[fc:to]) - mom[fc:to, ])^2) / (nrow(mom[fc:to, ]) * ncol(mom[fc:to, ])) # Moment cost
    } else if (criteria == "fmg") {
        cst <- sum((cbind(mom_tmp$FMhip[fc:to], mom_tmp$FMknee[fc:to], mom_tmp$FMankle[fc:to]) - fmg[fc:to, ])^2) / (nrow(fmg[fc:to, ]) * ncol(fmg[fc:to, ])) # Force cost
    } else if (critera == "both") {
        cst <- sum((cbind(
            (mom_tmp$Mhip - mean(mom_tmp$Mhip)) / sd(mom_tmp$Mhip),
            (mom_tmp$Mknee - mean(mom_tmp$Mknee)) / sd(mom_tmp$Mknee),
            (mom_tmp$Mankle - mean(mom_tmp$Mankle)) / sd(mom_tmp$Mankle)
        ) -
            cbind(
                (mom[, 1] - mean(mom[, 1])) / sd(mom[, 1]),
                (mom[, 2] - mean(mom[, 2])) / sd(mom[, 2]),
                (mom[, 3] - mean(mom[, 3])) / sd(mom[, 3])
            ))^2) / (nrow(mom) * ncol(mom)) +

            sum((cbind(
                (mom_tmp$FMhip - mean(mom_tmp$FMhip)) / sd(mom_tmp$FMhip),
                (mom_tmp$FMknee - mean(mom_tmp$FMknee)) / sd(mom_tmp$FMknee),
                (mom_tmp$FMankle - mean(mom_tmp$FMankle)) / sd(mom_tmp$FMankle)
            ) -
                cbind(
                    (fmg[, 1] - mean(fmg[, 1])) / sd(fmg[, 1]),
                    (fmg[, 2] - mean(fmg[, 2])) / sd(fmg[, 2]),
                    (fmg[, 3] - mean(fmg[, 3])) / sd(fmg[, 3])
                ))^2) / (nrow(fmg) * ncol(fmg))
    }
    #

    return(cst)
}


id_cutoff_optimal <- function(kin, frc, mom, fmg, fs, fc, to, criteria = "mom") {
    f <- function(fcstar) {
        return(opt_cost(fcstar, kin, frc, mom, fmg, fs, fc, to, criteria = criteria))
    }
    res <- optim(bounded_transform(rep(5, ncol(kin) + ncol(frc)), 0, fs / 2), f, method = "BFGS", control = list(trace = 0, maxit = 20000))
    fc_opt <- inv_bounded_transform(res$par, 0, fs / 2)
    resid <- res$value

    return(list(
        fc = fc_opt,
        resid = resid
    ))
}


# Optimal
post_opt_cost_mom <- function(fcut, y, cirt, fs, fc, to) {
    y_tmp <- filtbutter(y, butter(2, fcut / (fs / 2), type = "low"))
    cst <- rmse(y_tmp[fc:to], cirt[fc:to])
    return(cst)
}



id_cutoff_optimal_post_mom <- function(y, crit, fs, fc, to) {
    f <- function(fcstar) {
        return(post_opt_cost(fcstar, y, crit, fs, fc, to))
    }
    res <- optim(1, f, lower = 0, upper = fs / 2, method = "L-BFGS-B", control = list(trace = 3, maxit = 20000))
    fc_opt <- res$par
    resid <- res$value

    return(list(
        fc = fc_opt,
        resid = resid
    ))
}


post_opt_cost_fmg <- function(fcut, y, cirt, fs, fc, to) {
    y_tmp <- matrix(NA, nrow = nrow(y), ncol = 2)
    y_tmp[, 1] <- filtbutter(y[, 1], butter(2, fcut[1] / (fs / 2), type = "low"))
    y_tmp[, 2] <- filtbutter(y[, 2], butter(2, fcut[2] / (fs / 2), type = "low"))

    cst <- rmse(sqrt(y_tmp[fc:to, 1]^2 + y_tmp[fc:to, 2]^2), cirt[fc:to])
    return(cst)
}

id_cutoff_optimal_post_fmg <- function(y, crit, fs, fc, to) {
    f <- function(fcstar) {
        return(post_opt_cost_fmg(fcstar, y, crit, fs, fc, to))
    }
    res <- optim(c(10, 10), f, lower = c(0, 0), upper = c(fs / 2, fs / 2), method = "L-BFGS-B", control = list(trace = 0, maxit = 20000))
    fc_opt <- res$par
    resid <- res$value

    return(list(
        fc = fc_opt,
        resid = resid
    ))
}


dowling_cutoff_optimal <- function(x, xcriterion, t, fs) {
    cost <- function(fcstar, x, xcriterion, t, fs) {
        fc <- inv_bounded_transform((fcstar), 1, fs / 2)
        xhat <- filtbutter(x, butter(2, fc / (fs / 2), type = "low"))
        ddxhat <- finite_diff(xhat, t, 2)
        return(sqrt(mean((xcriterion - ddxhat)^2)))
    }
    f <- function(fcstar) {
        return(cost(fcstar, x, xcriterion, t, fs))
    }
    res <- optim(bounded_transform(10, 1, fs / 2), f, method = "BFGS", control = list(trace = 0, maxit = 20000))
    xf <- filtbutter(x, butter(2, inv_bounded_transform(res$par, 1, fs / 2) / (fs / 2), type = "low"))
    resid <- sqrt(sum((x - xf)^2))
    return(list(
        fc = inv_bounded_transform(res$par, 1, fs / 2),
        resid = resid
    ))
}


cutoff_power <- function(x, power, fs) {
    n <- length(x)
    n2 <- floor(n / 2)
    dft <- fft(x)
    sig_pow <- (abs(dft)^2) / n
    fp <- (seq_len(n) - 1) * fs / n

    sig_pow <- sig_pow[1:n2]
    fp <- fp[1:n2]
    cum_pow <- cumsum(sig_pow) / sum(sig_pow)

    fc <- fp[min(which(cum_pow > power))]
    xf <- filtbutter(x, butter(2, fc / (fs / 2), type = "low"))
    resid <- sqrt(sum((x - xf)^2))

    return(list(
        fc = fc,
        resid = resid
    ))
}
