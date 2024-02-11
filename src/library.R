# ==============================================================================
# library.R

# Contains key functions for the APSR package.

# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================

print("===== Welcome to AdPSplineR ===== ")
list.of.packages <- c(
    "here", "Rcpp", "posterior", "RColorBrewer", "coda", "bayesplot",
    "shiny", "DT", "shinycssloaders", "truncnorm"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages)) {
    print(cat(paste0("Installing the following packages: \n", new.packages, sep = "\n")))
    install.packages(new.packages)
}


library("RColorBrewer")
library("truncnorm")

pal <- brewer.pal(5, "Set1")
pal2 <- brewer.pal(6, "Dark2")
pal <- c(pal, pal2[c(4, 3)])

inv_scale <- function(x) {
    cent <- ifelse("scaled:center" %in% names(attributes(x)), attr(y_st, "scaled:center"), NA)
    scale <- ifelse("scaled:scale" %in% names(attributes(x)), attr(y_st, "scaled:scale"), NA)
    xout <- x
    if (sum(is.na(c(cent, scale))) == 2) {
        warning("x not a scaled object: returning x")
    }

    if (!is.na(scale)) {
        xout <- xout * scale
    }

    if (!is.na(cent)) {
        xout <- xout + cent
    }
    return(xout)
}

rmse <- function(x, y, avg = TRUE) {
    if (is.matrix(x) | is.data.frame(x)) {
        rmse_tmp <- rep(NA, ncol(x))
        for (i in 1:ncol(x)) {
            rmse_tmp[i] <- sqrt(mean((x[, i] - y[, i])^2))
        }
        if (avg) {
            return(mean(rmse_tmp))
        } else {
            return(rmse_tmp)
        }
    } else {
        return(sqrt(mean((x - y)^2)))
    }
}

tpower <- function(x, t, p) {
    # Constructs a truncated power basis of order p
    # over the domain t.
    return((x - t)^p * (x > t))
}


basis <- function(x, nK, deg = 3) {
    # TODO throw error if deg < -1
    xl <- min(x) # left boundary
    xr <- max(x) # right boundary

    dx <- (xr - xl) / nK # knot spacing

    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx) # expanded knot sequence (n+1+2*deg knots)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
    B <- (-1)^(deg + 1) * P %*% t(D)

    # Make B-splines exactly zero beyond their end knots
    nb <- ncol(B)
    sk <- knots[(1:nb) + deg + 1]
    Mask <- outer(x, sk, "<")
    B <- B * Mask
    return(list(
        matrix = B,
        deg = deg,
        knots = sk,
        x = x,
        nK = nK,
        deriative = 0,
        constraints = "none"
    ))
}

basis_mm <- function(P) {
    # Constructs fixed matrix X and random matrix Z of the
    # mixed model respresentation of a penalized spline with penalty
    # P

    S <- P$matrix
    S_rank <- qr(S)$rank
    eobj <- eigen(S)
    d <- eobj$values
    d[abs(eobj$values) < sqrt(.Machine$double.eps)] <- 0
    ord_tmp <- order(d)
    d <- d[ord_tmp]
    gam <- eobj$vectors[, ord_tmp]
    gam1 <- gam[, 1:(nrow(S) - S_rank)]
    gam2 <- gam[, -(1:(nrow(S) - S_rank))]
    X_tilde <- gam1
    Z_tilde <- gam2 %*% diag(d[-(1:(nrow(S) - S_rank))]^(-0.5))
    return(list(X_tilde = X_tilde, Z_tilde = Z_tilde))
}

dbasis <- function(B, deriative) {
    nK <- dim(B$matrix)[2] - B$deg

    dB <- basis(B$x, nK, deg = B$deg - deriative)
    h <- (max(B$x) - min(B$x)) / nK
    dB <- dB$matrix %*% diff(diag(nK + B$deg), diff = deriative) / h^deriative
    return(list(
        matrix = dB,
        deg = deg,
        knots = B$knots,
        x = B$x,
        nK = nK,
        deriative = deriative
    ))
}

penalty <- function(Q, pDeg) {
    # Construct a deriative penalty on the deg deriative for a B-Spline basis

    D <- diff(diag(rep(1, Q)), diff = pDeg)
    S <- t(D) %*% D
    return(list(
        D = D,
        matrix = S,
        deg = pDeg
    ))
}



MM_likelihood <- function(theta, X, Z, y) {
    # Constructs the loglikelihood of the mixed model representation of
    # a pspline with fixed effects X and random effects Z.
    # given theta with theta[1] being log transformed variance of random effects...

    # For use in maximisation of the profile likelihood of a MM representation.
    # From Wood 2017 pp. 79

    ## untransform parameters...
    sigma.b <- exp(theta[1])
    sigma <- exp(theta[2])
    ## extract dimensions...
    n <- length(y)
    pr <- ncol(Z)
    pf <- ncol(X)
    ## obtain \hat \beta, \hat b...
    X1 <- cbind(X, Z)
    ipsi <- c(rep(0, pf), rep(1 / sigma.b^2, pr))
    b1 <- solve(
        crossprod(X1) / sigma^2 + diag(ipsi),
        t(X1) %*% y / sigma^2
    )
    ## compute log|Z’Z/sigma^2 + I/sigma.b^2|...
    ldet <- sum(log(diag(chol(crossprod(Z) / sigma^2 +
        diag(ipsi[-(1:pf)])))))
    ## compute log profile likelihood...
    l <- (-sum((y - X1 %*% b1)^2) / sigma^2 - sum(b1^2 * ipsi) -
        n * log(sigma^2) - pr * log(sigma.b^2) - 2 * ldet - n * log(2 * pi)) / 2
    attr(l, "alp") <- as.numeric(b1)[1:pf] ## return \hat beta and \hat b
    attr(l, "gam") <- as.numeric(b1)[(pf + 1):length(b1)]
    -l
}

mle_pspl <- function(B, MM, y) {
    f <- function(x) MM_likelihood(x, B$matrix %*% MM$X_tilde, B$matrix %*% MM$Z_tilde, y)
    res <- optim(rep(0, 2), f, hessian = FALSE, control = list(maxit = 10000))
    alp_hat <- attr(MM_likelihood(res$par, B$matrix %*% MM$X_tilde, B$matrix %*% MM$Z_tilde, y), "alp")
    gam_hat <- attr(MM_likelihood(res$par, B$matrix %*% MM$X_tilde, B$matrix %*% MM$Z_tilde, y), "gam")
    theta_hat <- MM$X_tilde %*% alp_hat + MM$Z_tilde %*% gam_hat
    tau_hat <- exp(res$par[1])
    sigma_hat <- exp(res$par[2])

    return(list(
        theta_hat = theta_hat,
        tau_hat = tau_hat,
        sigma_hat = sigma_hat,
        alp_hat = alp_hat,
        gam_hat = gam_hat
    ))
}

# Save results files
save_rmse_results <- function(
    args,
    model,
    type,
    criterion,
    kin_hat,
    frc_hat,
    mom_hat,
    fc = 42,
    to = 83,
    ord = "pre",
    file = "./results/AdPSpline_RMSE_results.csv") {
    if (ord == "pre") {
        if (model == "Filter") {
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    rmse(criterion$kin[fc:to, ], kin_hat[fc:to, ]),
                    rmse(criterion$frc[fc:to, 1:2], frc_hat[fc:to, 1:2]),
                    rmse(criterion$frc[fc:to, 3], frc_hat[fc:to, 3]),
                    rmse(criterion$mom[fc:to, 1], mom_hat$Mhip[fc:to]),
                    rmse(criterion$mom[fc:to, 2], mom_hat$Mknee[fc:to]),
                    rmse(criterion$mom[fc:to, 3], mom_hat$Mankle[fc:to]),
                    rmse(criterion$fmg[fc:to, 1], mom_hat$FMhip[fc:to]),
                    rmse(criterion$fmg[fc:to, 2], mom_hat$FMknee[fc:to]),
                    rmse(criterion$fmg[fc:to, 3], mom_hat$FMankle[fc:to])
                ), collapse = ","),
                file = file,
                sep = "\n",
                append = TRUE
            )
        } else { # Bayes
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    rmse(criterion$kin[fc:to, ], kin_hat),
                    rmse(criterion$frc[fc:to, 1:2], frc_hat[, 1:2]),
                    rmse(criterion$frc[fc:to, 3], frc_hat[, 3]),
                    rmse(criterion$mom[fc:to, 1], mom_hat$Mhip),
                    rmse(criterion$mom[fc:to, 2], mom_hat$Mknee),
                    rmse(criterion$mom[fc:to, 3], mom_hat$Mankle),
                    rmse(criterion$fmg[fc:to, 1], mom_hat$FMhip),
                    rmse(criterion$fmg[fc:to, 2], mom_hat$FMknee),
                    rmse(criterion$fmg[fc:to, 3], mom_hat$FMankle)
                ), collapse = ","),
                file = file,
                sep = "\n",
                append = TRUE
            )
        }
    } else if (ord == "post") {
        if (model == "Filter") {
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    rmse(criterion$mom[fc:to, 1], mom_hat$Mhip[fc:to]),
                    rmse(criterion$mom[fc:to, 2], mom_hat$Mknee[fc:to]),
                    rmse(criterion$mom[fc:to, 3], mom_hat$Mankle[fc:to]),
                    rmse(criterion$fmg[fc:to, 1], mom_hat$FMhip[fc:to]),
                    rmse(criterion$fmg[fc:to, 2], mom_hat$FMknee[fc:to]),
                    rmse(criterion$fmg[fc:to, 3], mom_hat$FMankle[fc:to])
                ), collapse = ","),
                file = file,
                sep = "\n",
                append = TRUE
            )
        } else { # Bayes
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    rmse(criterion$mom[fc:to, 1], mom_hat$Mhip),
                    rmse(criterion$mom[fc:to, 2], mom_hat$Mknee),
                    rmse(criterion$mom[fc:to, 3], mom_hat$Mankle),
                    rmse(criterion$fmg[fc:to, 1], mom_hat$FMhip),
                    rmse(criterion$fmg[fc:to, 2], mom_hat$FMknee),
                    rmse(criterion$fmg[fc:to, 3], mom_hat$FMankle)
                ), collapse = ","),
                file = file,
                sep = "\n",
                append = TRUE
            )
        }
    }
}

save_ess_results <- function(
    args,
    model,
    type,
    time,
    ess = NA,
    ord = "pre",
    file = "./results/AdPSpline_ESS_results.csv") {
    if (ord == "pre") {
        if (model == "Filter") {
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model, type, ord,
                    time[3],
                    rep(NA, 8)
                ), collapse = ","),
                file = file,
                append = TRUE,
                sep = "\n"
            )
        } else { # Bayes
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model, type, ord,
                    time[3],
                    mean(ess$kin_ess) / time[3],
                    mean(ess$frc_ess) / time[3],
                    mean(ess$mom_ess$Mhip) / time[3],
                    mean(ess$mom_ess$Mknee) / time[3],
                    mean(ess$mom_ess$Mankle) / time[3],
                    mean(ess$mom_ess$FMhip) / time[3],
                    mean(ess$mom_ess$FMknee) / time[3],
                    mean(ess$mom_ess$FMankle) / time[3]
                ), collapse = ","),
                file = file,
                append = TRUE,
                sep = "\n"
            )
        }
    } else if (ord == "post") {
        if (model == "Filter") {
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    time[3],
                    rep(NA, 6)
                ), collapse = ","),
                file = file,
                append = TRUE,
                sep = "\n"
            )
        } else { # Bayes
            cat(
                paste0(c(
                    args[1],
                    args[2],
                    model,
                    type, ord,
                    time[3],
                    mean(ess$mom_ess$Mhip) / time[3],
                    mean(ess$mom_ess$Mknee) / time[3],
                    mean(ess$mom_ess$Mankle) / time[3],
                    mean(ess$mom_ess$FMhip) / time[3],
                    mean(ess$mom_ess$FMknee) / time[3],
                    mean(ess$mom_ess$FMankle) / time[3]
                ), collapse = ","),
                file = file,
                append = TRUE,
                sep = "\n"
            )
        }
    }
}


bayes_est <- function(fit, fc_star = 11, to_star = 52, ord = "pre") {
    if (ord == "pre") {
        kin_hat <- matrix(NA, nrow = fit$data$T, ncol = ncol(fit$data$y_kin))
        for (i in 1:ncol(fit$data$y_kin)) {
            kin_hat[, i] <- apply(fit$mom_samples$y_kin_samps[, , i], 2, median)
        }
        frc_hat <- matrix(NA, nrow = fit$data$T, ncol = 3)
        for (i in 1:3) {
            frc_hat[, i] <- apply(fit$mom_samples$y_frc_samps[, , i], 2, median)
        }

        mom_hat <- list(
            Mhip = apply(fit$mom_samples$Mhip_samps, 2, mean)[fc_star:to_star],
            Mknee = apply(fit$mom_samples$Mknee_samps, 2, mean)[fc_star:to_star],
            Mankle = apply(fit$mom_samples$Mankle_samps, 2, mean)[fc_star:to_star],
            FMhip = apply(fit$mom_samples$FMhip_samps, 2, mean)[fc_star:to_star],
            FMknee = apply(fit$mom_samples$FMknee_samps, 2, mean)[fc_star:to_star],
            FMankle = apply(fit$mom_samples$FMankle_samps, 2, mean)[fc_star:to_star]
        )
        return(list(
            kin_hat = kin_hat[fc_star:to_star, ],
            frc_hat = frc_hat[fc_star:to_star, ],
            mom_hat = mom_hat
        ))
    } else if (ord == "post") {
        mom_hat <- list(
            Mhip = apply(fit$mom_samples$Mhip_samps, 2, mean)[fc_star:to_star],
            Mknee = apply(fit$mom_samples$Mknee_samps, 2, mean)[fc_star:to_star],
            Mankle = apply(fit$mom_samples$Mankle_samps, 2, mean)[fc_star:to_star],
            FMhip = apply(fit$mom_samples$FMhip_samps, 2, mean)[fc_star:to_star],
            FMknee = apply(fit$mom_samples$FMknee_samps, 2, mean)[fc_star:to_star],
            FMankle = apply(fit$mom_samples$FMankle_samps, 2, mean)[fc_star:to_star]
        )
        return(list(mom_hat = mom_hat))
    }
}


ess_est <- function(fit, fc_star = 11, to_star = 52, ord = "pre") {
    if (ord == "pre") {
        y_kin_ess <- matrix(NA, nrow = fit$data$T, ncol = ncol(fit$data$y_kin))
        for (i in 1:ncol(fit$data$y_kin)) {
            ch <- list()
            for (j in 1:fit$pars$nchains) {
                ch[[j]] <- fit$mom_samples$y_kin_samps[((j - 1) * (fit$pars$nsampling / fit$pars$thin) + 1):(j * (fit$pars$nsampling / fit$pars$thin)), , i]
            }
            class(ch) <- "mcmc.list"
            y_kin_ess[, i] <- effectiveSize(ch)
        }

        y_frc_ess <- matrix(NA, nrow = fit$data$T, ncol = 3)
        for (i in 1:3) {
            ch <- list()
            for (j in 1:fit$pars$nchains) {
                ch[[j]] <- fit$mom_samples$y_frc_samps[((j - 1) * (fit$pars$nsampling / fit$pars$thin) + 1):(j * (fit$pars$nsampling / fit$pars$thin)), , i]
            }
            class(ch) <- "mcmc.list"
            y_frc_ess[, i] <- effectiveSize(ch)
        }

        mom_vars <- c("Mhip_samps", "Mknee_samps", "Mankle_samps", "FMhip_samps", "FMknee_samps", "FMankle_samps")
        mom_ess <- matrix(NA, nrow = fit$data$T, ncol = length(mom_vars))
        for (i in 1:length(mom_vars)) {
            ch <- list()
            for (j in 1:fit$pars$nchains) {
                ch[[j]] <- fit$mom_samples[[mom_vars[i]]][((j - 1) * (fit$pars$nsampling / fit$pars$thin) + 1):(j * (fit$pars$nsampling / fit$pars$thin)), ]
            }
            class(ch) <- "mcmc.list"
            mom_ess[, i] <- effectiveSize(ch)
        }

        return(list(
            kin_ess = y_kin_ess[fc_star:to_star, ],
            frc_ess = y_frc_ess[fc_star:to_star, ],
            mom_ess = list(
                Mhip = mom_ess[fc_star:to_star, 1],
                Mknee = mom_ess[fc_star:to_star, 2],
                Mankle = mom_ess[fc_star:to_star, 3],
                FMhip = mom_ess[fc_star:to_star, 4],
                FMknee = mom_ess[fc_star:to_star, 5],
                FMankle = mom_ess[fc_star:to_star, 6]
            )
        ))
    } else if (ord == "post") {
        mom_vars <- c("Mhip_samps", "Mknee_samps", "Mankle_samps", "FMhip_samps", "FMknee_samps", "FMankle_samps")
        mom_ess <- matrix(NA, nrow = fit$data$T, ncol = length(mom_vars))
        for (i in 1:length(mom_vars)) {
            ch <- list()
            for (j in 1:fit$pars$nchains) {
                ch[[j]] <- fit$mom_samples[[mom_vars[i]]][((j - 1) * (fit$pars$nsampling / fit$pars$thin) + 1):(j * (fit$pars$nsampling / fit$pars$thin)), ]
            }
            class(ch) <- "mcmc.list"
            mom_ess[, i] <- effectiveSize(ch)
        }

        return(list(
            mom_ess = list(
                Mhip = mom_ess[fc_star:to_star, 1],
                Mknee = mom_ess[fc_star:to_star, 2],
                Mankle = mom_ess[fc_star:to_star, 3],
                FMhip = mom_ess[fc_star:to_star, 4],
                FMknee = mom_ess[fc_star:to_star, 5],
                FMankle = mom_ess[fc_star:to_star, 6]
            )
        ))
    }
}



order_diff <- function(x, ord) {
    if (!(ord %in% c(2, 3, 4))) stop("Order must be either 2, 3, or 4")
    K <- length(x)
    out <- rep(0, K)
    if (ord == 2) {
        for (k in (ord + 1):K) {
            out[k] <- 2 * out[k - 1] - out[k - 2]
        }
    } else if (ord == 3) {
        for (k in (ord + 1):K) {
            out[k] <- 3 * out[k - 1] - 3 * out[k - 2] + out[k - 3]
        }
    } else if (ord == 4) {
        for (k in (ord + 1):K) {
            out[k] <- 4 * out[k - 1] - 6 * out[k - 2] + 4 * out[k - 3] - out[k - 4]
        }
    }
    return(out[(ord + 1):K])
}

get_prior_hyperparms_opt <- function(ddB, criterion, model = "PSpline", pDeg = 2, atol = 1, nsim = 1000) {
    K <- dim(ddB$matrix)[2]
    print(paste0("            Obtaining hyperparameters for ", model))
    # Cost functions
    cost_pspl <- function(x, ddB, criterion) {
        xi <- exp(x)
        max_derv <- rep(NA, nsim)
        for (i in 1:length(max_derv)) {
            # Generate tau2 ~ t_3+(0, xi^2)
            x <- rchisq(1, df = 3)
            v <- (3 * xi^2) / x # v ~ invchi2(nu, xi^2)
            tau2 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt(v)) # tau2 ~ t_3+(0, xi^2)
            tau <- sqrt(tau2)
            # theta | . is 2nd order RW
            theta <- rep(0, K)
            theta[(pDeg + 1):K] <- order_diff(theta, pDeg) + tau * rnorm(K - pDeg, 0, 1)
            # for (k in (pDeg + 1):K) {
            #     theta[k] <- 2 * theta[k - 1] - theta[k - 2] + tau * rnorm(1)
            # }
            # # Compute max of 2nd deriative
            max_derv[i] <- max(abs(ddB$matrix %*% theta))
        }
        # Compute 90th quantile of max(|y''(t)|)
        qderv <- as.numeric(quantile(max_derv, 0.90, na.rm = TRUE))
        return(sqrt((qderv - criterion)^2))
    }

    cost_aps_ind <- function(x, ddB, criterion) {
        dxi <- exp(x[1])
        da0 <- exp(x[2])
        nu <- 3
        max_derv <- rep(NA, nsim)
        for (i in 1:length(max_derv)) {
            # generate a0 from half-t prior with delta_a0 = 0.002
            x <- rchisq(1, df = 3)
            a0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * da0^2) / x)) # a0^2 ~ invchi2(nu, xi^2)
            # Generate xi ~ t_3+(0, dxi^2)
            x <- rchisq(1, df = 3)
            v_xi <- (3 * dxi^2) / x # tau2_k ~ invchi2(nu, xi^2)
            xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt(v_xi))

            # tau2_k ~ invchi2(nu, xi^2)
            ltau2_k <- rep(NA, K - pDeg)
            ltau2_k <- log(a0^2) + xi * rnorm(K - pDeg, 0, 1) # ltau2_k ~ N(log(a0^2), xi^2)
            tau2_k <- exp(ltau2_k)
            # theta | . is 2nd order RW
            theta <- rep(0, K)
            theta[(pDeg + 1):K] <- order_diff(theta, pDeg) + sqrt(tau2_k) * rnorm(K - pDeg, 0, 1)

            # Compute max of 2nd deriative
            max_derv[i] <- max(abs(ddB$matrix %*% theta))
        }
        qderv <- as.numeric(quantile(max_derv, 0.90, na.rm = TRUE))
        return(sqrt((qderv - criterion)^2))
    }

    cost_aps_ar <- function(x, ddB, criterion) {
        dxi <- exp(x[1])
        da0 <- exp(x[2])
        max_derv <- rep(NA, nsim)
        for (i in 1:length(max_derv)) {
            # generate \alpha_0 from half-t prior with delta_a0 = 0.002
            x <- rchisq(1, df = 3)
            a0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * da0^2) / x))
            # generate xi from half t priro with scale d_xi
            x <- rchisq(1, df = 3)
            xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * dxi^2) / x))
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
            theta[(pDeg + 1):K] <- order_diff(theta, pDeg) + sqrt(exp(ltau2_k)) * rnorm(K - pDeg, 0, 1)
            # Compute max of 2nd deriative
            max_derv[i] <- max(abs(ddB$matrix %*% theta))
        }
        # generate 90th quantile of deriative
        qderv <- quantile(max_derv, 0.90, na.rm = T)
        return(sqrt((qderv - criterion)^2))
    }

    cost_aps_spl <- function(x, ddB, criterion) {
        bDeg_sm <- 3
        C <- basis(seq(0, 1, length.out = dim(ddB$matrix)[2] - pDeg), nK_sm, deg = bDeg_sm)
        K_C <- dim(C$matrix)[2] # Define grid of prior scales
        dxi <- exp(x[1])
        da0 <- exp(x[2])
        max_derv <- rep(NA, nsim)
        for (i in 1:length(max_derv)) {
            # generate \alpha_0 from half-t prior with delta_a0 = 0.002
            x <- rchisq(1, df = 3)
            a0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * da0^2) / x)) # tau2_k ~ invchi2(nu, xi^2)
            # generate xi from half t priro with scale d_xi
            x <- rchisq(1, df = 3)
            xi <- rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = sqrt((3 * dxi^2) / x))
            # gamma ~ N(log(a0^2), xi^2)
            # ltau2 | phi, xi, a0 ~ Cubic Spline given basis matrix C
            gamma <- log(a0^2) + xi * rnorm(K_C, 1)
            ltau2_k <- C$matrix %*% gamma
            # theta | . is 2nd order RW

            theta <- rep(0, K)
            theta[(pDeg + 1):K] <- order_diff(theta, pDeg) + sqrt(exp(ltau2_k)) * rnorm(K - pDeg, 0, 1)
            max_derv[i] <- max(abs(ddB$matrix %*% theta), na.rm = TRUE)
        }
        # generate 90th quantile of deriative
        qderv <- quantile(max_derv, 0.95, na.rm = TRUE)
        return(sqrt((qderv - criterion)^2))
    }


    if (model == "PSpline") {
        out <- suppressWarnings(optim(log(0.01), cost_pspl, ddB = ddB, criterion = criterion, control = list(trace = 0, abstol = atol)))
        hyperparms <- list(xi = exp(out$par))
    } else if (model == "APS_ind") {
        out <- suppressWarnings(optim(log(c(0.01, 0.01)), cost_aps_ind, ddB = ddB, criterion = criterion, control = list(trace = 0, abstol = atol)))
        hyperparms <- list(xi = exp(out$par[1]), a0 = exp(out$par[2]))
    } else if (model == "APS_ar") {
        out <- suppressWarnings(optim(log(c(0.01, 0.01)), cost_aps_ar, ddB = ddB, criterion = criterion, control = list(trace = 0, abstol = atol)))
        hyperparms <- list(xi = exp(out$par[1]), a0 = exp(out$par[2]))
    } else if (model == "APS_spl") {
        out <- suppressWarnings(optim(log(c(0.01, 0.01)), cost_aps_spl, ddB = ddB, criterion = criterion, control = list(trace = 0, abstol = atol)))
        hyperparms <- list(xi = exp(out$par[1]), a0 = exp(out$par[2]))
    }

    return(hyperparms)
}
