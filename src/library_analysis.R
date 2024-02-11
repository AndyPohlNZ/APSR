# ==============================================================================
# library_analysis.R

# Contains key functions for the analysis of  APS model output.

# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================

if (!require("coda")) install.packages("coda")

library("coda")




rmse <- function(x, y) {
    # Compute the Root Mean Square Error (RMSE) between two numeric vectors.
    #
    # Args:
    #   x: A numeric vector containing observed or actual values.
    #   y: A numeric vector containing predicted or estimated values.
    #
    # Returns:
    #   A numeric value representing the RMSE between x and y.
    #
    # Details:
    #   The RMSE is a measure of the average deviation between the observed values (x)
    #   and the predicted values (y). A lower RMSE indicates a better model fit.
    #
    # Example:
    #   observed <- c(2.0, 3.0, 4.0, 5.0, 6.0)
    #   predicted <- c(1.8, 3.2, 4.2, 5.1, 5.8)
    #   rmse_value <- rmse(observed, predicted)
    #   cat("Root Mean Square Error:", rmse_value)
    #
    # References:
    #   - Wikipedia: https://en.wikipedia.org/wiki/Root-mean-square_deviation
    #
    # See Also:
    #   - mean: Compute the mean of numeric values.
    #   - sd: Compute the standard deviation of numeric values.
    #
    # Note:
    #   This function assumes that the lengths of x and y are the same.    sqrt(mean((x - y)^2))

    return(sqrt(mean((x - y)^2)))
}




bayes_est <- function(fit) {
    # Compute Bayesian estimates for y, ddy, and tau from a fitted model.
    #
    # Args:
    #   fit: A list or data structure containing Bayesian model output with the following elements:
    #        - y: A MCMC list  of sampled y-values from the model.
    #        - ddy: A MCMC list of sampled second derivatives.
    #        - tau: A MCMC list of sampled precision or noise estimates.
    #
    # Returns:
    #   A list containing the following Bayesian estimates and confidence intervals:
    #     - y_hat: Bayesian estimate of y (median).
    #     - y_ci: 90% Credible interval for y (5th and 95th percentiles).
    #     - ddy_hat: Bayesian estimate of ddy (median).
    #     - ddy_ci: 90% Credible interval for ddy (5th and 95th percentiles).
    #     - tau_hat: Bayesian estimate of tau (median).
    #     - tau_ci: 90% Credible interval for tau (5th and 95th percentiles).
    #
    # Details:
    #   This function takes a fitted Bayesian model (e.g., from a Markov Chain Monte Carlo simulation) and calculates Bayesian estimates and confidence intervals for the model parameters, including y (predicted values), ddy (second derivatives), and tau (precision).
    #
    # Example:
    #   fit <- list(
    #     y = matrix(rnorm(1000), ncol = 100),  # Simulated y values
    #     ddy = matrix(rnorm(1000), ncol = 100),  # Simulated second derivatives
    #     tau = matrix(runif(1000), ncol = 100)  # Simulated precision estimates
    #   )
    #   estimates <- bayes_est(fit)
    #

    return(list(
        y_hat = apply(fit$y, 2, quantile, 0.5),
        y_ci = cbind(apply(fit$y, 2, quantile, 0.05), apply(fit$y, 2, quantile, 0.95)),
        ddy_hat = apply(fit$ddy, 2, quantile, 0.5),
        ddy_ci = cbind(apply(fit$ddy, 2, quantile, 0.05), apply(fit$ddy, 2, quantile, 0.95)),
        tau_hat = apply(fit$tau, 2, quantile, 0.5),
        tau_ci = cbind(apply(fit$tau, 2, quantile, 0.05), apply(fit$tau, 2, quantile, 0.95))
    ))
}


ess_est <- function(fit, B, ddB) {
    # Compute Effective Sample Sizes (ESS) and Gelman-Rubin R-hat statistics for a Bayesian model.
    #
    # Args:
    #   fit: A fitted Bayesian model object.
    #   B: A BSpline Basis object containing B-Spline matrix B$matrix
    #   ddB: A BSpline basis object for the 2nd deriative of B
    # Returns:
    #   A list containing the following ESS and R-hat statistics:
    #     - ess_y: ESS for y (observed values).
    #     - ess_ddy: ESS for ddy (second derivatives).
    #     - ess_tau: ESS for tau (precision or noise estimates).
    #     - ess_sigma: ESS for sigma.
    #     - rhat_y: Gelman-Rubin R-hat statistic for y.
    #     - rhat_ddy: Gelman-Rubin R-hat statistic for ddy.
    #     - rhat_tau: Gelman-Rubin R-hat statistic for tau.
    #     - rhat_sigma: Gelman-Rubin R-hat statistic for sigma.
    #
    # Details:
    #   This function calculates Effective Sample Sizes (ESS) and Gelman-Rubin R-hat statistics for various model parameters, including y (observed values), ddy (second derivatives), tau (precision), and sigma. ESS measures the effective number of independent samples, and R-hat assesses the convergence of MCMC chains.
    #
    # Example:
    #   fit <- ...  # A fitted Bayesian model
    #   B <- ...    # B-spline basis matrix
    #   ddB <- ...  # Second derivatives of B-spline basis functions
    #   statistics <- ess_est(fit, B, ddB)
    #
    # References:
    #   - Gelman-Rubin diagnostic: https://en.wikipedia.org/wiki/Gelman%E2%80%93Rubin_diagnostic
    #

    nchains <- length(fit$samples)
    ysamps <- list(nchains)
    ddysamps <- list(nchains)
    tausamps <- list(nchains)
    for (i in 1:nchains) {
        theta_hat <- fit$samples[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
        yhat <- matrix(NA, nrow(theta_hat), nrow(B$matrix))
        ddyhat <- matrix(NA, nrow(theta_hat), nrow(B$matrix))
        tauhat <- matrix(NA, nrow(theta_hat), dim(fit$tau)[2])
        for (j in 1:nrow(theta_hat)) {
            yhat[j, ] <- B$matrix %*% theta_hat[j, ]
            ddyhat[j, ] <- ddB$matrix %*% theta_hat[j, ]
            if (dim(fit$tau)[2] > 1) {
                tauhat[j, ] <- fit$samples[[i]][j, paste0("tau[", 1:dim(fit$tau)[2], "]")]
            } else {
                tauhat[j, ] <- fit$samples[[i]][j, "tau"]
            }
        }
        colnames(yhat) <- paste0("y[", 1:dim(B$matrix)[1], "]")
        colnames(ddyhat) <- paste0("ddy[", 1:dim(B$matrix)[1], "]")
        if (dim(fit$tau)[2] > 1) {
            colnames(tauhat) <- paste0("tau[", 1:dim(fit$tau)[2], "]")
        } else {
            colnames(tauhat) <- "tau"
        }
        ysamps[[i]] <- yhat
        ddysamps[[i]] <- ddyhat
        tausamps[[i]] <- tauhat
    }
    class(ysamps) <- "mcmc.list"
    class(ddysamps) <- "mcmc.list"
    class(tausamps) <- "mcmc.list"

    ess_y <- effectiveSize(ysamps)
    ess_ddy <- effectiveSize(ddysamps)
    ess_tau <- effectiveSize(tausamps)
    ess_sigma <- effectiveSize(fit$samples[, "sigma"])

    rhat_y <- gelman.diag(ysamps, autoburnin = FALSE, multivariate = FALSE)
    rhat_ddy <- gelman.diag(ddysamps, autoburnin = FALSE, multivariate = FALSE)
    rhat_tau <- gelman.diag(tausamps, autoburnin = FALSE, multivariate = FALSE)

    rhat_sigma <- gelman.diag(fit$samples[, "sigma"], autoburnin = FALSE, multivariate = FALSE)

    return(list(
        ess_y = ess_y,
        ess_ddy = ess_ddy,
        ess_tau = ess_tau,
        ess_sigma = ess_sigma,
        rhat_y = rhat_y,
        rhat_ddy = rhat_ddy,
        rhat_tau = rhat_tau,
        rhat_sigma = rhat_sigma
    ))
}
