#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
if (!require("rjags")) install.packages("rjags")
if (!require("coda")) install.packages("coda")
if (!require("bayesplot")) install.packages("bayesplot")
if (!require("ggplot2")) install.packages("ggplot2")

library(shiny)
library(DT)
library(RColorBrewer)
library(shinycssloaders)
library(rjags)
library(posterior)
library(Rcpp)
Rcpp::sourceCpp("./R/mm_likelihood.cpp")

# Define server logic required to draw a histogram
function(input, output, session) {
  state <- reactiveValues(t_obs = NULL, y_obs = NULL, fit = NULL, output = NULL)

  # Set color palette
  pal <- brewer.pal(5, "Set1") # Set color pallette for plots

  # Hide tabs until data loaded
  hideTab(inputId = "mainPage", target = "basisTab")
  hideTab(inputId = "mainPage", target = "priorTab")
  hideTab(inputId = "mainPage", target = "modelTab")
  hideTab(inputId = "mainPage", target = "resultsTab")

  data <- reactive({
    if (is.null(input$file1)) {
      return(NULL)
    }

    df <- read.csv(input$file1$datapath, header = as.logical(input$header))

    updateSelectInput(session,
      inputId = "t_col", label = "Time Variable",
      choices = names(df), selected = names(df)[1]
    )
    updateSelectInput(session,
      inputId = "y_col", label = "Y Variable",
      choices = names(df), selected = names(df)[2]
    )
    return(df)
  })


  output$rawDataTable <- renderDT(
    {
      req(data()) # if user_data() is null than the reactive will not run
      data()
    },
    options = list(scrollX = TRUE, scrollY = TRUE, paging = FALSE, dom = "lt"),
  )

  output$rawDataPlot <- renderPlot({
    req(data()) # if data() is null than the reactive will not run
    dat <- data()
    ttmp <- as.numeric(dat[, input$t_col])
    ytmp <- as.numeric(dat[, input$y_col])

    # x <- data()[, c(input$t_col, input$y_col)]
    plot(
      x = ttmp, y = ytmp, xlab = "Time (s)", ylab = "y(t)", main = "Observed Data",
      pch = 16, col = adjustcolor("black", alpha = 0.2)
    )

    state$t_obs <- ttmp
    state$y_obs <- ytmp
  })


  # On click Basis
  observeEvent(input$upload2basis, {
    updateTabsetPanel(session, "mainPage",
      selected = "basisTab"
    )

    showTab("mainPage", "basisTab")

    x <- data()

    updateSliderInput(session,
      "nK",
      min = 5,
      max = floor(nrow(x) / 2),
      value = min(floor(nrow(x) / 4), 100)
    )
  })


  output$basisSizePlot <- renderPlot({
    deg <- input$deg

    # Load Data
    # dat <- data()
    # time_vals <- as.numeric(dat[, input$t_col])
    # y_vals <- as.numeric(dat[, input$y_col])

    time_vals <- state$t_obs
    y_vals <- state$y_obs
    y_std <- y_vals - mean(y_vals)
    t_std <- time_vals - mean(time_vals)

    nKs <- 10
    idx <- 2
    while (2 * nKs[idx - 1] < floor(length(t_std) / 1.5)) {
      nKs <- c(nKs, 2 * nKs[idx - 1])
      idx <- idx + 1
    }

    nsizes <- length(nKs)

    rss <- rep(NA, nsizes)
    # edf <- rep(NA, nsizes)
    withProgress(message = "Evaluating Basis Size", value = 0, {
      for (i in 1:nsizes) {
        Btmp <- basis(t_std, nKs[i], deg)
        Ptmp <- penalty(dim(Btmp$matrix)[2], 2)
        MMtmp <- basis_mm(Ptmp)
        mle <- mle_pspl(Btmp, MMtmp, y_std)
        y_hat <- Btmp$matrix %*% mle$theta_hat
        rss[i] <- sum((y_std - y_hat)^2)

        # D = Ptmp$D
        # S <- t(D) %*% diag(rep(mle$tau_hat^2, dim(D)[1])) %*% D
        # print(S[1:5, 1:5])
        # A <- Btmp$matrix %*% solve(t(Btmp$matrix) %*% Btmp$matrix + S) %*% t(Btmp$matrix)
        # print("=====")
        # print(A[1:5, 1:5])
        # edf[i] <- sum(diag(A))
        #

        print(paste0("i = ", i, " nK = ", nKs[i], " rss = ", rss[i]))
        incProgress(1 / nsizes, detail = paste("Processing Basis size", nKs[i]))
      }

      plot(nKs, log(rss),
        xlab = "Basis Size", ylab = "log(SSE)",
        type = "o", col = pal[2], lty = 2
      )
      # plot(nKs, edf, xlab = "Basis Size", ylab = "EDF",
      #      type = 'o', col = pal[2], lty = 2)
      # abline(a = deg, b = 1)
    })
  })

  observeEvent(input$basis2prior, {
    nK <- input$nK
    updateTabsetPanel(session, "mainPage",
      selected = "priorTab"
    )

    showTab("mainPage", "priorTab")
    updateSliderInput(session, "aps_spl_Kc",
      min = 3, max = floor(nK / 2),
      value = floor(nK / 5),
      step = 1
    )
  })

  # Model Descriptions
  output$pSplineModel <- renderUI({
    withMathJax(helpText("Observations \\(y_t\\) over discrete time values \\(t = 1, \\ldots, T\\) are modelled via a \\(d\\)th order B-Spline:
    $$y_t = \\left[ B\\theta \\right]_t + \\varepsilon_t \\quad \\quad \\varepsilon_t \\sim \\text{Normal}\\left(0, \\sigma^2\\right)$$
    where \\(B\\) is a \\(T \\times (K + d)\\) matrix of B-Splines with associated  coefficients \\(\\theta\\) and \\(\\varepsilon_t\\) is Gaussian white noise with variance \\(\\sigma^2\\).
    A 2nd order random walk is used for spline coefficients:
                         $$ \\theta_k = 2\\theta_{k-1} - \\theta_{k-2} + \\nu_k, \\quad k = 3, \\ldots, K$$
    where \\(\\nu_k\\) is a random variable with variance \\(\\tau^2\\):
                         $$\\nu_k \\sim \\text{Normal}\\left(0, \\tau^2\\right)$$"))
  })
  output$pSplinePrior <- renderUI({
    withMathJax(helpText("A centred, truncated \\(t\\) distribution with 3 degrees of freedom and scale \\(\\xi\\) provides a weekly regularizing prior for unknown smoothing parameter \\(\\tau^2\\):
                          $$ \\tau \\sim t^+_3\\left(0, \\xi^2\\right)$$
                          Remaining paramters are given uninformative priors:
                         $$ \\theta_1, \\theta_2 \\sim \\text{Normal}\\left(0, 100^2\\right)$$
                         $$ \\sigma^2 \\sim \\text{Inv-Gamma}\\left(0.001, 0.001\\right)$$"))
  })

  output$APSindModel <- renderUI({
    withMathJax(helpText("Observations \\(y_t\\) over discrete time values \\(t = 1, \\ldots, T\\) are modelled via a \\(d\\)th order B-Spline:
                        $$y_t = \\left[ B\\theta \\right]_t + \\varepsilon_t \\quad \\quad \\varepsilon_t \\sim \\text{Normal}\\left(0, \\sigma^2\\right)$$
    where \\(B\\) is a \\(T \\times (K + d)\\) matrix of B-Splines with associated  coefficients \\(\\theta\\) and \\(\\varepsilon_t\\) is Gaussian white noise with variance \\(\\sigma^2\\).
    A 2nd order random walk is used for spline coefficients:
                         $$ \\theta_k = 2\\theta_{k-1} - \\theta_{k-2} + \\nu_k, \\quad k = 3, \\ldots, K$$
  Temporal variation in curve smoothness (\\(\\tau_k^2\\)) is modelled as a series of independent random normal variables with mean \\(\\log\\left(\\log(\\alpha_0^2)\\right)\\) and variance \\(\\xi^2\\):
                         $$\\log(\\tau_k^2) \\sim \\text{Normal} \\left(\\log(\\alpha_0^2), \\xi^2 \\right).$$"))
  })
  output$APSindPrior <- renderUI({
    withMathJax(helpText("The prior for \\(\\xi\\) is half \\(t\\) with 3 degrees of freedom and scale parameter \\(\\delta_{\\xi}\\):
                         $$\\xi \\sim t^{+}_3 \\left(0, \\delta_{\\xi}^2\\right)$$
                         A half-t prior with 3 degrees of freedom and scale \\(\\delta_{\\alpha_0}\\) is specified for average smoothness:
                        $$\\alpha_0 \\sim t^{+}_3 \\left(0, \\delta_{\\alpha_0}^2\\right)$$
                         Uninformative priors are used for remaining parameters:
                         $$ \\theta_1, \\theta_2 \\sim \\text{Normal}\\left(0, 100^2\\right)$$
                         $$ \\sigma^2 \\sim \\text{Inv-Gamma}\\left(0.001, 0.001\\right)$$"))
  })
  output$APSarModel <- renderUI({
    withMathJax(helpText("Observations \\(y_t\\) over discrete time values \\(t = 1, \\ldots, T\\) are modelled via a \\(d\\)th order B-Spline:
                        $$y_t = \\left[ B\\theta \\right]_t + \\varepsilon_t \\quad \\quad \\varepsilon_t \\sim \\text{Normal}\\left(0, \\sigma^2\\right)$$
    where \\(B\\) is a \\(T \\times (K + d)\\) matrix of B-Splines with associated  coefficients \\(\\theta\\) and \\(\\varepsilon_t\\) is Gaussian white noise with variance \\(\\sigma^2\\).
    A 2nd order random walk is used for spline coefficients:
                         $$ \\theta_k = 2\\theta_{k-1} - \\theta_{k-2} + \\nu_k, \\quad k = 3, \\ldots, K$$
    Parameters \\(\\nu_k\\) are assumed normally distributed with differing variance \\(\\tau^2_k\\):
                         $$\\nu_k \\sim \\text{Normal}\\left(0, \\tau_k^2\\right)$$
    Temporal variation in curve smoothness is modelled on the log scale with variation about average smoothness \\(\\alpha_0^2\\) described via a first order autoregressive process:
                         $$\\log(\\tau_k^2) \\sim \\text{Normal} \\left(\\log(\\alpha_0^2) + \\phi\\left(\\log(\\tau_{k-1}^2) - \\log(\\alpha_0^2)\\right), \\xi^2 \\right)$$"))
  })

  output$APSarPrior <- renderUI({
    withMathJax(helpText("The prior for \\(\\xi\\) is half \\(t\\) with 3 degrees of freedom and scale parameter \\(\\delta_{\\xi}\\):
                         $$\\xi \\sim t^{+}_3 \\left(0, \\delta_{\\xi}^2\\right)$$
                         A half-t prior with 3 degrees of freedom and scale \\(\\delta_{\\alpha_0}\\) is specified for average smoothness:
                        $$\\alpha_0 \\sim t^{+}_3 \\left(0, \\delta_{\\alpha_0}^2\\right)$$
                         Correleation between adjacent smoothness terms is assumed to be positive and an uninformative prior is specified:
                         $$\\phi \\sim \\text{Uniform}(0, 1)$$
                         Uninformative priors are used for remaining parameters:
                         $$ \\theta_1, \\theta_2 \\sim \\text{Normal}\\left(0, 100^2\\right)$$
                         $$ \\sigma^2 \\sim \\text{Inv-Gamma}\\left(0.001, 0.001\\right)$$"))
  })


  output$APSsplModel <- renderUI({
    withMathJax(helpText("Observations \\(y_t\\) over discrete time values \\(t = 1, \\ldots, T\\) are modelled via a \\(d\\)th order B-Spline:
                        $$y_t = \\left[ B\\theta \\right]_t + \\varepsilon_t \\quad \\quad \\varepsilon_t \\sim \\text{Normal}\\left(0, \\sigma^2\\right)$$
    where \\(B\\) is a \\(T \\times (K + d)\\) matrix of B-Splines with associated  coefficients \\(\\theta\\) and \\(\\varepsilon_t\\) is Gaussian white noise with variance \\(\\sigma^2\\).
    A 2nd order random walk is used for spline coefficients:
                         $$ \\theta_k = 2\\theta_{k-1} - \\theta_{k-2} + \\nu_k, \\quad k = 3, \\ldots, K$$
    Parameters \\(\\nu_k\\) are assumed normally distributed with differing variance \\(\\tau^2_k\\):
                         $$\\nu_k \\sim \\text{Normal}\\left(0, \\tau_k^2\\right)$$
    Temporal variation in curve smoothness is determined through a cubic B-Spline on the log scale:
                         $$\\log(\\tau_k^2) = C\\gamma$$
                         where \\(C\\) is a \\(K \\times \\left(K_C+3\\right)\\) cubic B-Spline basis with \\(K_C << K\\) equidistant knots and corresponding parameters \\(\\gamma\\)
                         The distribution of \\(\\gamma\\) for \\(l\\) in (\\(1, \\ldots, K_C\\)) is assumed normal with mean \\(\\log(\\alpha_0^2)\\) and variance \\(\\xi^2\\)."))
  })

  output$APSsplPrior <- renderUI({
    withMathJax(helpText("The prior for \\(\\xi\\) is half \\(t\\) with 3 degrees of freedom and scale parameter \\(\\delta_{\\xi}\\):
                         $$\\xi \\sim t^{+}_3 \\left(0, \\delta_{\\xi}^2\\right)$$
                         A half-t prior with 3 degrees of freedom and scale \\(\\delta_{\\alpha_0}\\) is specified for average smoothness:
                        $$\\alpha_0 \\sim t^{+}_3 \\left(0, \\delta_{\\alpha_0}^2\\right)$$
                         Uninformative priors are used for remaining parameters:
                         $$ \\theta_1, \\theta_2 \\sim \\text{Normal}\\left(0, 100^2\\right)$$
                         $$ \\sigma^2 \\sim \\text{Inv-Gamma}\\left(0.001, 0.001\\right)$$"))
  })


  # Switch to MCMC tab
  observeEvent(input$prior2fit, {
    nK <- input$nK
    updateTabsetPanel(session, "mainPage",
      selected = "modelTab"
    )

    showTab("mainPage", "modelTab")
  })
  # Initialize MCMC Pane
  output$MCMCMainPane <-
    renderUI({
      fluidRow(
        align = "center",
        helpText("No model output ...")
      )
    })

  # Initilize MCMC buttons
  output$MCMCPaneButtons <-
    renderUI({
      fluidRow(
        align = "center",
        actionButton("mcmc2prior", label = div(icon("chevron-left"), "Specify Prior")),
        actionButton("run_sampler", label = "Fit Model"),
        actionButton("mcmc2results", label = div("Check Results", icon("chevron-right"))),
      )
    })

  # Run sampler
  observeEvent(
    input$run_sampler,
    {
      # Fit Jags
      model_class <- input$model_class
      bDeg <- input$deg
      nK <- input$nK

      # Get Data
      # dat <- data()
      time_vals <- state$t_obs
      y_vals <- state$y_obs

      y_mu <- mean(y_vals)
      y_sd <- sd(y_vals)
      y_std <- (y_vals - y_mu) / y_sd
      t_std <- time_vals - mean(time_vals)

      # Construct basis
      B <- basis(t_std, nK, bDeg)
      dB <- dbasis(B, 1)
      ddB <- dbasis(B, 2)
      P <- penalty(dim(B$matrix)[2], 2)
      MM <- basis_mm(P)

      if (model_class == "APS_spl") {
        nKC <- input$aps_spl_Kc
        t_c <- seq(0, 1, length.out = dim(B$matrix)[2] - P$deg)
        C <- basis(t_c, nKC, deg = 3)
      }


      # Get mle for standard PSpline (used for inilization)
      mle <- mle_pspl(B, MM, y_std)

      out <- NULL
      parameter_names <- NULL

      # MCMC parameters
      nchains <- input$mcmc_nchains
      parallel <- input$mcmc_parallel
      nwarmup <- input$mcmc_nwarmup
      nsampling <- input$mcmc_nsampling
      thin <- input$mcmc_thin



      if (model_class == "P-Spline") {
        xi <- input$pspl_xi
        print("Sampling P-Spline Model")
        init_fn <- function() {
          mle <- mle_pspl(B, MM, y_std)
          theta <- as.numeric(mle$theta_hat)
          list(
            theta = theta
          )
        }

        withProgress(message = "Compiling JAGS model...", value = 0, {
          incProgress(0.5)

          mod <- jags.model("./mdl/PSpline.jags",
            data = list(
              T = length(t_std),
              Q = dim(B$matrix)[2],
              B = B$matrix,
              y = y_std,
              xi = xi
            ),
            inits = init_fn,
            n.chains = nchains,
            n.adapt = nwarmup
          )
        })
        withProgress(message = "Sampling from JAGS model...", value = 0, {
          incProgress(0.25)
          out <- coda.samples(mod,
            c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
          )
          incProgress(0.75)
          nchains <- length(out)
          nsamps <- dim(out[[1]])[1]
          nsamps_total <- nsamps * nchains
          theta_samps <- matrix(NA, nrow = nsamps_total, ncol = dim(B$matrix)[2])
          for (i in 1:nchains) {
            theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
          }


          y_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          dy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          ddy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nsamps_total) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
          }

          tau_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          sigma_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nchains) {
            tau_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- rep(out[[i]][, "tau"], length(time_vals))
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- rep(out[[i]][, "sigma"], length(time_vals))
          }

          fit <- list(
            y = y_samps,
            dy = dy_samps,
            ddy = ddy_samps,
            tau = tau_samps,
            sigma = sigma_samps
          )
        })
      } else if (model_class == "APS_ind") {
        dxi <- input$aps_ind_dxi
        da0 <- input$aps_ind_da0
        print("Sampling APS_ind Model")
        init_fn <- function() {
          mle <- mle_pspl(B, MM, y_std)
          theta <- as.numeric(mle$theta_hat)
          list(
            theta = theta
          )
        }

        withProgress(message = "Compiling JAGS model...", value = 0, {
          incProgress(0.5)

          mod <- jags.model("./mdl/APS_ind.jags",
            data = list(
              T = length(t_std),
              Q = dim(B$matrix)[2],
              B = B$matrix,
              y = y_std,
              d_xi0 = dxi,
              d_a0 = da0
            ),
            inits = init_fn,
            n.chains = nchains,
            n.adapt = nwarmup
          )
        })
        withProgress(message = "Sampling from JAGS model...", value = 0, {
          incProgress(0.25)
          out <- coda.samples(mod,
            c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
          )
          incProgress(0.75)
          nchains <- length(out)
          nsamps <- dim(out[[1]])[1]
          nsamps_total <- nsamps * nchains
          theta_samps <- matrix(NA, nrow = nsamps_total, ncol = dim(B$matrix)[2])
          for (i in 1:nchains) {
            theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
          }
          y_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          dy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          ddy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nsamps_total) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
          }

          tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = length(t_std))
          t_tau <- seq(from = min(t_std), to = max(t_std), length.out = (dim(B$matrix)[2] - 2))
          idx <- 1
          for (i in 1:nchains) {
            for (j in 1:nsamps) {
              tau_tmp <- out[[i]][j, paste0("tau[", 1:(dim(B$matrix)[2] - 2), "]")]
              tau_samps[idx, ] <- approx(t_tau, tau_tmp, xout = t_std)$y
              idx <- idx + 1
            }
          }

          sigma_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- rep(out[[i]][, "sigma"], length(time_vals))
          }

          fit <- list(
            y = y_samps,
            dy = dy_samps,
            ddy = ddy_samps,
            tau = tau_samps,
            sigma = sigma_samps
          )
        })
      } else if (model_class == "APS_ar") {
        dxi <- input$aps_ar_dxi
        da0 <- input$aps_ar_da0
        print("Sampling APS_ind Model")
        init_fn <- function() {
          mle <- mle_pspl(B, MM, y_std)
          theta <- as.numeric(mle$theta_hat)
          list(
            theta = theta
          )
        }

        withProgress(message = "Compiling JAGS model...", value = 0, {
          incProgress(0.5)

          mod <- jags.model("./mdl/APS_ar.jags",
            data = list(
              T = length(t_std),
              Q = dim(B$matrix)[2],
              B = B$matrix,
              y = y_std,
              d_xi0 = dxi,
              d_a0 = da0
            ),
            inits = init_fn,
            n.chains = nchains,
            n.adapt = nwarmup
          )
        })
        withProgress(message = "Sampling from JAGS model...", value = 0, {
          incProgress(0.25)
          out <- coda.samples(mod,
            c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
          )
          incProgress(0.75)
          nchains <- length(out)
          nsamps <- dim(out[[1]])[1]
          nsamps_total <- nsamps * nchains
          theta_samps <- matrix(NA, nrow = nsamps_total, ncol = dim(B$matrix)[2])
          for (i in 1:nchains) {
            theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
          }
          y_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          dy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          ddy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nsamps_total) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
          }

          tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = length(t_std))
          t_tau <- seq(from = min(t_std), to = max(t_std), length.out = (dim(B$matrix)[2] - 2))
          idx <- 1
          for (i in 1:nchains) {
            for (j in 1:nsamps) {
              tau_tmp <- out[[i]][j, paste0("tau[", 1:(dim(B$matrix)[2] - 2), "]")]
              tau_samps[idx, ] <- approx(t_tau, tau_tmp, xout = t_std)$y
              idx <- idx + 1
            }
          }

          sigma_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- rep(out[[i]][, "sigma"], length(time_vals))
          }

          fit <- list(
            y = y_samps,
            dy = dy_samps,
            ddy = ddy_samps,
            tau = tau_samps,
            sigma = sigma_samps
          )
        })
      } else if (model_class == "APS_spl") {
        dxi <- input$aps_spl_dxi
        da0 <- input$aps_spl_da0
        print("Sampling APS_spl Model")
        init_fn <- function() {
          mle <- mle_pspl(B, MM, y_std)
          theta <- as.numeric(mle$theta_hat)
          list(
            theta = theta
          )
        }

        withProgress(message = "Compiling JAGS model...", value = 0, {
          incProgress(0.5)

          mod <- jags.model("./mdl/APS_spl.jags",
            data = list(
              T = length(t_std),
              Q = dim(B$matrix)[2],
              Q_sm = dim(C$matrix)[2],
              B = B$matrix,
              C = C$matrix,
              y = y_std,
              d_xi0 = dxi,
              d_a0 = da0
            ),
            inits = init_fn,
            n.chains = nchains,
            n.adapt = nwarmup
          )
        })
        withProgress(message = "Sampling from JAGS model...", value = 0, {
          incProgress(0.25)
          out <- coda.samples(mod,
            c("theta", "tau", "sigma"),
            n.iter = nsampling,
            thin = thin
          )
          incProgress(0.75)
          nchains <- length(out)
          nsamps <- dim(out[[1]])[1]
          nsamps_total <- nsamps * nchains
          theta_samps <- matrix(NA, nrow = nsamps_total, ncol = dim(B$matrix)[2])
          for (i in 1:nchains) {
            theta_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- out[[i]][, paste0("theta[", 1:dim(B$matrix)[2], "]")]
          }
          y_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          dy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          ddy_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nsamps_total) {
            y_samps[i, ] <- B$matrix %*% theta_samps[i, ] * y_sd + y_mu
            dy_samps[i, ] <- dB$matrix %*% theta_samps[i, ] * y_sd
            ddy_samps[i, ] <- ddB$matrix %*% theta_samps[i, ] * y_sd
          }

          tau_samps <- matrix(NA, nrow = nsamps * nchains, ncol = length(t_std))
          t_tau <- seq(from = min(t_std), to = max(t_std), length.out = (dim(B$matrix)[2] - 2))
          idx <- 1
          for (i in 1:nchains) {
            for (j in 1:nsamps) {
              tau_tmp <- out[[i]][j, paste0("tau[", 1:(dim(B$matrix)[2] - 2), "]")]
              tau_samps[idx, ] <- approx(t_tau, tau_tmp, xout = t_std)$y
              idx <- idx + 1
            }
          }

          sigma_samps <- matrix(NA, nrow = nsamps_total, ncol = length(time_vals))
          for (i in 1:nchains) {
            sigma_samps[((i - 1) * nsamps + 1):(i * nsamps), ] <- rep(out[[i]][, "sigma"], length(time_vals))
          }

          fit <- list(
            y = y_samps,
            dy = dy_samps,
            ddy = ddy_samps,
            tau = tau_samps,
            sigma = sigma_samps
          )
        })
      }

      parameter_names <- dimnames(out[[1]])[[2]]

      # Show Traceplot
      output$MCMCMainPane <- renderUI({
        fluidPage(
          h3("Check Traceplots:"),
          fluidRow(
            selectInput(
              "trace_selparm",
              "Select parameter:",
              parameter_names,
              selected = "sigma",
              selectize = FALSE
            )
          ),
          fluidRow(
            renderPlot({
              mcmc_combo(out,
                combo = c("trace", "dens_overlay"),
                pars = input$trace_selparm,
                gg_theme = ggplot2::theme_linedraw() + legend_none()
              )
            })
          ),
        )
      })

      # Update Fit object
      state$fit <- fit
    }
  )

  # Switch to results
  observeEvent(input$mcmc2results, {
    updateTabsetPanel(session, "mainPage",
      selected = "resultsTab"
    )
    showTab("mainPage", "resultsTab")
  })

  gen_output <- function() {
    output <- data.frame(t = state$t_obs)
    # output$t <- state$t_obs
    output$y_obs <- state$y_obs


    if ("y" %in% input$exportOpt) {
      output$y_hat <- apply(state$fit$y, 2, median)
    }
    if ("dy" %in% input$exportOpt) {
      output$dy_hat <- apply(state$fit$dy, 2, median)
    }
    if ("ddy" %in% input$exportOpt) {
      output$ddy_hat <- apply(state$fit$ddy, 2, median)
    }

    if ("tau" %in% input$exportOpt) {
      output$tau_hat <- apply(state$fit$tau, 2, median)
    }

    if ("sigma" %in% input$exportOpt) {
      output$sigma_hat <- apply(state$fit$sigma, 2, median)
    }
    # Add Intervals
    if (input$exportCI == TRUE) {
      interval_probs <- (1 - (input$exportCIAlpha / 100)) / 2

      if ("y" %in% input$exportOpt) {
        output$y_lwr <- apply(state$fit$y, 2, quantile, probs = interval_probs)
        output$y_upr <- apply(state$fit$y, 2, quantile, probs = 1 - interval_probs)
      }

      if ("dy" %in% input$exportOpt) {
        output$dy_lwr <- apply(state$fit$dy, 2, quantile, probs = interval_probs)
        output$dy_upr <- apply(state$fit$dy, 2, quantile, probs = 1 - interval_probs)
      }

      if ("ddy" %in% input$exportOpt) {
        output$ddy_lwr <- apply(state$fit$ddy, 2, quantile, probs = interval_probs)
        output$ddy_upr <- apply(state$fit$ddy, 2, quantile, probs = 1 - interval_probs)
      }
      if ("tau" %in% input$exportOpt) {
        output$tau_lwr <- apply(state$fit$tau, 2, quantile, probs = interval_probs)
        output$tau_upr <- apply(state$fit$tau, 2, quantile, probs = 1 - interval_probs)
      }

      if ("sigma" %in% input$exportOpt) {
        output$sigma_lwr <- apply(state$fit$sigma, 2, quantile, probs = interval_probs)
        output$sigma_upr <- apply(state$fit$sigma, 2, quantile, probs = 1 - interval_probs)
      }
    }

    return(output)
  }

  # Generate Figures
  output$ResultsMainPane <- renderUI({
    withProgress(message = "Generating Results...", value = 0, {
      incProgress(0.25)

      state$output <- gen_output()
      y_obs <- state$output$y_obs
      t_obs <- state$output$t
      idx <- sample(1:nrow(state$fit$y), 100)

      incProgress(0.5)
      fluidPage(
        if ("y" %in% input$exportOpt) {
          fluidRow(
            renderPlot({
              plot(state$t_obs, state$y_obs,
                pch = 16,
                ylab = "y(t)", xlab = "Time",
                col = adjustcolor("#4e4040", alpha = 0.1)
              )
              for (i in 1:length(idx)) {
                lines(state$t_obs, state$fit$y[idx[i], ], type = "l", col = adjustcolor(pal[1], alpha = 0.1))
              }
              lines(state$t_obs, state$output$y_hat, type = "l", lwd = 2, col = pal[1])
              lines(state$t_obs, state$output$y_lwr, type = "l", lwd = 1, col = pal[1], lty = 3)
              lines(state$t_obs, state$output$y_upr, type = "l", lwd = 1, col = pal[1], lty = 3)
            })
          )
        },
        if ("dy" %in% input$exportOpt) {
          fluidRow(
            renderPlot({
              plot(state$t_obs, state$output$dy_hat,
                type = "l", ylab = "y'(t)", xlab = "Time",
                col = pal[2]
              )
              for (i in 1:length(idx)) {
                lines(state$t_obs, state$fit$dy[idx[i], ], type = "l", col = adjustcolor(pal[2], alpha = 0.1))
              }
              lines(state$t_obs, state$output$dy_hat, lwd = 2, col = pal[2])
              lines(state$t_obs, state$output$dy_lwr, type = "l", lwd = 1, col = pal[2], lty = 3)
              lines(state$t_obs, state$output$dy_upr, type = "l", lwd = 1, col = pal[2], lty = 3)
            })
          )
        },
        if ("ddy" %in% input$exportOpt) {
          fluidRow(
            renderPlot({
              plot(state$t_obs, state$output$ddy_hat,
                type = "l", ylab = "y''(t)", xlab = "Time",
                col = pal[3]
              )
              for (i in 1:length(idx)) {
                lines(state$t_obs, state$fit$ddy[idx[i], ], type = "l", col = adjustcolor(pal[3], alpha = 0.1))
              }
              lines(state$t_obs, state$output$ddy_hat, lwd = 2, col = pal[3])
              lines(state$t_obs, state$output$ddy_lwr, type = "l", lwd = 1, col = pal[3], lty = 3)
              lines(state$t_obs, state$output$ddy_upr, type = "l", lwd = 1, col = pal[3], lty = 3)
            })
          )
        },
        if ("tau" %in% input$exportOpt) {
          fluidRow(
            renderPlot({
              plot(state$t_obs, state$output$tau_hat,
                type = "l", ylab = "tau(t)", xlab = "Time",
                col = pal[3]
              )
              for (i in 1:length(idx)) {
                lines(state$t_obs, state$fit$tau[idx[i], ], type = "l", col = adjustcolor(pal[4], alpha = 0.1))
              }
              lines(state$t_obs, state$output$tau_hat, lwd = 2, col = pal[4])
              lines(state$t_obs, state$output$tau_lwr, type = "l", lwd = 1, col = pal[4], lty = 3)
              lines(state$t_obs, state$output$tau_upr, type = "l", lwd = 1, col = pal[4], lty = 3)
            })
          )
        },
        if ("sigma" %in% input$exportOpt) {
          fluidRow(
            renderPlot({
              plot(state$t_obs, state$output$sigma_hat,
                type = "l", ylab = "sigma(t)", xlab = "Time",
                col = pal[3]
              )
              for (i in 1:length(idx)) {
                lines(state$t_obs, state$fit$sigma[idx[i], ], type = "l", col = adjustcolor(pal[5], alpha = 0.1))
              }
              lines(state$t_obs, state$output$sigma_hat, lwd = 2, col = pal[5])
              lines(state$t_obs, state$output$sigma_lwr, type = "l", lwd = 1, col = pal[5], lty = 3)
              lines(state$t_obs, state$output$sigma_upr, type = "l", lwd = 1, col = pal[5], lty = 3)
            })
          )
        },
      )
    })
  })

  # Export Data button
  output$exportResults <- downloadHandler(
    filename = function() {
      paste0("APS_results", ".csv")
    },
    content = function(file) {
      write.csv(state$output, file)
    }
  )

  output$exportFit <- downloadHandler(
    filename = function() {
      paste0("APS_fit", ".RDS")
    },
    content = function(file) {
      saveRDS(state$fit, file)
    }
  )
}
