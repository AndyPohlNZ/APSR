#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
# Define UI for application that draws a histogram
fluidPage(

  # Application title
  navbarPage("Adaptive P-Splines for Biomechanics",
    id = "mainPage",
    footer = column(12,
      align = "center",
      hr(),
      "(c) 2023 Andy Pohl | Licensed under ",
      tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.en.html", target = "_blank", "GNU GPLv3"),
      "| download the source",
      tags$a(href = "www.github.com/AndyPohlNZ", target = "_blank", "here"), "." # TODO update github
    ),

    # Alignment inline field tag (use tags$div(id = "inline", ....))
    # Upload Data tab
    tabPanel(
      title = "Data", value = "uploadTab",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            h3("Upload Data"),
            fileInput("file1", "Choose CSV File",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
              )
            ),
            radioButtons("header", "File Header:",
              c(Yes = TRUE, No = FALSE),
              inline = TRUE
            ),
            radioButtons("sep", "Separator:",
              c(
                Comma = ",",
                Semicolon = ";",
                Tab = "\t"
              ),
              ",",
              inline = TRUE
            ),
            selectInput("t_col", "Time Variable", ""),
            selectInput("y_col", "Signal Variable", "", selected = ""),
          ),
          fluidRow(
            align = "center",
            tags$hr(style = "border-color: #1f1e1c;"),
            actionButton("upload2basis", label = div("Define Basis", icon("chevron-right"))),
          )
        ),
        mainPanel(
          fluidRow(column(
            12,
            plotOutput("rawDataPlot")
          )),
          fluidRow(column(
            12,
            DT::dataTableOutput("rawDataTable", height = 500)
          ))
        )
      )
    ),

    # Basis Setup tab
    tabPanel(
      title = "Basis", value = "basisTab",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            h3("Basis Parameters"),
            numericInput("deg", "Basis Degree (d)", 5, min = 1, max = 7, step = 1),
            sliderInput("nK", "Basis Size (K)", min = 2, max = 500, value = 100, step = 1),
          ),
          fluidRow(
            align = "center",
            tags$hr(style = "border-color: #1f1e1c;"),
            actionButton("basis2upload", label = div(icon("chevron-left"), "Upload Data")),
            actionButton("basis2prior", label = div("Specify Priors", icon("chevron-right"))),
          )
        ),
        mainPanel(
          plotOutput("basisSizePlot")
        )
      )
    ),
    # Prior Tab
    tabPanel(
      title = "Priors", value = "priorTab",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            h3("Model Class"),
            selectInput("model_class", NULL, c("P-Spline", "APS_ind", "APS_ar", "APS_spl"), selected = "P-Spline", width = "800px"),
            tags$hr(style = "border-color: #1f1e1c;"),
            h3("Prior Parmaeters"),
            # PSpline
            conditionalPanel(
              condition = "input.model_class == 'P-Spline'",
              splitLayout(withMathJax(helpText("\\(\\xi:\\)")),
                numericInput("pspl_xi", NULL, value = 0.002, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              )
            ),
            # APSind
            conditionalPanel(
              condition = "input.model_class == 'APS_ind'",
              splitLayout(withMathJax(helpText("\\(\\delta_{\\xi}\\):")),
                numericInput("aps_ind_dxi", NULL, value = 1.65, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              ),
              splitLayout(withMathJax(helpText("\\(\\delta_{\\alpha_{0}}\\):")),
                numericInput("aps_ind_da0", NULL, value = 0.002, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              )
            ),

            # APSar
            conditionalPanel(
              condition = "input.model_class == 'APS_ar'",
              splitLayout(withMathJax(helpText("\\(\\delta_{\\xi}\\):")),
                numericInput("aps_ar_dxi", NULL, value = 1.250, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              ),
              splitLayout(withMathJax(helpText("\\(\\delta_{\\alpha_{0}}\\):")),
                numericInput("aps_ar_da0", NULL, value = 0.002, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              )
            ),
            # APSspl
            conditionalPanel(
              condition = "input.model_class == 'APS_spl'",
              splitLayout(withMathJax(helpText("\\(\\delta_{\\xi}\\):")),
                numericInput("aps_spl_dxi", NULL, value = 1.000, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              ),
              splitLayout(withMathJax(helpText("\\(\\delta_{\\alpha_{0}}\\):")),
                numericInput("aps_spl_da0", NULL, value = 0.002, min = 0.0, max = 10.0, step = 0.001),
                cellWidths = c("20%", "80%")
              ),
              splitLayout(withMathJax(helpText("\\(K_C\\):")),
                sliderInput("aps_spl_Kc", NULL, value = 10, min = 3, max = 20, step = 1),
                cellWidths = c("20%", "80%")
              )
            ),
          ),
          fluidRow(
            align = "center",
            tags$hr(style = "border-color: #1f1e1c;"),
            actionButton("prior2basis", label = div(icon("chevron-left"), "Upload Data")),
            actionButton("prior2fit", label = div("Fit Model", icon("chevron-right"))),
          )
        ),
        mainPanel(
          fluidRow(
            h3("Model Description:"),
            conditionalPanel(
              condition = "input.model_class == 'P-Spline'",
              uiOutput("pSplineModel")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_ind'",
              uiOutput("APSindModel")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_ar'",
              uiOutput("APSarModel")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_spl'",
              uiOutput("APSsplModel")
            ),
            h4("Priors:"),
            conditionalPanel(
              condition = "input.model_class == 'P-Spline'",
              uiOutput("pSplinePrior")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_ind'",
              uiOutput("APSindPrior")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_ar'",
              uiOutput("APSarPrior")
            ),
            conditionalPanel(
              condition = "input.model_class == 'APS_spl'",
              uiOutput("APSsplPrior")
            ),
          )
        )
      )
    ),

    # Model Fitting Tab
    tabPanel(
      title = "Model", value = "modelTab",
      sidebarLayout(
        sidebarPanel(
          h3("MCMC Parmeters"),
          splitLayout(helpText("nChains:"),
            numericInput("mcmc_nchains", NULL, value = 4, min = 2, max = 100, step = 1),
            cellWidths = c("35%", "65%")
          ),
          # splitLayout(helpText("Parallel:"),
          #   div(
          #     style = "margin-top : -5px;",
          #     checkboxInput("mcmc_parallel", NULL, value = TRUE, width = "100%")
          #   ),
          #   cellWidths = c("35%", "65%")
          # ),
          splitLayout(helpText("nWarmup:"),
            numericInput("mcmc_nwarmup", NULL, value = 5000, min = 100, max = 1e6, step = 1),
            cellWidths = c("35%", "65%")
          ),
          splitLayout(helpText("nSampling:"),
            numericInput("mcmc_nsampling", NULL, value = 10000, min = 100, max = 1e6, step = 1),
            cellWidths = c("35%", "65%")
          ),
          splitLayout(helpText("thin:"),
            numericInput("mcmc_thin", NULL, value = 1, min = 1, max = 100, step = 1),
            cellWidths = c("35%", "65%")
          ),
          fluidRow(
            align = "center",
            tags$hr(style = "border-color: #1f1e1c;"),
            actionButton("mcmc2prior", label = div(icon("chevron-left"), "Specify Prior")),
            actionButton("run_sampler", label = "Sample Model"),
            actionButton("mcmc2results", label = div("Check Results", icon("chevron-right"))),
          )
        ),
        mainPanel(
          uiOutput("MCMCMainPane")
        )
      )
    ),
    tabPanel(
      title = "Results", value = "resultsTab",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            h3("Results"),
            withMathJax(),
            checkboxGroupInput("exportOpt", "Include",
              selected = c("y", "dy", "ddy", "tau", "sigma"),
              choiceNames = list(
                "\\(y(t)\\)", "\\(\\dot{y}(t)\\)", "\\(\\ddot{y}(t)\\)",
                "Smoothing: \\(\\tau(t)\\)",
                "Measurement Noise: \\(\\sigma\\)"
              ),
              choiceValues = list("y", "dy", "ddy", "tau", "sigma")
            )
          ),
          fluidRow(
            radioButtons("exportCI", "Include Credible Interval:",
              c(Yes = TRUE, No = FALSE),
              inline = TRUE
            ),
            splitLayout(helpText("% CI"),
              numericInput("exportCIAlpha", NULL, value = 90, min = 1, max = 99, step = 1),
              cellWidths = c("15%", "85%")
            ),
          ),
          fluidRow(
            align = "center",
            tags$hr(style = "border-color: #1f1e1c;"),
            downloadButton("exportResults", label = "Export Results", style = "width:400px"),
          ),
          fluidRow(
            align = "center",
            downloadButton("exportFit", label = "Export Fit Object", style = "width:400px")
          ),
          fluidRow(
            align = "center",
            actionButton("results2mcmc", label = div(icon("chevron-left"), "Fit Model"), style = "width:400px"),
          )
        ),
        mainPanel(
          uiOutput("ResultsMainPane")
        )
      )
    )
  )
)
