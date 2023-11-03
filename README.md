# APSR

An R package for implementation of adaptive P-splines for challenging signal processing problems in biomechanics. This code currently fits a selection of splines with adaptive penalty terms to noisy signals with rapidally changing deriatives.  Traditional signal processing techniques such as recursive filters or splines with a fixed penalty can struggle to accurately capture the underlying signal in these situations.  The adaptive penalty terms allow the spline to adapt to the local signal characteristics and provide a more accurate representation of the underlying signal and its derivatives.

The majority of  methods considered are Bayesian in nature and rely on Markov Chain Monte Carlo methods for inference.  A combination of Gibbs and slice sampling is used via the [JAGS](https://mcmc-jags.sourceforge.io/)   and its acompanying [R interface](https://cran.r-project.org/web/packages/rjags/index.html).  `JAGS` should be installed prior to using the code in this repository.

A walkthrough of model fitting is provided in the `APS_Fit.R` file, with `APS_Analysis.R` containing additonal frequently used visualisations which form part of the acompanying manuscript.  In addition a GUI has been developed to allow for interactive model fitting and visualisation of results.  This is provided in the `/app` directory and can be run using the following command:

```r
if.packages("shiny") install.packages("shiny"); library("shiny")
shiny::runApp("./app/")
```

This package is currently under development and is not yet registered.  While the package is designed for use for common data-processing problems in biomechancis (noteably inverse kinematics and inverse dynamics applications) its use is not limited to this application.

AdPSpline was developed by:

[Andrew J. Pohl](https://andypohlnz.github.io/)
andrew.pohl@ucalgary.ca

If you wish to contribute to the development of this package or identify any issues please raise these via an issue or pull request on the github repository.

## Installation
`APSR` requires a working installation of [R](https://www.r-project.org/) and was tested on R version 4.2.1.  Please see the [R](https://www.r-project.org/) documetation for installation instructions.

Importantly `APSR` relies on a working instation of [JAGS](https://mcmc-jags.sourceforge.io/) and its acompanying [R interface](https://cran.r-project.org/web/packages/rjags/index.html).  JAGS version 4.3.1 was used for development and testing.

## Publications
A peer reviewed manuscript detailing this work is currently under review and will be liked once accepted for publication.

## Acknowledgements
This work was supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) and an Alberta Graduate Scholarship. The author would like to thank the following people for their contributions to this work:

Dr. Matthew Schofield
Department of Mathematics and Statistics, Otago University

Dr. Reed Ferber
Faculty of Kinesiology, University of Calgary

Dr W. Brent Edwards
Faculty of Kinesiology, University of Calgary