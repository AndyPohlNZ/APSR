## ==============================================================================
# run_impactSim.R
#
# Simulation of a point mass impacting a compliant surface

# The model is a simple point mass with a spring and damper in parallel
# with the mass. The mass is dropped from a height of 1m and the impact
# is simulated. The model is solved using the deSolve package.

# Created by:   Andy Pohl
#               Feb 2024
#               andrew.pohl@ucalgary.ca
# ==============================================================================


library("deSolve")
gen_impact_data <- function(
    times = seq(0, 5, by = 0.01),
    parms = list(
        g = 9.81,
        m = 100.0, k = 25000, c = 1500
    ),
    inits = c(y = 1.0, dy = 5)) {
    # ODE function
    dydt <- function(t, y, parms) {
        if (y[1] > 0) {
            dy <- y[2]
            dy2 <- -parms$g
            res <- list(c(dy, dy2))
        } else {
            dy <- y[2]
            dy2 <- -(parms$k / parms$m) * y[1] - (parms$c / parms$m) * y[2] - parms$g
            res <- list(c(dy, dy2))
        }
        return(res)
    }
    # Force
    force <- function(t, y, dy, parms) {
        res <- rep(NA, length(y))
        for (i in 1:length(y)) {
            if (y[i] >= 0 + 1e-5) {
                res[i] <- -parms$m * parms$g
            } else {
                res[i] <- -parms$k * y[i] - parms$c * dy[i] - parms$m * parms$g
            }
        }
        return(res)
    }
    # Sovle ode
    out <- ode(y = inits, times = times, func = dydt, parms = parms, method = "adams")
    f <- force(out[, "time"], out[, "y"], out[, "dy"], parms)
    return(list(
        time = times,
        parms = parms,
        y = out[, "y"],
        dy = out[, "dy"],
        ddy = f / parms$m
    ))
}


times <- seq(0, 3, by = 0.002)
parms <- list(g = 9.81, m = 100.0, k = 20000, c = 1500)
init <- c(y = 1.0, dy = 5) # inital conditions start from a height of 1m with upward velcoity of 5m/s

out <- gen_impact_data(times, parms, init)


saveRDS(out, file = "./data/SimulatedImpact.rds")

par(mfrow = c(1, 3))
plot(out[, "time"], out[, "y"], type = "l", ylab = "Height (m)")
abline(h = 0, col = "red")
plot(out[, "time"], out[, "dy"], type = "l", ylab = "Velocity (m/s)")
plot(out[, "time"], finite_diff(out[, "ddy"], out[, "time"], 2),
    type = "l",
    ylab = "Acceleration (m/s^2)"
)
