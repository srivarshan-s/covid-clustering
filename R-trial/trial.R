library("fda")

set.seed(999)
n_obs <- 67
time_span <- 100
time <- sort(runif(n_obs, 0, time_span))
wiener <- cumsum(rnorm(n_obs)) / sqrt(n_obs)
y_obs <- wiener + rnorm(n_obs, 0, .05)

times_basis <- seq(0, time_span, 1)
knots <- c(seq(0, time_span, 5)) # Location of knots
n_knots <- length(knots) # Number of knots
n_order <- 4 # order of basis functions: cubic bspline: order = 3 + 1
n_basis <- length(knots) + n_order - 2
basis <- create.bspline.basis(
     c(min(times_basis), max(times_basis)),
     n_basis,
     n_order,
     knots
)

wiener_obj <- smooth.basis(argvals = time, y = y_obs, fdParobj = basis)
print(wiener_obj)

plot(time, wiener,
     type = "l", xlab = "time", ylab = "f(time)",
     main = "Comparison of fda package and naive smoothing estimates",
     col = "grey"
)
lines(wiener_obj, lwd = 1, col = "blue")
