# Adapted from recurr function in tseriesChaos
CLrecurr <- function (series, m, d, start.time = stats::start(series), end.time = stats::end(series),
                       maincolor="forestgreen", invcol=FALSE, pos.name, ...)
{
    xyz <- tseriesChaos::embedd(stats::window(series, start = start.time, end = end.time),
        m = m, d = d)
    D <- stats::dist(xyz)
    D <- as.matrix(D)/max(D)
    grid.x <- stats::time(series)[1:nrow(D)]
    # filled.contour(grid.x, grid.x, D, color.palette = function(n) gray(0:(n -
    #     1)/(n - 1)), xlab = "time", ylab = "time", main = "Recurrence plot",
    graphics::filled.contour(grid.x,
                   grid.x,
                   D,
                   color.palette = function(n) {
                       if (invcol) {grDevices::colorRampPalette(c("white",maincolor),space = "rgb")(n)
                       } else { grDevices::colorRampPalette(c(maincolor, "white"),space = "rgb")(n)
                       }
                     },
                   #color.palette = function(n) cm.colors(n),
                   #color.palette = function(n) gray((n - 1):0 /(n - 1)),
                   #color.palette = function(n) gray(log(1/(1+ seq(n - 1)/(n - 1)):0),
                   xlab = "nt position / AUG",
                   ylab = "nt position / AUG",
                   main = pos.name, #"Recurrence plot",
        ...)
}


