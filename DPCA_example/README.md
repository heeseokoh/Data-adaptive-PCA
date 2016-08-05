
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **DPCA_example** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : DPCA_example

Published in : 'Lim, Y., Oh, H. S. (2016): A data-adaptive principal component analysis: Use of
composite asymmetric Huber function. Journal of Computational and Graphical Statistics'

Description : 'Examples of using DPCA with daily maximum precipitation data in August during the
period of year 1997 to 2008.'

Keywords : principal component, pca, asymmetric, bimodal, Huber norm, weighted, percipitation

See also : DPCA_source_code

Author : OH H. S., Lim Y.

Submitted : 20160728 Burdejova Petra

Datafile : max_prec.RData

Example : 'Examples of using DPCA with daily maximum precipitation data in August during the period
of year 1997 to 2008. The data are available with loading max_prec.RData file. First from
grid_finder function we find the indices of East Asia region, perform and plot loading matrices of
first 4 PCs. We may produce some results, which is dpca_loading2.pdf. Finally, we plot 4 PCs as a
time series, which is dpca_pcs.pdf.'

```


### R Code:
```r
grid_finder = function(lon, lat, full.lon, full.lat, plot = TRUE) {
    
    library(fields)
    log.grid = seq.default(0, 360, , full.lon)
    lat.grid = seq.default(-90, 90, , full.lat)
    xy.grid  = expand.grid(log.grid, lat.grid)
    
    
    start = which.min(abs(xy.grid[, 1] - lon[1]) + abs(xy.grid[, 2] - lat[1]))
    mid1  = which.min(abs(xy.grid[, 1] - lon[2]) + abs(xy.grid[, 2] - lat[1]))
    mid2  = which.min(abs(xy.grid[, 1] - lon[1]) + abs(xy.grid[, 2] - lat[2]))
    final = which.min(abs(xy.grid[, 1] - lon[2]) + abs(xy.grid[, 2] - lat[2]))
    
    
    st = ASIA = c(start:mid1)
    for (m in 1:((mid2 - start)/full.lon)) {
        ASIA = cbind(ASIA, st + full.lon * m)
    }
    
    if (plot == TRUE) {
        A <- matrix(5, nrow = 1, ncol = full.lon * full.lat)
        A[1, ASIA] = 1
        image.fit = as.image(A, x = xy.grid, grid = list(x = log.grid, y = lat.grid))
        image.plot(image.fit, main = "Region", zlim = c(1:2))
        map("world2", ylim = c(-90, 90), xlim = c(0, 360), add = TRUE)
    }
    
    print(paste("number of lon : ", length(which(xy.grid[ASIA, ] == xy.grid[ASIA[1], 
        1]))))
    print(paste("number of lat : ", length(which(xy.grid[ASIA, ] == xy.grid[ASIA, 2][1]))))
    return(as.vector(ASIA))
    
}  #end of function grin_finder

ASIA <- grid_finder(c(120, 140), c(30, 50), 360, 180, plot = FALSE)

y = max_prec[, ASIA]

fit0 = comp_Pseudo(y, 4)

log.grid = seq(0, 360, length = 360)
lat.grid = seq(-90, 90, length = 180)
xy.grid  = expand.grid(log.grid, lat.grid)
xy.grid  = xy.grid[ASIA, ]

par(mfrow = c(1, 4))
for (i in 1:4) {
    fit = fastTps(xy.grid, fit0[[1]][i, ], lambda = 1, theta = 10)
    surface.Krig(fit, type = "C", col = rev(terrain.colors(64, alpha = 0.5)), xlab = "", 
        ylab = "", zlim = c(-0.2, 0.2))
    world(shift = T, add = T)
    
}

pc = y %*% t(fit0[[1]])
plot(pc[, 1], type = "l", col = 1, ylim = c(-320, 370))
lines(pc[, 2], type = "l", col = 2)
lines(pc[, 3], type = "l", col = 3)
lines(pc[, 4], type="l", col=4)

```
