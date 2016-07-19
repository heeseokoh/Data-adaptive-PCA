 grid_finder<-function(lon, lat, full.lon, full.lat,plot=TRUE){
	
	
	library(fields)
	log.grid <- seq.default(0,360,, full.lon)
	lat.grid <- seq.default(-90,90,, full.lat)
	xy.grid <- expand.grid(log.grid,lat.grid)
	
	
	start = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[1] )    )
	mid1 = which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[1] )    )
	mid2 = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[2] )    )
	final=  which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[2] )    )
	
	
	st=ASIA= c(start:mid1 )
	for(m in 1:( (mid2-start)/full.lon)){
		ASIA<-  cbind( ASIA, st+ full.lon*m   )}
	
	
	if( plot==TRUE){
		A<-matrix(5,nrow=1, ncol= full.lon* full.lat)
		A[1,ASIA]<- 1
		image.fit <- as.image(A,x=xy.grid, grid=list(x=log.grid, y=lat.grid))
		image.plot(image.fit, main='Region', zlim=c(1:2))
		map("world2", ylim=c(-90,90), xlim = c(0,360), add = TRUE)

	}
	
	print( paste('number of lon : '  ,length(which(xy.grid[ASIA,]==xy.grid[ASIA[1],1])) ) )
	print(paste('number of lat : '  , length(which(xy.grid[ASIA,]==xy.grid[ASIA,2][1]))))
	return (as.vector(ASIA))
	
}



ASIA<-grid_finder( c(120, 140) , c(30, 50)  ,  360,180,plot=FALSE)

y= max_prec[,ASIA]; 

fit0 = comp_Pseudo(y, 4); 

log.grid <- seq(0,360,length=360)
lat.grid <- seq(-90,90,length=180)
xy.grid <- expand.grid(log.grid,lat.grid)
xy.grid<-xy.grid[ASIA,]

par(mfrow=c(1,4))
for(i in 1:4){
	fit= fastTps(xy.grid, fit0[[1]][i, ], lambda=1, theta=10)
  surface.Krig(fit, type="C", col =rev(terrain.colors(64, alpha = .5)) ,xlab='',ylab='', zlim=c(-0.2, 0.2))
  world(shift = T,add=T)

}

pc <- y %*% t(fit0[[1]])
plot(pc[,1], type="l", col=1, ylim=c(-320, 370))
lines(pc[,2], type="l", col=2)
lines(pc[,3], type="l", col=3)
lines(pc[,4], type="l", col=4)
