library(ncvreg)
library(lars)
library(MASS)
library(fastICA)
library(geoR)
library(fields)
library(glmnet)


log.grid <- seq(0,360,,144)
lat.grid <- seq(-90,90,,73)
xy.grid <- expand.grid(log.grid,lat.grid)

num.full.grid =p=10512
reduced.grid=2500
yr=20
y1.hat_en<-matrix(nrow=yr, ncol=p)
result<-matrix(nrow=yr, ncol=reduced.grid)



#########################################################################
#########################################################################
#########################################################################
#############################				Functions		   		  		###############
#########################################################################
#########################################################################
#########################################################################
#########################################################################


cv.scadglm<-function(x,y,K=5,method=c('ridge','lasso','scad'),type=c('ica','pca', 'dpca')){
			
			

			
			n=nrow(y)	
			all.folds <- cv.folds(n-1, K)

			index = c(1:10)
			residmat <-matrix(0, length(index), K)
			y1.hat_en <-y
			
		for(r in index){
			 		cat(paste('###########', r, 'th: components Start'), '\n')
			for (i in seq(K)) {
	
		omit <- all.folds[[i]]
									
			alpha.t<-newA <-list()
 		
 		for (mm in 1:9){
 			
 			
 			XX<-x[[mm]][[n]][-omit,]
 
 		if(type=='ica'){
 			A<-fastICA((XX),r)
 			alpha.t[[mm]]<-(A$S)	
 			newA[[mm]]<-matrix(nrow=length(omit), ncol=r)
 			newA[[mm]]<-x[[mm]][[n]][omit,] %*%(A$K%*%A$W)  }
 			
 			
 		if(type=='pca'){
 			A<-svd(cov(XX))
 			pc=XX %*%  A$v[,1:r] 
 			alpha.t[[mm]]<-pc
 			newA[[mm]]<-matrix(nrow=length(omit), ncol=r)
 			newA[[mm]]<-  x[[mm]][[n]][omit,] %*% A$v[,1:r]   }
 			
 		if(type=='dpca'){
 			A<-comp_Pseudo(XX, r)[[1]]
 			pc=XX %*% t(A)
			alpha.t[[mm]]<-pc
 			newA[[mm]]<-matrix(nrow=length(omit), ncol=r)
 			newA[[mm]]<-  x[[mm]][[n]][omit,] %*% t(A)  }
 		
 		
 			
 		}
 		
 
 
 			A <-cbind((alpha.t[[1]]))
 			for(m in 2:9){ A = cbind(A , cbind((alpha.t[[m]])))}
 			XX= A 


 	    
 	    for(grid in 1:ncol(y)){
 	    	

if(method=='lasso'){	
			lasso_pred<-lars(XX,y[-omit, grid],type='lasso')
			cv_result<-cv.lars(XX, y[-omit, grid])	
			opt_frac<-cv_result$index[c(which(cv_result$cv==min(cv_result$cv)))]
			ridge_pred <- predict(lasso_pred, s=opt_frac, type="coef", mode="fraction")$coefficients}

if(method=='scad'){				
	scad_pred<-ncvreg(XX,y[-omit, grid], family="gaussian", penalty="SCAD")
	cv_result<-cv.ncvreg(XX, y[-omit, grid], family="gaussian", penalty="SCAD", nfolds=5)
	ridge_pred<-scad_pred$beta[, cv_result$min][-1]}
if(method=='ridge'){		
    	 cv_fit=cv.glmnet(XX, y[-omit, grid],lambda=seq(0.01,3,by=0.01), nfolds=5)
 				ridge_pred <- glmnet(XX, y[-omit, grid], family='gaussian',alpha=0,lambda= cv_fit $lambda.min)$beta
 				ridge_pred<-as.vector(ridge_pred) }
    
    
    	#fit <-lm.ridge(y[-omit, grid]~ XX-1,  lambda = seq(0.01,3,by=0.01))
  	   #opt_lambda_ridge <-fit$lambda[ which(fit $GCV==min(fit $GCV))[[1]]] 		
  	   #ridge_pred <-lm.ridge(y[-omit, grid]~ XX-1, lambda = opt_lambda_ridge)$coef}
 		
 				
 			A <-cbind(newA[[1]])
 			for(m in 2:9){ A = cbind(A , cbind(newA[[m]]))}
 			newX= A 
 				
 		y1.hat_en[omit, grid] = newX %*%(ridge_pred)
 		}
 		residmat[r, i] <- mean(apply((y1.hat_en[omit,]- y[omit,])^2, 2, mean))}}
 		
 		
 		cv <- apply(residmat, 1, mean)
 		cv.error <- sqrt(apply(residmat, 1, var)/K)
		object <- list(index = index, cv = cv, cv.error = cv.error)
		invisible(object)
}


#type='dpca'
#reg_type='ridge'
#reg_type='lasso'
#reg_type='scad'



ICARR<-function(mean_obs_prec,model,type,reg_type,ASIA){


	yr=20
	p= 10512
	reduced.grid =2500
	y1.hat_en<-matrix(nrow=yr, ncol=p)
	result<-matrix(nrow=yr, ncol= reduced.grid)
	
	for(n in 1:20){ #Change to 20
 	
 	
 		cat(paste('###########', n, 'th: Iteration Start'), '\n')
 		
 		Y<-(mean_obs_prec[[n]][-n,ASIA])
 		X <-list()
 		for(m in 1:9){ 	
 			X[[m]]<-list()
 			for(year in 1:(yr-1)){
 				yy=c(1:yr)[-n]
 			X[[m]][[year]]<-(model[[m]][[yy[year]]][-n,ASIA])}}

 		cv.fit=cv.scadglm(X,Y,K=5, reg_type, type)
 		
 		print("cv.scadglm is done")
 	   
		index = c(1:10)
 	   r = index[which.min(cv.fit$cv)]
 		print(r)
 		
 	
	DD<-cbind(xy.grid, c(1: num.full.grid))
	y1 <- sample.geodata(as.geodata(DD), reduced.grid)
	sp_sample=  as.numeric(rownames(y1$coords))




	sam<-as.numeric(row.names(y1$coords))
	
	sam.obs<-list()
	for(year in 1:yr){
		sam.obs[[year]]<-matrix(nrow=yr, ncol=length(sam))
		for(yyear in 1:yr){
		sam.obs[[year]][yyear,] = mean_obs_prec[[year]][yyear,sam]}
	}
	
	sam.model<-list()
for(K in 1:9){
	sam.model[[K]]<-list()
	for(year in 1:yr){
		sam.model[[K]][[year]]<-matrix(nrow=yr, ncol=length(sam))
		for(yyear in 1:yr){
		sam.model[[K]][[year]][yyear,] = model[[K]][[year]][yyear,sam]}
	}}
	

 		
 		Y<-(sam.obs[[n]][-n,])
 
 		alpha.t<-newA <-list()

 		
 		for (K in 1:9){
 			
 			
 		 	XX=sam.model[[K]][[n]][-n,]	
 		 	
 		 	if(type=='ica'){
 				A<-fastICA((XX),r)
 				alpha.t[[K]]<-(A$S)	
 				newA[[K]]<-matrix(nrow=1, ncol=r)
 				newA[[K]]<-(sam.model[[K]][[n]][n,]%*% A$K%*% A$W )  }
 			
 			if(type=='pca'){
 				A<-svd(cov(XX))
 				pc=XX %*%  A$v[,1:r] 
 				alpha.t[[K]]<-pc
 				newA[[K]]<-matrix(nrow=1, ncol=r)
 		   	 newA[[K]]<-t(t( A$v[,1:r])    %*%  sam.model[[K]][[n]][n,] )  } 
 		
 		
 			if(type=='dpca'){
 				A<-comp_Pseudo(XX, r)[[1]]
 				pc=XX %*% t(A)
 				alpha.t[[K]]<-pc
 				newA[[K]]<-matrix(nrow=1, ncol=r)
 		   	 newA[[K]]<-t(A    %*%  sam.model[[K]][[n]][n,] )  } 
 		}
 		
 
 			A <-cbind((alpha.t[[1]]))
 			for(m in 2:9){ A = cbind(A , cbind((alpha.t[[m]])))}
 			XX= A 
 				
 	
 	    
 	    for(i in 1:length(sam)){
 	    	
 			if(reg_type=='ridge'){
   				 cv_fit=cv.glmnet(XX, Y[,i],lambda=seq(0.01,3,by=0.01))
				ridge_pred <- glmnet(XX, Y[,i], family='gaussian',alpha=0,lambda= cv_fit $lambda.min)$beta
				ridge_pred<-as.vector(ridge_pred)}
   				 #fit <-lm.ridge(Y[,i]~ XX-1,  lambda = seq(0.01,3,by=0.01))
  	 			 #opt_lambda_ridge <-fit$lambda[ which(fit $GCV==min(fit $GCV))[[1]]] 		
  	  			 #ridge_pred <-lm.ridge(Y[,i]~ XX-1, lambda = opt_lambda_ridge)$coef}
 
  			if(reg_type=='scad'){
				scad_pred<-ncvreg(XX, Y[,i], family="gaussian", penalty="SCAD")
				cv_result<-cv.ncvreg(XX, Y[,i], family="gaussian", penalty="SCAD", nfolds=5)
				ridge_pred<-scad_pred$beta[, cv_result$min][-1]}

	 		if(reg_type=='lasso'){
	 			lasso_pred<-lars(XX,Y[,i],type='lasso')
				cv_result<-cv.lars(XX, Y[,i])	
				opt_frac<-cv_result$index[c(which(cv_result$cv==min(cv_result$cv)))]
				ridge_pred <- predict(lasso_pred, s=opt_frac, type="coef", mode="fraction")$coefficients}

		
 				
 			A <-cbind(newA[[1]])
 			for(m in 2:9){ A = cbind(A , cbind(newA[[m]]))}
 			newX= A 
 				
 				result[n,i] = newX %*%(ridge_pred)
 		}
 

	pred.grid<-y1$coords
	krig.fit<-Krig(pred.grid,(result[n,]), theta=10)
	y1.hat_en[n,]<- predict(krig.fit,xy.grid)
}
	return(y1.hat_en )
}







ica_ridge=ICARR(mean_obs_prec,model,type='ica',reg_type='ridge',ASIA)
ica_scad=ICARR(mean_obs_prec,model,type='ica',reg_type='scad',ASIA)
pca_ridge=ICARR(mean_obs_prec,model,type='pca',reg_type='ridge',ASIA)
pca_scad=ICARR(mean_obs_prec,model,type='pca',reg_type='scad',ASIA)
dpca_ridge=ICARR(mean_obs_prec,model,type='dpca',reg_type='ridge',ASIA)
dpca_scad=ICARR(mean_obs_prec,model,type='dpca',reg_type='scad',ASIA)


w_asia_mean<-function(X){

	up.asia<-0
	down.asia<-0
	p=10512
	
	#ASIA<-c(7814:7834,7670:7690,7526:7545,7382:7400,7238:7255,7094:7110,6950:6964,6806:6816)
 ASIA<-sort(c(6806:6826, 6950:6970,7094:7114,7238:7258,7382:7402,7526:7546,7670:7690,7814:7834,7958: 7978))
 
	temp<-vector(length=p)
	temp[ASIA]<-X
	
	for(i in ASIA){		
		up.asia<-up.asia+(temp[i]*cos((2.5*i%/%144-90)*pi/180))
		down.asia<-down.asia+cos((2.5*i%/%144-90)*pi/180)
		
		}
	
	avg.ASIA<- up.asia/down.asia
	
	return(avg.ASIA)

}


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#######                            										ASIA			   		  		###############
#########################################################################
#########################################################################
#########################################################################
#########################################################################

 
ASIA<-sort(c(6806:6826, 6950:6970,7094:7114,7238:7258,7382:7402,7526:7546,7670:7690,7814:7834,7958: 7978))

avg.ASIA.obs <-avg.ASIA.ica.ridge <-avg.ASIA.ica.lasso<-avg.ASIA.ica.scad <-vector()
 for(year in 1:20){	avg.ASIA.ica.ridge[year]<-w_asia_mean(ica_ridge[year,ASIA]) }
  for(year in 1:20){avg.ASIA.ica.lasso[year]<-w_asia_mean(ica_lasso[year,ASIA]) }
 for(year in 1:20){	avg.ASIA.ica.scad[year]<-w_asia_mean(ica_scad[year,ASIA]) }
 for(year in 1:20){	avg.ASIA.obs[year]<-w_asia_mean(mean_obs_prec[[year]][year,ASIA]) }



avg.ASIA.pca.ridge <-avg.ASIA.pca.lasso<-avg.ASIA.pca.scad <-vector()
 for(year in 1:20){	avg.ASIA.pca.ridge[year]<-w_asia_mean(pca_ridge[year, ASIA]) }
  for(year in 1:20){avg.ASIA.pca.lasso[year]<-w_asia_mean(pca_lasso[year, ASIA]) }
 for(year in 1:20){	avg.ASIA.pca.scad[year]<-w_asia_mean(pca_scad[year, ASIA]) }


avg.ASIA.model1<-avg.ASIA.model2<-avg.ASIA.model3<-avg.ASIA.model4<-avg.ASIA.model5<-avg.ASIA.model6<-avg.ASIA.model7<-avg.ASIA.model8<-avg.ASIA.model9<-vector(length=20)
for(year in 1:20){
 
	avg.ASIA.model1[year]<-w_asia_mean(model[[1]][[year]][year,][ASIA])
 	avg.ASIA.model2[year]<-w_asia_mean(model[[2]][[year]][year,][ASIA])
 	avg.ASIA.model3[year]<-w_asia_mean(model[[3]][[year]][year,][ASIA])
 	avg.ASIA.model4[year]<-w_asia_mean(model[[4]][[year]][year,][ASIA]) 
 	avg.ASIA.model5[year]<-w_asia_mean(model[[5]][[year]][year,][ASIA])
 	avg.ASIA.model6[year]<-w_asia_mean(model[[6]][[year]][year,][ASIA])
 	avg.ASIA.model7[year]<-w_asia_mean(model[[7]][[year]][year,][ASIA])
 	avg.ASIA.model8[year]<-w_asia_mean(model[[8]][[year]][year,][ASIA]) 
 	avg.ASIA.model9[year]<-w_asia_mean(model[[9]][[year]][year,][ASIA])
 }
#load('/Volumes/Macintosh HD2/my code/multimodel/all_parameter/MMA_temp.RData')		
load('~/Desktop/Papers/climate_code/multimodel/all_parameter/MMA.RData')	
	
avg.ASIA.MMA <-vector(length=20)

for(year in 1:20){

	avg.ASIA.MMA[year]<-w_asia_mean(MMA[[year]][year,ASIA]) 
}

  avg.ASIA.dpca.ridge <-avg.ASIA.dpca.scad <-vector()
 for(year in 1:20){	avg.ASIA.dpca.ridge[year]<-w_asia_mean(dpca_ridge[year,ASIA]) }
 for(year in 1:20){	avg.ASIA.dpca.scad[year]<-w_asia_mean(dpca_scad[year,ASIA]) }
  
par(cex.main=1.2, cex.axis=1.2,omi=c(0,0.5,0.5,0),lwd=2.1,cex.lab=1.2)
A<-c(cor(avg.ASIA.MMA,avg.ASIA.obs),cor(avg.ASIA.ica.ridge,avg.ASIA.obs),cor(avg.ASIA.ica.lasso,avg.ASIA.obs),cor(avg.ASIA.ica.scad,avg.ASIA.obs),cor(avg.ASIA.pca.ridge,avg.ASIA.obs),cor(avg.ASIA.pca.lasso,avg.ASIA.obs),cor(avg.ASIA.pca.scad,avg.ASIA.obs), cor(avg.ASIA.dpca.ridge, avg.ASIA.obs), cor(avg.ASIA.dpca.scad, avg.ASIA.obs))
A<-round(A,3)
plot(1983:2002, avg.ASIA.obs, type='l',xlab='time(1983~2002)', ylab='Stand.Anom.',ylim=c(-.8,1))
#points(1983:2002, avg.ASIA.model2,type='l',col=2)
points(1983:2002, avg.ASIA.MMA,type='l',col=2)
points(1983:2002, avg.ASIA.pca.ridge,type='l',col=3)
#points(1983:2002, avg.ASIA.pca.lasso,type='l',col=4)
points(1983:2002, avg.ASIA.pca.scad,type='l',col=4)
points(1983:2002, avg.ASIA.ica.scad,type='l',col=5)
points(1983:2002, avg.ASIA.ica.ridge,type='l',col=10)
points(1983:2002, avg.ASIA.dpca.ridge, type='l', col=6)
points(1983:2002, avg.ASIA.dpca.scad,type='l',col='orange')


legend('topleft',c("obs" ,paste("MMA : ",A[1],sep=''),paste("PCA/ridge : ",A[5],sep=''),paste("PCA/SCAD : ",A[7],sep=''),paste("ICA/ridge : ",A[2],sep=''),paste("ICA/SCAD : ",A[4],sep=''),paste("DPCA/ridge : ",A[8],sep=''), paste("DPCA/SCAD : ",A[9],sep='')),lty=1,col=c(1,2,3,4,7,5,6,'orange'),cex=1.0)
title(main='Temporal correlation in East Asia', line=1)

#######ridge
plot(1983:2002, avg.ASIA.obs, type='l',xlab='time(1983~2002)', ylab='Stand.Anom.',ylim=c(-.8,1), cex=1)
#points(1983:2002, avg.ASIA.model2,type='l',col=2)
points(1983:2002, avg.ASIA.MMA,type='l',col=2)
points(1983:2002, avg.ASIA.pca.ridge,type='l',col=3)
#points(1983:2002, avg.ASIA.pca.lasso,type='l',col=4)
#points(1983:2002, avg.ASIA.ica.ridge,type='l',col=7)
points(1983:2002, avg.ASIA.dpca.ridge, type='l', col=6)


legend('topleft',c("obs" ,paste("MMA : ",A[1],sep=''),paste("PCA/ridge : ",A[5],sep=''),paste("DPCA/ridge : ",A[8],sep='')) ,lty=1,col=c(1,2,3,6),cex=1.0)
title(main='Temporal correlation in East Asia', line=1)

#######SCAD
plot(1983:2002, avg.ASIA.obs, type='l',xlab='time(1983~2002)', ylab='Stand.Anom.',ylim=c(-.8,1))
#points(1983:2002, avg.ASIA.model2,type='l',col=2)
points(1983:2002, avg.ASIA.MMA,type='l',col=2)
points(1983:2002, avg.ASIA.pca.scad,type='l',col=3)
#points(1983:2002, avg.ASIA.pca.lasso,type='l',col=4)
#points(1983:2002, avg.ASIA.ica.scad,type='l',col=7)
points(1983:2002, avg.ASIA.dpca.scad, type='l', col=6)


legend('topleft',c("obs" ,paste("MMA : ",A[1],sep=''),paste("PCA/SCAD : ",A[7],sep=''),paste("DPCA/SCAD : ",A[9],sep='')) ,lty=1,col=c(1,2,3,6),cex=1.0)
title(main='Temporal correlation in East Asia', line=1)



#########################################################################
#######                            					  RMSE	 (Table 2)	- ASIA		     		###############
#########################################################################

elnino_yr<-c(1987,1991,1992,1994,1997,2002)
elyr<-c(9,13,14,16,19,24)-4
lanina_yr<-c(1985,1988,1998,1999)
layr<-c(7,10,20,21)-4
normal_yr<- c(1983,1984,1986,1989:1990,1993,1995:1996,2000:2001)
nryr<-c(1:20)[-c(layr,elyr)]
monsoon <- c(6,7,8,12,16,19,21,23,24)-4
monsoon_below <-c(1:20)[-monsoon]
allyr<-c(1:20)


TIME<-list(allyr, elyr, layr, nryr, monsoon, monsoon_below)
names(TIME)<-c("all year", "elnino", "lanina", "normal", "postEASMI", "negEASMI")

for(iter in 1:6){
	time=TIME[[iter]]

MSEE<-list()
for(m in 1:18){MSEE[[m]]<-vector()}

for(i in time){
	MSEE[[1]][i]<-sqrt(mean((avg.ASIA.model1[i]-avg.ASIA.obs[i])^2))
	MSEE[[2]][i]<-sqrt(mean((avg.ASIA.model2[i]-avg.ASIA.obs[i])^2))
	MSEE[[3]][i]<-sqrt(mean((avg.ASIA.model3[i]-avg.ASIA.obs[i])^2))
	MSEE[[4]][i]<-sqrt(mean((avg.ASIA.model4[i]-avg.ASIA.obs[i])^2))
	MSEE[[5]][i]<-sqrt(mean((avg.ASIA.model5[i]-avg.ASIA.obs[i])^2))
	MSEE[[6]][i]<-sqrt(mean((avg.ASIA.model6[i]-avg.ASIA.obs[i])^2))
	MSEE[[7]][i]<-sqrt(mean((avg.ASIA.model7[i]-avg.ASIA.obs[i])^2))
	MSEE[[8]][i]<-sqrt(mean((avg.ASIA.model8[i]-avg.ASIA.obs[i])^2))
	MSEE[[9]][i]<-sqrt(mean((avg.ASIA.model9[i]-avg.ASIA.obs[i])^2))
	MSEE[[10]][i]<-sqrt(mean((avg.ASIA.MMA[i]-avg.ASIA.obs[i])^2))

	MSEE[[11]][i]<-sqrt(mean((avg.ASIA.pca.ridge[i]-avg.ASIA.obs[i])^2))
	MSEE[[12]][i]<-sqrt(mean((avg.ASIA.pca.lasso[i]-avg.ASIA.obs[i])^2))
	MSEE[[13]][i]<-sqrt(mean((avg.ASIA.pca.scad[i]-avg.ASIA.obs[i])^2))
		
	MSEE[[14]][i]<-sqrt(mean((avg.ASIA.ica.ridge[i]-avg.ASIA.obs[i])^2))
	MSEE[[15]][i]<-sqrt(mean((avg.ASIA.ica.lasso[i]-avg.ASIA.obs[i])^2))
	MSEE[[16]][i]<-sqrt(mean((avg.ASIA.ica.scad[i]-avg.ASIA.obs[i])^2))
	
	MSEE[[17]][i]<-sqrt(mean((avg.ASIA.dpca.scad[i]-avg.ASIA.obs[i])^2))
	
	}
rm<-rm_sd<-vector()


for(m in 1:18){
rm[m]<-round(mean(MSEE[[m]][time]),3)
rm_sd[m]<-round(sd(MSEE[[m]][time]),3)}

print(names(TIME)[iter])
print(c("MMA", "PCAridge", "PCAscad", "ICAridge", "ICAscad", "DPCAridge", "DPCAscad"))
 print(paste( rm[10],'(',rm_sd[10],') &',rm[11],'(',rm_sd[11],') &' ,rm[13],'(',rm_sd[13],') &' ,rm[14],'(',rm_sd[14],') &'  ,rm[16],'(',rm_sd[16],')&'  ,rm[17],'(',rm_sd[17],')&', rm[18],'(',rm_sd[18],')',sep=''))

}

#########################################################################
#########################################################################
#######                            			  Pattern	 (Fig  5)				     				###############
#######                  Figure data are saved for the other graphical tools				###############
#########################################################################
#########################################################################

library(fields)
library(spam)
library(QuantPsyc)



## ASIA ##

ASIA<-sort(c(6806:6826, 6950:6970,7094:7114,7238:7258,7382:7402,7526:7546,7670:7690,7814:7834,7958: 7978))
xy.grid <- xy.grid[ASIA, ] 
log.grid <- seq(93.14685,143.49650,by=2.5)
lat.grid <- seq(27.5,47.5,by=2.5)





##  Input year
timelist<-list(elyr, layr, nryr, monsoon, monsoon_below, allyr)
for(iter in 1:6){
	dev.new()
time= timelist[[iter]]

mean_hat_y_ica.lasso<-apply(ica_lasso[time,ASIA],2,mean)
mean_hat_y_ica.scad<-apply(ica_scad[time,ASIA],2,mean)
mean_hat_y_ica.ridge<-apply(ica_ridge[time,ASIA],2,mean)

mean_hat_y_pca.lasso<-apply(pca_lasso[time,ASIA],2,mean)
mean_hat_y_pca.scad<-apply(pca_scad[time,ASIA],2,mean)
mean_hat_y_pca.ridge<-apply(pca_ridge[time,ASIA],2,mean)

#mean_hat_y_dpca.lasso<-apply(dpca_lasso[time,ASIA],2,mean)
mean_hat_y_dpca.scad<-apply(dpca_scad[time,ASIA],2,mean)
mean_hat_y_dpca.ridge<-apply(dpca_ridge[time,ASIA],2,mean)
mean_model <-list()

MO<-list()
for(m in 1:9){
	MO[[m]]<-matrix(nrow=20,ncol=length(ASIA))
	for(i in time){
	MO[[m]][i,]<-model[[m]][[i]][i,ASIA]}
	mean_model[[m]]<-apply(MO[[m]][time,],2,mean)
}
	MO_MME <-matrix(nrow=20,ncol=length(ASIA))
	for(i in time){
	MO_MME[i,]<-MMA[[i]][i,ASIA]}

mean_mma_cate <-apply(MO_MME[time,],2,mean)
asia_Y<-list()
for(k in 1:20){
asia_Y[[k]]<-mean_obs_prec[[k]][,ASIA]}

A<-matrix(nrow=20, ncol=length(ASIA))
	for(j in 1:20){
		A[j,]<-asia_Y[[j]][j,]
}
mean_obs_asia<-apply(A[time,],2,mean)


 

  Result<-matrix(nrow=189, ncol=8)
 par(mfrow=c(3,2), mar=c(2,2,1,1), cex.main=1, cex.axis=1,omi=c(0,0,0.5,0))
 zlim=c(-3,3)
  AA<-Make.Z(mean_obs_asia)
 AA[which(Make.Z(mean_obs_asia)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_obs_asia)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_obs_asia)> zlim[2])  )))
 fit <- fastTps(xy.grid, AA, lambda=1, theta=10)
 surface(fit, type="I",main='Observation', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,1]=fit$fitted.values


AA<-Make.Z(mean_mma_cate)
 AA[which(Make.Z(mean_mma_cate)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_mma_cate)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_mma_cate)> zlim[2])  )))
 fit <- fastTps(xy.grid,AA, lambda=1, theta=10)
 surface(fit, type="I",main='MMA', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,2]=fit$fitted.values
 
   AA<-Make.Z(mean_hat_y_pca.ridge)
 AA[which(Make.Z(mean_hat_y_pca.ridge)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_hat_y_pca.ridge)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_pca.ridge)> zlim[2])  )))
 fit <- fastTps(xy.grid, AA, lambda=1, theta=10)
 surface(fit, type="I",main='PCA/ridge', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,3]=fit$fitted.values
 
   AA<-Make.Z(mean_hat_y_pca.scad)
 AA[which(Make.Z(mean_hat_y_pca.scad)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_hat_y_pca.scad)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_pca.scad)> zlim[2])  )))
 fit <- fastTps(xy.grid,AA, lambda=1, theta=10)
 surface(fit, type="I",main='PCA/SCAD', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,4]=fit$fitted.values
 
   # AA<-Make.Z(mean_hat_y_ica.ridge)
 # AA[which(Make.Z(mean_hat_y_ica.ridge)< zlim[1])]=zlim[1]
 # AA[which(Make.Z(mean_hat_y_ica.ridge)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_ica.ridge)> zlim[2])  )))
 # fit <- fastTps(xy.grid, AA, lambda=1, theta=10)
 # surface(fit, type="I",main='ICA/ridge', zlim= zlim,col= terrain.colors(50))
 # world(shift = T,add=T); Result[,5]=fit$fitted.values
 
   # AA<-Make.Z(mean_hat_y_ica.scad)
 # AA[which(Make.Z(mean_hat_y_ica.scad)< zlim[1])]=zlim[1]
 # AA[which(Make.Z(mean_hat_y_ica.scad)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_ica.scad)> zlim[2])  )))
 # fit <- fastTps(xy.grid,AA, lambda=1, theta=10)
 # surface(fit, type="I",main='ICA/SCAD', zlim= zlim,col= terrain.colors(50))
 # world(shift = T,add=T); Result[,6]=fit$fitted.values


  AA<-Make.Z(mean_hat_y_dpca.ridge)
 AA[which(Make.Z(mean_hat_y_dpca.ridge)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_hat_y_dpca.ridge)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_dpca.ridge)> zlim[2])  )))
 fit <- fastTps(xy.grid, AA, lambda=1, theta=10)
 surface(fit, type="I",main='DPCA/ridge', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,7]=fit$fitted.values
 
   AA<-Make.Z(mean_hat_y_dpca.scad)
 AA[which(Make.Z(mean_hat_y_dpca.scad)< zlim[1])]=zlim[1]
 AA[which(Make.Z(mean_hat_y_dpca.scad)> zlim[2])]=zlim[2]-abs(rnorm(length(which(Make.Z(mean_hat_y_dpca.scad)> zlim[2])  )))
 fit <- fastTps(xy.grid,AA, lambda=1, theta=10)
 surface(fit, type="I",main='DPCA/SCAD', zlim= zlim,col= terrain.colors(50))
 world(shift = T,add=T); Result[,8]=fit$fitted.values

}
 write(t(Result),file='~/Desktop/b.txt', ncolumns=8)
 
 

