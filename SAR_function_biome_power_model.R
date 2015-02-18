###############################################################################################
###
### Title: Function for global SAR Model (power model)
###
### Authors: Cornelius Senf  
###
### Date: Nov 2014
###
###############################################################################################

SARmodelBiome <- function(SAR_database_output, modelrun_name, dataFolder="L:/global_lubio/chapter_2/CHAPTER_II", kOptStrategy="AIC"){
    
  ### Load workspace
  modelingdata <- read.csv(SAR_database_output)
  
  f <- dataFolder%+%"/data/lib/resampling_approach/models/"
    
  ###############################################################################################
  ###############################################################################################
  
  ### Packages
  library(reshape)
  library(sp)
  library(spdep)
  library(maptools)
  library(plyr)
  library(ncf)
  library(stringr)
  library(caret)
  
  ### Packages needed for parallellization
  library(doParallel)
  library(foreach)
  
  ### Global functions
  "%+%" <- function(x,y)paste(x,y,sep="")
  
  ###############################################################################################
  ### Model calibration
  
  #krange <- seq(100,10000,100)
  krange <- 1:20
  
  parameterization <- function(k, modelingdata){
    coords <- as.matrix(modelingdata[,c("x","y")])
    knn <- knearneigh(coords, k=k, longlat=TRUE)
    nb <- knn2nb(knn, row.names=modelingdata$id)
    listw <- nb2listw(nb, style="W", zero.policy=TRUE)
    fit <- spautolm(log(srTotal)~log(biome_area_1+1)+log(biome_area_2+1)+log(biome_area_3+1)+
                      log(biome_area_4+1)+log(biome_area_5+1)+log(biome_area_6+1)+log(biome_area_7+1)+
                      log(biome_area_8+1)+log(biome_area_9+1)+log(biome_area_10+1)+log(biome_area_11+1)+
                      log(biome_area_12+1)+log(biome_area_13+1)+log(biome_area_14+1),data=modelingdata, listw=listw)  
    result <- c(AIC(fit), 
                summary(fit,Nagelkerke=TRUE)$NK,
                sum(abs(correlog(x=modelingdata[,c("x")], y=modelingdata[,c("y")], z=as.numeric(residuals(fit)), increment=10, resamp=1, latlon=TRUE)$correlation)))
    return(result)
  }
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  parameterization.results <- foreach(i=krange, .combine='cbind', .packages=c("spdep","ncf")) %dopar% {
    parameterization(i, modelingdata)
  } 
  
  stopCluster(cl)
  
  kMinAIC <- krange[which(parameterization.results[1,]==min(parameterization.results[1,]))]
  kMinNK <- krange[which(parameterization.results[2,]==min(parameterization.results[2,]))]
  kMinRSA <- krange[which(parameterization.results[3,]==min(parameterization.results[3,]))]
  
  # Set optimal k based on AIC or RSA
  if(kOptStrategy=="AIC"){k.opt <- kMinAIC}else if(kOptStrategy=="RSA"){k.opt <- kMinRSA}
  
  # Check if 20 iteration were enough
  #if(k.opt==max(krange)){print("Warning: k.opt occured at maximum of krange: "%+%max(krange))}
  
  parameterization.plot <- melt(parameterization.results)
  parameterization.plot$k <- as.integer(substr(as.character(parameterization.plot[,2]),8,9))
  parameterization.plot$Measure <- c("AIC","R²","RSA")
  
  # Write out table with parameterization results
  parameterization.results.df <- as.data.frame(parameterization.results)
  colnames(parameterization.results.df) <- paste("K",krange,sep="=")
  rownames(parameterization.results.df) <- c("AIC","RSquare","RSA")
  write.csv(parameterization.results.df, file=f%+%modelrun_name%+%"_parameterization.csv")
  
  ###############################################################################################
  ### Model validation (10-fold cross validation)
  
  biomeAreaCols <- which(names(modelingdata)=="biome_area_1"):which(names(modelingdata)=="biome_area_14")
  biomeCols <- which(names(modelingdata)=="biome_1"):which(names(modelingdata)=="biome_14")
  
  folds <- createFolds(modelingdata$srTotal, k=10, list=TRUE, returnTrain=FALSE)
  
  crossvalidation <- function(i, folds){
    test <- modelingdata[folds[[i]],]
    train <- modelingdata[unlist(folds[-i]),]
    
    nb <- knn2nb(knearneigh(as.matrix(train[,c("x","y")]), k=k.opt, RANN=TRUE), row.names=train$id)
    listw <- nb2listw(nb, style="W", zero.policy=TRUE)
    fit <- spautolm(log(srTotal)~log(biome_area_1+1)+log(biome_area_2+1)+log(biome_area_3+1)+
                      log(biome_area_4+1)+log(biome_area_5+1)+log(biome_area_6+1)+log(biome_area_7+1)+
                      log(biome_area_8+1)+log(biome_area_9+1)+log(biome_area_10+1)+log(biome_area_11+1)+
                      log(biome_area_12+1)+log(biome_area_13+1)+log(biome_area_14+1), data=train, listw=listw)  
    fit.summary <- summary(fit)
    prediction <- apply(log(test[,biomeAreaCols]+1), 1, function(x){sum(fit.summary$Coef[2:15]*x+fit.summary$Coef[1])})
    results <- c(cor(prediction,log(test$srTotal))^2,
                 sqrt(mean((prediction-log(test$srTotal))^2)))
    return(results)
  }
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  crossvalidation.results <- foreach(i=1:10, .combine='cbind', .packages=c("spdep","ncf")) %dopar% {
    crossvalidation(i, folds)
  } 
  
  stopCluster(cl)
  
  crossvalidation.results.df <- data.frame(crossvalidation.results, Mean=apply(crossvalidation.results,1,mean), SD=apply(crossvalidation.results,1,sd))
  rownames(crossvalidation.results.df) <- c("Rsquare","RMSE")
  write.csv(crossvalidation.results.df, f%+%modelrun_name%+%"_crossvalidation.csv")
  
  ###############################################################################################
  ### Final Model
  
  nb <- knn2nb(knearneigh(as.matrix(modelingdata[,c("x","y")]), k=k.opt, longlat=TRUE), row.names=modelingdata$id)
  #nb <- dnearneigh(as.matrix(modelingdata[,c("x","y")]), d1=0, d2=dist.opt, row.names=modelingdata$id, longlat=TRUE)
  listw <- nb2listw(nb, style="W", zero.policy=TRUE)
  fit.final <- spautolm(log(srTotal)~log(biome_area_1+1)+log(biome_area_2+1)+log(biome_area_3+1)+
                          log(biome_area_4+1)+log(biome_area_5+1)+log(biome_area_6+1)+log(biome_area_7+1)+
                          log(biome_area_8+1)+log(biome_area_9+1)+log(biome_area_10+1)+log(biome_area_11+1)+
                          log(biome_area_12+1)+log(biome_area_13+1)+log(biome_area_14+1), data=modelingdata, listw=listw)  
  summary(fit.final, Nagelkerke=TRUE)
  
  save(fit.final, file=f%+%modelrun_name%+%"_final_model.RData")
}

