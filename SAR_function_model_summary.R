###############################################################################################
###
### Title: Function for extracting resampled model information
###
### Authors: Cornelius Senf  
###
### Date: Feb 2015
###
###############################################################################################

SARmodelSummary <- function(databases, modelname, dataFolder="L:/global_lubio/chapter_2/CHAPTER_II"){
  
  "%+%" <- function(x,y)paste(x,y,sep="")
  
  ### Get cross-validation results
  
  f <- dataFolder%+%"/data/lib/resampling_approach/models/"
  
  rsquare <- c()
  rmse <- c()
  
  for(i in databases){
    modelrun_name <- paste(modelname,substr(basename(i),1,nchar(basename(i))-4),sep="_")
    x <- read.csv(f%+%modelrun_name%+%"_crossvalidation.csv")
    rsquare <- c(rsquare,as.matrix(x[,2:11])[1,])
    rmse <- c(rmse,as.matrix(x[,2:11])[2,])
  }
  
  crossvalidation <- data.frame(MEAN=c(mean(rsquare), mean(rmse)), 
                                SD=c(sd(rsquare), sd(rmse)),
                                SE=c(sd(rsquare)/sqrt(length(rsquare)), mean(rmse)/sqrt(length(rmse))))
  rownames(crossvalidation) <- c("Rsquare","RMSE")
  
  coefficients <- matrix(NA,ncol=20,nrow=length(databases)) #ncol must be larger than the maximum numbers of predictors in a model
  a <- 0
  for(i in databases){
    a <- a+1
    modelrun_name <- paste(modelname,substr(basename(i),1,nchar(basename(i))-4),sep="_")
    load(f%+%modelrun_name%+%"_final_model.RData")
    coefficients[a,1] <- fit.final$fit$coefficients[1]
    coefficients[a,2:(1+length(fit.final$fit$coefficients[-1]))] <- fit.final$fit$coefficients[-1]
  }
  
  coefficients <- as.data.frame(coefficients)
  coefficients <- coefficients[,colSums(is.na(coefficients)) != nrow(coefficients)]
  colnames(coefficients) <- names(fit.final$fit$coefficients)
  
  return(list(crossvalidation,coefficients))
  
}




