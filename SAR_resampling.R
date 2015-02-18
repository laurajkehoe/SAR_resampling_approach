### Packages needed for parallellization
library(doParallel)
library(foreach)

### Set number of iterations
n <- 50

# Set path to data folder
dataFolder

### Load functions
source(dataFolder%+%"/data/lib/resampling_approach/SAR_function.R")
source(dataFolder%+%"/data/lib/resampling_approach/SAR_function_global_power_model.R")
source(dataFolder%+%"/data/lib/resampling_approach/SAR_function_biome_power_model.R")
source(dataFolder%+%"/data/lib/resampling_approach/SAR_function_model_summary.R")

### Create unique database names
databases <- paste(dataFolder%+%"/data/lib/resampling_approach/SAR_databases/SAR_",1:n,".csv",sep="")

###############################################################################################
### Create SAR databases

cl <- makeCluster(30)
registerDoParallel(cl)

foreach(SAR_database_output=databases, 
        .packages=c("ggplot2","reshape","raster","foreign",
                    "sp","spdep","maptools","plyr","ncf",
                    "stringr","caret","data.table","classInt")) %dopar% {
  createSAR(SAR_database_output, dataFolder=dataFolder)
} 

stopCluster(cl)

###############################################################################################
### Models

### Calculate models and save outputs (parameterization/crossvalidation) with uniqe extension

for(SAR_database_output in databases){
  # Global power
#   SARmodelGlobal(SAR_database_output, 
#                  modelrun_name=paste("global_power_model",substr(basename(SAR_database_output),1,nchar(basename(SAR_database_output))-4),sep="_"), 
#                  kOptStrategy="AIC",
#                  dataFolder=dataFolder)
  # Biome power
  SARmodelBiome(SAR_database_output, 
                modelrun_name=paste("biome_power_model",substr(basename(SAR_database_output),1,nchar(basename(SAR_database_output))-4),sep="_"), 
                kOptStrategy="AIC",
                dataFolder=dataFolder)
}

### Calculate biome models and save outputs (parameterization/crossvalidation) with uniqe extension

for(SAR_database_output in databases){
    
}

###############################################################################################
### Get summary information from models

results_global <- SARmodelSummary(databases, modelname="global_power_model", dataFolder)
