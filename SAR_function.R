###############################################################################################
###
### Title: Function for creating SAR database
###
### Authors: Cornelius Senf  
###
### Date: Feb 2015
###
###############################################################################################


createSAR <- function(SAR_database_output, dataFolder="L:/global_lubio/chapter_2/CHAPTER_II"){
  
  ###############################################################################################
  ### Load packages and functions
  
  "%+%" <- function(x,y)paste(x,y,sep="")
  
  source(dataFolder%+%"/data/ode_to_r/helper_functions.R") 
  
  library(ggplot2)
  library(reshape)
  library(raster)
  library(foreign)
  library(maps)
  library(sp)
  library(spdep)
  library(maptools)
  library(plyr)
  library(ncf)
  library(stringr)
  library(caret)
  library(data.table)
  library(classInt)
  
  ### Define pathes 
  root_dir <- dataFolder%+%"/"
  
  # Input
  worldGridPath <- root_dir%+%"360x114global/360x114global.dbf"
  worldGridAttributesPath <- root_dir%+%"ode_to_r/360x114global/er_mammals_update_oct.dbf"
  
  ###############################################################################################
  ### Create samples
  
  ### Load and prepare data
  gridXY <- read.dbf(worldGridPath)
  gridXY <- gridXY[with(gridXY, order(HBWID)),]
  
  # Create grid
  gridHBWID <- raster(ncol=360, nrow=113)
  values(gridHBWID) <- gridXY$HBWID
  
  # Create mask using land/water threshold and mask grids
  gridLandWater <- raster(ncol=360, nrow=113)
  vals <- gridXY$PROP0_0062
  vals[which(gridXY$HBWID>=39355)] <- 0
  values(gridLandWater) <- vals
  landMask <- cut(gridLandWater, breaks=c(0.5,1))
  gridHBWID <- mask(gridHBWID, landMask)
  #plot(gridHBWID)
  ### Get Biome information
  gridXYattributes <- read.dbf(worldGridAttributesPath)
  biomes <- gridXYattributes[,c("HBWID","BIOME")]
  biomes$BIOME <- str_replace_all(biomes$BIOME, "[^[:alnum:]]", "") # Get rid of special characters!
  biomes$BIOME <- gsub(" ","",biomes$BIOME, fixed=TRUE)
  gridXYbiome <- merge(gridXY, biomes, by="HBWID")
  biomenames_area <- data.frame(biome=as.character(unique(biomes[,"BIOME"])),
                                code=c("nothing",paste("biome_area",1:(length(as.character(unique(biomes[,"BIOME"])))-1),sep="_")))
  
  ### Get LC information and recode land cover 
  lc <- read.dbf(root_dir%+%"/modis_landcover/documents-export-2015-01-08/grid_landcov_modis_dis.dbf")
  DT <- as.data.table(lc)
  lc_reclass <- DT[ , GRIDCODE[which.max(area_lc)], by = "HBWID"]
  gridXYlc <- merge(gridXY, lc_reclass, by="HBWID")
  
  ### Sample from HBWID grid, recording XY/coordinates, central HBWID, focal size, valid pixels within focal window, area and neighborhoud HBWID
  
  nSamples <- 100 # Number of samples
  #maxFocal <- 9 # Define maximum window size for focal filter, must be odd
  maxFocal <- sample(1:15,1)
  colnames <- c("id","x","y","HBWID","focal","ncell","area","neighbours", as.character(biomenames_area[,"biome"]), paste("LC", sort(unique(lc_reclass$V1)), sep="_"))
  
  i <- 0
  j <- 1
  samples <- as.data.frame(matrix(0, nSamples, length(colnames)))
  colnames(samples) <- colnames
  rowcolTracker <- c()
  
  while(i < nSamples){
    
    row <- sample(1:nrow(landMask),1,replace=FALSE)
    col <- sample(1:ncol(landMask),1,replace=FALSE)
    focal <- sample(seq(1,maxFocal,2),1)
    neighbours <- getValuesBlock(gridHBWID, row-floor(focal/2), focal, col-floor(focal/2), focal)
    if(length(neighbours)==1){center<-neighbours}else{center <- neighbours[ceiling(length(neighbours)/2)]}
    
    
    if(!is.na(center)){
      
      if(!center%in%rowcolTracker){
        
        rowcolTracker <- c(rowcolTracker,center)
        
        i<-i+1
        
        samples[i,"id"] <- i
        samples[i,"x"] <- gridXY[center,"X_COORD"]
        samples[i,"y"] <- gridXY[center,"Y_COORD"]
        samples[i,"HBWID"] <- center
        samples[i,"focal"] <- focal
        samples[i,"ncell"] <- length(na.omit(neighbours))
        samples[i,"area"] <- sum(gridXY[na.omit(neighbours),"AREA_KM2"])
        biome_area <- aggregate(gridXYbiome[na.omit(neighbours),"AREA_KM2"], by=list(gridXYbiome[na.omit(neighbours),"BIOME"]), FUN=sum)
        targetcol <- colnames(samples)%in%biome_area$Group.1
        samples[i,targetcol] <- biome_area$x
        
        lc_area <- aggregate(gridXYlc[na.omit(neighbours),"AREA_KM2"], by=list(gridXYlc[na.omit(neighbours),"V1"]), FUN=sum)
        targetcol <- colnames(samples)%in%paste("LC",lc_area$Group.1,sep="_")
        samples[i,targetcol] <- lc_area$x
        
        samples[i,"neighbours"] <- gsub(", ",",",toString(neighbours))
      }
    }  
  }
  
  colnames(samples)[colnames(samples)%in%biomenames_area$biome] <- as.character(biomenames_area$code)
  
  ###############################################################################################
  ### Prepare species richness data
  
  # Load data
  load(root_dir%+%"species_richness/RangeMaps_Grid_IUCNWR.RData")
  load(root_dir%+%"species_richness/TraitData_IUCNWR.Rdata")
  traits_kissling <- read.table(root_dir%+%"traits/MammalDIET_v1.0.txt", sep="\t", header=TRUE)
  
  # Change scientific names in traits kissling to lower case to match elton traits names
  traits_kissling$Scientific <- paste(tolower(as.character(traits_kissling$Genus)), traits_kissling$Species, sep=" ")
  
  # Check for species missing either in kissling or elton
  in_kissling_but_not_in_elton <- traits_kissling[which(!traits_kissling$Scientific%in%TraitData_IUCNWR$Scientific),]$Scientific
  in_elton_but_not_in_kissling <- TraitData_IUCNWR[which(!TraitData_IUCNWR$Scientific%in%traits_kissling$Scientific),]$Scientific
  
  # Combine kissling and elton by scientific names
  traits_all <- merge(TraitData_IUCNWR, traits_kissling, by="Scientific", all=TRUE)
  
  # Combine data sets
  names(RangeMaps_Grid_IUCNWR) <- c("Scientific", "HBWID")
  dataSR <- merge(RangeMaps_Grid_IUCNWR, traits_all, by="Scientific", all=TRUE)
  
  
  ###############################################################################################
  ### Calculate species richness for sample locations
  
  trophicLevel <- aggregate(dataSR$TrophicLevel, by=list(dataSR$Scientific), FUN=unique)
  
  threatLevel <- aggregate(dataSR$RedListStatusBroad, by=list(dataSR$Scientific), FUN=unique)
  
  activity <- aggregate(dataSR$ActivityOrdinal, by=list(dataSR$Scientific), FUN=unique)
  
  ### Body mass
  
  # First create the new column
  dataSR$BM_rc <- NA
  
  # Then copy the data from the existing column into the new one.
  dataSR$BM_rc <- dataSR$BodyMass.Value
  
  # Recode to lt10kg, gt10kg etc from helper functions
  dataSR$BM_rc <- sapply(dataSR$BM_rc, FUN=recodeBodyMass)
  body_mass <- aggregate(dataSR$BM_rc, by=list(dataSR$Scientific), FUN=unique)
  
  ### Range size
  
  # First create the new column
  dataSR$RS_rc <- NA
  
  # Then copy the data from the existing column into the new one.
  dataSR$RS_rc <- dataSR$RangeSize
  
  # Rrecode to RS1,2 AND 3 etc from helper functions
  dataSR$RS_rc <- sapply(dataSR$RS_rc, FUN=recodeRangeSize_T)
  range_size <- aggregate(dataSR$RS_rc, by=list(dataSR$Scientific), FUN=unique)
  
  ### Combine everything
  
  gridXYattributes <- read.dbf(worldGridAttributesPath)
  
  # Add biome information to SR data!
  biomes <- gridXYattributes[,c("HBWID","BIOME")]
  biomes$BIOME <- str_replace_all(biomes$BIOME, "[^[:alnum:]]", "") # Get rid of special characters!
  biomes$BIOME <- gsub(" ","",biomes$BIOME, fixed=TRUE)
  dataSR <- merge(dataSR, biomes, by="HBWID", all.x=TRUE) 
  
  # Add LAND COVER information to SR data!
  lc <- gridXYlc[,c("HBWID","V1")]
  colnames(lc) <- c("HBWID","landcover")
  dataSR <- merge(dataSR, lc, by="HBWID", all.x=TRUE) 
  
  # Add traits information to SR data!
  traits <- gridXYattributes[,c("HBWID","gt100kg","lt100kg","lt10kg","rs1stQ","rs2ndQ","rs3rdQ","rs4thQ","CR","DD","EN","EW","LC","NT","VU")]
  dataSR <- merge(dataSR, traits, by="HBWID", all.x=TRUE) 
    
  # Create empty vectors for SR and empty DF for SR in biomes
  srTotal <- c()
  srCarnivore <- c()
  srOmnivore <- c()
  srHerbivore <- c()
  
  srNonthreatened <- c()
  srThreatened <- c()
  
  gt100kg <- c()
  lt100kg <- c()
  lt10kg <- c()
  
  rs1 <- c()
  rs2 <- c()
  rs3 <- c()
    
  ##"ActivityOrdinal" estimates the degree of diurnality (activity during light times): I converted the thee dummy
  #variables describing the activity periods into an ordinal scale (1=nocturnal only; 2=nocturnal and crepuscular; 3=crepuscular only; 4=nocturnal, crepuscular and diurnal; 5=crepuscular and diurnal; 6=diurnal only; the combination nocturnal-diurnal doesn't exist in the Elton Trait dataset).
  
  nocturnal <- c()
  noct_crep <- c()
  crepuscular <- c()
  noct_crep_diurn <- c()
  crep_diurn <- c()
  diurnal <- c()
  
  srBiomes <- as.data.frame(matrix(0, nrow(samples), length(unique(biomes[,"BIOME"]))))
  biomenames <- data.frame(biome=as.character(unique(biomes[,"BIOME"])),
                           code=c("nothing",paste("biome",1:(length(as.character(unique(biomes[,"BIOME"])))-1),sep="_")))
  colnames(srBiomes) <- biomenames[,"biome"]
  
  srlc <- as.data.frame(matrix(0, nrow(samples), length(unique(lc[,"landcover"]))))
  lcnames <- paste("LC_sp",sort(unique(lc_reclass$V1)),sep="_")
  colnames(srlc) <- lcnames
  
  for(i in 1:nrow(samples)){
    # Get neighbour list
    neighbours <- unlist(strsplit(samples[i,"neighbours"],","))
    neighbours <- as.integer(neighbours[neighbours!="NA"])
    # Select unqie species within neighbourhood
    species <- unique(dataSR[which(dataSR$HBWID%in%neighbours),"Scientific"])
    # Calculate total number of sapecies
    srTotal <- c(srTotal, length(species))
    
    # Calculate number of species within biome
    if(length(species)!=0){
      srB <- aggregate(dataSR[which(dataSR$HBWID%in%neighbours),"Scientific"], by=list(dataSR[which(dataSR$HBWID%in%neighbours),"BIOME"]), FUN=function(x){length(unique(x))})
      targetcol <- biomenames[,"biome"]%in%as.character(as.matrix(srB["Group.1"]))
      srBiomes[i,targetcol] <- srB[match(biomenames[,"biome"][targetcol],srB$Group.1),"x"]
    }else{
      srBiomes[i,] <- 0
    }
    
    # Calculate number of species within lc class in each sample
    if(length(species)!=0){
      srLC <- aggregate(dataSR[which(dataSR$HBWID%in%neighbours),"Scientific"], by=list(dataSR[which(dataSR$HBWID%in%neighbours),"landcover"]), FUN=function(x){length(unique(x))})
      targetcol <- lcnames%in%paste("LC_sp",srLC$Group.1,sep="_")
      srlc[i,targetcol] <- srLC[match(lcnames[targetcol],paste("LC_sp",srLC$Group.1,sep="_")),"x"]
      #srlc[i,targetcol] <- srLC[match(lcnames[targetcol],srLC$Group.1),"x"]
    }else{
      srlc[i,] <- 0
    }
    
    # Calculate number of species within trophic levels
    speciesTrophicLevles <- as.character(trophicLevel[which(trophicLevel$Group.1%in%species),"x"])
    srCarnivore <- c(srCarnivore, length(species[which(speciesTrophicLevles=="Carnivore")]))
    srOmnivore <- c(srOmnivore, length(species[which(speciesTrophicLevles=="Omnivore")]))
    srHerbivore <- c(srHerbivore, length(species[which(speciesTrophicLevles=="Herbivore")]))
        
    # Calculate number of species threatened/nonthreatened
    speciesThreatenedLevles <- as.character(threatLevel[which(threatLevel$Group.1%in%species),"x"])
    srNonthreatened <- c(srNonthreatened, length(species[which(speciesThreatenedLevles=="NonThreatened")]))
    srThreatened <- c(srThreatened, length(species[which(speciesThreatenedLevles=="Threatened")]))
    
    #calculate number of small med. and large bodied species per sample 
    species_mass <- as.character(body_mass[which(body_mass$Group.1%in%species),"x"])
    lt10kg <- c(lt10kg, length(species[which(species_mass =="lt10kg")]))
    lt100kg <- c(lt100kg, length(species[which(species_mass =="lt100kg")]))
    gt100kg <- c(gt100kg, length(species[which(species_mass =="gt100kg")]))
    
    #calculate number of small med. and large ranged species per sample 
    species_range <- as.character(range_size[which(range_size$Group.1%in%species),"x"])
    rs1 <- c(rs1, length(species[which(species_range =="rs1")]))
    rs2 <- c(rs2, length(species[which(species_range =="rs2")]))
    rs3 <- c(rs3, length(species[which(species_range =="rs3")]))
    
    #calculate number of diurnal/nocturnal species
    species_activity <- as.character(activity[which(activity$Group.1%in%species),"x"])
        
    nocturnal <- c(nocturnal, length(species[which(species_activity =="1")]))
    noct_crep <- c(noct_crep, length(species[which(species_activity =="2")]))
    crepuscular <- c(crepuscular, length(species[which(species_activity =="3")]))
    noct_crep_diurn <- c(noct_crep_diurn, length(species[which(species_activity =="4")]))
    crep_diurn <- c(crep_diurn, length(species[which(species_activity =="5")]))
    diurnal <- c(diurnal, length(species[which(species_activity =="6")]))
    
  }
  
  colnames(srBiomes) <- as.character(biomenames[,"code"])
  
  samples$srTotal <- srTotal
  samples$srCarnivore <- srCarnivore
  samples$srOmnivore <- srOmnivore
  samples$srHerbivore <- srHerbivore
  samples$srNonthreatened <- srNonthreatened
  samples$srThreatened <- srThreatened
  samples$lt10kg <- lt10kg
  samples$lt100kg <- lt100kg
  samples$gt100kg <- gt100kg
  samples$rs1 <- rs1
  samples$rs2 <- rs2
  samples$rs3 <- rs3
  samples$nocturnal <- nocturnal
  samples$crepuscular <- crepuscular
  samples$diurnal <- diurnal
  samples$crep_diurn <- crep_diurn
  samples$noct_crep_diurn <- noct_crep_diurn
  samples$noct_crep <- noct_crep
  
  samples <- cbind(samples,srBiomes)
  samples <- cbind(samples,srlc)
    
  ###############################################################################################
  ### Extract predictors for samples
  
  gridXYattributes <- read.dbf(worldGridAttributesPath)
  gridXYattributes[gridXYattributes==-9999] <- NA
  gridXYattributes[which(gridXYattributes$h_mean<0),"h_mean"] <- -0.001
  
  h_mean <- c()
  h_pc0 <- c()
  fer_mean <- c()
  fer_pc <- c()
  irr_area_p <- c()
  rice_ha <- c()
  maize_ha <- c()
  wheat_ha <- c()
  cereal_ha <- c()
  oil_p_ha <- c()
  soy_ha <- c()
  
  for(i in 1:nrow(samples)){
    neighbours <- unlist(strsplit(samples[i,"neighbours"],","))
    neighbours <- as.integer(neighbours[neighbours!="NA"])
    h_mean <- c(h_mean, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"h_mean"], na.rm=TRUE))
    h_pc0 <- c(h_pc0, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"h_pc0"], na.rm=TRUE))
    fer_mean <- c(fer_mean, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"fer_mean"], na.rm=TRUE))
    fer_pc <- c(fer_pc, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"fer_pc"], na.rm=TRUE))
    irr_area_p <- c(irr_area_p, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"irr_area_p"], na.rm=TRUE))
    rice_ha <- c(rice_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"rice_ha"], na.rm=TRUE))
    cereal_ha <- c(cereal_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"cereal_ha"], na.rm=TRUE))
    wheat_ha <- c(wheat_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"wheat_ha"], na.rm=TRUE))
    maize_ha <- c(maize_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"maize_ha"], na.rm=TRUE))
    oil_p_ha <- c(oil_p_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"oil_p_ha"], na.rm=TRUE))
    soy_ha <- c(soy_ha, mean(gridXYattributes[which(gridXYattributes$HBWID%in%neighbours),"soy_ha"], na.rm=TRUE))  
  }
  
  samples$h_mean <- h_mean
  samples$h_pc0 <- h_pc0
  samples$fer_mean <- fer_mean
  samples$fer_pc <- fer_pc
  samples$irr_area_p <- irr_area_p
  samples$rice_ha <- rice_ha
  samples$cereal_ha <- cereal_ha
  samples$wheat_ha <- wheat_ha
  samples$maize_ha <- maize_ha
  samples$oil_p_ha <- oil_p_ha
  samples$soy_ha <- soy_ha
  
  #samples$rice_ha[samples$rice_ha == 0] <- 0.0001
  
  # Cut LUI into categories
  #samples$h_mean_cut <- cut(samples$h_mean, quantile(samples$h_mean, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$h_mean_cut <- cut(samples$h_mean, breaks=classIntervals(samples$h_mean, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$h_mean_excut <- cut(samples$h_mean, quantile(samples$h_mean, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
    
  #samples$fer_mean_cut <- cut(samples$fer_mean, quantile(samples$fer_mean, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$fer_mean_cut <- cut(samples$fer_mean, breaks=classIntervals(samples$fer_mean, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$fer_mean_excut <- cut(samples$fer_mean, quantile(samples$fer_mean, probs=c(0,0.1,0.9,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$fer_pc_cut <- cut(samples$fer_pc, quantile(samples$fer_pc, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$fer_pc_cut <- cut(samples$fer_pc, breaks=classIntervals(samples$fer_pc, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$fer_pc_excut <- cut(samples$fer_pc, quantile(samples$fer_pc, probs=c(0,0.1,0.9,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$irr_area_p_cut <- cut(samples$irr_area_p, quantile(samples$irr_area_p, probs=c(0,0.33,0.66000003,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$irr_area_p_cut <- cut(samples$irr_area_p, breaks=classIntervals(samples$irr_area_p, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$irr_area_p_excut <- cut(samples$irr_area_p, quantile(samples$irr_area_p, probs=c(0,0.1,0.9,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$rice_ha_cut <- cut(samples$rice_ha, quantile(samples$rice_ha, probs=c(0,0.3332,0.6678,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$rice_ha_cut <- cut(samples$rice_ha, breaks=classIntervals(samples$rice_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$rice_ha_excut <- cut(samples$rice_ha, quantile(samples$rice_ha, probs=c(0,0.1,0.9,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$cereal_ha_cut <- cut(samples$cereal_ha, quantile(samples$cereal_ha, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$cereal_ha_cut <- cut(samples$cereal_ha, breaks=classIntervals(samples$cereal_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$cereal_ha_excut <- cut(samples$cereal_ha, quantile(samples$cereal_ha, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$wheat_ha_cut <- cut(samples$wheat_ha, quantile(samples$wheat_ha, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$wheat_ha_cut <- cut(samples$wheat_ha, breaks=classIntervals(samples$wheat_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$wheat_ha_excut <- cut(samples$wheat_ha, quantile(samples$wheat_ha, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$maize_ha_cut <- cut(samples$maize_ha, quantile(samples$maize_ha, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$maize_ha_cut <- cut(samples$maize_ha, breaks=classIntervals(samples$maize_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$maize_ha_excut <- cut(samples$maize_ha, quantile(samples$maize_ha, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$oil_p_ha_cut <- cut(samples$oil_p_ha, quantile(samples$oil_p_ha, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$oil_p_ha_cut <- cut(samples$oil_p_ha, breaks=classIntervals(samples$oil_p_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$oil_p_ha_excut <- cut(samples$oil_p_ha, quantile(samples$oil_p_ha, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  #samples$soy_ha_cut <- cut(samples$soy_ha, quantile(samples$soy_ha, probs=c(0,0.33,0.66,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  samples$soy_ha_cut <- cut(samples$soy_ha, breaks=classIntervals(samples$soy_ha, 3, style="jenks")$brks, labels=c("Low","Medium","High"))
  #samples$soy_ha_excut <- cut(samples$soy_ha, quantile(samples$soy_ha, probs=c(0,0.10,0.90,1), na.rm=TRUE), labels=c("Low","Medium","High"))
  
  ###############################################################################################
  ### Save everything
  
  #save.image(file=SAR_database_output)
  write.csv(samples,SAR_database_output)
  
}

