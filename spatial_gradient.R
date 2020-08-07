setwd("~/Dropbox/Research/Trait work/Trait niches/Final R scripts") #!!!!Change filepaths for GitHub

##################################
## Function to sample spatial gradients in taxonomic and trait structure using gradient forests

##Arguments:
load("trait.niches.spatial.gf.data.v7.2.RData")

## Names for the taxa to be included in the analysis. These names must match filenames (without file extension) of the SDM rasters contained in file.path
taxa.names <- names(genera.mods)[1:3]

## Trait data (standardised if fuzzy - see Chevenet et al. 1994) with taxa as rows and trait modalities as columns. The first column should contain the name for each taxon, which should match exactly the names in taxa.names. Trait category and modality names should be separated by a period in the column names, i.e. "category.modality"
traits <- read.csv("traits.standardised.csv", row.names=1)[1:3,]
row.names(traits) <- taxa.names

## Altitude (raster) with the same resolution, extent and origin of SDM rasters for each taxon
alt <- euro.env$altitude
rm(euro.env)

## File path to SDM rasters for each taxon. Must all be of the same resolution, extent and origin of alt
file.path <- "~/Dropbox/Research/Trait work/Trait niches/v7/SDM rasters v3/"

## Number of grid cells to select per sample. Suggest n>=1000 e.g. at 100 m resolution and continental extent
n=10

## Number of times to sample n grid cells. Suggest k>=1000
k=10

## (Optional) number of cores to use for parallel processing. Defaults to 1
ncores=2

## Value:
## A list of six dataframes, each providing mean and standard deviations for the turnover functions along eastings, northings and altitude for both taxonomic structure and trait structure, the overall goodness-of-fit for taxonomic and trait models, and the goodness-of-fit per taxon and per trait
##################################

spatial_gradient <- function(taxa.names, traits, alt, file.path, n, k, ncores){
  library(raster)
  library(doParallel)
  library(foreach)
  require(matrixStats)
  library(abind)
  library(gradientForest)
  
  if(missing(ncores)){
    ncores=1
  }
  
  ## Sample n cells k times
  get.sample.cell <- function(x){
    repeat{
      sample.cell.no <- sample(1:ncell(alt), size=1)
      sample.cell.preds <- extract(alt, sample.cell.no)
      if(!is.na(sample.cell.preds)) break
    }
    return(sample.cell.no)
  }
  sample.cells <- array(dim=c(n,k))
  for(i in 1:n){
    for(j in 1:k){
      sample.cells[i,j] <- get.sample.cell() 
    }
  } #Rows (i) are cells sampled [n], columns (j) are samples of n cells (m)
  
  #Get predicted occurrence probability in each sampled cell for each genus
  get.taxon.pred <- function(taxon){
    taxon.raster <- raster(paste(file.path, taxon, ".tif", sep=""))
    pred.cells <- matrix(data=NA, nrow=n, ncol=k)
    for(i in 1:k){
      pred.cells[,i] <- taxon.raster[sample.cells[,i]]
    }
    pred.cells
  }
  taxa.preds.samples <- lapply(taxa.names, FUN=get.taxon.pred) #Gives list of matrices, one matrix per species with rows as sampled cells and columns as samples of n cells
  taxa.preds.samples <- array(as.numeric(unlist(taxa.preds.samples)), dim=c(n, k, length(taxa.preds.samples)))
  taxa.preds.samples <- aperm(taxa.preds.samples, c(1,3,2)) #Each slice is now sampled cells (rows) and species (columns); third dimensions are samples of n cells
  colnames(taxa.preds.samples) <- taxa.names
  
  ##GF
  #Prepare traits
  get.categories <- function(x) { sub("(.*?)[\\.|:].*", "\\1", colnames(x)) } #Functions to extract category names on the fly
  categories <- unique(get.categories(traits)) #Category names
  
  #get CWMs
  traits.preds.samples <- array(dim=c(n, ncol(traits), k))
  for(i in 1:k){
    for(j in 1:n){
      rel.abundance <- taxa.preds.samples[j,,i]/sum(taxa.preds.samples[j,,i])
      weights <- apply(traits, MARGIN=2, FUN=function(x) rel.abundance*x)
      traits.preds.samples[j,,i] <- colSums(weights)
    }
  }
  row.names(traits.preds.samples) <- seq(1:n)
  colnames(traits.preds.samples) <- colnames(traits)
  
  #Reprojecting coordinates of sampled cells
  get.coords <- function(x){
    coords.m <- as.data.frame(coordinates(alt)[x,])
    coordinates(coords.m) <- c("x", "y")
    proj4string(coords.m) <- CRS(proj4string(alt))
    coords.m <- spTransform(coords.m, "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") #Transforms To Lambert Azimuthal Equal Area projection centred on Europe - see https://gis.stackexchange.com/questions/182417/what-is-the-proj4string-of-the-lambert-azimuthal-equal-area-projection
  }
  coords.m <- apply(sample.cells, MARGIN=2, FUN=get.coords) #List of spatial points objects, one for each set of 100 sampled cells
  
  #Getting altitude data for sampled cells
  altitude.samples <- apply(sample.cells, MARGIN=2, FUN=function(x) extract(alt, x)) #Matrix of altitude in sampled cells*sets of sampled cells
  
  #Running GF for taxonomic and trait matrices
  gf.data.taxonomic <- list()
  gf.data.functional <- list()
  for(k in 1:k){
    gf.data.taxonomic[[k]] <- as.data.frame(cbind(coords.m[[k]]@coords/1000000, altitude.samples[,k], taxa.preds.samples[,,k]))
    colnames(gf.data.taxonomic[[k]])[1:3] <- c("easting", "northing", "altitude")
    gf.data.functional[[k]] <- as.data.frame(cbind(coords.m[[k]]@coords/1000000, altitude.samples[,k], traits.preds.samples[,,k]))
    colnames(gf.data.functional[[k]])[1:3] <- c("easting", "northing", "altitude")
  }
  
  get.gf <- function(x){
    message('Getting gradientForest results for sample ', x, ' of ', k, sep="")
    #Getting gf objects
    gf.mod.taxonomic <- gradientForest(data=gf.data.taxonomic[[x]], predictor.vars=c("easting", "northing", "altitude"),
                                       response.vars=colnames(taxa.preds.samples[,,x]), ntree=500, trace=TRUE)
    gf.mod.functional <- gradientForest(data=gf.data.functional[[x]], predictor.vars=c("easting", "northing", "altitude"),
                                        response.vars=colnames(traits.preds.samples[,,x]), ntree=500, trace=TRUE)
    #Function to extract R2 per species/trait and predictor
    get.gf.R2 <- function(model, mod.names){
      gf.R2.vars <- model$res.u[,c("spec", "var", "rsq.var")]
      gf.R2.vars <- lapply(c("easting", "northing", "altitude"), FUN=function(x) gf.R2.vars[gf.R2.vars$var==x,])
      gf.R2.vars <- do.call(cbind, gf.R2.vars)[,c(1,3,6,9)]
      row.names(gf.R2.vars) <- gf.R2.vars$spec
      gf.R2.vars <- gf.R2.vars[,-1]
      colnames(gf.R2.vars) <- c("easting", "northing", "altitude")
      '%ni%' <- Negate('%in%')
      names.no.fit <- mod.names[which(mod.names %ni% row.names(gf.R2.vars))]
      no.fit <- matrix(data=0, nrow=length(names.no.fit), ncol=3)
      row.names(no.fit) <- names.no.fit; colnames(no.fit) <- c("easting", "northing", "altitude")
      rbind(gf.R2.vars, no.fit)
    }
    #Function to extract turnover function per predictor
    get.turnover <- function(variable, model){
      turnover <- lapply(unique(model$res$spec), FUN=function(x) subset(model$res, spec==x & var==variable))
      turnover <- do.call(rbind, lapply(turnover, FUN=function(x) as.data.frame(x[order(x$split),c("split", "improve")])))
      colnames(turnover) <- c('split', 'improve')
      turnover <- turnover[order(turnover$split),]
      turnover <- as.data.frame(cbind(turnover$split, cumsum(turnover$improve))[seq(1, nrow(turnover), 100),]) #Only taking every 100th row to reduce data volume; essentially no detail is lost from the plot
      colnames(turnover) <- c('split', 'cum.imp')
      imp <- importance(gf.mod.taxonomic)
      imp <- imp[names(imp)==variable]
      turnover$cum.imp <- imp*(turnover$cum.imp/max(turnover$cum.imp))
      turnover
    }
    #Preparing results, including assigning NA in any cases where forest was empty
    if(!is.null(gf.mod.taxonomic)){
      gf.taxonomic.imp <- importance(gf.mod.taxonomic)
      gf.R2.per.taxon <- get.gf.R2(gf.mod.taxonomic, colnames(taxa.preds.samples[,,x]))
      gf.turnover.taxonomic <- lapply(c("easting", "northing", "altitude"), FUN=function(x) get.turnover(x, gf.mod.taxonomic)) #List of 3 turnover functions for the taxonomic model, each corresponding to a predictor variable
      names(gf.turnover.taxonomic) <- c("easting", "northing", "altitude")
    } else {
      gf.taxonomic.imp <- c(easting=NA, northing=NA, altitude=NA)
      gf.R2.per.taxon <- t(replicate(length(colnames(taxa.preds.samples[,,x])), c(easting=NA, northing=NA, altitude=NA)))
      row.names(gf.R2.per.taxon) <- colnames(taxa.preds.samples[,,x])
      gf.turnover.taxonomic <- list(easting=c(split=NA, cum.imp=NA), northing=c(split=NA, cum.imp=NA), altitude=c(split=NA, cum.imp=NA))
    }
    if(!is.null(gf.mod.functional)){
      gf.functional.imp <- importance(gf.mod.functional)
      gf.R2.per.trait <- get.gf.R2(gf.mod.functional, colnames(traits.preds.samples[,,x]))
      gf.turnover.functional <- lapply(c("easting", "northing", "altitude"), FUN=function(x) get.turnover(x, gf.mod.functional)) #List of 3 turnover functions for the functional model, each corresponding to a predictor variable
      names(gf.turnover.functional) <- c("easting", "northing", "altitude")
    } else {
      gf.functional.imp <- c(easting=NA, northing=NA, altitude=NA)
      gf.R2.per.trait <- t(replicate(length(colnames(traits.preds.samples[,,x])), c(easting=NA, northing=NA, altitude=NA)))
      row.names(gf.R2.per.trait) <- colnames(traits.preds.samples[,,x])
      gf.turnover.functional <- list(easting=c(split=NA, cum.imp=NA), northing=c(split=NA, cum.imp=NA), altitude=c(split=NA, cum.imp=NA))
    }
    
    #Listing results
    results <- list(taxonomic.imp=gf.taxonomic.imp, functional.imp=gf.functional.imp, taxonomic.R2=ifelse(is.null(gf.mod.taxonomic$result), NA, gf.mod.taxonomic$result), functional.R2=ifelse(is.null(gf.mod.functional$result), NA, gf.mod.functional$result), R2.per.taxon=gf.R2.per.taxon, R2.per.trait=gf.R2.per.trait, turnover.taxonomic=gf.turnover.taxonomic, turnover.functional=gf.turnover.functional)
    rm(gf.mod.taxonomic, gf.mod.functional)
    results
  }
  
  cl <- parallel::makeCluster(ncores, setup_timeout = 0.5)
  registerDoParallel(cl)
  gf.results <- foreach(i=1:k, .packages='gradientForest', .verbose=T) %dopar% get.gf(i)
  stopCluster(cl)
  
  ##Extracting gf results
  taxonomic.imp <- t(rbind(colMeans(do.call(rbind, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$taxonomic.imp))),
                           colSds(do.call(rbind, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$taxonomic.imp)))))
  colnames(taxonomic.imp) <- c("mean", "sd")
  
  functional.imp <- t(rbind(colMeans(do.call(rbind, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$functional.imp))),
                            colSds(do.call(rbind, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$functional.imp)))))
  colnames(functional.imp) <- c("mean", "sd")
  
  taxonomic.R2 <- c(do.call(mean, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$taxonomic.R2)),
                    sd(do.call(c, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$taxonomic.R2))))
  functional.R2 <- c(do.call(mean, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$functional.R2)),
                     sd(do.call(c, lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$functional.R2))))
  
  R2.per.taxon <- abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.taxon), along=3)
  R2.per.taxon <- aperm(R2.per.taxon, c(3,1,2))
  R2.per.taxon[,,"easting"] #For example; colums are taxa, rows are samples
  
  R2.per.trait <- abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.trait), along=3)
  R2.per.trait <- aperm(R2.per.trait, c(3,1,2))
  R2.per.trait[,,"easting"] #For example; columns are traits, rows are samples
  
  R2.per.taxon.mean <- as.data.frame(apply(abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.taxon), along=3), c(1,2), mean))
  R2.per.taxon.sd <- as.data.frame(apply(abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.taxon), along=3), c(1,2), sd))
  
  R2.per.trait.mean <- as.data.frame(apply(abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.trait), along=3), c(1,2), mean))
  R2.per.trait.sd <- as.data.frame(apply(abind(lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$R2.per.trait), along=3), c(1,2), sd))
  
  sample.easting.ranges <- do.call(rbind, lapply(seq(1:length(gf.results)), function(x) c(min(gf.results[[x]]$turnover.taxonomic$easting[,1]), max(gf.results[[x]]$turnover.taxonomic$easting[,1]))))
  sample.northing.ranges <- do.call(rbind, lapply(seq(1:length(gf.results)), function(x) c(min(gf.results[[x]]$turnover.taxonomic$northing[,1]), max(gf.results[[x]]$turnover.taxonomic$northing[,1]))))
  sample.altitude.ranges <- do.call(rbind, lapply(seq(1:length(gf.results)), function(x) c(min(gf.results[[x]]$turnover.taxonomic$altitude[,1]), max(gf.results[[x]]$turnover.taxonomic$altitude[,1]))))
  sample.ranges <- as.data.frame(rbind(c(max(sample.easting.ranges[,1]), min(sample.easting.ranges[,2])),
                                       c(max(sample.northing.ranges[,1]), min(sample.northing.ranges[,2])),
                                       c(max(sample.altitude.ranges[,1]), min(sample.altitude.ranges[,2]))))
  colnames(sample.ranges) <- c("common.start", "common.end"); row.names(sample.ranges) <- c("easting", "northing", "altitude")
  #sample.ranges gives the common starts and ends from which to extract turnover functions for both taxo and func models (but need to extract them first, then take only values in this range)
  
  turnover.taxonomic.easting <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.taxonomic$easting)
  turnover.taxonomic.easting <- lapply(turnover.taxonomic.easting, function(x) x[x$split>=sample.ranges$common.start[1] & x$split<=sample.ranges$common.end[1],2])
  min(do.call(c, lapply(turnover.taxonomic.easting, function(x) length(x)))) #29582 is the minimum vector length, so take every 29th value up to 29000
  turnover.taxonomic.easting <- lapply(turnover.taxonomic.easting, function(x) x[seq.int(1L, 29000, 29L)])
  turnover.taxonomic.easting <- as.data.frame(cbind(seq(from=sample.ranges[1,1], to=sample.ranges[1,2], length.out=1000), 
                                                    rowMeans(do.call(cbind, turnover.taxonomic.easting)),
                                                    rowSds(do.call(cbind, turnover.taxonomic.easting))))
  colnames(turnover.taxonomic.easting) <- c("easting", "mean", "sd")
  
  turnover.taxonomic.northing <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.taxonomic$northing)
  turnover.taxonomic.northing <- lapply(turnover.taxonomic.northing, function(x) x[x$split>=sample.ranges$common.start[2] & x$split<=sample.ranges$common.end[2],2])
  min(do.call(c, lapply(turnover.taxonomic.northing, function(x) length(x)))) #29586 is the minimum vector length, so take every 29th value up to 29000
  turnover.taxonomic.northing <- lapply(turnover.taxonomic.northing, function(x) x[seq.int(1L, 29000, 29L)])
  turnover.taxonomic.northing <- as.data.frame(cbind(seq(from=sample.ranges[2,1], to=sample.ranges[2,2], length.out=1000), 
                                                     rowMeans(do.call(cbind, turnover.taxonomic.northing)),
                                                     rowSds(do.call(cbind, turnover.taxonomic.northing))))
  colnames(turnover.taxonomic.northing) <- c("northing", "mean", "sd")
  
  turnover.taxonomic.altitude <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.taxonomic$altitude)
  turnover.taxonomic.altitude <- lapply(turnover.taxonomic.altitude, function(x) x[x$split>=sample.ranges$common.start[3] & x$split<=sample.ranges$common.end[3],2])
  min(do.call(c, lapply(turnover.taxonomic.altitude, function(x) length(x)))) #29582 is the minimum vector length, so take every 29th value up to 29000
  turnover.taxonomic.altitude <- lapply(turnover.taxonomic.altitude, function(x) x[seq.int(1L, 29000, 29L)])
  turnover.taxonomic.altitude <- as.data.frame(cbind(seq(from=sample.ranges[3,1], to=sample.ranges[3,2], length.out=1000), 
                                                     rowMeans(do.call(cbind, turnover.taxonomic.altitude)),
                                                     rowSds(do.call(cbind, turnover.taxonomic.altitude))))
  colnames(turnover.taxonomic.altitude) <- c("altitude", "mean", "sd")
  
  turnover.functional.easting <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.functional$easting)
  turnover.functional.easting <- lapply(turnover.functional.easting, function(x) x[x$split>=sample.ranges$common.start[1] & x$split<=sample.ranges$common.end[1],2])
  min(do.call(c, lapply(turnover.functional.easting, function(x) length(x)))) #6781 is the minimum vector length, so take every 7th value for a gradient 968 values long
  turnover.functional.easting <- lapply(turnover.functional.easting, function(x) x[seq.int(1L, 968*7, 7L)])
  turnover.functional.easting <- as.data.frame(cbind(seq(from=sample.ranges[1,1], to=sample.ranges[1,2], length.out=968), 
                                                     rowMeans(do.call(cbind, turnover.functional.easting)),
                                                     rowSds(do.call(cbind, turnover.functional.easting))))
  colnames(turnover.functional.easting) <- c("easting", "mean", "sd")
  
  turnover.functional.northing <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.functional$northing)
  turnover.functional.northing <- lapply(turnover.functional.northing, function(x) x[x$split>=sample.ranges$common.start[2] & x$split<=sample.ranges$common.end[2],2])
  min(do.call(c, lapply(turnover.functional.northing, function(x) length(x)))) #6781 is the minimum vector length, so take every 7th value for a gradient 968 values long
  turnover.functional.northing <- lapply(turnover.functional.northing, function(x) x[seq.int(1L, 968*7, 7L)])
  turnover.functional.northing <- as.data.frame(cbind(seq(from=sample.ranges[2,1], to=sample.ranges[2,2], length.out=968), 
                                                      rowMeans(do.call(cbind, turnover.functional.northing)),
                                                      rowSds(do.call(cbind, turnover.functional.northing))))
  colnames(turnover.functional.northing) <- c("northing", "mean", "sd")
  
  turnover.functional.altitude <- lapply(seq(1:length(gf.results)), function(x) gf.results[[x]]$turnover.functional$altitude)
  turnover.functional.altitude <- lapply(turnover.functional.altitude, function(x) x[x$split>=sample.ranges$common.start[3] & x$split<=sample.ranges$common.end[3],2])
  min(do.call(c, lapply(turnover.functional.altitude, function(x) length(x)))) #6781 is the minimum vector length, so take every 7th value for a gradient 968 values long
  turnover.functional.altitude <- lapply(turnover.functional.altitude, function(x) x[seq.int(1L, 968*7, 7L)])
  turnover.functional.altitude <- as.data.frame(cbind(seq(from=sample.ranges[3,1], to=sample.ranges[3,2], length.out=968), 
                                                      rowMeans(do.call(cbind, turnover.functional.altitude)),
                                                      rowSds(do.call(cbind, turnover.functional.altitude))))
  colnames(turnover.functional.altitude) <- c("altitude", "mean", "sd")
  
  list(turnover.taxonomic.easting=turnover.taxonomic.easting, turnover.taxonomic.northing=turnover.taxonomic.northing,
       turnover.taxonomic.altitude=turnover.taxonomic.altitude, taxonomic.R2, R2.per.taxon,
       turnover.functional.easting=turnover.functional.easting, turnover.functional.northing=turnover.functional.northing,
       turnover.functional.altitude=turnover.functional.altitude, functional.R2, R2.per.trait)
}

