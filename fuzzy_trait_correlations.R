fuzzy_trait_correlations <- function(tr, n.samples, n.cors){
  library(plyr)
  
  ## Extract category names
  get.categories <- function(x) { sub("(.*?)[\\.|:].*", "\\1", colnames(x)) } #Functions to extract category names on the fly
  categories <- unique(get.categories(tr)) #Category names
  
  ## Fuzzy trait scoring system
  fuzzy.row.max <- vector()
  fuzzy.row.min <- vector()
  fuzzy.max <- vector()
  for(i in 1:length(categories)){
    y <- tr[, which(get.categories(tr) %in% categories[i])]
    fuzzy.row.max[i] <- max(rowSums(y))
    fuzzy.row.min[i] <- min(rowSums(y))
    fuzzy.max[i] <- max(y)
  }
  names(fuzzy.row.max) <- categories #Gives maximum possible total fuzzy score for each trait category. Shows that using raw fuzzy codes would introduce bias as some categories have >double the maximum than others
  names(fuzzy.row.min) <- categories #Gives minimum total fuzzy score in the data for each trait category
  names(fuzzy.max) <- categories #Gives maximum possible total fuzzy score for individual traits in a category
  
  ## Resampling traits given the fuzzy trait scoring system
  resample.traits <- function(trait.data, trait.categories, n, trait.maxs, trait.row.mins, trait.row.maxs){
    gen.traits <- function (x, max, rowMin, rowMax){
      gen <- function(x, max, rowMin, rowMax){
        while(TRUE){
          z <- sample(seq(1:(max+1))-1, size=length(x), replace=T)
          if(sum(z)<=rowMax & sum(z)>=rowMin) return(z) }
      }
      y <- array(dim=c(n,length(x)))
      for(i in 1:n){
        y[i,] <- gen(x, max, rowMin, rowMax)
      }
      y
    } #Function to resample trait distributions within trait categories
    resampled.traits <- list()
    for(i in 1:length(trait.categories)){
      y <- (trait.data[1, ])[,which(get.categories(tr) %in% trait.categories[i])]
      resampled.traits[[i]] <- gen.traits(y, trait.maxs[i], trait.row.mins[i], trait.row.maxs[i])
    }
    resampled.traits <- do.call(cbind, resampled.traits)
    colnames(resampled.traits) <- colnames(trait.data)
    resampled.traits
  } #Function to resample trait distributions across all trait categories
  
  resampled.traits <- resample.traits(trait.data=tr, trait.categories=categories, n=n.samples, trait.maxs=fuzzy.max, trait.row.mins=fuzzy.row.min, trait.row.maxs=fuzzy.row.max) #Matrix of n hypothetical species (rows) and their hypothethical traits (columns) given the fuzzy scoring system used
  
  ## Unique trait combinations in resampled trait matrix
  max.utcs <- vector()
  for(i in 1:length(categories)){
    max.utcs[i] <- nrow(count(resampled.traits[, which(get.categories(tr) %in% categories[i])]))
  } #Number of theoretical combinations by trait category
  names(max.utcs) <- categories
  utcs.theoretical <- prod(max.utcs) #Number of theoretical combinations as product of individual category combinations
  
  ## Observed unique trait profiles
  singular.taxa <- nrow(unique(tr))
  
  ## Null correlations
  get.cor <- function(x){
    samplei <- x[c(sample(1:nrow(tr), nrow(tr), replace=F)),]
    samplei <- as.data.frame(scale(samplei, center=T, scale=F))
    cor(samplei, method="spearman")
  } #Function samples ntaxa (number of taxa in trait database) from resampled trait distributions and gives pairwise trait correlation matrix
  
  cors.samples <- list()
  for(i in 1:n.cors){
    cors.samples[[i]] <- get.cor(resampled.traits)
  } #Null correlation matrix for comparison with observed correlation matrix
  cors.samples
}
