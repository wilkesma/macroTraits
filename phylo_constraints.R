setwd("~/Dropbox/Research/Trait work/Trait niches/Final R scripts") #!!!!Change filepaths for GitHub

##################################
## Function to quantify phylogenetic signals and produce phylogenetic correlograms
## Assumes that each taxon in the trait data may correspond zero, one or multiple tips in the phylogenetic tree

## Arguments:
## Trait data (standardised if fuzzy - see Chevenet et al. 1994) with taxa as rows and trait modalities as columns. The first column should contain the short name for each taxon. Trait category and modality names should be separated by a period in the column names, i.e. "category.modality"
traits <- read.csv("traits.standardised.csv", row.names=1)

## Phylogenetic tree in NEXUS format
phylo.tree <- read.nexus("trait.niches.mrbayes.nex.con.tre")

## Data frame relating tip labels to unique species short name in the trait database (NB: As a minumum this should include columns named "short.name" and "tip.label", which rows corresponding to all tip labels in the phylogenetic tree)
label.info <- read.csv("tip_label_names.csv")

## Number of replicate samples to take from the full phylogenetic tree. Suggest nreps>=100 for trait data which includes taxa assigned at a coarse taxonomic resolution (e.g. order)
nreps <- 4

## (optional) Return plots of p values and correlograms. Defaults to TRUE. If TRUE, function outputs two pdfs to the working directory: ("Phylosignal p values.pdf") a boxplot of pvalues per trait category; ("Phylosignal correlations.pdf") correlograms of cumulative phylogenetic correlations per trait category
plots <- T

## (optional) Number of cores for computing phylogenetic signal. Defaults to 1
cores <- 4

## Value:
## A list containing two data frames: (1) sample means, standard deviations, minima and maxima of p values for each trait category from phyloSignal function (phylosignal package); (2) means and 95% confidence intervals of phylogenetic signals at 100 regular intervals of phylogenetic distance, as well as the null correlation (Moran's I)
##################################

phylo_constraints <- function(traits, phylo.tree, label.info, nreps, plots, cores){
  library(ape)
  library(phylosignal)
  library(phylobase)
  library(fields)
  
  if(missing(cores)){
    cores=1
  }
  if(cores>1){
    library(parallel)
  }
  if(missing(plots)){
    plots=T
  }
  
  ## Extract category names
  get.categories <- function(x) { sub("(.*?)[\\.|:].*", "\\1", colnames(x)) } #Functions to extract category names on the fly
  categories <- unique(get.categories(traits)) #Category names
  
  ## Sample tip labels if required (NB: unneccessary if each taxon in the trait database corresponds to just a single tip in the tree)
  sample.species <- function(k){
    sampled.species <- list()
    for(i in 1:k){ 
      sampled.species[[i]] <- sapply(unique(label.info$short.name),
                                     function(x) sample(label.info$tip.label[label.info$short.name==x], 1 ))
    }
    sampled.species
  }
  
  sampled.species <-  sample.species(k=nreps) #List of k matrices. Each matrix is a sample of of n tips from tree (where n is number of taxa in trait database for which phylogenetic data were available)
  
  ## Sample phylogenetic tree
  sample.tips <- function(samples){
    tips.to.drop <- phylo.tree$tip.label[!(phylo.tree$tip.label %in% gsub(" ", "_", samples))]
    tree.subset <- drop.tip(phylo.tree, tips.to.drop)
    tree.subset$tip.label <- names(samples)
    tree.subset
  }
  phylo.trees.subsets <- lapply(sampled.species, sample.tips) #List of k trees as selections for unique taxa from the trait database for which phylogenetic data were available
  
  get.phylostats <- function(tree.sample){
    traits.tree <- traits[tree.sample$tip.label,]
    phylo4d.subset <- phylo4d(tree.sample, tip.data=traits.tree)
    phylosig.subset <- phyloSignal(phylo4d.subset)
    phylosig.traits.subset <- row.names(phylosig.subset$pvalue)[phylosig.subset$pvalue$K.star<0.05] #Gives list of traits significant at 0.05 level under k.star statistic (Blomberg S.P., Garland Jr T. & Ives A.R. (2003) Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57, 717-745)
    phylosig.signal.subset <- data.frame(k.star=phylosig.subset$stat$K.star, p=phylosig.subset$pvalue$K.star)
    row.names(phylosig.signal.subset) <- colnames(traits.tree)
    
    get.phylocor <- function(x){
      traits.tree.subset <- traits.tree[,get.categories(traits.tree)==x]
      phylocor.tree.subset <- phyloCorrelogram(phylo4d.subset, trait=colnames(traits.tree.subset))
      h0 <- -1/(phylocor.tree.subset$n - 1)
      phylocor.tree.res.i <- cbind(phylocor.tree.subset$res, rep(h0, nrow(phylocor.tree.subset$res)))
      colnames(phylocor.tree.res.i) <- c("phylo.dist", "lower", "upper", "mean", "null.I")
      phylocor.tree.res.i
    }
    phylocor.tree.res <- lapply(categories, get.phylocor)
    names(phylocor.tree.res) <- categories
    list(phylosig=phylosig.signal.subset, phylocor=phylocor.tree.res)
  }
  
  phylostats <- mclapply(phylo.trees.subsets, FUN=get.phylostats, mc.cores=cores) #List of k lists of (1) K* and p valuesfrom phyloSignal and (2) data from phyloCorrelogram for plotting
  
  ## Calculate p values from sample
  phylosig.p <- t(do.call(cbind, lapply(phylostats, FUN=function(x) x$phylosig$p)))
  colnames(phylosig.p) <- row.names(phylostats[[1]]$phylosig)
  phylosig.p.summary <- data.frame(mean=colMeans(phylosig.p), sd=apply(phylosig.p, 2, sd), min=apply(phylosig.p, 2, min), max=apply(phylosig.p, 2, max))
  row.names(phylosig.p.summary) <- row.names(phylostats[[1]]$phylosig)
  
  if(plots==T){
    pdf("Phylosignal p values.pdf", paper='a4')
    par(mfrow=c(1,1), mgp=c(1.6, 0.6, 0), mar=c(4.1, 8.1, 1.1, 8.1))
    boxplot(phylosig.p, las=2, horizontal=T, xaxt='n', xlab="Phylosignal p-value")
    axis(1)
    abline(v=0.05, lty=2)
    dev.off()
  }
  
  ## Calculate data for correlograms from sample
  phylosig.cor <- list()
  for(i in 1:length(categories)){
    phylosig.trait.list <- lapply(phylostats, FUN=function(x) x$phylocor[names(x$phylocor)==categories[[i]]])
    phylosig.trait.array <- array(as.numeric(unlist(phylosig.trait.list)), dim=c(100, 5, nreps))
    colnames(phylosig.trait.array) <- c("dist", "lower", "upper", "mean", "null.I")
    phylosig.cor[[i]] <- phylosig.trait.array
  }
  names(phylosig.cor) <- categories
  
  phylosig.cor.means <- lapply(phylosig.cor, FUN=function(x) apply(x, c(1,2), mean))
  
  if(plots==T){
    plot.phylosig.cor <- function(x){
      trait.to.plot <- as.array(phylosig.cor[names(phylosig.cor)==x])[[1]]
      plot(trait.to.plot[,"dist",1], trait.to.plot[,"mean",1], ylim=range(-0.4,1), type='n', main=x, xlab="Phylogenetic distance", ylab="Trait category distance")
      for(i in 1:nreps){
        polygon(c(trait.to.plot[,"dist",i], rev(trait.to.plot[,"dist",i])), c(trait.to.plot[,"upper",i], rev(trait.to.plot[,"lower",i])),
                col=rgb(0,0,0, max=255,alpha=255/nreps), border=NA)
      }
      abline(h=trait.to.plot[,"null.I",1], lty=2)
    }
    
    pdf("Phylosignal correlations.pdf", paper='a4r')
    par(mgp=c(1.6, 0.6, 0), mar=c(3.1, 3.1, 1.1, 1.1))
    for(i in 1:length(categories)) { plot.phylosig.cor(categories[i]) }
    plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    image.plot(matrix(seq(1:10), seq(1:10)), nlevel = nreps, col=rev(gray((nreps:(nreps*2))/(nreps*2))), legend.only=T)
    dev.off()
  }
  
  list(phylosig=phylosig.p.summary, phylocor=phylosig.cor)
}
