# macroTraits
A set of functions to test trait correlations, phylogenetic constraints and spatial variability at large scales.

# fuzzy_trait_correlations

Description:
Function to produce null correlation matrix for fuzzy traits given the trait categories and fuzzy scoring rules used

Usage:
fuzzy_trait_correlations(tr, n.samples, n.cors)

Arguments:<br/>
tr         A data frame containing fuzzy traits (numeric). Trait category and modality names should be separated by a period in the column names e.g. category.modality. Taxon names as row names<br/>
n.samples The number of times to resample to fuzzy trait matrix (integer). Should be >>nrow(tr) and preferably a very large number, e.g. 1 million)<br/>
n.cors    The number of times to draw a sample of nrow(tr) from the resampled trait distributions for calculating pairwise null trait correlations

Value:
A list of null correlation matrices with length equal to n.cors for comparison with observed pairwise trait correlations

# phylo_constraints

Description:
Function to quantify phylogenetic signals and produce phylogenetic correlograms.

Usage:
phylo_constraints(traits, phylo.tree, label.info, nreps, plots, cores)

Arguments:<br/>
traits      Trait data (standardised if fuzzy - see Chevenet et al. 1994) with taxa as rows and trait modalities as columns. The first column should contain the short name for each taxon. Trait category and modality names should be separated by a period in the column names, i.e. "category.modality"<br/>
phylo.tree  Phylogenetic tree in NEXUS format. Function assumes that each taxon in traits may correspond zero, one or multiple tips in phylo.tree<br/>
label.info  Data frame relating tip labels to unique species short names in traits (NB: As a minumum this should include columns named "short.name" and "tip.label", which rows corresponding to all tip labels in the phylogenetic tree)<br/>
n.reps      Number of replicate samples to take from the full phylogenetic tree (integer). Suggest nreps>=100 for trait data which includes taxa assigned at a coarse taxonomic resolution (e.g. order, family)<br/>
plots       (optional) Return plots of p values and correlograms (logical). Defaults to TRUE. If TRUE, function outputs two pdfs to the working directory: ("Phylosignal p values.pdf") a boxplot of pvalues per trait category; ("Phylosignal correlations.pdf") correlograms of cumulative phylogenetic correlations per trait category<br/>
cores       (optional) Number of cores for computing phylogenetic signal (integer). Defaults to 1<br/>

Value:
A list containing two data frames: (1) sample means, standard deviations, minima and maxima of p values for each trait category from phyloSignal function (phylosignal package); (2) means and 95% confidence intervals of phylogenetic signals at 100 regular intervals of phylogenetic distance, as well as the null correlation (Moran's I)

# spatial_gradient

Description:
Function to sample spatial gradients in taxonomic and trait structure using gradient forests

Usage:
spatial_gradient(taxa.names, traits, alt, file.path, n, k, ncores)

Arguments:<br/>
taxa.names    Names for the taxa to be included in the analysis. These names must match filenames (without file extension) of the SDM rasters contained in file.path<br/>
traits        Trait data (standardised if fuzzy - see Chevenet et al. 1994) with taxa as rows and trait modalities as columns. The first column should contain the name for each taxon, which should match exactly the names in taxa.names. Trait category and modality names should be separated by a period in the column names, i.e. "category.modality"<br/>
alt           Altitude (raster) with the same resolution, extent and origin of SDM rasters for each taxon<br/>
file.path     File path to SDM rasters for each taxon. Must all be of the same resolution, extent and origin of alt<br/>
n             Number of grid cells to select per sample. Suggest n>=1000 e.g. at 100 m resolution and continental extent<br/>
k             Number of times to sample n grid cells. Suggest k>=1000<br/>
ncores        (optional) Number of cores to use for parallel processing. Defaults to 1

Value:
A list of dataframes: means and standard deviations for the turnover functions along eastings, northings and altitude for both taxonomic structure and trait structure, the overall goodness-of-fit for taxonomic and trait models, and the goodness-of-fit per taxon and per trait.
