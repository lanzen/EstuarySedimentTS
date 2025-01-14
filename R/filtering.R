############################################################
## Methods for filtering rare OTUs and cross contaminants ## 
## Adapted for Metabridge, Anders Lanzen 2024-06-26       ##
############################################################



# Compares abs. read abundance to total OTU reads/k as lower limit -> 
# FNs in undersampled. Similar to UNCROSS (Edgar et al.)
# Based on relative abundance which is better for uneven sample sizes
# If method is "mean" we use average abundances. Otherwise max. abundances
filterCrossContaminants2 = function(otus, minRAQuote = 200, taxaAsRows = FALSE,
                                    method="max"){
  if (taxaAsRows) otus = as.data.frame(t(otus))
  otusM = as.matrix(otus)
  require(vegan)
  ra = decostand(otus,method="total")
  raM = as.matrix(ra)
  if (method=="mean") meanTaxonAbundance = as.vector(colMeans(ra))
  else if (method=="max") {
    meanTaxonAbundance = as.vector(apply(ra,2,max))
  }
  otus.clean = matrix(nrow=dim(otus)[1],ncol=dim(otus)[2])
  #Iterate over each column (OTU)
  for (i in c(1:dim(otus)[2])) {
    otus.clean[,i] = ifelse(raM[,i] >= meanTaxonAbundance[i] / minRAQuote, otusM[,i], 0)
  }
  otus.clean = as.data.frame(otus.clean)
  row.names(otus.clean) = row.names(otus)
  names(otus.clean) = names(otus)
  if (taxaAsRows) otus.clean = as.data.frame(t(otus.clean))
  return(otus.clean)
}


# dropRareByMaxAbundance Removes whole taxa if maximum relative abundance is below cutoff
# Good for uneven sampling and uneven distribution, to avoid losing relatively
# rare but discrimnant taxa.
# Recommended minRelativeAbundance is 10/(the number of reads in the least sequeneced sample)

dropRareByMaxAbundance = function(abundances, minRelativeAbundance, 
                                  taxaAsRows = FALSE) {
  # Transpose if taxa are rows
  if (taxaAsRows) abundances = as.data.frame(t(abundances))
  if (rowSums(abundances[1,])>1) ra=decostand(abundances,method="total")
  else ra = abundances
  
  # Filter and return
  maxTaxonAbundances = sapply(ra, max)
  filtered = abundances[,maxTaxonAbundances >= minRelativeAbundance]
  if (taxaAsRows) filtered = (as.data.frame(t(filtered)))
  return(filtered)
}

# dropRareByAvgAbundance Removes whole taxa if average relative abundance is below cutoff
# Good for uneven sampling and similar samples (even distribution)
dropRareByAvgAbundance = function(relativeAbundances, minRelativeAbundance, 
                                  taxaAsRows = FALSE) {
  require(vegan)
  # Transpose if taxa are rows
  if (taxaAsRows) relativeAbundances = as.data.frame(t(relativeAbundances))
  
  # Filter and return
  # avgTaxonAbundances = colMeans(relativeAbundances)
  filtered = relativeAbundances[,avgTaxonAbundances >= minRelativeAbundance]
  if (taxaAsRows) filtered = (as.data.frame(t(filtered)))
  return(filtered)
}

# dropRareTaxa only filters individual occurences (bad with uneven sampling -> FNs in samples with many reads)
# Accept and returns either rel. or abs. abundance df (2017-01-24)
dropRareTaxa <- function(df, minRelativeAbundance=0, maxRelativeAbundance=Inf){
  dfM = as.matrix(df)
  ra = decostand(df,method="total")
  raM = as.matrix(ra)
  
  # Filter below min relative abundances
  fdist=as.data.frame(ifelse(raM>minRelativeAbundance | raM<maxRelativeAbundance, 
                             raM,dfM,0))
  
  names(fdist)=names(df)
  row.names(fdist)=row.names(df)
  nonzero = fdist[,colSums(fdist)>0]
  return(nonzero)
}


## Remove all OTUs with exactly one read
dropSingletons <- function(df){
  require(vegan)
  good = df[,colSums(df)>1]
  return(good)
}

