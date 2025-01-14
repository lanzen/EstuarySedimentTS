seasonOfMonth = function(month){
  season = NA
  if (month == 12 | month<=2) season = 1
  if (month >= 3 & month<6) season = 2
  if (month >= 6 & month<9) season = 3
  if (month >= 9 & month<12) season = 4
  return(season)
}

# Finds season using IndVal approach incl Spring+Autumn combined
# metadata must have Month variable with numeric for month of sampling
findSeasonality = function(abundances, metadata, indAlpha = 0.01){

  require(labdsv)
  
  ## Define and handle seasons
  
  seasons = c("Winter","Spring","Summer","Autumn", "Spring + Autumn")
  
  metadata$Season = sapply(metadata$Month, seasonOfMonth)
  
  ## Use IndVal approach to find indicators
  seasonInd = indval(abundances, metadata$Season)
  inds = seasonInd$pval<=indAlpha
  print(paste(sum(inds),"indicator taxa out of",dim(abundances)[2]))

  # Both autumn and spring - override if already existing
  metadata$SeasonAlt = metadata$Season
  metadata$SeasonAlt[metadata$SeasonAlt==4] <- 2
  seasonIndAlt = indval(abundances, metadata$SeasonAlt)
  indsAlt = seasonIndAlt$pval<=indAlpha
  
  indsAll = data.frame(row.names=names(abundances),
                       p.m1 =  seasonInd$pval,
                       p.m2 = seasonIndAlt$pval,
                       season.m1 = seasonInd$maxcls,
                       season.m2 = seasonIndAlt$maxcls,
                       p = seasonInd$pval, season = seasonInd$maxcls
  )
  
  for(ti in c(1:dim(indsAll)[1])){
    if (indsAll$season.m2[ti]==2 & indsAll$p.m2[ti] <= indsAll$p.m1[ti]){
      indsAll$p[ti] = indsAll$p.m2[ti]
      indsAll$season[ti] = 5
    }
  }
  
  indsAll = indsAll[indsAll$p<=indAlpha,]
  indsAll$season = as.factor(indsAll$season)
  
  ## Make final data frame without data on the two models to return
  indsFinal = data.frame(row.names=row.names(indsAll),
                         p=indsAll$p, season = indsAll$season,
                         season_name = seasons[indsAll$season])
  
  return(indsFinal)
}


getCol = function(n){
  return (c(rainbow(n-25),"#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
            "#0075DC", "#F0A3FF", "#993F00","#4C005C","#2BCE48","#FFCC99",
            "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
            "#FFA8BB","#426600","#5EF1F2","#00998F",
            "black","grey"))
}
