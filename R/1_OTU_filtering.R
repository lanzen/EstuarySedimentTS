# This script filters the classified OTU tables from the Oct 2022 clustering (sed only)
# 1 Comparison of Spatial heterogeneity and export for network analysis (EOK20 and bf and after dredge)
# 2 test of new vs old PowerSoil kits
#
# First rare OTUs below detection are removed, then sediment samples are filtered for plankton
# Extraction and sampling replicates are compared and then pooled. 
# Then data is exported for CREST re-run into CREST_final dirs

# SilvaMod138PR2 used for classification

## CHANGE TO WORKING DIR
## setwd("~/projects/IndiRed/SWARM_Oct2022_S/")

source('R/filtering.R')
source('R/taxaplot.R')
source('R/diversity.r')
source('R/mergeOTUTable.R')
source("RseasonalityIndVal.R")
source("RtaxFunctions.R")

require(vegan)

## Read and check metadata
md.all = read.table("metadata_all.csv",header=T,row.names=1,sep="\t")
md.all$X16SName <- gsub("\\-","\\.",md.all$X16SName) #for 16S
md.16S = md.all[md.all$X16S,]
dim(md.16S) #827 samples
# md.16S.s = md.16S[md.16S$Type == "S",]
# dim(md.16S.s) #510 samples

## Read and check annotated SV (SWARM) table
otus.all.16S = read.delim("SWARM_table_curated.tsv.gz",
                          sep="\t",header=T,row.names=1)
dim(otus.all.16S) #204,225 SVs, 346 samples

# Remove Classification column
otus.all.16S = otus.all.16S[,-dim(otus.all.16S)[2]]

# Parse CREST4 taxonomy results
tax.16S=makeTaxonomy(crest4_assignment_file = "assignments.txt.gz")
summary(tax.16S)
table(row.names(tax.16S) %in% row.names(otus.all.16S)) # Different order

## Check that all included samples are in metadata or should be removed if not
names(otus.all.16S)[!(names(otus.all.16S) %in% md.16S$X16SName)]

# Include only those sampels that are in the OTU table
md.16S = md.16S[md.16S$X16SName %in% names(otus.all.16S),]

## Transpond OTU table and put in order of Metadata
otus.t.16S = as.data.frame(t(otus.all.16S[,md.16S$X16SName]))

## Check that SV and metadata names are identical and set names
table(row.names(otus.t.16S) == md.16S$X16SName) 
row.names(otus.t.16S) = row.names(md.16S)
tax.16S = tax.16S[names(otus.t.16S),]

sum(otus.t.16S) # 26,773,625

write.csv(file = "16S_unfiltered_div.csv",data.frame(reads=rowSums(otus.t.16S),
                                                       richness=specnumber(otus.t.16S)),
            quote=F)

# ----------- Filtering of potential contaminants ------------

## Inspect OTUs present in blanks

otus_in_blanks = otus.t.16S[,colSums(otus.t.16S[md.16S$Type=="C",])>0]
dim(otus_in_blanks) # 2572 OTUs in controls with no overlap

otus.ra.16S = decostand(otus.t.16S,method="total")
otus.blank.ra = otus.ra.16S[,colSums(otus.t.16S[md.16S$Type=="C",]) > 0 &
                              colSums(otus.t.16S[md.16S$Type!="C",]) > 0]

md.16S$RunPlot = md.16S$Run16S
md.16S[md.16S$Type=="C",]$RunPlot = "Control"
otus.blank.rap = mergeOTUTable(otus.blank.ra, by="RunPlot",metadata=md.16S)
otus.blank.rap = decostand(otus.blank.rap,method="total")
names(otus.blank.rap) = make.names(paste(names(otus.blank.rap), tax.16S[names(otus.blank.rap),]$bestTx,sep="_"))

pdf("img/16S/OTUs_in_controls_RA.pdf", height=8, width=10)
taxaplot(30, data.frame(row.names=row.names(otus.blank.rap), 
                        row.names(otus.blank.rap)),
         otus.blank.rap)
dev.off()
# SWARM_4556_Vibrio looks like contaminant

# Categories

otus.byType = mergeOTUTable(otus.ra.16S, metadata = md.16S, behaviour = "mean", by="Type")
md.byType = mergeMetadata(md.16S, by="Type")
cat.byType = mergeOTUsByTaxa(otus.byType, taxon_table = tax.16S, rank="category")

pdf("img/16S/categories_bf_filtering.pdf", height=8, width=8)
taxaplot(12, data.frame(row.names=row.names(cat.byType), 
                        row.names(cat.byType)),
         cat.byType)
dev.off()

## Filter using decontam
require(decontam)
om = as.matrix(otus.t.16S)
predicted_contaminants = isContaminant(om, neg=(md.16S$Type=="C"), batch=md.16S$Run16S)
                                       #conc=md.16S$DNA.16S, method="combined")
summary(predicted_contaminants) #114 predicted contaminants at p<.1
sc = (!is.na(predicted_contaminants$p) & predicted_contaminants$p<=.05) 
sum(sc) #47 contaminants w p<0.05

sigCont = predicted_contaminants[sc,]
sigCont$Tax = tax.16S$bestTx[row.names(tax.16S) %in% row.names(sigCont)]
sigCont$Cat = droplevels(tax.16S$category[row.names(tax.16S) %in% row.names(sigCont)])


summary(as.factor(sigCont$Tax)) 
# 4xGammaproteobacteria and then a little of each

summary(sigCont$Cat)
# Chemoautotrophic aerobic prok Heterotrophic aerobic prokaryotes      Heterotrophic anaerobic prok              Host associated prok 
# 3                                30                                 6                                 5 
# Mitochondria  Prok-prok symbionts or parasites                              NA's 
#                                 1                                 1                                 1 

## Remove contaminants
otus.16S.clean = otus.t.16S[,!(names(otus.t.16S) %in% row.names(sigCont))]

## Remove blank samples
otus.16S.clean = otus.16S.clean[md.16S$Type!="C",]
otus.16S.clean = otus.16S.clean[,colSums(otus.16S.clean)>0]
dim(otus.16S.clean) #204,106 OTUs
sum(otus.16S.clean) #26,748,297
sum(otus.t.16S) - sum(otus.16S.clean) #25,328 removed
md.16S.s = droplevels(md.16S[md.16S$Type!="C",])

# ---- Filter probable cross-contamination reads ----

otus.16S.clean = filterCrossContaminants2(otus.16S.clean)
sum(otus.16S.clean) #  26,658,865 (1.1k removed)

# ---- Filtering all unclassified and eukaryotes  -----

for (taxonDel in c("No hits","Eukaryota")){
  toDel = row.names(tax.16S)[grep(taxonDel,tax.16S$classification)]
  if (length(toDel)>0){
    otus.16S.clean = otus.16S.clean[,!(names(otus.16S.clean) %in% toDel)]
    print(paste("Removing ",taxonDel))
    print(dim(otus.16S.clean)[2])
    print(sum(otus.16S.clean))
  }
}
# 1] "Removing  No hits"
# [1] 201492
# [1] 26657162
# [1] "Removing  Eukaryota"
# [1] 198455
# [1] 23667288

# 3M reads removed, mostly chloroplasts

# ---- Abundance filtering ----

table(row.names(md.16S.s) == row.names(otus.16S.clean))

summary(rowSums(otus.16S.clean))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#483   43514   64304   71290   91447  430874 

## Remove all samples with coverage below 5000 reads
row.names(otus.16S.clean)[rowSums(otus.16S.clean)<5000]
# "EOK20-Spatial21"  "EUBalenciaga-A03"

otus.16S.clean = otus.16S.clean[rowSums(otus.16S.clean)>5000,]
md.16S.s = md.16S.s[row.names(otus.16S.clean),]


## Abundance filtering:
## Remove rare OTUs below detection limit of ~5 reads in at least 1 dataset 

otus.s.16S = dropRareByMaxAbundance(otus.16S.clean,1E-3)
dim(otus.s.16S) # 3474 OTUs (out of 197k)
sum(otus.s.16S) # 16,063,855 reads (7M removed)

# Write diversity stats
writeDivStats("16SDiversity_Sed_abFiltered.csv",otus.s.16S)
div.s.16S = read.csv("16SDiversity_Sed_abFiltered.csv", header=T,row.names = 1)
table(row.names(md.16S.s) == row.names(div.s.16S))


# ------------- Merge replicates and export final OTU table  ------------

require(lubridate)
md.16S.s$Month = month(md.16S.s$Date)

md.16S.final = mergeMetadata(md.16S.s, by = "Proper_name")
otus.s.16S.final = mergeOTUTable(otus.s.16S, metadata = md.16S.s, by = "Proper_name")

table(row.names(md.16S.final) == row.names(otus.s.16S.final)) #309 merged

# Write cleaned up OTU table
write.table(as.data.frame(t(otus.s.16S.final)), file = "16S/SWARM_table_final_S.tsv",
            sep="\t",quote=F,row.names=T, col.names=NA)

md.16S.final$Date = as.Date(md.16S.final$Date)
stationN = as.numeric(as.factor(md.16S.final$Station))
md.16S.final$StationDate = paste(md.16S.final$Station, md.16S.final$Date)
md.16S.final$dayFromFirst = md.16S.final$Date - min(md.16S.final$Date)
md.16S.final$dayOfYear = as.numeric(strftime(md.16S.final$Date, format="%j"))
md.16S.final$winterness = abs(md.16S.final$dayOfYear - 182.5) # Days away from midsummer


writeDivStats("16SDiversity_merged.csv",otus.s.16S.final)
div.16S.final = read.csv("16SDiversity_merged.csv", header=T, row.names=1)[,c("Reads","Richness","Rarefied.richness",
                                                                              "H","J")]


table(row.names(div.16S.final) == row.names(md.16S.final))

div.16S.final$Station = md.16S.final$Station
div.16S.final$dayOfYear = md.16S.final$dayOfYear
div.16S.final$Month = as.factor(md.16S.final$Month)
printANOVA(md.16S.final[,c(9,16,18:21,24,25)], div.16S.final,.01) # nothing interesting

# taxa and categories
tax.16S.final = tax.16S[names(otus.s.16S.final),]
table(row.names(tax.16S.final) == names(otus.s.16S.final))
cat.16S.final = decostand(mergeOTUsByTaxa(otus.s.16S.final, tax.16S.final, rank="category"), method="total")


# ---- Boxplot of div per month and station -----

require(ggplot2)

## Rarefied richness
pdf("img/16S/RRichness_by_Station.pdf",width=10,height=5)
ggplot(div.16S.final, aes(x=Station, y=Rarefied.richness)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                                          outlier.size=2) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/16S/RRichness_by_Month.pdf",width=12,height=5)
ggplot(div.16S.final, aes(x=Month, y=Rarefied.richness)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                                        outlier.size=2) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/16S/RRichness_time.pdf",width=7,height=5)
plot(Rarefied.richness~dayOfYear, data=div.16S.final, col = stationN, pch= stationN, cex=.8)

# Add lowess lines
for(st in unique(stationN)){
  div.st = div.16S.final[stationN==st,]
  lines(lowess(div.st$dayOfYear,div.st$Rarefied.richness), col=st, lwd=2, lt=2)
}
legend("bottomright",pch=unique(stationN),col=unique(stationN),
       legend=unique(md.16S.final$Station),ncol=3,cex=.7)
dev.off()

## Shannon
pdf("img/16S/Shannon_by_Station.pdf",width=10,height=5)
ggplot(div.16S.final, aes(x=Station, y=H)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                          outlier.size=2) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/16S/Shannon_by_Month.pdf",width=12,height=5)
ggplot(div.16S.final, aes(x=Month, y=H)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                        outlier.size=2) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/16S/Shannon_time.pdf",width=12,height=5)
plot(H~dayOfYear, data=div.16S.final, col = stationN, pch= stationN, cex=.8)

# Add lowess lines
for(st in unique(stationN)){
  div.st = div.16S.final[stationN==st,]
  lines(lowess(div.st$dayOfYear,div.st$H), col=st, lwd=2, lt=2)
}
legend("bottomright",pch=unique(stationN),col=unique(stationN),
       legend=unique(md.16S.final$Station),ncol=3,cex=.7)
dev.off()

# ------ Total div per site -------

mdf = md.16S.final[md.16S.final$Station %in% c("EBI20","ED10","EO10","EOI15","EOK05","EU08") & md.16S.final$TimeSeries=="S",]
otusf = otus.s.16S.final[md.16S.final$Station %in% c("EBI20","ED10","EO10","EOI15","EOK05","EU08") & md.16S.final$TimeSeries=="S",]
otus.perSite = mergeOTUTable(otusf, mdf, by="Station")
writeDivStats("Diversity_per_site.tsv", otus.perSite)


# -------- Make subsets, export and draw graphs of each site ---------

otus.s.16S.final.t = as.data.frame(t(otus.s.16S.final))
at = mergeOTUTable(otus.s.16S.final.t, metadata = tax.16S.final, by="bestTx")
assignments.16S.s.t.ra = decostand(as.data.frame(t(at)),method="total")
dim(assignments.16S.s.t.ra) # 830 taxa
assignments.16S.s.t.ra = dropRareByMaxAbundance(assignments.16S.s.t.ra,1E-3)
dim(assignments.16S.s.t.ra) # 827 taxa

library(tidyverse)
library(ggplot2)
INC_LIMIT_OTU = .005
INC_LIMIT_TAXA = .005

table(row.names(md.16S.final) == row.names(assignments.16S.s.t.ra)) #309

#for (st in unique(md.16S.final$Station)){
# Station EBI20 with 402 days and 30 observations.
# Station ED10 with 315 days and 22 observations."
# Station EO10 with 408 days and 31 observations
# Station EOI15 with 538 days and 22 observations
# Station EOK05 with 353 days and 22 observations
# Station EOK20 with 375 days and 21 observations
# Station EU08 with 433 days and 31 observations

md.st = md.16S.final[md.16S.final$Station==st,]
obs = dim(md.st)[1]
obs
ts_length = as.numeric(max(md.st$Date) - min(md.st$Date))


## Select and export taxon table
tt = assignments.16S.s.t.ra[row.names(md.st),]
tt = tt[,colSums(tt)>0]

tc = cat.16S.final[row.names(md.st),]
tc = tc[,colSums(tc)>0]

## Graph of assigned including SX

tt.ab = tt[,colMeans(tt)> INC_LIMIT_TAXA]
noT = dim(tt.ab)[2] 

tt.ab$Date = md.st$Date

pdf(paste("img/16S/", st, "_stacked_ntaxa.pdf", sep=""), height=7, width=ts_length/100)

colours <- c(rainbow(noT-23),"#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
             "#0075DC", "#F0A3FF", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
             "#FFA8BB","#426600","#5EF1F2","#00998F");
             

tt.ab %>%
  gather(variable, value, 1:noT) %>%
  ggplot(aes(x = Date, y = value, fill = variable)) +
  scale_fill_manual(values=colours) +
  geom_area(color="black",lwd=.1)
dev.off()

## Category plot

pdf(paste("img/16S/", st, "_stacked_categories.pdf", sep=""), height=7, width=ts_length/100)
tc$Date = md.st$Date
tc %>%
  gather(variable, value, 1:dim(tc)[2]-1) %>%
  ggplot(aes(x = Date, y = value, fill = variable)) +
  scale_fill_manual(values=colours) +
  geom_area(color="black",lwd=.1)
dev.off()

## Remove extra samples
md.st2 = md.st[md.st$TimeSeries=="S",]
tt = assignments.16S.s.t.ra[row.names(md.st2),]
write.table(cbind(names(tt), as.data.frame(t(tt))), row.names=FALSE, col.names=c("taxon",row.names(tt)),
            file=paste("16S/",st,"_final_taxa.csv",sep=""),quote=F,sep=",")

ts_length = as.numeric(max(md.st2$Date) - min(md.st2$Date))
print(paste("Station",st,"with",ts_length,"days and",dim(md.st2)[1],"observations."))

tt.ab = tt[,colMeans(tt)> INC_LIMIT_TAXA]
noT = dim(tt.ab)[2] 
noT # 39 in EBI20, 45 for ED10, 41 for EO10, 36 in EOI15, 36 for EOK5, 41 in EOK20, 38 in EU08
tt.ab$Date = md.st2$Date

pdf(paste("img/16S/", st, "_stacked_taxa.pdf", sep=""), height=7, width=ts_length/25)

colours <- c(rainbow(noT-23),"#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
             "#0075DC", "#F0A3FF", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
             "#FFA8BB","#426600","#5EF1F2","#00998F");
             
tt.ab %>%
  gather(variable, value, 1:noT) %>%
  ggplot(aes(x = Date, y = value, fill = variable)) +
  scale_fill_manual(values=colours) +
  geom_area(color="black",lwd=.1)
dev.off()

# Categories
tc = cat.16S.final[row.names(md.st2),]

pdf(paste("img/16S/", st, "_stacked_categories_IR.pdf", sep=""), height=7, width=ts_length/25)
tc$Date = md.st2$Date
tc %>%
  gather(variable, value, 1:dim(tc)[2]-1) %>%
  ggplot(aes(x = Date, y = value, fill = variable)) +
  scale_fill_manual(values=colours) +
  geom_area(color="black",lwd=.1)
dev.off()

#}
