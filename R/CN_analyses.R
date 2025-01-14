# This script is for building the Consensus Networks (CNs) from the individual networks of those estuaries with same disturbance level, and analyze the network properties.
# Networks are built for each of the six estuaries using the OTU tables and after applying the prevalence filter of 50%
# The CNs are built and exported to Cytoscape for visualization and first analyses. 
# The networks are analyzed in Cytoscape and exported (graphml format) for further analyses in R.

setwd("/home/leire/Documentos/IndiRed/wTO/")

library(wTO)
library(igraph)
library(NetIndices)
source("../R/utils/filterCoOccuringwTO.R")

# Individual networks of slightly disturbed estuaries (EBI20 and EO10) 

# OTU tables is read and a prevalence filter is applied before network calculation 
# EBI20 
ebi16 <- read.csv("../SWARM_Oct2022_S/16S/EBI20_final_taxa.csv", row.names = 1)
ebi16.pa<- decostand(ebi16, method = "pa")
ebi16.pa.pr<- ebi16.pa[rowSums(ebi16.pa)>=15,]
ebi16.prev<- ebi16[rownames(ebi16) %in% rownames(ebi16.pa.pr),]
ebi16.prev<-as.data.frame(ebi16.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(ebi16.prev)){
  acf(t(ebi16.prev[i,]))
} 

#Network calculation
ebi16.complete = wTO.Complete(k = 1, n = 250, Data = ebi16.prev, Overlap = row.names(ebi16.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
dim(ebi16.complete$wTO) # 85905

#EO10
eor16 <- read.csv("../SWARM_Oct2022_S/16S/EO10_final_taxa.csv", row.names = 1)
eor16.pa<- decostand(eor16, method = "pa")
eor16.pa.pr<- eor16.pa[rowSums(eor16.pa)>=15,]
eor16.prev<- eor16[rownames(eor16) %in% rownames(eor16.pa.pr),]
eor16.prev<-as.data.frame(eor16.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(eor16.prev)){
  acf(t(eor16.prev[i,]))
}

# Network calculation
eor16.complete = wTO.Complete(k = 1, n = 250, Data = eor16.prev, Overlap = row.names(eor16.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
dim(eor16.complete$wTO) # 142845

# CN building and analyses of slightly disturbed estuaries (EBI20 and EO10) 

# CN calculation
ebi_eor_16 = wTO.Consensus(data = list (ebi16.complete = data.frame
                                        (Node.1 = ebi16.complete$wTO$Node.1, 
                                          Node.2 = ebi16.complete$wTO$Node.2, 
                                          wTO = ebi16.complete$wTO$wTO_sign,
                                          pval = ebi16.complete$wTO$pval_sig), 
                                        eor16.complete = data.frame
                                        (Node.1 = eor16.complete$wTO$Node.1, 
                                          Node.2 = eor16.complete$wTO$Node.2, 
                                          wTO = eor16.complete$wTO$wTO_sign,
                                          pval = eor16.complete$wTO$pval_sig)))

# Total common nodes: 398
# Filtering correlations with CN value > 0.4 and p values > 0.001
ebi_eor_16_filt = subset(ebi_eor_16, abs(ebi_eor_16$CN)> 0.4 & ebi_eor_16$pval.fisher < 0.001)

# Exporting CN for visualization and analyzing in Cytoscape
write.table(ebi_eor_16_filt, file = "consensus_nets/EBI20-EO10_CN.txt", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

# Network properties 
ebi_eor_net <- read_graph("consensus_nets/EBI20-EO10_CN.graphml", format = c("graphml"))
plot.igraph(ebi_eor_net, vertex.size=8, vertex.label.cex=.2, edge.arrow.size=.2)

vcount(ebi_eor_net) # 243 nodes
ecount(ebi_eor_net) # 1433 edges

ebi_eor_net.adjacency<-get.adjacency(ebi_eor_net, sparse=F)
ebi_eor_net.properties <- GenInd(ebi_eor_net.adjacency) 
ebi_eor_net.properties # Link Density 5.897119 ;  Connectance 0.02436826

ebi_eor_net_comm <- walktrap.community(ebi_eor_net)
modularity(ebi_eor_net, membership(ebi_eor_net_comm)) # 0.3316436


# Individual networks of moderately disturbed estuaries (EU08 and ED10) 

# OTU tables is read and a prevalence filter is applied before network calculation 
# EU08
eur16 <- read.csv("../SWARM_Oct2022_S/16S/EU08_final_taxa.csv", row.names = 1)
eur16.pa<- decostand(eur16, method = "pa")
eur16.pa.pr<- eur16.pa[rowSums(eur16.pa)>=15,]
eur16.prev<- eur16[rownames(eur16) %in% rownames(eur16.pa.pr),]
eur16.prev<-as.data.frame(eur16.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(eur16.prev)){
  acf(t(eur16.prev[i,]))
}

# Network calculation
eur16.complete = wTO.Complete(k = 1, n = 250, Data = eur16.prev, Overlap = row.names(eur16.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
dim(eur16.complete$wTO) # 103740

# ED10
ede16 <- read.csv("../SWARM_Oct2022_S/16S/ED10_final_taxa.csv", row.names = 1)
ede16.pa<- decostand(ede16, method = "pa")
ede16.pa.pr<- ede16.pa[rowSums(ede16.pa)>=11,]
ede16.prev<- ede16[rownames(ede16) %in% rownames(ede16.pa.pr),]
ede16.prev<-as.data.frame(ede16.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(ede16.prev)){
  acf(t(ede16.prev[i,]))
}

# Network calculation
ede16.complete = wTO.Complete(k = 1, n = 250, Data = ede16.prev, Overlap = row.names(ede16.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
dim(ede16.complete$wTO) # 138601

# CN building and analyses of moderately disturbed estuaries (EU08 and ED10) 

# CN calculation
ede_eur = wTO.Consensus(data = list (ede16.complete = data.frame
                                     (Node.1 = ede16.complete$wTO$Node.1, 
                                       Node.2 = ede16.complete$wTO$Node.2, 
                                       wTO = ede16.complete$wTO$wTO_sign,
                                       pval = ede16.complete$wTO$pval_sig), 
                                     eur16.complete = data.frame
                                     (Node.1 = eur16.complete$wTO$Node.1, 
                                       Node.2 = eur16.complete$wTO$Node.2, 
                                       wTO = eur16.complete$wTO$wTO_sign,
                                       pval = eur16.complete$wTO$pval_sig)))

#Total common nodes: 408

# Filtering correlations with CN value > 0.4 and p values > 0.001
ede_eur_filt = subset(ede_eur, abs(ede_eur$CN)> 0.4 & ede_eur$pval.fisher < 0.001)

# Exporting CN for visualization and analyzing in Cytoscape
write.table(ede_eur_filt, file = "consensus_nets/ED10_EU08_CN.txt", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

# Network properties 
eur_ede_net <- read_graph("consensus_nets/ED10_EU08_CN.graphml", format = c("graphml"))
plot.igraph(eur_ede_net, vertex.size=8, vertex.label.cex=.2, edge.arrow.size=.2)

vcount(eur_ede_net) # 230 nodes
ecount(eur_ede_net) # 1195 edges

eur_ede_net.adjacency<-get.adjacency(eur_ede_net, sparse=F)
eur_ede_net.properties <- GenInd(eur_ede_net.adjacency) 
eur_ede_net.properties # Link Density 5.195652 ;  Connectance 0.02268844

eur_ede_net_comm <- walktrap.community(eur_ede_net)
modularity(eur_ede_net, membership(eur_ede_net_comm)) # 0.5006103


# Individual networks of heavily disturbed estuaries (EOK05 and EOI15) 

# OTU tables is read and a prevalence filter is applied before network calculation 
# EOK05
eok516 <- read.csv("../SWARM_Oct2022_S/16S/EOK05_final_taxa.csv", row.names = 1)
eok516.pa<- decostand(eok516, method = "pa")
eok516.pa.pr<- eok516.pa[rowSums(eok516.pa)>=11,]
eok516.prev<- eok516[rownames(eok516) %in% rownames(eok516.pa.pr),]
eok516.prev<-as.data.frame(eok516.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(eok516.prev)){
  acf(t(eok516.prev[i,]))
}

#Network calculation
eok516.complete = wTO.Complete(k = 1, n = 250, Data = eok516.prev, Overlap = row.names(eok516.prev), method = 's' , 
                               method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                               expected.diff = 0.2, plot = T)
dim(eok516.complete$wTO) # 92235

#EOI15
eoi16 <- read.csv("../SWARM_Oct2022_S/16S/EOI15_final_taxa.csv", row.names = 1)
eoi16.pa<- decostand(eoi16, method = "pa")
eoi16.pa.pr<- eoi16.pa[rowSums(eoi16.pa)>=11,]
eoi16.prev<- eoi16[rownames(eoi16) %in% rownames(eoi16.pa.pr),]
eoi16.prev<-as.data.frame(eoi16.prev)

# Autocorrelation checked for selecting the appropriate lag in the wTO.Complete function
par(mfrow = c(3,3))
for ( i in 1:nrow(eoi16.prev)){
  acf(t(eoi16.prev[i,]))
}

#Network calculation
eoi16.complete = wTO.Complete(k = 1, n = 250, Data = eoi16.prev, Overlap = row.names(eoi16.prev), method = 's' , 
                              method_resampling = 'BlockBootstrap', lag = 1, pvalmethod = 'BH', savecor = TRUE, 
                              expected.diff = 0.2, plot = T)
dim(eoi16.complete$wTO) #56616

# CN building and analyses of heavily disturbed estuaries (EOK05 and EOI15)

eoi_eok5_16 = wTO.Consensus(data = list (eoi16.complete = data.frame
                                         (Node.1 = eoi16.complete$wTO$Node.1, 
                                           Node.2 = eoi16.complete$wTO$Node.2, 
                                           wTO = eoi16.complete$wTO$wTO_sign,
                                           pval = eoi16.complete$wTO$pval_sig), 
                                         eok516.complete = data.frame
                                         (Node.1 = eok516.complete$wTO$Node.1, 
                                           Node.2 = eok516.complete$wTO$Node.2, 
                                           wTO = eok516.complete$wTO$wTO_sign,
                                           pval = eok516.complete$wTO$pval_sig)))

#Total common nodes: 273

# Filtering correlations with CN value > 0.4 and p values > 0.001
eoi_eok5_16_filt = subset(eoi_eok5_16, abs(eoi_eok5_16$CN)> 0.4 & eoi_eok5_16$pval.fisher < 0.001)

# Exporting CN for visualization and analyzing in Cytoscape
write.table(eoi_eok5_16_filt, file = "consensus_nets/EOI15_EOK05_CN.txt", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

# Network properties
eoi_eok05_net <- read_graph("consensus_nets/EOI15_EOK05_CN.graphml", format = c("graphml"))
plot.igraph(eoi_eok05_net, vertex.size=8, vertex.label.cex=.2, edge.arrow.size=.2) 

vcount(eoi_eok05_net) # 135 nodes
ecount(eoi_eok05_net) # 340 edges

eoi_eok05_net.adjacency<-get.adjacency(eoi_eok05_net, sparse=F)
eoi_eok05_net.properties <- GenInd(eoi_eok05_net.adjacency) 
eoi_eok05_net.properties # Link Density 2.518519 ;  Connectance 0.01879491

eoi_eok05_net_comm <- walktrap.community(eoi_eok05_net) 
modularity(eoi_eok05_net, membership(eoi_eok05_net_comm)) # 0.5939144
