library(pheatmap)
library(dplyr)
Bac.counts <- read.csv(file = "Normalized_data.txt", header = TRUE, sep= "\t", row.names=1)
Bac2 = Bac.counts[, 4:ncol(Bac.counts)]
Bac.factors <- read.csv(file = "metadata.txt", header = TRUE, sep= "\t", row.names=1)
Bac.factorsDS <- select(Bac.factors, group)
pdf('Heatmap_NanoString.pdf') 
pheatmap(Bac2, annotation_col = Bac.factorsDS)
dev.off()

