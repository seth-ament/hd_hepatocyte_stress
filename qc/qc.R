# check for outlier samples by MDS and hierarchical clustering

library( edgeR )
source("/proj/price1/sament/resources/WGCNA/R/Functions.R")

setwd("/proj/price4/hd/hepatocyte_rnaseq/qc")

counts0 = read.csv("/proj/price4/hd/hepatocyte_rnaseq/20160725/counts.csv")
meta0 = read.csv("/proj/price4/hd/hepatocyte_rnaseq/20160725/metadata.csv")

counts = counts0[,-1]
rownames(counts) = counts0[,1]

meta = meta0
meta$Sample = paste( "X" , meta$Sample , sep="" )
all( colnames(counts) == meta$Sample )
# TRUE

group = paste( meta$Genotype , meta$Group , meta$Time , sep = "_" )
design = model.matrix( ~ 0 + group )
colnames(design) = gsub("group","",colnames(design))

y = round( counts )
y[ is.na(y) ] = 0
y = DGEList( counts = y , group = group )
y = calcNormFactors( y )

pdf("mds.pdf")
plotMDS( y , labels = group )
plotMDS( y , labels = meta$Batch )
dev.off()

cpm = cpm( y )
sampleTree = hclust(dist(t(cpm)), method = "average")
traitColors = data.frame(
   genotype = as.numeric( meta$Genotype ) ,
   treatment = as.numeric( meta$Group ) , 
   time = as.numeric( meta$Time ) , 
   batch = as.numeric( meta$Batch ) )

pdf("hclust.pdf", height = 7 , width = 10 )
plotDendroAndColors( sampleTree, traitColors ,
     main = "Sample clustering to detect outliers" )
dev.off()

# both MDS and hclust reveal that the 24-hour samples are
# very different from the other timepoints.
# Jeff Carroll confirms that this is consistent
# with what they saw by qPCR for a handful of genes.

# hclust also reveals a handful of outliers that
# do not match any of the biological variables.

# remove the outliers not explained by time
# retain the 24-hour samples but analyze them separately
# from the other time points.

clusters = cutree( sampleTree , h = 80000 )
counts2 = counts[ , clusters != 2 ]
meta2 = meta[ clusters != 2 , ]

setwd("/proj/price4/hd/hepatocyte_rnaseq/20160802")
write.csv( counts2 , file="counts.csv" )
write.csv( meta2 , file="metadata.csv" )










