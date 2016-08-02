# find DEGs with voom/limma

library( limma )
library( edgeR )
library( splines )

setwd("/proj/price4/hd/hepatocyte_rnaseq/20160802")

counts = read.csv("counts.csv",row.names=1)
meta = read.csv("metadata.csv",row.names=1)

setwd("/proj/price4/hd/hepatocyte_rnaseq/degs")

group = paste( meta$Genotype , meta$Group , meta$Time , sep="_" )

# normalize the counts with voom

counts[ is.na(counts) ] = 0
y = DGEList( counts , group = group )
drop = apply( y , 1 , median ) < 10
y = y[ drop == F , ]
y = voom( y )

design = model.matrix( ~ 0 + group )
colnames(design) = gsub("group","",colnames(design) )

# model fitting
# strategy 1: test for a genotype x treatment effect at each time point
# test for significant separately for each time point and simultaneously across all time points
# focus on contrasts related to interactions

contrasts = makeContrasts(
  genoXtrt.Etop.0.5 = (Het_Etop_0.5-Het_Cont_0.5) - (WT_Etop_0.5-WT_Cont_0.5) ,
  genoXtrt.Etop.1 = (Het_Etop_1-Het_Cont_1) - (WT_Etop_1-WT_Cont_1) ,
  genoXtrt.Etop.2 = (Het_Etop_2-Het_Cont_2) - (WT_Etop_2-WT_Cont_2) ,
  genoXtrt.Etop.3 = (Het_Etop_3-Het_Cont_3) - (WT_Etop_3-WT_Cont_3) ,
  genoXtrt.Etop.6 = (Het_Etop_6-Het_Cont_6) - (WT_Etop_6-WT_Cont_6) ,
  genoXtrt.Etop.24 = (Het_Etop_24-Het_Cont_24) - (WT_Etop_24-WT_Cont_24) ,
  genoXtrt.Etmx.0.5 = (Het_Etmx_0.5-Het_Cont_0.5) - (WT_Etmx_0.5-WT_Cont_0.5) ,
  genoXtrt.Etmx.1 = (Het_Etmx_1-Het_Cont_1) - (WT_Etmx_1-WT_Cont_1) ,
  genoXtrt.Etmx.2 = (Het_Etmx_2-Het_Cont_2) - (WT_Etmx_2-WT_Cont_2) ,
  genoXtrt.Etmx.3 = (Het_Etmx_3-Het_Cont_3) - (WT_Etmx_3-WT_Cont_3) ,
  genoXtrt.Etmx.6 = (Het_Etmx_6-Het_Cont_6) - (WT_Etmx_6-WT_Cont_6) ,
  genoXtrt.Etmx.24 = (Het_Etmx_24-Het_Cont_24) - (WT_Etmx_24-WT_Cont_24) ,
  levels = design )

fit = lmFit( y , design )
fit = contrasts.fit( fit , contrasts )
fit = eBayes( fit )

topTable( fit , coef = 2 )
topTable( fit , coef = 12 )

# strategy 2: test for a single interaction effect across all time points
# use a spline to correct for a main effect of time
# compute a single genotype x treatment effect across all time points

# fit a spline to the time data
X = ns(meta$Time,df=3)

group = paste( meta$Genotype , meta$Group , sep = "_" )
design = model.matrix( ~ 0 + group + X )
colnames(design) = gsub("group","",colnames(design))

contrasts = makeContrasts(
  genoXtrt.Etop = (Het_Etop-Het_Cont)-(WT_Etop-WT_Cont) ,
  genoXtrt.Etmx = (Het_Etmx-Het_Cont)-(WT_Etmx-WT_Cont) ,
  WT.EtopVsCont = WT_Etop-WT_Cont ,
  Het.EtopVsCont = Het_Etop-Het_Cont ,
  WT.EtmxVsCont = WT_Etmx-WT_Cont ,
  Het.EtmxVsCont = Het_Etmx-Het_Cont ,
  HetVsWT.Cont = Het_Cont-WT_Cont ,
  HetVsWT.Etop = Het_Etop-WT_Etop ,
  HetVsWT.Etmx = Het_Etmx-WT_Etmx ,
  levels = design )

fit2 = lmFit( y , design )
fit2 = contrasts.fit( fit2 , contrasts )
fit2 = eBayes( fit2 )

genoXtrt.Etop = topTable( fit2 , coef = 1 , number = Inf , sort.by = "none" )
genoXtrt.Etmx = topTable( fit2 , coef = 2 , number = Inf , sort.by = "none" )
WT.EtopVsCont = topTable( fit2 , coef = 3 , number = Inf , sort.by = "none" )
Het.EtopVsCont = topTable( fit2 , coef = 4 , number = Inf , sort.by = "none" )
WT.EtmxVsCont = topTable( fit2 , coef = 5 , number = Inf , sort.by = "none" )
Het.EtmxVsCont = topTable( fit2 , coef = 6 , number = Inf , sort.by = "none" )
HetVsWT = topTable( fit2 , coef = 7:9 , number = Inf , sort.by = "none" )
Etop = topTable( fit2 , coef = 3:4 , number = Inf , sort.by = "none" )
Etmx = topTable( fit2 , coef = 5:6 , number = Inf , sort.by = "none" )

genoXtrt.Etmx[ order( genoXtrt.Etmx$P.Value )[1:10] , ]
topTable( fit2 , coef = 3:4 )
topTable( fit2 , coef = 5:6 )
topTable( fit2 , coef = 7:9 )

etop.gxt = rownames( genoXtrt.Etop[ genoXtrt.Etop$P.Value < 0.05 , ] )
etmx.gxt = rownames( genoXtrt.Etmx[ genoXtrt.Etmx$P.Value < 0.05 , ] )
etop = rownames( Etop[ Etop$adj.P.Val < 0.05 , ] )
etmx = rownames( Etmx[ Etmx$adj.P.Val < 0.05 , ] )
geno = rownames( HetVsWT[ HetVsWT$adj.P.Val < 0.05 , ] )

genes = rownames(Etop)

degs = data.frame(
  Etop = genes %in% etop ,
  Etmx = genes %in% etmx ,
  Genotype = genes %in% geno ,
  EtopXGeno = genes %in% etop.gxt ,
  EtmxXGeno = genes %in% etmx.gxt )

vennCounts( degs[,c(1,3,4) ] )
vennCounts( degs[,c(2,3,5) ] )

fisher.test( degs[,2] , degs[,1] )
#p-value = 0.004409
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.111911 1.831797
#sample estimates:
#odds ratio 
#  1.434112 


fisher.test( degs[,2] , degs[,3] )
#data:  degs[, 2] and degs[, 3]
#p-value = 1.408e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.1298233 0.5523278
#sample estimates:
#odds ratio 
#  0.286883 


fisher.test( degs[,1] , degs[,3] )
#data:  degs[, 1] and degs[, 3]
#p-value = 6.592e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.4396637 0.7749731
#sample estimates:
#odds ratio 
# 0.5885222 

pdf( "venn.pdf" )
vennDiagram( degs[,1:3] )
dev.off()





 









