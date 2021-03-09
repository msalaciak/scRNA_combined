library(SCENIC)
library(SCopeLoomR)

cd8.subset.cc.1 <-subset(cd8.subset.cc, Timepoint==1)
cellInfo.tp1 <- data.frame(cd8.subset.cc.1=Idents(cd8.subset.cc.1))
rm(cd8.subset.cc.1)

cd8.subset.cc.6 <-subset(cd8.subset.cc, Timepoint==6)
cellInfo.tp6 <- data.frame(cd8.subset.cc.6=Idents(cd8.subset.cc.6))
rm(cd8.subset.cc.6)


pyScenicDir <- "/home/matthew/Research/scenic"

cd8tp1loom <- file.path(pyScenicDir, "/cd8TP1.loom_SCENIC.loom")
cd8tp1loom <- open_loom(cd8tp1loom, mode="r")

regulons_incidMat <- get_regulons(cd8tp1loom,column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(cd8tp1loom,column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(cd8tp1loom)

GRNBoost_linkList.tp1 <- importArboreto(file.path(pyScenicDir,  "/adj.tsv"))
head(GRNBoost_linkList.tp1)

motifsDf.tp1 <- data.table::fread(file.path(pyScenicDir, "reg-fix.csv"), header = T, sep=",")

maxRows <- 20 # (low value only for the tutorial)

# Visualize
tableSubset <- motifsDf.tp1[TF=="STAT6"]
tableSubset <- tableSubset[1:maxRows,] 
colsToShow <- colnames(motifsDf.tp1)[-c(2,9:11)]
viewMotifs(motifsDf.tp1,options=list(pageLength=5),colsToShow=colsToShow)

#Regulators for known cell types or clusters
regulonActivity_byCellType <- sapply(split(rownames(cellInfo.tp1), cellInfo.tp1$cd8.subset.cc.1),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

# rss
rss <- calcRSS(AUC=getAUC(regulonsAUC), cellAnnotation=cellInfo.tp1[colnames(regulonsAUC), "cd8.subset.cc.1"])
plotRSS_oneSet(rss, setName = c("CD8 TEM-1")) + plotRSS_oneSet(rss, setName = c("CD8 TEFF"))

##### TIME POINT 6

cd8tp6loom <- file.path(pyScenicDir, "/cd8TP6.loom_SCENIC.loom")
cd8tp6loom <- open_loom(cd8tp6loom, mode="r")

regulons_incidMat.tp6 <- get_regulons(cd8tp6loom,column.attr.name='Regulons')
regulons.tp6 <- regulonsToGeneLists(regulons_incidMat.tp6)
regulonsAUC.tp6 <- get_regulons_AUC(cd8tp6loom,column.attr.name='RegulonsAUC')
regulonsAucThresholds.tp6 <- get_regulon_thresholds(cd8tp6loom)

GRNBoost_linkList.tp6 <- importArboreto(file.path(pyScenicDir,  "/adj-tp6.tsv"))
head(GRNBoost_linkList.tp6)

motifsDf.tp6 <- data.table::fread(file.path(pyScenicDir, "reg-tp6-fix.csv"), header = T, sep=",")

viewMotifs(motifsDf.tp6,options=list(pageLength=5))

# Visualize
 
colsToShow <- colnames(motifsDf.tp6)[-c(2,9:11)]
viewMotifs(motifsDf.tp6,options=list(pageLength=5),colsToShow=colsToShow)

#Regulators for known cell types or clusters
regulonActivity_byCellType.tp6 <- sapply(split(rownames(cellInfo.tp6), cellInfo.tp6$cd8.subset.cc.6),
                                     function(cells) rowMeans(getAUC(regulonsAUC.tp6)[,cells]))
regulonActivity_byCellType_Scaled.tp6 <- t(scale(t(regulonActivity_byCellType.tp6), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled.tp6, name="Regulon activity")

topRegulators.tp6 <- reshape2::melt(regulonActivity_byCellType_Scaled.tp6)
colnames(topRegulators.tp6) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators.tp6 <- topRegulators.tp6[which(topRegulators.tp6$RelativeActivity>0),]
viewTable(topRegulators.tp6)

# rss
rss.tp6 <- calcRSS(AUC=getAUC(regulonsAUC.tp6), cellAnnotation=cellInfo.tp6[colnames(regulonsAUC.tp6), "cd8.subset.cc.6"])
plotRSS_oneSet(rss.tp6, setName = c("CD8 TEM-1")) + plotRSS_oneSet(rss.tp6, setName = c("CD8 TEFF"))



auc_mtx1 <- read.csv("/home/matthew/Research/auc_mtx1.csv", header = TRUE)
auc_mtx6 <- read.csv("/home/matthew/Research/auc_mtx6.csv", header = TRUE)

auc_mtx6 <- as.data.frame(auc_mtx6)
auc_mtx1 <- as.data.frame(auc_mtx1)


plot(density(auc_mtx1$TCF7_...))
abline(v = names(regulonsAucThresholds)[regulonsAucThresholds == "TCF7(+)"])


dr_coords.tp6 <- Embeddings(subset(cd8.subset.cc, Timepoint==6), reduction="umap")

tfs <- c("EOMES(+)","E2F2(+)","TCF7(+)")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords.tp6, cellsAUC=selectRegulons(regulonsAUC.tp6, tfs), plots = "AUC")


