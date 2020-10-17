library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(RCurl)
library(cowplot)
library(pheatmap)
library(clusterProfiler)
library("org.Dm.eg.db",character.only = TRUE)
library(DOSE)
library(plotly)
theme_set(theme_cowplot())

#load datasets of each time point
timepoint_1.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/151231/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_2.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/160411/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_3.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/161962/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_4.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/171094/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_5.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/171642/outs/filtered_gene_bc_matrices/GRCh38")
timepoint_6.data <- Read10X(data.dir ="/home/matthew/datatransfer/mercier/180251/outs/filtered_gene_bc_matrices/GRCh38")


#create seurat objects for each
timepoint_1 <-CreateSeuratObject(counts = timepoint_1.data, project = "pbmc timepoint 1", min.cells =3,min.features=200)
timepoint_2 <-CreateSeuratObject(counts = timepoint_2.data, project = "pbmc timepoint 2", min.cells =3,min.features=200)
timepoint_3 <-CreateSeuratObject(counts = timepoint_3.data, project = "pbmc timepoint 3", min.cells =3,min.features=200)
timepoint_4 <-CreateSeuratObject(counts = timepoint_4.data, project = "pbmc timepoint 4", min.cells =3,min.features=200)
timepoint_5 <-CreateSeuratObject(counts = timepoint_5.data, project = "pbmc timepoint 5", min.cells =3,min.features=200)
timepoint_6 <-CreateSeuratObject(counts = timepoint_6.data, project = "pbmc timepoint 6", min.cells =3,min.features=200)

#add metadata for each object
timepoint_1@meta.data$Timepoint <- "1"
timepoint_1@meta.data$Disease_status <- "HL"
timepoint_1@meta.data$Pembrolizumab <- "No"
timepoint_1@meta.data$iRAE <- "No"

timepoint_2@meta.data$Timepoint <- "2"
timepoint_2@meta.data$Disease_status <- "Remission"
timepoint_2@meta.data$Pembrolizumab <- "Yes"
timepoint_2@meta.data$iRAE <- "Yes"

timepoint_3@meta.data$Timepoint <- "3"
timepoint_3@meta.data$Disease_status <- "Remission"
timepoint_3@meta.data$Pembrolizumab <- "Yes"
timepoint_3@meta.data$iRAE <- "No"

timepoint_4@meta.data$Timepoint <- "4"
timepoint_4@meta.data$Disease_status <- "Relapse"
timepoint_4@meta.data$Pembrolizumab <- "No"
timepoint_4@meta.data$iRAE <- "No"

timepoint_5@meta.data$Timepoint <- "5"
timepoint_5@meta.data$Disease_status <- "Relapse"
timepoint_5@meta.data$Pembrolizumab <- "Yes"
timepoint_5@meta.data$iRAE <- "Yes"

timepoint_6@meta.data$Timepoint <- "6"
timepoint_6@meta.data$Disease_status <- "Relapse"
timepoint_6@meta.data$Pembrolizumab <- "No"
timepoint_6@meta.data$iRAE <- "No"

#QC

#Add MT to each time point
timepoint_1[["percent.mt"]] <- PercentageFeatureSet(timepoint_1, pattern = "^MT-")
timepoint_2[["percent.mt"]] <- PercentageFeatureSet(timepoint_2, pattern = "^MT-")
timepoint_3[["percent.mt"]] <- PercentageFeatureSet(timepoint_3, pattern = "^MT-")
timepoint_4[["percent.mt"]] <- PercentageFeatureSet(timepoint_4, pattern = "^MT-")
timepoint_5[["percent.mt"]] <- PercentageFeatureSet(timepoint_5, pattern = "^MT-")
timepoint_6[["percent.mt"]] <- PercentageFeatureSet(timepoint_6, pattern = "^MT-")

#combine seurat objects into a list
seuratObjectList <- list(timepoint_1,timepoint_2,timepoint_3,timepoint_4,timepoint_5,timepoint_6)
  
for (i in 1:length(seuratObjectList)){
print(FeatureScatter(seuratObjectList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(seuratObjectList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

}



#using information from plots, we can filter based on these values
#t ime point 1 -6
# mt= 17, 17, 17, 17, 19, 15
# feature rna= 4700, 4000, 4700, 3800, 3600, 5100




#subset and filter cells based on MT/Counts
seuratObjectList[[1]] <- subset(seuratObjectList[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 17)
seuratObjectList[[2]] <- subset(seuratObjectList[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 17)
seuratObjectList[[3]] <- subset(seuratObjectList[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 17)
seuratObjectList[[4]] <- subset(seuratObjectList[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 3800 & percent.mt < 17)
seuratObjectList[[5]] <- subset(seuratObjectList[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 3600 & percent.mt < 19)
seuratObjectList[[6]] <- subset(seuratObjectList[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 5100 & percent.mt < 15)



#normalize data IF USING FOR LOOP FOR LIST
# for (i in 1:length(seuratObjectList)){
#   seuratObjectList[[i]] <- NormalizeData(seuratObjectList[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
#   seuratObjectList[[i]] <- FindVariableFeatures(seuratObjectList[[i]], selection.method = "vst",
#                                                 dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 2000)
# }

seuratObjectList <- lapply(X = seuratObjectList, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)
})



# 
#integrate data
pbmc.anchors <- FindIntegrationAnchors(object.list = seuratObjectList, dims = 1:20)
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors, dims = 1:20)

DefaultAssay(pbmc.combined) <- "integrated"


#MERGE OPTION IF NEEDED
# pbmc.combined <- merge(seuratObjectList[[1]],list(seuratObjectList[[2]],seuratObjectList[[3]],seuratObjectList[[4]],
#                                                   seuratObjectList[[5]],seuratObjectList[[6]]), project = "combined_pbmc")
# pbmc.combined <- NormalizeData(pbmc.combined, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst",dispersion.cutoff = c(0.5, Inf),
#                                       mean.cutoff = c(0.0125, 3),nfeatures = 2000)
# pbmc.combined <- ScaleData(pbmc.combined)



pbmc.combined <- ScaleData(pbmc.combined, verbose = TRUE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = TRUE)

pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined,resolution = 0.8)

pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:20)

# Visualization
umap_timepoint <- DimPlot(pbmc.combined, reduction = "umap", split.by = "Timepoint", label=TRUE)
oldClusterID <- DimPlot(pbmc.combined, reduction = "umap", label = TRUE)
print(umap_timepoint)
print(oldClusterID)

#integrated count = 26739
Reduce("+",table ( Idents(pbmc.combined) ) )

pbmc.combined[["old.ident"]] <- Idents(object = pbmc.combined)

new.cluster.ids <- c("CD14 Mono 1", "CD14 Mono 2", "CD4 TCM", "CD4 Naive","CD14 Mono 3"," CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3", "CD14 Mono 4"
                     ,"NK" , "CD16 Mono " ,"B Cell" ,"CD8 TEM","CD14 Mono 5" ,"cpDC","CD14 Mono 6", "CD4 Treg" ,"CD8 Exhausted" ,"CD8 Teff" ,"pDC")

names(new.cluster.ids) <- levels(pbmc.combined)
pbmc.combined <- RenameIdents(pbmc.combined, new.cluster.ids)


DimPlot(pbmc.combined,reduction = "umap", split.by = "iRAE")
DimPlot(pbmc.combined,reduction = "umap", split.by = "Disease_status")

DimPlot(pbmc.combined,reduction = "umap",label = TRUE)

timepoint_1.subset <- subset(x = pbmc.combined, subset = Timepoint == "1")
timepoint_2.subset <- subset(x = pbmc.combined, subset = Timepoint == "2")
timepoint_3.subset <- subset(x = pbmc.combined, subset = Timepoint == "3")
timepoint_4.subset <- subset(x = pbmc.combined, subset = Timepoint == "4")
timepoint_5.subset <- subset(x = pbmc.combined, subset = Timepoint == "5")
timepoint_6.subset <- subset(x = pbmc.combined, subset = Timepoint == "6")


t1 <- DimPlot(timepoint_1.subset,reduction = "umap", label = TRUE ,repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 1") + NoLegend()
t2 <-DimPlot(timepoint_2.subset,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 2") +  NoLegend()
t3 <- DimPlot(timepoint_3.subset,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 3") +  NoLegend()
t4 <- DimPlot(timepoint_4.subset,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 4") +  NoLegend()
t5 <- DimPlot(timepoint_5.subset,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 5") +  NoLegend()
t6 <-DimPlot(timepoint_6.subset,reduction = "umap" , label = TRUE,  repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 6") +  NoLegend()



plotlists <-list(t1,t2,t3,t4,t5,t6)


#switch back pbmc.combined assay to RNA from integrated for find markers but switch back to integration for visual
DefaultAssay(pbmc.combined) <- "RNA"


pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

for (i in 0:19){
  print("clustering ")
  print(i)
  markers <- FindMarkers(pbmc.combined, ident.1 = i, only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

  markers <- cbind(gene = rownames(markers), markers)
  rownames(markers) <- 1:nrow(markers)

  markers$cluster <- i
  markers<- markers[, c(2,3,4,5,6,1,7)]

  nam <- paste("cluster", i,".markers", sep = "")
  assign(nam, markers)

}

pbmc_full.markers <-rbind(cluster0.markers, cluster1.markers, cluster2.markers,cluster3.markers,cluster4.markers,cluster5.markers
                          ,cluster6.markers,cluster7.markers,cluster8.markers,cluster9.markers,cluster10.markers,cluster11.markers
                          ,cluster12.markers,cluster13.markers,cluster14.markers,cluster15.markers,cluster16.markers,cluster17.markers
                          ,cluster18.markers,cluster19.markers)

#how many cells per timepoint / cluster
table(Idents(pbmc.combined),pbmc.combined$Timepoint)



# saveRDS(pbmc.combined, file = "pbmc_combined.rds")


cluster16_logfc <- cluster16.markers$avg_logFC


names(cluster16_logfc) <- cluster16.markers$gene
cluster16_logfc <- cluster16_logfc[!is.na(cluster16_logfc)]
cluster16_logfc = sort(cluster16_logfc, decreasing = TRUE)


gse_cluster16 <- gseGO(
  cluster16_logfc,
  ont = "ALL",
  OrgDb =org.Hs.eg.db,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

cluster2_UNIPROT <- mapIds(org.Hs.eg.db, cluster2.markers$gene, 'UNIPROT', 'SYMBOL')
cluster2_logfc <- cluster2.markers$avg_logFC
names(cluster2_logfc) <- cluster2_UNIPROT
cluster2_logfc <- cluster2_logfc[!is.na(names(cluster2_logfc))]
cluster2_logfc <-sort(cluster2_logfc, decreasing = TRUE)


kegg_cluster2 <- gseKEGG(geneList= cluster2_logfc,
                keyType = "uniprot",
               organism     = "hsa",
               minGSSize    = 20,
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               eps=0)


FeaturePlot(object = pbmc.combined, features = c("CD3D","CD3G","CD3E"), cols = c("grey", "red"), reduction = "umap", label=TRUE)
FeaturePlot(object = pbmc.combined, features = c("CD4"), cols = c("grey", "red"), reduction = "umap", label=TRUE)
FeaturePlot(object = pbmc.combined, features = c("CD8A","CD8B"), cols = c("grey", "red"), reduction = "umap", label=TRUE)
FeaturePlot(object = pbmc.combined, features = c("CD8A","CD8B"), cols = c("grey", "red"), reduction = "umap", label=TRUE,split.by="Timepoint")
FeaturePlot(object = pbmc.combined, features = c("CD8A","CD8B"), cols = c("grey", "red"), reduction = "umap", label=TRUE,split.by="iRAE")
FeaturePlot(object = pbmc.combined, features = c("PDCD1"), cols = c("grey", "red"), reduction = "umap", label=TRUE,split.by="Timepoint")
FeaturePlot(object = pbmc.combined, features = c("CD4"), cols = c("grey", "red"), reduction = "umap", label=TRUE,split.by="iRAE")
FeaturePlot(object = timepoint_2.subset, features = c("CD38"), reduction = "umap", label=TRUE)

cluster11_compare <- subset(pbmc.combined, idents = "11")
Idents(cluster11_compare) <- "Timepoint"
cluster11_compare_avg <- log1p(AverageExpression(cluster11_compare, verbose = FALSE)$RNA)
cluster11_compare_avg$gene <- rownames(cluster11_compare_avg)



cluster11_1v2 <-ggplot(cluster11_compare_avg, aes_(as.name(1), as.name(2), 
                      text = paste("Gene: ", cluster11_compare_avg$gene, "\n",
                                   sep = ""))) + 
                                    geom_point() + ggtitle("Cluster 11 Timepoint 1 vs. 2")

ggplot(cluster11_compare_avg, aes_(as.name(2), as.name(3))) + geom_point() + ggtitle("Cluster 11 Timepoint 2 vs. 3")
ggplot(cluster11_compare_avg, aes_(as.name(3), as.name(4))) + geom_point() + ggtitle("Cluster 11 Timepoint 3 vs. 4")

cluster11_4v5 <-ggplot(cluster11_compare_avg, aes_(as.name(4), as.name(5),
                                   text = paste("Gene: ", cluster11_compare_avg$gene, "\n",
                                                sep = ""))) +
                                    geom_point() + ggtitle("Cluster 11 Timepoint 4 vs. 5")


ggplot(cluster11_compare_avg, aes_(as.name(5), as.name(6))) + geom_point() + ggtitle("Cluster 11 Timepoint 5 vs. 6")

print(cluster11_1v2)
ggplotly(cluster11_1v2, tooltip = "text")
ggplotly(cluster11_4v5, tooltip = "text")


#for loop to subset timepoint avglogfoldchange DEG data
# for (i in 0:19){
#   print("cluster ")
#   print(i)
#   cluster <- subset(pbmc.combined, idents = i)
#   Idents(cluster) <- "Timepoint"
#   
#   cluster_compare_avg <- log1p(AverageExpression(cluster, verbose = FALSE)$RNA)
#   cluster_compare_avg$gene <- rownames(cluster_compare_avg)
#   
#   
#   nam <- paste("cluster", i,".avgLOG_timepoint", sep = "")
#   assign(nam, cluster_compare_avg)
#   
#   rm(cluster)
#   rm(cluster_compare_avg)
#   
# }

fileNames = list(cluster0.avgLOG_timepoint,
                 cluster1.avgLOG_timepoint,
                 cluster2.avgLOG_timepoint,
                 cluster3.avgLOG_timepoint,
                 cluster4.avgLOG_timepoint,
                 cluster5.avgLOG_timepoint,
                 cluster6.avgLOG_timepoint,
                 cluster7.avgLOG_timepoint,
                 cluster8.avgLOG_timepoint,
                 cluster9.avgLOG_timepoint,
                 cluster10.avgLOG_timepoint,
                 cluster11.avgLOG_timepoint,
                 cluster12.avgLOG_timepoint,
                 cluster13.avgLOG_timepoint,
                 cluster14.avgLOG_timepoint,
                 cluster15.avgLOG_timepoint,
                 cluster16.avgLOG_timepoint,
                 cluster17.avgLOG_timepoint,
                 cluster18.avgLOG_timepoint,
                 cluster19.avgLOG_timepoint)


save(plotlists,file="plotlists.RData")

top20_cluster18 <- cluster18.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


