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
p1 <- DimPlot(pbmc.combined, reduction = "umap", group.by = "Timepoint")
p2 <- DimPlot(pbmc.combined, reduction = "umap", label = TRUE)
print(p1)
print(p2)

#integrated count = 26739
Reduce("+",table ( Idents(pbmc.combined) ) )



DimPlot(pbmc.combined,reduction = "umap", split.by = "iRAE")
DimPlot(pbmc.combined,reduction = "umap", split.by = "Disease_status")

timepoint_1.subset <- subset(x = pbmc.combined, subset = Timepoint == "1")
timepoint_2.subset <- subset(x = pbmc.combined, subset = Timepoint == "2")
timepoint_3.subset <- subset(x = pbmc.combined, subset = Timepoint == "3")
timepoint_4.subset <- subset(x = pbmc.combined, subset = Timepoint == "4")
timepoint_5.subset <- subset(x = pbmc.combined, subset = Timepoint == "5")
timepoint_6.subset <- subset(x = pbmc.combined, subset = Timepoint == "6")


DimPlot(timepoint_1.subset,reduction = "umap") +ggtitle("Timepoint 1")
DimPlot(timepoint_2.subset,reduction = "umap") +ggtitle("Timepoint 2")
DimPlot(timepoint_3.subset,reduction = "umap") +ggtitle("Timepoint 3")
DimPlot(timepoint_4.subset,reduction = "umap") +ggtitle("Timepoint 4")
DimPlot(timepoint_5.subset,reduction = "umap") +ggtitle("Timepoint 5")
DimPlot(timepoint_6.subset,reduction = "umap") +ggtitle("Timepoint 6")




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

#CONTROL 
# control.data <- Read10X(data.dir ="/home/matthew/Research/pbmc_controls/filtered_gene_bc_matrices/GRCh38")

# control <-CreateSeuratObject(counts = control.data, project = "Controls", min.cells =3,min.features=200)
# 
# control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")
# print(FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt"))
# print(FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
# 
# control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 5300 & percent.mt < 16)
# 
# control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
# control<- FindVariableFeatures(control, selection.method = "vst",
#                                                 dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)
# control <- ScaleData(control, verbose = TRUE)
# control <- RunPCA(control, npcs = 30, verbose = TRUE)
# 
# control <- FindNeighbors(control, reduction = "pca", dims = 1:20)
# control <- FindClusters(control,resolution = 0.8)
# 
# control <- RunUMAP(control, dims = 1:20)

# Visualization

DimPlot(control, reduction = "umap", label = TRUE)

# saveRDS(pbmc.combined, file = "pbmc_combined.rds")


cluster2_logfc <- cluster2.markers$avg_logFC


cluster2_entrez <- mapIds(org.Hs.eg.db, cluster2.markers$gene, 'ENTREZID', 'SYMBOL')

names(cluster2_logfc) <- cluster2.markers$gene
cluster2_logfc <- cluster2_logfc[!is.na(cluster2_logfc)]
cluster2_logfc = sort(cluster2_logfc, decreasing = TRUE)


gse_cluster2 <- gseGO(
  cluster2_logfc,
  ont = "ALL",
  OrgDb =org.Hs.eg.db,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
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



