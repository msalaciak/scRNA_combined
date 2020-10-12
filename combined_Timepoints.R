library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(RCurl)
library(cowplot)
library(pheatmap)

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
seuratObjectList <- c(timepoint_1,timepoint_2,timepoint_3,timepoint_4,timepoint_5,timepoint_6)
  
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

#combined seurat objects
combinedSeuratObjects <- merge(seuratObjectList[[1]],list(seuratObjectList[[2]],seuratObjectList[[3]],seuratObjectList[[4]],
                                                       seuratObjectList[[5]],seuratObjectList[[6]]), project = "combined_pbmc")
#normalize, find variable genes and scale data
combinedSeuratObjects <- NormalizeData(combinedSeuratObjects, normalization.method = "LogNormalize", scale.factor = 10000)
combinedSeuratObjects <- FindVariableFeatures(combinedSeuratObjects, selection.method = "vst",
                                                dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 2000)
combinedSeuratObjects <- ScaleData(combinedSeuratObjects)

#run PCA

combinedSeuratObjects <- RunPCA(combinedSeuratObjects, features = VariableFeatures(object = combinedSeuratObjects))
print(combinedSeuratObjects[["pca"]], dims = 1:5, nfeatures = 5)
print(VizDimLoadings(combinedSeuratObjects, dims = 1:2, reduction = "pca"))
print(DimPlot(combinedSeuratObjects, reduction = "pca"))
print(DimHeatmap(combinedSeuratObjects, dims = 1, cells = 500, balanced = TRUE))
print(DimHeatmap(combinedSeuratObjects, dims = 1:15, cells = 500, balanced = TRUE))
combinedSeuratObjects <- JackStraw(combinedSeuratObjects, num.replicate = 100)
combinedSeuratObjects <- ScoreJackStraw(combinedSeuratObjects, dims = 1:20)
print(JackStrawPlot(combinedSeuratObjects, dims = 1:20))
print(ElbowPlot(combinedSeuratObjects))



combinedSeuratObjects <- FindNeighbors(combinedSeuratObjects, dims = 1:15)
combinedSeuratObjects <- FindClusters(combinedSeuratObjects, resolution = 1.2, algorithm=4)


combinedSeuratObjects.umap <- RunUMAP(combinedSeuratObjects, dims = 1:15, label=TRUE)
combinedSeuratObjects <- RunTSNE(combinedSeuratObjects, dims = 1:15,label = TRUE)

print(DimPlot(combinedSeuratObjects, reduction = "tsne",pt.size=.75))
print(DimPlot(combinedSeuratObjects.umap, reduction = "umap",pt.size=.75))






