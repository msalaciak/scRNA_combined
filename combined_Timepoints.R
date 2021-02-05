library(dplyr)	
library(Seurat)	
library(patchwork)	
library(ggplot2)	
library(Matrix)	
library(tidyverse)
library(plyr)
library(RCurl)	
library(cowplot)
library(clusterProfiler)
library("org.Dm.eg.db",character.only = TRUE)
library(org.Dm.eg.db)
library(DOSE)
library(plotly)
theme_set(theme_cowplot())
library(enrichplot)
library(reshape2)
library(pheatmap)
library(scRepertoire)

#proportion csv
cell_proportion <- read.csv(file = 'pbmc_combined_percent_r.csv')
cell_proportion <-cell_proportion[ , -c(1)]
cell_proportion <-cell_proportion[ -c(21), ]

cell_proportion_per <- read.csv(file = 'pbmc_combined_percent_format.csv')


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

new.cluster.ids <- c("CD14 Mono 1", "CD14 Mono 2", "CD4 TCM", "CD4 Naive","CD14 Mono 3","CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3", "CD14 Mono 4"
                     ,"NK" , "CD16 Mono" ,"B Cell" ,"CD8 TEM","CD14 Mono 5" ,"cpDC","CD14 Mono 6", "CD4 Treg" ,"CD8 Exhausted" ,"CD8 Teff" ,"pDC")

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
cell.ident.pro<-table(Idents(pbmc.combined),pbmc.combined$Timepoint)

cell.ident.pro<-as.data.frame.matrix(cell.ident.pro)
cell.ident.pro <- cbind(cell.ident.pro = rownames(cell.ident.pro), cell.ident.pro)
rownames(cell.ident.pro) <- 1:nrow(cell.ident.pro)

cell.ident.pro <- cell.ident.pro %>%
  rename(
    Identity = cell.ident.pro)



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
  nPerm  = 1000,
  exponent = 1,
  minGSSize = 100,
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


FeaturePlot(object = pbmc.combined, features = c("CD3D","CD3G","CD3E"), cols = c("grey", "red"), reduction = "umap")
FeaturePlot(object = pbmc.combined, features = c("CD4"), cols = c("grey", "red"), reduction = "umap")
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

#sum rows for cell identiy
cell.ident.pro %>% summarize_if(is.numeric, sum, na.rm=TRUE)

#use copy of pbmc.combined to make it seperated by clusters/timepoints
pbmc.combined.timepoint <-pbmc.combined
pbmc.combined.timepoint$celltype.timepoint <- paste(Idents(pbmc.combined.timepoint), pbmc.combined.timepoint$Timepoint, sep = "_")
pbmc.combined.timepoint$celltype <- Idents(pbmc.combined.timepoint)
Idents(pbmc.combined.timepoint) <- "celltype.timepoint"

# test.test <- FindMarkers(pbmc.combined.timepoint, ident.1 = "NK_1", ident.2 = "NK_2", only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)
# 
# test.test <- cbind(gene = rownames(test.test), test.test)
# rownames(test.test) <- 1:nrow(test.test)
# 
# test.test$cluster <- 9
# test.test<- test.test[, c(2,3,4,5,6,1,7)]


#arrays of timepoint cluster names for for loop
cd14mono_1 <-c("CD14 Mono 1_1","CD14 Mono 1_2","CD14 Mono 1_3","CD14 Mono 1_4","CD14 Mono 1_5","CD14 Mono 1_6")
cd14mono_2 <-c("CD14 Mono 2_1","CD14 Mono 2_2","CD14 Mono 2_3","CD14 Mono 2_4","CD14 Mono 2_5","CD14 Mono 2_6")
cd14mono_3 <-c("CD14 Mono 3_1","CD14 Mono 3_2","CD14 Mono 3_3","CD14 Mono 3_4","CD14 Mono 3_5","CD14 Mono 3_6")
cd14mono_4 <-c("CD14 Mono 4_1","CD14 Mono 4_2","CD14 Mono 4_3","CD14 Mono 4_4","CD14 Mono 4_5","CD14 Mono 4_6")
cd14mono_5 <-c("CD14 Mono 5_1","CD14 Mono 5_2","CD14 Mono 5_3","CD14 Mono 5_4","CD14 Mono 5_5","CD14 Mono 5_6")
cd14mono_6 <-c("CD14 Mono 6_1","CD14 Mono 6_2","CD14 Mono 6_3","CD14 Mono 6_4","CD14 Mono 6_5","CD14 Mono 6_6")

cd4_tcm <-c("CD4 TCM_1","CD4 TCM_2","CD4 TCM_3","CD4 TCM_4","CD4 TCM_5","CD4 TCM_6")
cd4_naive <-c("CD4 Naive_1","CD4 Naive_2","CD4 Naive_3","CD4 Naive_4","CD4 Naive_5","CD4 Naive_6")

cd8_1 <- c("CD8 t-cell 1_1","CD8 t-cell 1_2","CD8 t-cell 1_3","CD8 t-cell 1_4","CD8 t-cell 1_5","CD8 t-cell 1_6")
cd8_2 <- c("CD8 t-cell 2_1","CD8 t-cell 2_2","CD8 t-cell 2_3","CD8 t-cell 2_4","CD8 t-cell 2_5","CD8 t-cell 2_6")
cd8_3 <- c("CD8 t-cell 3_1","CD8 t-cell 3_2","CD8 t-cell 3_3","CD8 t-cell 3_4","CD8 t-cell 3_5","CD8 t-cell 3_6")

nk<-c("NK_1","NK_2","NK_3","NK_4","NK_5","NK_6")

cd16_mono<-c("CD16 Mono _1","CD16 Mono _2","CD16 Mono _3","CD16 Mono _4","CD16 Mono _5","CD16 Mono _6")

bcell<-c("B Cell_1","B Cell_2","B Cell_3","B Cell_4","B Cell_5","B Cell_6")

cd8tem <-c("CD8 TEM_1","CD8 TEM_2","CD8 TEM_3","CD8 TEM_4","CD8 TEM_5","CD8 TEM_6")

cpdc<-c("cpDC_1","cpDC_2","cpDC_3","cpDC_4","cpDC_5","cpDC_6")

cd4treg<-c("CD4 Treg_1","CD4 Treg_2","CD4 Treg_3","CD4 Treg_4","CD4 Treg_5","CD4 Treg_6")

cd8ex<-c("CD8 Exhausted_1","CD8 Exhausted_2","CD8 Exhausted_3","CD8 Exhausted_4","CD8 Exhausted_5","CD8 Exhausted_6")

cd8ef<-c("CD8 Teff_1","CD8 Teff_2","CD8 Teff_3","CD8 Teff_4","CD8 Teff_5","CD8 Teff_6")

pdc<-c("pDC_1","pDC_2","pDC_3","pDC_4","pDC_5","pDC_6")
new.cluster.ids <- c("CD14 Mono 1", "CD14 Mono 2", "CD4 TCM", "CD4 Naive","CD14 Mono 3","CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3", "CD14 Mono 4"
                     ,"NK" , "CD16 Mono" ,"B Cell" ,"CD8 TEM","CD14 Mono 5" ,"cpDC","CD14 Mono 6", "CD4 Treg" ,"CD8 Exhausted" ,"CD8 Teff" ,"pDC")

timepoint.all.clusters <- array(c(cd14mono_1,cd14mono_2,
                            cd4_tcm,cd4_naive,cd14mono_3,cd8_1,cd8_2,cd8_3,cd14mono_4,nk,
                           cd16_mono, bcell,cd8tem,cd14mono_5,cpdc,cd14mono_6,cd4treg,cd8ex,cd8ef,pdc),dim=c(6,20,1))

#first number is timepoint, 2nd is cell ident, 3 is just set to 1
timepoint.all.clusters[1,1,1]


for (i in 1:20){
  print("clustering ")
  
  for(j in 1:5){
    timepoint_name_1 <-timepoint.all.clusters[j,i,1]
    timepoint_name_2 <-timepoint.all.clusters[j+1,i,1]
    print("timepoint ")
    print(j)
    print("cells")
    print(timepoint_name_1)
    print(timepoint_name_2)
    
    markers <- FindMarkers(pbmc.combined.timepoint, ident.1 = timepoint_name_1,ident.2=timepoint_name_2,
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)
    
    markers <- cbind(gene = rownames(markers), markers)
    rownames(markers) <- 1:nrow(markers)
    
    markers$cluster <- i
    markers<- markers[, c(2,3,4,5,6,1,7)]
    
    nam <- paste(timepoint_name_1, "  ", timepoint_name_2,".markers", sep = "")
    assign(nam, markers)
    
  }
}





#organize by top 20 avg_logfc
cd_ex_top_10_1v2_up <- head(arrange(eval(as.name("CD8 Exhausted_1  CD8 Exhausted_2.markers")), desc(avg_logFC)), n = 20L)
cd_ex_top_10_1v2_down <- head(arrange(eval(as.name("CD8 Exhausted_1  CD8 Exhausted_2.markers")), avg_logFC), n = 20L)


cdex_cluster <- subset(pbmc.combined, idents = "CD8 Exhausted")
cdex_cluster1 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 1")
cdex_cluster2 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 2")
cdex_cluster3 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 3")
cdex_cluster4 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 4")
cdex_cluster5 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 5")
cdex_cluster6 <-WhichCells(object = cdex_cluster, expression = orig.ident == "pbmc timepoint 6")

DimPlot(pbmc.combined,split.by="Timepoint",cells.highlight = c(cdex_cluster1,cdex_cluster2,cdex_cluster3,cdex_cluster4,cdex_cluster5,cdex_cluster6)  ,
        cols.highlight = "black")

FeaturePlot(object = pbmc.combined, features = c("EOMES"), cols = c("grey", "red"), reduction = "umap", label=TRUE
                      ,cells= WhichCells(object = pbmc.combined, expression = orig.ident == "pbmc timepoint 1"),pt.size = 2,order=TRUE,
            label.size = 3,repel=TRUE) / 
            FeaturePlot(object = pbmc.combined, features = c("EOMES"), cols = c("grey", "red"), reduction = "umap", label=TRUE
                      ,cells= WhichCells(object = pbmc.combined, expression = orig.ident == "pbmc timepoint 2"),pt.size = 2,order=TRUE,
                      label.size = 3,repel=TRUE)



cd_ex_1v2_logfc <- eval(as.name("CD8 Exhausted_1  CD8 Exhausted_2.markers"))$avg_logFC


names(cd_ex_1v2_logfc) <- eval(as.name("CD8 Exhausted_1  CD8 Exhausted_2.markers"))$gene
cd_ex_1v2_logfc <- cd_ex_1v2_logfc[!is.na(cd_ex_1v2_logfc)]
cd_ex_1v2_logfc = sort(cd_ex_1v2_logfc, decreasing = TRUE)


gse_test <- gseGO(
  cd_ex_1v2_logfc,
  ont = "bp",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  nPerm  = 1000,
  minGSSize = 100,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

upsetplot(gse_test)



# cluster2_UNIPROT <- mapIds(org.Hs.eg.db, cluster2.markers$gene, 'UNIPROT', 'SYMBOL')
# cluster2_logfc <- cluster2.markers$avg_logFC
# names(cluster2_logfc) <- cluster2_UNIPROT
# cluster2_logfc <- cluster2_logfc[!is.na(names(cluster2_logfc))]
# cluster2_logfc <-sort(cluster2_logfc, decreasing = TRUE)
# 
# 
# kegg_cluster2 <- gseKEGG(geneList= cluster2_logfc,
#                          keyType = "uniprot",
#                          organism     = "hsa",
#                          minGSSize    = 20,
#                          pvalueCutoff = 0.05,
#                          verbose      = TRUE,
#                          eps=0)


# exhausted cd8 tcell plot of genes
cluster17_test_plot <- cluster17.avgLOG_timepoint[c("DUSP2" ,    "GZMK"  ,    "TUBB"    ,  "HLA-DRB5",  "TRBV5-1" ,  "TRBV11-2" , "TRBV7-9"  , "MTRNR2L8" , "TUBA1B"  ,  "TRAV3"    , "TRAV14DV4",
                         "TRDV3"  ,   "DUT"    ,   "MTRNR2L1" , "CCL5"    ,  "NOSIP"  ,  
                         "GZMK"    ,  "HIST1H3B" , "HIST1H4C" , "TUBB"   ,   "TRBV7-2" ,  "TRBV11-2" , "TRBV7-9" ,  "TRBV19" ,   "TUBA1B"  ,  "TRAV14DV4", "CCL4" ,    
                         "GZMK"    ,"TRBV6-1", "TRBV30" , "TRDV1"  ,
                         "JUN"      ,"GNLY"  ,   "DUSP1" ,   "TRBV3-1" , "TRBV7-9" , "TRBV20-1", "MAP3K8" ,  "KLRF1"  ,  "NFKBIA"  , "FOS"  ,    "CCL4"   ,  "CST3" ,   
                         "IER2"     ,"ZFP36" ,  
                         "GNLY"     ,    "TRBV7-9"  ,    "RP11-290D2.6", "IFI27"  ,      "CST3"      ,   "IGLV2-14"  ,   "IGLL5"   ),]
cluster17_test_plot<-cluster17_test_plot[!duplicated(cluster17_test_plot), ]

pheatmap(cluster17_test_plot[,1:6],cluster_cols=FALSE)

cluster17_test_plot <-melt(cluster17_test_plot, id.vars=c("gene"))
cluster17_test_plot$timepoint <- as.numeric(as.character(cluster17_test_plot$variable))

ggplot(cluster17_test_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none") + ggtitle("CD8 Exhausted Across 6 Time ")




ggplot(data = cluster17_test_plot, mapping = aes(x = timepoint,
                                                       y = reorder(gene,value),
                                                       fill = value)) +
  geom_tile() +
  xlab(label = "Timepoints") + theme_gray()   + scale_fill_distiller(palette = "Spectral")




##cluster 18 t-eff

cluster18_test_plot <- cluster18.avgLOG_timepoint[c("DENND2D"    ,   "FAM212B"    ,   "TXNIP"    ,     "S100A9"   ,     "AQP10"  ,       "PYCR2"     ,    "ACP1"     ,     "UBXN2A"     ,  
                                                    "ZFP36L2"    ,   "GNLY"    ,      "PDK1"      ,    "ATP5G3"  ,      "SSFA2"    ,     "UBE2F"     ,    "SEPT2"    ,     "CRBN"       ,  
                                                    "VGLL4"      ,  "TMEM40"   ,     "CMTM6"      ,   "FAM212A"   ,    "MITF"      ,    "MTRNR2L12"  ,   "TACC3"      ,   "LINC00989"    ,
                                                    "RP11-473L15.3", "NDUFS4"  ,      "CSNK1A1"   ,    "HIST1H2BN"  ,   "SRSF3"    ,     "PTCRA"     ,    "CNPY3"    ,     "GTPBP2"     ,  
                                                    "SNX3"      ,    "MAN1A1"  ,      "TBPL1"     ,    "MTPN"    ,      "TRBV11-2"  ,    "INSIG1"   ,     "DNAJB6"    ,    "SEPT6"       , 
                                                    "ANK1"     ,     "CA2"    ,       "TMEM64"    ,    "VPS28"   ,      "ZFAND5"    ,    "SYK"      ,     "MTRNR2L8"  ,    "RARRES3"      ,
                                                    "MALAT1"   ,     "CD3D"   ,       "ARHGAP21"  ,    "BAMBI"   ,      "KLRB1"     ,    "CLEC1B"    ,    "TUBA1B"    ,    "NFE2"        , 
                                                    "RHOF"     ,     "USP12"  ,       "TPP2"    ,      "ANKRD10"  ,     "TRDC"     ,     "FAM177A1"  ,    "MAP4K5"  ,      "MAX"         , 
                                                    "ETFA"    ,      "SMG1"     ,     "LAT"     ,      "UBE2G1"  ,      "RPL26"    ,     "MTRNR2L1"  ,    "ICAM2"    ,     "SEC14L1"      ,
                                                    "FAM110A" ,      "TBXA2R"    ,    "GMFG"    ,      "NUCB1"   ,      "ATF4"     ,     "DYRK1A"   ,    
                                                    "CMPK1"   ,  "RNF11"  ,   "JUN"   ,    "S100A8" ,   "GNLY"    ,  "CXCR4"   ,  "RPL29" ,    "SNCA"   ,   "MMRN1"   ,  "DUSP1"     ,"HIST1H3H" ,
                                                    "PRKAR1B" ,  "TSPAN13" ,  "PTPN12" ,   "PRKAR2B" ,  "LINC-PINT" ,"MTRNR2L10", "RPL7A"  ,   "GZMH"  ,    "RPS2"    ,  "RPL23A"  ,  "FKBP1A",   
                                                    "MYL9"   ,   "YIF1B"    , "NKG7"    ,  "SF3A1"  ,   "MT-ND2"  ,  "IGKV1-5"  , "IGHV5-51" ,
                                                    "S100A9"  ,  "S100A12"  , "S100A8"  ,  "GNLY" ,     "MTRNR2L12", "SNCA"    ,  "VCAN"    ,  "CD74"  ,    "TBPL1"   ,  "MTRNR2L8" , "LYZ"   ,   
                                                    "IFI27"   ,  "MTRNR2L1" , "JUNB"    ,  "IGKV1-5" ,  "IGHV5-51" ,
                                                    "NEXN"   ,   "S100A9"   , "S100A8" ,   "FCER1G" ,   "MPC2"   ,   "GNLY"  ,    "MTRNR2L12" ,"FGFBP2" ,   "SNCA"  ,    "GZMA"   ,   "TIMP1" ,   
                                                    "FAM127A"  , "IFITM3"   , "MTRNR2L8" , "ANXA7"  ,   "TRAV8-4" , "GZMH" ,     "FOS" ,      "MTRNR2L1" , "ABCC3"   ,  "NKG7"    ,  "PACSIN2" , 
                                                    "FHL3"    ,      "NEXN"      ,    "IFI44L"   ,     "FAM212B"  ,     "KIFAP3"    ,    "ADIPOR1" ,      "CNST"   ,       "TRIM58"   ,    
                                                    "CMPK2"  ,       "MEIS1"   ,      "PDK1"     ,     "C2orf88" ,      "FRMD4B"    ,    "MTRNR2L12" ,    "PPBP"   ,       "SEPT11"   ,    
                                                    "CASP3"   ,      "RP11-367G6.3" , "HIST1H2BJ"  ,   "PTCRA"    ,     "TBPL1"      ,   "AHR"      ,     "NT5C3A" ,       "GATA1"     ,   
                                                    "MTRNR2L10" ,    "CA2"     ,      "LCN2"     ,     "EGFL7"   ,      "MTRNR2L8"   ,   "ESAM"    ,      "RBM17"    ,     "NCOA4"      ,  
                                                    "NFE2"     ,     "C12orf76"   ,   "ELF1"      ,    "IFI27"   ,      "GTF2A2"     ,   "CD68"   ,       "MTRNR2L1"  ,    "HEXIM1"       ,
                                                    "MMD"       ,    "ICAM2"     ,    "PSTPIP2"    ,   "FAM110A"  ,     "SMOX"       ,   "TOP1"    ,      "TBXA2R"     ,   "RP11-678G14.3",
                                                    "PTGIR"       ,  "PARVB"    ,     "MT-ND2"    ,    "MT-CO1"   ,     "MT-ATP6"     ,  "MT-CO3"    ),]
cluster18_test_plot<-cluster18_test_plot[!duplicated(cluster18_test_plot), ]
pheatmap(cluster18_test_plot[,1:6],cluster_cols=FALSE)

cluster18_test_plot <-melt(cluster18_test_plot, id.vars=c("gene"))
cluster18_test_plot$timepoint <- as.numeric(as.character(cluster18_test_plot$variable))

ggplot(cluster18_test_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none") + ggtitle("CD8 Teff ")

# b cell plots

cluster11_test_plot <- cluster11.avgLOG_timepoint[c("IGKV4-1"  , "IGLV4-69" , "IGLV5-52" , "IGLV2-23" , "IGKV1-5"  , "IGKV3D-15" ,"IGHV3-7"  , "IGHV1-18",  "IGHV3-21" , "IGHV1-24"  ,"IGHV3-30" ,
                                                    "IGHV3-33" , "IGLV1-51" , "IGLV1-47" , "IGLV1-44" , "IGLV3-21",  "IGLV2-11" ,
                                                    "TSC22D3" , "FOS"    ,  "IGHV4-39", "IGLV2-8" , "IGHV1-18" ,"IGKV1-27",
                                                    "JUN"     , "IGKV4-1"  ,"IGKV3-20", "CD83"    , "FOS"     , "IGHV4-39" ,"IL32"   ,  "JUNB"   ,  "IER2"  ,   "IGLV2-8" , "IGHV1-18", "IGLV1-51",
                                                    "IGKV1-27",
                                                    "TXNIP"  ,  "IGKV3-20" ,"FOS"   ,   "PPP1R15A" ,"IGHV3-53" ,"IGLV3-21", "IGKV1-27",
                                                    "JUN"    ,  "GNLY"   ,  "IGKV3-11", "DUSP2"   , "CD69"  ,   "NFKBIA"  , "FOS"   ,   "IER2"  ,   "IGLV2-23" ,"IGHV1-24" ,"IGLV5-45" ,"IGKV1-27"
),]

cluster11.subset.markers<-subset(cluster11.markers, gene %in% c("IGKV4-1"  , "IGLV4-69" , "IGLV5-52" , "IGLV2-23" , "IGKV1-5"  , "IGKV3D-15" ,"IGHV3-7"  , "IGHV1-18",  "IGHV3-21" , "IGHV1-24"  ,"IGHV3-30" ,
                            "IGHV3-33" , "IGLV1-51" , "IGLV1-47" , "IGLV1-44" , "IGLV3-21",  "IGLV2-11" ,
                            "TSC22D3" , "FOS"    ,  "IGHV4-39", "IGLV2-8" , "IGHV1-18" ,"IGKV1-27",
                            "JUN"     , "IGKV4-1"  ,"IGKV3-20", "CD83"    , "FOS"     , "IGHV4-39" ,"IL32"   ,  "JUNB"   ,  "IER2"  ,   "IGLV2-8" , "IGHV1-18", "IGLV1-51",
                            "IGKV1-27",
                            "TXNIP"  ,  "IGKV3-20" ,"FOS"   ,   "PPP1R15A" ,"IGHV3-53" ,"IGLV3-21", "IGKV1-27",
                            "JUN"    ,  "GNLY"   ,  "IGKV3-11", "DUSP2"   , "CD69"  ,   "NFKBIA"  , "FOS"   ,   "IER2"  ,   "IGLV2-23" ,"IGHV1-24" ,"IGLV5-45" ,"IGKV1-27"))

cluster11.subset.markers %>% filter(pct.1 > .1)

cluster11_test_plot <- cluster11.avgLOG_timepoint[c("CD83","CD69","JUN","TSC22D3","DUSP2","FOS"),]


cluster11_test_plot<-cluster11_test_plot[!duplicated(cluster11_test_plot), ]

cluster11_test_plot <-melt(cluster11_test_plot, id.vars=c("gene"))
cluster11_test_plot$timepoint <- as.numeric(as.character(cluster11_test_plot$variable))

cluster11_test_plot <- ggplot(cluster11_test_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none",plot.title = element_text(face = "plain"),) + ggtitle("B Cell Across 6 Timepoints")

ggsave("cluster11_test_plot.png",width = 16 ,height = 5,dpi = 300)

# NK plots

cluster9_test_plot <- cluster9.avgLOG_timepoint[c("S100A8" ,
                                                    "JUN"  ,   "DUSP1" ,  "TSC22D3", "FOS"    ,
                                                    "JUN"   ,   "S100A8" ,  "NFKBIA" ,  "FOS"  ,    "PPP1R15A",
                                                    "JUN"   ,   "NR4A2"  ,  "DUSP1"  ,  "TNFAIP3" , "CD69"  ,   "NFKBIA"  , "FOS"  ,    "IER2" ,    "JUND"  ,   "PPP1R15A",
                                                    "IFI44L"
),]
cluster9_test_plot<-cluster9_test_plot[!duplicated(cluster9_test_plot), ]

cluster9_test_plot <-melt(cluster9_test_plot, id.vars=c("gene"))
cluster9_test_plot$timepoint <- as.numeric(as.character(cluster9_test_plot$variable))

ggplot(cluster9_test_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none") + ggtitle("NK")


#cd8 exhaust markers only
cluster17_ex_mark_plot <- cluster17.avgLOG_timepoint[c("PDCD1","EOMES","LAG3","HAVCR2","IFNG","CD244"),]
cluster17_ex_mark_plot<-cluster17_ex_mark_plot[!duplicated(cluster17_ex_mark_plot), ]

cluster17_ex_mark_plot[,7] <- c("PDCD1","EOMES","LAG3","HAVCR2 (TIM3)","IFNG","CD244 (2B4)")

cluster17_ex_mark_plot <-melt(cluster17_ex_mark_plot, id.vars=c("gene"))
cluster17_ex_mark_plot$timepoint <- as.numeric(as.character(cluster17_ex_mark_plot$variable))

cd8explot <- ggplot(cluster17_ex_mark_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line()  +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none",plot.title = element_text(face = "plain")) + ggtitle("CD8 Exhaustion Markers")

ggsave("cd8explot.png",width = 16 ,height = 5,dpi = 300)

#cd8 teff markers

cluster18_ex_mark_plot <- cluster18.avgLOG_timepoint[c("IL7R","CD44","SELL","GZMB","IFNG","PRF1"),]
cluster18_ex_mark_plot<-cluster18_ex_mark_plot[!duplicated(cluster18_ex_mark_plot), ]

cluster18_ex_mark_plot <-melt(cluster18_ex_mark_plot, id.vars=c("gene"))
cluster18_ex_mark_plot$timepoint <- as.numeric(as.character(cluster18_ex_mark_plot$variable))

ggplot(cluster18_ex_mark_plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none") + ggtitle("CD8 Teff ")






#proportion graph
cell_proportion<-melt(cell.ident.pro, id.vars=c("Identity"))


prop.ex <- ggplot(cell_proportion,aes(x=value, y=Identity)) + 
  geom_bar(aes(fill = as.double(variable)),colour="black", stat="identity" , position = "fill") +
  scale_fill_continuous(type = "viridis", name = "Time points") +   labs(y='Cell Identity', x='% of Cluster Relative to Sample')




ggsave("prop.ex.png",width = 8, height = 5,dpi = 300)


# geneNames <-data.frame(geneNames)
testtest<-filter(cluster17.markers, gene == c("DUSP2" ,    "GZMK"  ,    "TUBB"    ,  "HLA-DRB5",  "TRBV5-1" ,  "TRBV11-2" , "TRBV7-9"  , "MTRNR2L8" , "TUBA1B"  ,  "TRAV3"    , "TRAV14DV4",
                                                 "TRDV3"  ,   "DUT"    ,   "MTRNR2L1" , "CCL5"    ,  "NOSIP"  ,  
                                                 "GZMK"    ,  "HIST1H3B" , "HIST1H4C" , "TUBB"   ,   "TRBV7-2" ,  "TRBV11-2" , "TRBV7-9" ,  "TRBV19" ,   "TUBA1B"  ,  "TRAV14DV4", "CCL4" ,    
                                                 "GZMK"    ,"TRBV6-1", "TRBV30" , "TRDV1"  ,
                                                 "JUN"      ,"GNLY"  ,   "DUSP1" ,   "TRBV3-1" , "TRBV7-9" , "TRBV20-1", "MAP3K8" ,  "KLRF1"  ,  "NFKBIA"  , "FOS"  ,    "CCL4"   ,  "CST3" ,   
                                                 "IER2"     ,"ZFP36" ,  
                                                 "GNLY"     ,    "TRBV7-9"  ,    "RP11-290D2.6", "IFI27"  ,      "CST3"      ,   "IGLV2-14"  ,   "IGLL5"   ))

testtest<-subset(cluster17.markers, gene %in% c("DUSP2" ,    "GZMK"  ,    "TUBB"    ,  "HLA-DRB5",  "TRBV5-1" ,  "TRBV11-2" , "TRBV7-9"  , "MTRNR2L8" , "TUBA1B"  ,  "TRAV3"    , "TRAV14DV4",
                                                "TRDV3"  ,   "DUT"    ,   "MTRNR2L1" , "CCL5"    ,  "NOSIP"  ,  
                                                "GZMK"    ,  "HIST1H3B" , "HIST1H4C" , "TUBB"   ,   "TRBV7-2" ,  "TRBV11-2" , "TRBV7-9" ,  "TRBV19" ,   "TUBA1B"  ,  "TRAV14DV4", "CCL4" ,    
                                                "GZMK"    ,"TRBV6-1", "TRBV30" , "TRDV1"  ,
                                                "JUN"      ,"GNLY"  ,   "DUSP1" ,   "TRBV3-1" , "TRBV7-9" , "TRBV20-1", "MAP3K8" ,  "KLRF1"  ,  "NFKBIA"  , "FOS"  ,    "CCL4"   ,  "CST3" ,   
                                                "IER2"     ,"ZFP36" ,  
                                                "GNLY"     ,    "TRBV7-9"  ,    "RP11-290D2.6", "IFI27"  ,      "CST3"      ,   "IGLV2-14"  ,   "IGLL5"   ))





## DE for pembro-naive relapse (in timepoint 1) and pembro-treated relapse (time point 4)

bcell_1v4 <- FindMarkers(pbmc.combined, ident.1 = "1", ident.2="4",group.by="Timepoint",
                            only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "B Cell")
bcell_1v4 <- cbind(gene = rownames(bcell_1v4), bcell_1v4)
rownames(bcell_1v4) <- 1:nrow(bcell_1v4)
bcell_1v4<- bcell_1v4[, c(2,3,4,5,6,1)]


tcellex_1v4 <- FindMarkers(pbmc.combined, ident.1 = "1", ident.2="4",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Exhausted")
tcellex_1v4 <- cbind(gene = rownames(tcellex_1v4), tcellex_1v4)
rownames(tcellex_1v4) <- 1:nrow(tcellex_1v4)
tcellex_1v4<- tcellex_1v4[, c(2,3,4,5,6,1)]


tcellef_1v4 <- FindMarkers(pbmc.combined, ident.1 = "1", ident.2="4",group.by="Timepoint",
                             only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Teff")
tcellef_1v4 <- cbind(gene = rownames(tcellef_1v4), tcellef_1v4)
rownames(tcellef_1v4) <- 1:nrow(tcellef_1v4)
tcellef_1v4<- tcellef_1v4[, c(2,3,4,5,6,1)]

nk_1v4 <- FindMarkers(pbmc.combined, ident.1 = "1", ident.2="4",group.by="Timepoint",
                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "NK")
nk_1v4 <- cbind(gene = rownames(nk_1v4), nk_1v4)
rownames(nk_1v4) <- 1:nrow(nk_1v4)
nk_1v4<- nk_1v4[, c(2,3,4,5,6,1)]

#compare 2 and 3 in remission and compare 5 and 6 in presence of disease, and IRAE 2 and 5

bcell_2v3 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="3",group.by="Timepoint",
                         only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "B Cell")
bcell_2v3 <- cbind(gene = rownames(bcell_2v3), bcell_2v3)
rownames(bcell_2v3) <- 1:nrow(bcell_2v3)
bcell_2v3<- bcell_2v3[, c(2,3,4,5,6,1)]


tcellex_2v3 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="3",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Exhausted")
tcellex_2v3 <- cbind(gene = rownames(tcellex_2v3), tcellex_2v3)
rownames(tcellex_2v3) <- 1:nrow(tcellex_2v3)
tcellex_2v3<- tcellex_2v3[, c(2,3,4,5,6,1)]


tcellef_2v3 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="3",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Teff")
tcellef_2v3 <- cbind(gene = rownames(tcellef_2v3), tcellef_2v3)
rownames(tcellef_2v3) <- 1:nrow(tcellef_2v3)
tcellef_2v3<- tcellef_2v3[, c(2,3,4,5,6,1)]

nk_2v3 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="3",group.by="Timepoint",
                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "NK")
nk_2v3 <- cbind(gene = rownames(nk_2v3), nk_2v3)
rownames(nk_2v3) <- 1:nrow(nk_2v3)
nk_2v3<- nk_2v3[, c(2,3,4,5,6,1)]

# 5 and 6

bcell_5v6 <- FindMarkers(pbmc.combined, ident.1 = "5", ident.2="6",group.by="Timepoint",
                         only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "B Cell")
bcell_5v6 <- cbind(gene = rownames(bcell_5v6), bcell_5v6)
rownames(bcell_5v6) <- 1:nrow(bcell_5v6)
bcell_5v6<- bcell_5v6[, c(2,3,4,5,6,1)]


tcellex_5v6 <- FindMarkers(pbmc.combined, ident.1 = "5", ident.2="6",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Exhausted")
tcellex_5v6 <- cbind(gene = rownames(tcellex_5v6), tcellex_5v6)
rownames(tcellex_5v6) <- 1:nrow(tcellex_5v6)
tcellex_5v6<- tcellex_5v6[, c(2,3,4,5,6,1)]


tcellef_5v6 <- FindMarkers(pbmc.combined, ident.1 = "5", ident.2="6",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Teff")
tcellef_5v6 <- cbind(gene = rownames(tcellef_5v6), tcellef_5v6)
rownames(tcellef_5v6) <- 1:nrow(tcellef_5v6)
tcellef_5v6<- tcellef_5v6[, c(2,3,4,5,6,1)]

nk_5v6 <- FindMarkers(pbmc.combined, ident.1 = "5", ident.2="6",group.by="Timepoint",
                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "NK")
nk_5v6 <- cbind(gene = rownames(nk_5v6), nk_5v6)
rownames(nk_5v6) <- 1:nrow(nk_5v6)
nk_5v6<- nk_5v6[, c(2,3,4,5,6,1)]


# 2 and 5

bcell_2v5 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="5",group.by="Timepoint",
                         only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "B Cell")
bcell_2v5 <- cbind(gene = rownames(bcell_2v5), bcell_2v5)
rownames(bcell_2v5) <- 1:nrow(bcell_2v5)
bcell_2v5<- bcell_2v5[, c(2,3,4,5,6,1)]


tcellex_2v5 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="5",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Exhausted")
tcellex_2v5 <- cbind(gene = rownames(tcellex_2v5), tcellex_2v5)
rownames(tcellex_2v5) <- 1:nrow(tcellex_2v5)
tcellex_2v5<- tcellex_2v5[, c(2,3,4,5,6,1)]


tcellef_2v5 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="5",group.by="Timepoint",
                           only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 Teff")
tcellef_2v5 <- cbind(gene = rownames(tcellef_2v5), tcellef_2v5)
rownames(tcellef_2v5) <- 1:nrow(tcellef_2v5)
tcellef_2v5<- tcellef_2v5[, c(2,3,4,5,6,1)]

nk_2v5 <- FindMarkers(pbmc.combined, ident.1 = "2", ident.2="5",group.by="Timepoint",
                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "NK")
nk_2v5 <- cbind(gene = rownames(nk_2v5), nk_2v5)
rownames(nk_2v5) <- 1:nrow(nk_2v5)
nk_2v5<- nk_2v5[, c(2,3,4,5,6,1)]




#organize by top/bottom 10 avg_logfc for timepoint 1v4
 tcellef_1v4_top10<- head(arrange(tcellef_1v4, desc(avg_logFC)), n = 20L)
 tcellef_1v4_bottom10<- head(arrange(tcellef_1v4, avg_logFC), n = 20L)
 
 tcellex_1v4_top10<- head(arrange(tcellex_1v4, desc(avg_logFC)), n = 20L)
 tcellex_1v4_bottom10<- head(arrange(tcellex_1v4, avg_logFC), n = 20L)
 
 #organize by top/bottom 10 avg_logfc for timepoint 2v3
 
 tcellef_2v3_top10<- head(arrange(tcellef_2v3, desc(avg_logFC)), n = 20L)
 tcellef_2v3_bottom10<- head(arrange(tcellef_2v3, avg_logFC), n = 20L)
 
 tcellex_2v3_top10<- head(arrange(tcellex_2v3, desc(avg_logFC)), n = 20L)
 tcellex_2v3_bottom10<- head(arrange(tcellex_2v3, avg_logFC), n = 20L)

 #organize by top/bottom 10 avg_logfc for timepoint 2v3
 
 tcellef_5v6_top10<- head(arrange(tcellef_5v6, desc(avg_logFC)), n = 20L)
 tcellef_5v6_bottom10<- head(arrange(tcellef_5v6, avg_logFC), n = 20L)
 
 tcellex_5v6_top10<- head(arrange(tcellex_5v6, desc(avg_logFC)), n = 20L)
 tcellex_5v6_bottom10<- head(arrange(tcellex_5v6, avg_logFC), n = 20L)
 
 # tox 
 View(filter(cluster12.avgLOG_timepoint, gene %in% c("TCF7","TIGIT","TOX","PDCD1","LAG3","HAVCR2","CD244","EOMES","IFNG","KLRG1")))
 View(filter(cluster17.avgLOG_timepoint, gene %in% c("TCF7","TIGIT","TOX","PDCD1","LAG3","HAVCR2","CD244","EOMES","IFNG","KLRG1")))
 View(filter(cluster18.avgLOG_timepoint, gene %in% c("TCF7","TIGIT","TOX","PDCD1","LAG3","HAVCR2","CD244","EOMES","IFNG","KLRG1")))
 
 
 
 # to test and see if the findMarkers function works
 cd8.tcell1.1v4 <- FindMarkers(pbmc.combined, ident.1 = "1", ident.2="4",group.by="Timepoint",
                               only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "CD8 t-cell 1")
 cd8.tcell1.1v4 <- cbind(gene = rownames(cd8.tcell1.1v4), cd8.tcell1.1v4)
 rownames(cd8.tcell1.1v4) <- 1:nrow(cd8.tcell1.1v4)
 cd8.tcell1.1v4<- cd8.tcell1.1v4[, c(2,3,4,5,6,1)]
 
 
 
FeaturePlot(object = pbmc.combined, features = c("CD68","FCGR2A","CSF1R"), cols = c("grey", "red"), reduction = "umap")


#cluster 17 line plots

tox.cluster17.plot <- filter(cluster17.avgLOG_timepoint, gene %in% c("NFATC1","TOX","PDCD1"))
tox.cluster17.plot <-melt(tox.cluster17.plot, id.vars=c("gene"))
tox.cluster17.plot$timepoint <- as.numeric(as.character(tox.cluster17.plot$variable))

ggplot(tox.cluster17.plot, aes(x = timepoint, y = value, color = gene, group = gene)) + 
  geom_point()+geom_line() +facet_wrap(vars(gene),scales = "free_y") +theme(legend.position="none",plot.title = element_text(face = "plain"),) + ggtitle("CD8 T ex TOX/TCF7/PDCD1 ")

# heatmap between cluster 17/18

cluster17.18.top30 <-pbmc.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)


DefaultAssay(pbmc.combined) <- "RNA"
DoHeatmap(subset(pbmc.combined,idents = c("CD8 Exhausted","CD8 TEM"),downsample=100), feature=filter(cluster17.18.top30, cluster == 17| cluster == 12)$gene)


DoHeatmap(pbmc.combined, feature=filter(cluster17.18.top30, cluster == 17| cluster == 12)$gene)



