#SCTransform integration test

# seuratObjectList.sct <- list(timepoint_1,timepoint_2,timepoint_3,timepoint_4,timepoint_5,timepoint_6)
# 
# 
# #subset and filter cells based on MT/Counts
# seuratObjectList.sct[[1]] <- subset(seuratObjectList.sct[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 17)
# seuratObjectList.sct[[2]] <- subset(seuratObjectList.sct[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 17)
# seuratObjectList.sct[[3]] <- subset(seuratObjectList.sct[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 17)
# seuratObjectList.sct[[4]] <- subset(seuratObjectList.sct[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 3800 & percent.mt < 17)
# seuratObjectList.sct[[5]] <- subset(seuratObjectList.sct[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 3600 & percent.mt < 19)
# seuratObjectList.sct[[6]] <- subset(seuratObjectList.sct[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 5100 & percent.mt < 15)
# 
# 
# 
# 
# 
# for (i in 1:length(seuratObjectList.sct)) {
#   seuratObjectList.sct[[i]] <- SCTransform(seuratObjectList.sct[[i]],vars.to.regress = "percent.mt", verbose = TRUE)
# }
# 
# sct.features <- SelectIntegrationFeatures(object.list = seuratObjectList.sct, nfeatures = 3000)
# seuratObjectList.sct <- PrepSCTIntegration(object.list = seuratObjectList.sct, anchor.features = sct.features, 
#                                     verbose = TRUE)
# 
# 
# sct.anchors <- FindIntegrationAnchors(object.list = seuratObjectList.sct, normalization.method = "SCT", 
#                                            anchor.features = sct.features, verbose = TRUE)
# 
# 
# 
# sct.integrated <- IntegrateData(anchorset = sct.anchors, normalization.method = "SCT", 
#                                      verbose = TRUE)



DefaultAssay(sct.integrated) <- "integrated"
sct.integrated <- RunPCA(sct.integrated, verbose = FALSE)
sct.integrated <- RunUMAP(sct.integrated, dims = 1:15)




sct.integrated <- FindNeighbors(sct.integrated, dims = 1:15, verbose = FALSE)
sct.integrated <- FindClusters(sct.integrated, verbose = FALSE,resolution = 0.8)

DimPlot(sct.integrated, label=T,split.by="Timepoint")

DefaultAssay(sct.integrated) <- "RNA"
FeaturePlot(sct.integrated,features = c("CD3D","CD8A","CD4","FCGR3A","CD14","CD19"))
FeaturePlot(sct.integrated,features = c("PDCD1"))

sct.integrated.markers <- FindAllMarkers(sct.integrated, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.25)
