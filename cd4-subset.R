s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cd4.subset.cc <- subset(x = pbmc.combined.cellrep, idents = c("CD4 TCM","CD4 Treg","CD4 Naive"))


cd4.subset.list <- SplitObject(cd4.subset.cc, split.by = "Timepoint")

cd4.subset.list <- lapply(X = cd4.subset.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)
})


# 
#integrate data
cd4.subset.anchors <- FindIntegrationAnchors(object.list = cd4.subset.list, dims = 1:30)
cd4.subset.cc <- IntegrateData(anchorset = cd4.subset.anchors, dims = 1:30)

rm(cd4.subset.list)
rm(cd4.subset.anchors)

DefaultAssay(cd4.subset.cc) <- "integrated"

cd4.subset.cc <- ScaleData(cd4.subset.cc, verbose = TRUE)

# cd4.subset.cc <- CellCycleScoring(cd4.subset.cc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# RidgePlot(cd4.subset.cc, features = c("PCNA", "TOP2A", "MCM6", "MKI67"))
# 
# cd4.subset.cc <- RunPCA(cd4.subset.cc, features = c(s.genes, g2m.genes))
# DimPlot(cd4.subset.cc)
# 
# 
# cd4.subset.cc$CC.Difference <- cd4.subset.cc$S.Score - cd4.subset.cc$G2M.Score
# cd4.subset.cc <- ScaleData(cd4.subset.cc, vars.to.regress = "CC.Difference", features = rownames(cd4.subset.cc))




cd4.subset.cc <- RunPCA(cd4.subset.cc, npcs = 30, verbose = TRUE)

cd4.subset.cc <- FindNeighbors(cd4.subset.cc, reduction = "pca", dims = 1:12)
cd4.subset.cc <- FindClusters(cd4.subset.cc,resolution = 0.5)

cd4.subset.cc <- RunUMAP(cd4.subset.cc, dims = 1:12)

# look at the UMAP dimplot now
DimPlot(cd4.subset.cc, label= TRUE)
DimPlot(cd4.subset.cc, label= TRUE,group.by="old.ident")

DefaultAssay(cd4.subset.cc) <- "RNA"
cd4.subset.cc <- NormalizeData(cd4.subset.cc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cd4.subset.cc)
cd4.subset.cc <- ScaleData(cd4.subset.cc, features = all.genes)


cd4.subset.cc.full.markers <- FindAllMarkers(cd4.subset.cc, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.05)



plot_density(pbmc.combined.cellrep, features= c("CD8A","CD4"), joint=TRUE)
FeaturePlot(cd4.subset.cc, c("CD8A"))
