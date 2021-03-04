

# cell cycle scoring



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cd8.subset.cc <- subset(x = pbmc.combined.cellrep, idents = c("CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3","CD8 TEM","CD8 Exhausted" ,"CD8 Teff"))


cd8.subset.list <- SplitObject(cd8.subset.cc, split.by = "Timepoint")

cd8.subset.list <- lapply(X = cd8.subset.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)
})



# 
#integrate data
cd8.subset.anchors <- FindIntegrationAnchors(object.list = cd8.subset.list, dims = 1:30)
cd8.subset.cc <- IntegrateData(anchorset = cd8.subset.anchors, dims = 1:30)

rm(cd8.subset.list)
rm(cd8.subset.anchors)

DefaultAssay(cd8.subset.cc) <- "integrated"

cd8.subset.cc <- ScaleData(cd8.subset.cc, verbose = TRUE)

cd8.subset.cc <- CellCycleScoring(cd8.subset.cc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(cd8.subset.cc, features = c("PCNA", "TOP2A", "MCM6", "MKI67"))

cd8.subset.cc <- RunPCA(cd8.subset.cc, features = c(s.genes, g2m.genes))
DimPlot(cd8.subset.cc)


cd8.subset.cc$CC.Difference <- cd8.subset.cc$S.Score - cd8.subset.cc$G2M.Score
cd8.subset.cc <- ScaleData(cd8.subset.cc, vars.to.regress = "CC.Difference", features = rownames(cd8.subset.cc))

cd8.subset.cc <- ScaleData(cd8.subset.cc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cd8.subset.cc))



cd8.subset.cc <- RunPCA(cd8.subset.cc, npcs = 30, verbose = TRUE)

cd8.subset.cc <- FindNeighbors(cd8.subset.cc, reduction = "pca", dims = 1:12)
cd8.subset.cc <- FindClusters(cd8.subset.cc,resolution = 0.8)

cd8.subset.cc <- RunUMAP(cd8.subset.cc, dims = 1:12)

# look at the UMAP dimplot now
DimPlot(cd8.subset.cc, label= TRUE)
DimPlot(cd8.subset.cc, label= TRUE,group.by="old.ident")

DefaultAssay(cd8.subset.cc) <- "RNA"
cd8.subset.cc <- NormalizeData(cd8.subset.cc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cd8.subset.cc)
cd8.subset.cc <- ScaleData(cd8.subset.cc, features = all.genes)

cd8.subset.cc.full.markers.test <- FindAllMarkers(cd8.subset.cc, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.00)

cd8.subset.cc.full.cluster0v6.4<- FindMarkers(cd8.subset.cc, ident.1 = "0",ident.2= c("6","4"),
                                      only.pos = FALSE, min.pct = 0.10,logfc.threshold = 0.05)

cd8.subset.cc.full.cluster0v6.4 <- cbind(gene = rownames(cd8.subset.cc.full.cluster0v6.4), cd8.subset.cc.full.cluster0v6.4)
cd8.subset.cc.full.cluster0v6.4<- cd8.subset.cc.full.cluster0v6.4[, c(2,3,4,5,6,1)]


## quick code to double check pct.1/pct.2
length(WhichCells(object = cd8.subset.cc, idents = "3"))
length(WhichCells(object = cd8.subset.cc, idents = "3", expression = JUN > 0))

new.cluster.ids <- c("CD8 TEM-like", "CD8 TCM-like", "CD8 TEM-1", "CD8 TEFF", "CD8 TEM-2", "CD8 NAIVE", 
                     "CD8 TEM-3", "CD8 TEM-4", "CD8 TEX-1","CD8 TEX-2","CD8 Unknown-1","CD8 Unknown-2")
names(new.cluster.ids) <- levels(cd8.subset.cc)
cd8.subset.cc <- RenameIdents(cd8.subset.cc, new.cluster.ids)

cd8.subset.cc$new.idents <- Idents(cd8.subset.cc)


