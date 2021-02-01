# subcluster cd8+ t-cells

cd8.subset <- subset(x = pbmc.combined, idents = c("CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3","CD8 TEM","CD8 Exhausted" ,"CD8 Teff"))

cd8.subset <- subset(x = pbmc.combined, idents = c("CD8 Exhausted" ,"CD8 Teff"))

cd8.subset <- subset(x = pbmc.combined, idents = c("CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3","CD8 Exhausted"))

cd8.subset <- FindVariableFeatures(cd8.subset, selection.method = "vst",
                                   dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)


# scale the data
all.genes <- rownames(cd8.subset)
cd8.subset <- ScaleData(cd8.subset, features = all.genes)




# perform PCA
cd8.subset <- RunPCA(cd8.subset, features = VariableFeatures(object = cd8.subset))



# now we do a elbow plot to determine how many components we should include in our analysis 
# Where the bend (elbow) occurs gives us a good idea of which components hold the most variability (the important data!)
ElbowPlot(cd8.subset)

# we choose the dims based on the elbow plot, the resolution paramter sometimes requires tuning to give the best clustering results
# this part is always tricky and sometimes you may need to change the dims a few times, but from looking at the elbow plot lets go with 12 dims


# default algorithm is louvain (which is what we were using in the previous analysis!)
# also by looking at their results in the params section of the reduction table within the seurat object, they used a resolution of 1
# we will do the same and see how it looks
cd8.subset <- FindNeighbors(cd8.subset, dims = 1:10)
cd8.subset <- FindClusters(cd8.subset, resolution = 1)

# now we can run a umap and tsne, keep the dims parameter the same value from FindNeighbours
cd8.subset <- RunUMAP(cd8.subset, dims = 1:10)



# look at the UMAP dimplot now
DimPlot(cd8.subset, label= TRUE)
DimPlot(cd8.subset, label= TRUE,split.by="old.ident")
DimPlot(cd8.subset, label= TRUE,group.by="old.ident")
DimPlot(cd8.subset, label= TRUE,group.by="old.ident",split.by="Timepoint")

FeaturePlot(object = cd8.subset, features = c("PDCD1"), cols = c("grey", "blue"))



table(Idents(cd8.subset), cd8.subset$Timepoint)

prop.table(table(Idents(cd8.subset), cd8.subset$Timepoint), margin = 2)



cd8.subset.markers <- FindAllMarkers(cd8.subset, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.05)


cd8.subset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



# cd8.subset.markers10v9 <- FindMarkers(cd8.subset, ident.1 = "10", ident.2="9",group.by="Timepoint",
#                               only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05,subset.ident = "4")

cd8.subset.markers9v10.mast <- FindMarkers(cd8.subset, ident.1 = "9", ident.2="10",group.by="seurat_clusters",
                                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05,test.use="MAST")

cd8.subset.markers9v10.mast <- cbind(gene = rownames(cd8.subset.markers9v10.mast), cd8.subset.markers9v10.mast)
rownames(cd8.subset.markers9v10.mast) <- 1:nrow(cd8.subset.markers9v10.mast)
cd8.subset.markers9v10.mast<- cd8.subset.markers9v10.mast[, c(2,3,4,5,6,1)]

cd8.subset.markers9v10.strict <- FindMarkers(cd8.subset, ident.1 = "9", ident.2="10",group.by="seurat_clusters",
                                      only.pos = FALSE, min.pct = .10,logfc.threshold = 1)

cd8.subset.cluster.top10 <-cd8.subset.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
View(filter(cd8.subset.cluster.top10, cluster == 9| cluster == 10))


DoHeatmap(subset(cd8.subset,idents = c("9","10")), feature=filter(cd8.subset.cluster.top10, cluster == 9| cluster == 10)$gene)

DefaultAssay(pbmc.combined) <- "RNA"
DoHeatmap(subset(pbmc.combined,idents = c("CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3","CD8 TEM","CD8 Exhausted" ,"CD8 Teff")), features=c("IL7R","CCR7","EOMES","GZMK","CXCR3","FGFBP2","KLRD1","TNF","IFNG","FOS","JUN","LAG3","HAVCR2","PDCD1","GZMB","ENTPD1","ITGAE")) +  scale_fill_gradientn(colors = c("blue", "white", "red"))

DotPlot(subset(pbmc.combined,idents = c("CD8 Exhausted")),
        features=c("IL7R","CCR7","EOMES","GZMK","CXCR3","FGFBP2","KLRD1","TNF","IFNG","FOS","JUN","LAG3","HAVCR2","PDCD1","GZMB","ENTPD1","ITGAE"),
        split.by = "Timepoint",cols = c("blue", "red","orange","yellow","pink","purple"), dot.scale = 8) +   RotatedAxis()

FeaturePlot(subset(cd8.subset,idents = c("9","10")), features=c("TOX","PDCD1","NFATC1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 


FeaturePlot(subset(cd8.subset,idents = c("9","10")), features=c("PRF1","NME1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 


cd8.subset <- AddModuleScore(object = cd8.subset, features = proliferation.set, name = "proliferation.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = resting.set, name = "resting.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = ifn.response.set, name = "ifn.response.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = cd8.cytotoxic.set, name = "cd8.cytotoxic.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = cd8.cytokine.set, name = "cd8.cytokine.set")


FeaturePlot(object = cd8.subset, features = "proliferation.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "resting.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "ifn.response.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "cd8.cytotoxic.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "cd8.cytokine.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")

FeaturePlot(object = cd8.subset, features = "proliferation.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "resting.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "ifn.response.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "cd8.cytotoxic.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "cd8.cytokine.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)



