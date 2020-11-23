DimPlot(pbmc.combined,reduction = "umap",label = TRUE)

table(Idents(pbmc.combined))
table(Idents(pbmc.combined),pbmc.combined$Timepoint)

pbmc.combined.cellrep <-pbmc.combined

## loading tcr data



## loading bcr data