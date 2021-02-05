# subcluster cd8+ t-cells

cd8.subset <- subset(x = pbmc.combined.cellrep, idents = c("CD8 t-cell 1","CD8 t-cell 2", "CD8 t-cell 3","CD8 TEM","CD8 Exhausted" ,"CD8 Teff"))


cd8.subset.list <- SplitObject(cd8.subset, split.by = "Timepoint")

cd8.subset.list <- lapply(X = cd8.subset.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            dispersion.cutoff = c(0.5, Inf),mean.cutoff = c(0.0125, 3),nfeatures = 3000)
})



# 
#integrate data
cd8.subset.anchors <- FindIntegrationAnchors(object.list = cd8.subset.list, dims = 1:30)
cd8.subset <- IntegrateData(anchorset = cd8.subset.anchors, dims = 1:30)

rm(cd8.subset.list)
rm(cd8.subset.anchors)

DefaultAssay(cd8.subset) <- "integrated"

cd8.subset <- ScaleData(cd8.subset, verbose = TRUE)
cd8.subset <- RunPCA(cd8.subset, npcs = 30, verbose = TRUE)

cd8.subset <- FindNeighbors(cd8.subset, reduction = "pca", dims = 1:12)
cd8.subset <- FindClusters(cd8.subset,resolution = 0.8)

cd8.subset <- RunUMAP(cd8.subset, dims = 1:12)

# look at the UMAP dimplot now
DimPlot(cd8.subset, label= TRUE)
DimPlot(cd8.subset, label= TRUE,group.by="old.ident")
DimPlot(cd8.subset, label= TRUE,split.by="old.ident")
DimPlot(cd8.subset, label= TRUE,group.by="old.ident",split.by="Timepoint")

DefaultAssay(cd8.subset) <- "RNA"

FeaturePlot(object = cd8.subset, features = c("PDCD1"), cols = c("grey", "blue"),min.cutoff = "q10", max.cutoff = "q90")

table(Idents(cd8.subset), cd8.subset$Timepoint)

prop.table(table(Idents(cd8.subset), cd8.subset$Timepoint), margin = 2)


cd8.subset.markers <- FindAllMarkers(cd8.subset, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.05)



# cd8.subset.markers10v9 <- FindMarkers(cd8.subset, ident.1 = "10", ident.2="9",group.by="Timepoint",
#                               only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05,subset.ident = "4")

cd8.subset.markers9v8 <- FindMarkers(cd8.subset, ident.1 = "9", ident.2="8",group.by="seurat_clusters",
                                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

cd8.subset.markers9v8 <- cbind(gene = rownames(cd8.subset.markers9v8), cd8.subset.markers9v8)
rownames(cd8.subset.markers9v8) <- 1:nrow(cd8.subset.markers9v8)
cd8.subset.markers9v8<- cd8.subset.markers9v8[, c(2,3,4,5,6,1)]



cd8.subset.cluster.top10 <-cd8.subset.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(filter(cd8.subset.cluster.top10, cluster == 9| cluster == 10))

DefaultAssay(cd8.subset) <- "integrated"
DoHeatmap(subset(cd8.subset,idents = c("8","9")), feature=filter(cd8.subset.cluster.top10, cluster == 8| cluster == 9)$gene)

DoHeatmap(subset(cd8.subset,idents = c("10","11")), feature=filter(cd8.subset.cluster.top10, cluster == 10| cluster == 11)$gene)


FeaturePlot(subset(cd8.subset,idents = c("9","8")), features=c("GZMK","PDCD1","NFATC1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 


FeaturePlot(subset(cd8.subset,idents = c("9","10")), features=c("PRF1","NME1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 


DefaultAssay(cd8.subset) <- "RNA"

resting.set<- list(c("WDR86","IL7R","TSPAN2","LRRC2","SNAI3","CD52","TIMP1","YPEL4","CALHM2","S100A4","TSC22D3","EPHA4","ZNF831","ZCCHC18","GLIPR1","AIM1L","KLRB1","CNPY4","GDPD5","TC2N","AHNAK","PBXIP1","HIST1H2BD","BCL9L","EEPD1","ADAM23","THRA","KCTD7","CTSF","CAMK2N1","DCP1B","SEPT4","FXYD2","CECR1","DPYSL2","CCDC65","CDC25B","FAM229A","RP11-111M22.2","BCO2","HHAT","TGFB3","ANTXR2","AQP3","CRIP2","MYO1F","MPP7","NMT2","UTRN","NLRP3","BTD","KLF2","ZFP36L2","SUN2","FXYD1","BEST4","IGFBP6","SFXN3","CEP128","PLCD1","CYB561","ANKMY2","NBPF11","BAZ2B","IL11RA","ITGB7","ACSF2","ATXN7L1","PINK1","PDCD4"))
ifn.response.set <-list(c("IFIT3","CMPK2","IFIT2","IFI44L","GBP1","LGALS9","GBP7","PARP9","XAF1","STAT1","MX1","SAMD9L","GBP4","GBP5","HAPLN3","OAS3","LGALS3BP","DDX60","SECTM1","ETV7","NEXN","OAS1","IRF7","IFI6","TRAFD1","ISG15","IFI35","UBE2L6","IFIT5","IFI30","EPSTI1","DTX3L","IRF1","BCL2L14","LPIN2","IFI44","APOL6","STAT2","PARP14","TRIM22","IFIH1","IFITM3","OAS2","FBXO6","PARP12","ERAP2","HELZ2","JAK2","ISG20","ALPK1","SLC37A3","APOL1","RALB","TNFSF10","TAP2","SOCS1","PARP10","USP18","PLPP1","TRIM21","TAF4B","BST2","CES4A","HERC6","EIF2AK2","RTP4","APOL2","LPP","SP110","FAM19A2"))
proliferation.set <-list(c("PYCR1","NPW","LIF","IL2","SLC29A1","GINS2","ODC1","CENPV","HPDL","NME1","FABP5","CD3EAP","NLN","POLR3G","GJB6","NOP16","DLL3","CCDC86","CDC20","METTL1","ATAD3B","UCK2","SRM","TOMM40","EIF4EBP1","RRP9","ANKLE1","EBNA1BP2","TMEM97","ORC6","MRTO4","BOP1","G0S2","F12","WDR4","FOSL1","C17orf96","PUS7","ECE2","GCK","PDSS1","SMKR1","FKBP4","MRPL12","C16orf59","ATAD3A","DGAT2","RPP25","CENPN","POLR3H","MFSD2A","TLCD1","CHEK1","NOLC1","IFRD2","CYP27B1","ANKRD13B","CKS2","DDIAS","CTPS1","TTLL12","HSPE1","SH2D4A","YRDC","C10orf2","TRAP1","TIMM8A"))
cd8.cytotoxic.set <-list(c("CCL5","GZMK","GNLY","TRGC2","FGFBP2","C1orf21","KLRF1","FCGR3A","PTGDR","KLRC2","EOMES","S1PR5","CLIC3","AOAH","CADM1","TRGC1","DTHD1","LILRB1","SAMD3","ZNF683","KLRD1","NCR1","FAM49A","KLRG1","CTSW","CD244","CMC1","APOBEC3H","CST7","CX3CR1","FCRL6","TMCC3","PLA2G16","TYROBP","TPRG1","C12orf75","PLCG2","PLEK","RCAN2","DKK3","ADRB2","FCRL3","NKG7","PPP2R2B","SYNGR1","KLRC4","HLA-DPB1","DAPK2","F2R","KIR3DL2","B3GAT1","CD8B","TTC16","GALNT3","SCD5","PDGFD","ABCB1","MXRA7","CTBP2","CD8A","ZEB2","SYTL2","CHN2","FGR","TGFBR3","SETBP1","COLGALT2","KIR2DL4","FKBP1B","ADGRG1"))
cd8.cytokine.set <-list(c("CCL3","CCL3L3","CCL4L2","IFNG","CCL4","GZMB","XCL1","XCL2","CSF2","CCR1","IL10","HOPX","FASLG","GZMH","BATF3","TIMD4","FEZ1","HAVCR2","ZBTB32","SLC27A2","LAG3","PRF1","IL18RAP","GZMA","PTMS","PHLDA1","NCS1","IQCG","JAKMIP1","KLRC1","CRTAM","CCR5","EGR2","TNFRSF9","TNF","SLC4A10","ADAM19","MSC","METRNL","SLAMF7","NKG7","CXCR6","CD70","ATP8B4","SDC4","ELOVL6","ADAP1","IL26","RDH10","LTA","MICAL2","TNIP3","CCL20","ULBP2","CEBPD","RCAN2","ZEB2","CD226","SEMA4A","SLAMF1","SEMA7A","FAM3C","SPAG1","PRDM1","KLRD1","ST8SIA4","IER3","PLEK"))


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

FeaturePlot(object = cd8.subset, features = c("MKI67"),min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")


# clones to track
clones.track.cd8 <-c("CAVNNHFNKFYF","CASSQALATDTQYF","CASSLGAASTDTQYF","CASSEGHISSSTDTQYF","CASSPAEMNTEAFF","CSAALARYNEQFF","CATAGGLGNTIYF","CAMREGGYQKVTF","CASSLTTGNEQFF","CASSHRLAGGYNEQFF","CASSLFISTFPDGELFF","CASSYGQGVNEQFF","CASTWSGANVLTF","CASSHSATGESYEQYF","CASSLGVDEQFF","CASSYGQGVNEQFF","CASSVVSEGEAFF","CASSDGRLGNTEAFF","CASSLGLQNEQFF","CASRRARGGAYNEQFF","CASSVAGEHEQFF")


clones.cd8 <-filter(cd8.subset@meta.data, grepl(paste(clones.track.cd8, collapse="|"), cd8.subset$t_cdr3s_aa))

DimPlot(cd8.subset, label=F,  cells.highlight= filter(cd8.subset@meta.data, t_cdr3s_aa == "TRB:CASSVAGEHEQFF")$t_barcode,split.by="Timepoint")

CASSHSATGESYEQYF
CASSYGQGVNEQFF