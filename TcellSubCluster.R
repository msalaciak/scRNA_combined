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

cd8.per<-table(Idents(cd8.subset.cc), cd8.subset.cc$Timepoint)

prop.table(table(Idents(cd8.subset.cc), cd8.subset.cc$Timepoint), margin = 2)


cd8.subset.markers <- FindAllMarkers(cd8.subset, only.pos = FALSE, min.pct = 0.1,logfc.threshold = 0.05)



cd8.subset.markers.8.t4<- FindMarkers(subset(cd8.subset,Timepoint ==4), ident.1 = "8",group.by="seurat_clusters",
                              only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

cd8.subset.markers9v8.t4 <- FindMarkers(subset(cd8.subset,Timepoint ==4), ident.1 = "9", ident.2="8",group.by="seurat_clusters",
                                      only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

cd8.subset.markers9v8.t4 <- cbind(gene = rownames(cd8.subset.markers9v8.t4), cd8.subset.markers9v8.t4)
cd8.subset.markers9v8.t4<- cd8.subset.markers9v8.t4[, c(2,3,4,5,6,1)]

#rownames(cd8.subset.markers9v8) <- cd8.subset.markers9v8$gene IF YOU WANT TO RENAME INDEX WITH GENES FOR PLOTTING.

Idents(cd8.subset.cc) <- "seurat_clusters"
# for loop to subset timepoint avglogfoldchange DEG data
for (i in 0:11){
  print("cluster ")
  print(i)
  cluster <- subset(cd8.subset.cc, idents = i)
  Idents(cluster) <- "Timepoint"

  cluster_compare_avg <- as.data.frame(log1p(AverageExpression(cluster, verbose = FALSE)$RNA))
  cluster_compare_avg$gene <- rownames(cluster_compare_avg)


  nam <- paste("cd8.subset.cc.cluster", i,".avgLOG_timepoint", sep = "")
  assign(nam, cluster_compare_avg)

  rm(cluster)
  rm(cluster_compare_avg)
  rm(nam)


}

test.1<-log1p(AverageExpression(subset(cd8.subset.cc,Timepoint==6,idents="CD8 TEFF"), verbose = FALSE)$RNA)
test.1 <-as.data.frame(test.1)
test.1$gene <- rownames(test.1)




cd8.subset.cluster.top10 <-cd8.subset.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(filter(cd8.subset.cluster.top10, cluster == 9| cluster == 10))

DefaultAssay(cd8.subset) <- "integrated"
DoHeatmap(subset(cd8.subset,idents = c("8","9")), feature=filter(cd8.subset.cluster.top10, cluster == 8| cluster == 9)$gene)

DoHeatmap(subset(cd8.subset,idents = c("10","11")), feature=filter(cd8.subset.cluster.top10, cluster == 10| cluster == 11)$gene)


FeaturePlot(subset(cd8.subset,idents = c("9","8")), features=c("GZMK","PDCD1","NFATC1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 


FeaturePlot(subset(cd8.subset,idents = c("9","10")), features=c("PRF1","NME1"),min.cutoff = "q10", max.cutoff = "q90",split.by="Timepoint") 

DoHeatmap(subset(cd8.subset,Timepoint==4,downsample = 50), feature=c('ABCB7','ACAA1','ACAA2','ACADM','ACADSB','ACADVL','ACAT1','ACO2','AFG3L2','AIFM1','ALAS1','ALDH6A1','ATP1B1','ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5MC1','ATP5MC2','ATP5MC3','ATP5ME','ATP5MF','ATP5MG','ATP5PB','ATP5PD','ATP5PF','ATP5PO','ATP6AP1','ATP6V0B','ATP6V0C','ATP6V0E1','ATP6V1C1','ATP6V1D','ATP6V1E1','ATP6V1F','ATP6V1G1','ATP6V1H','BAX','BCKDHA','BDH2','CASP7','COX10','COX11','COX15','COX17','COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A2','COX7A2L','COX7B','COX7C','COX8A','CPT1A','CS','CYB5A','CYB5R3','CYC1','CYCS','DECR1','DLAT','DLD','DLST','ECH1','ECHS1','ECI1','ETFA','ETFB','ETFDH','FDX1','FH','FXN','GLUD1','GOT2','GPI','GPX4','GRPEL1','HADHA','HADHB','HCCS','HSD17B10','HSPA9','HTRA2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','IMMT','ISCA1','ISCU','LDHA','LDHB','LRPPRC','MAOB','MDH1','MDH2','MFN2','MGST3','MPC1','MRPL11','MRPL15','MRPL34','MRPL35','MRPS11','MRPS12','MRPS15','MRPS22','MRPS30','MTRF1','MTRR','MTX2','NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFAB1','NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFC1','NDUFC2','NDUFS1','NDUFS2','NDUFS3','NDUFS4','NDUFS6','NDUFS7','NDUFS8','NDUFV1','NDUFV2','NNT','NQO2','OAT','OGDH','OPA1','OXA1L','PDHA1','PDHB','PDHX','PDK4','PDP1','PHB2','PHYH','PMPCA','POLR2F','POR','PRDX3','RETSAT','RHOT1','RHOT2','SDHA','SDHB','SDHC','SDHD','SLC25A11','SLC25A12','SLC25A20','SLC25A3','SLC25A4','SLC25A5','SLC25A6','SUCLA2','SUCLG1','SUPV3L1','SURF1','TCIRG1','TIMM10','TIMM13','TIMM17A','TIMM50','TIMM8B','TIMM9','TOMM22','TOMM70','UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ','VDAC1','VDAC2','VDAC3'))

  
DoHeatmap(subset(cd8.subset,Timepoint==4,downsample = 50), feature=c('ABCB6','ADORA2B','AGL','AGRN','AK3','AK4','AKR1A1','ALDH7A1','ALDH9A1','ALDOA','ALDOB','ALG1','ANG','ANGPTL4','ANKZF1','ARPP19','ARTN','AURKA','B3GALT6','B3GAT1','B3GAT3','B3GNT3','B4GALT1','B4GALT2','B4GALT4','B4GALT7','BIK','BPNT1','CACNA1H','CAPN5','CASP6','CD44','CDK1','CENPA','CHPF','CHPF2','CHST1','CHST12','CHST2','CHST4','CHST6','CITED2','CLDN3','CLDN9','CLN6','COG2','COL5A1','COPB2','CTH','CXCR4','CYB5A','DCN','DDIT4','DEPDC1','DLD','DPYSL4','DSC2','ECD','EFNA3','EGFR','EGLN3','ELF3','ENO1','ENO2','ERO1A','EXT1','EXT2','FAM162A','FBP2','FKBP4','FUT8','G6PD','GAL3ST1','GALE','GALK1','GALK2','GAPDHS','GCLC','GFPT1','GLCE','GLRX','GMPPA','GMPPB','GNE','GNPDA1','GOT1','GOT2','GPC1','GPC3','GPC4','GPR87','GUSB','GYS1','GYS2','HAX1','HDLBP','HK2','HMMR','HOMER1','HS2ST1','HS6ST2','HSPA5','IDH1','IDUA','IER3','IGFBP3','IL13RA1','IRS2','ISG20','KDELR3','KIF20A','KIF2A','LCT','LDHA','LDHC','LHPP','LHX9','MDH1','MDH2','ME1','ME2','MED24','MERTK','MET','MIF','MIOX','MPI','MXI1','NANP','NASP','NDST3','NDUFV3','NOL3','NSDHL','NT5E','P4HA1','P4HA2','PAM','PAXIP1','PC','PDK3','PFKFB1','PFKP','PGAM1','PGAM2','PGK1','PGLS','PGM2','PHKA2','PKM','PKP2','PLOD1','PLOD2','PMM2','POLR3K','PPFIA4','PPIA','PPP2CB','PRPS1','PSMC4','PYGB','PYGL','QSOX1','RARS1','RBCK1','RPE','RRAGD','SAP30','SDC1','SDC2','SDC3','SDHC','SLC16A3','SLC25A10','SLC25A13','SLC35A3','SLC37A4','SOD1','SOX9','SPAG4','SRD5A3','STC1','STC2','STMN1','TALDO1','TFF3','TGFA','TGFBI','TKTL1','TPBG','TPI1','TPST1','TSTA3','TXN','UGP2','VCAN','VEGFA','VLDLR','XYLT2','ZNF292'))
DefaultAssay(cd8.subset) <- "RNA"

resting.set<- list(c("WDR86","IL7R","TSPAN2","LRRC2","SNAI3","CD52","TIMP1","YPEL4","CALHM2","S100A4","TSC22D3","EPHA4","ZNF831","ZCCHC18","GLIPR1","AIM1L","KLRB1","CNPY4","GDPD5","TC2N","AHNAK","PBXIP1","HIST1H2BD","BCL9L","EEPD1","ADAM23","THRA","KCTD7","CTSF","CAMK2N1","DCP1B","SEPT4","FXYD2","CECR1","DPYSL2","CCDC65","CDC25B","FAM229A","RP11-111M22.2","BCO2","HHAT","TGFB3","ANTXR2","AQP3","CRIP2","MYO1F","MPP7","NMT2","UTRN","NLRP3","BTD","KLF2","ZFP36L2","SUN2","FXYD1","BEST4","IGFBP6","SFXN3","CEP128","PLCD1","CYB561","ANKMY2","NBPF11","BAZ2B","IL11RA","ITGB7","ACSF2","ATXN7L1","PINK1","PDCD4"))
ifn.response.set <-list(c("IFIT3","CMPK2","IFIT2","IFI44L","GBP1","LGALS9","GBP7","PARP9","XAF1","STAT1","MX1","SAMD9L","GBP4","GBP5","HAPLN3","OAS3","LGALS3BP","DDX60","SECTM1","ETV7","NEXN","OAS1","IRF7","IFI6","TRAFD1","ISG15","IFI35","UBE2L6","IFIT5","IFI30","EPSTI1","DTX3L","IRF1","BCL2L14","LPIN2","IFI44","APOL6","STAT2","PARP14","TRIM22","IFIH1","IFITM3","OAS2","FBXO6","PARP12","ERAP2","HELZ2","JAK2","ISG20","ALPK1","SLC37A3","APOL1","RALB","TNFSF10","TAP2","SOCS1","PARP10","USP18","PLPP1","TRIM21","TAF4B","BST2","CES4A","HERC6","EIF2AK2","RTP4","APOL2","LPP","SP110","FAM19A2"))
proliferation.set <-list(c("PYCR1","NPW","LIF","IL2","SLC29A1","GINS2","ODC1","CENPV","HPDL","NME1","FABP5","CD3EAP","NLN","POLR3G","GJB6","NOP16","DLL3","CCDC86","CDC20","METTL1","ATAD3B","UCK2","SRM","TOMM40","EIF4EBP1","RRP9","ANKLE1","EBNA1BP2","TMEM97","ORC6","MRTO4","BOP1","G0S2","F12","WDR4","FOSL1","C17orf96","PUS7","ECE2","GCK","PDSS1","SMKR1","FKBP4","MRPL12","C16orf59","ATAD3A","DGAT2","RPP25","CENPN","POLR3H","MFSD2A","TLCD1","CHEK1","NOLC1","IFRD2","CYP27B1","ANKRD13B","CKS2","DDIAS","CTPS1","TTLL12","HSPE1","SH2D4A","YRDC","C10orf2","TRAP1","TIMM8A"))
cd8.cytotoxic.set <-list(c("CCL5","GZMK","GNLY","TRGC2","FGFBP2","C1orf21","KLRF1","FCGR3A","PTGDR","KLRC2","EOMES","S1PR5","CLIC3","AOAH","CADM1","TRGC1","DTHD1","LILRB1","SAMD3","ZNF683","KLRD1","NCR1","FAM49A","KLRG1","CTSW","CD244","CMC1","APOBEC3H","CST7","CX3CR1","FCRL6","TMCC3","PLA2G16","TYROBP","TPRG1","C12orf75","PLCG2","PLEK","RCAN2","DKK3","ADRB2","FCRL3","NKG7","PPP2R2B","SYNGR1","KLRC4","HLA-DPB1","DAPK2","F2R","KIR3DL2","B3GAT1","CD8B","TTC16","GALNT3","SCD5","PDGFD","ABCB1","MXRA7","CTBP2","CD8A","ZEB2","SYTL2","CHN2","FGR","TGFBR3","SETBP1","COLGALT2","KIR2DL4","FKBP1B","ADGRG1"))
cd8.cytokine.set <-list(c("CCL3","CCL3L3","CCL4L2","IFNG","CCL4","GZMB","XCL1","XCL2","CSF2","CCR1","IL10","HOPX","FASLG","GZMH","BATF3","TIMD4","FEZ1","HAVCR2","ZBTB32","SLC27A2","LAG3","PRF1","IL18RAP","GZMA","PTMS","PHLDA1","NCS1","IQCG","JAKMIP1","KLRC1","CRTAM","CCR5","EGR2","TNFRSF9","TNF","SLC4A10","ADAM19","MSC","METRNL","SLAMF7","NKG7","CXCR6","CD70","ATP8B4","SDC4","ELOVL6","ADAP1","IL26","RDH10","LTA","MICAL2","TNIP3","CCL20","ULBP2","CEBPD","RCAN2","ZEB2","CD226","SEMA4A","SLAMF1","SEMA7A","FAM3C","SPAG1","PRDM1","KLRD1","ST8SIA4","IER3","PLEK"))
OXIDATIVE_PHOSPHORYLATION<-list(c('ABCB7','ACAA1','ACAA2','ACADM','ACADSB','ACADVL','ACAT1','ACO2','AFG3L2','AIFM1','ALAS1','ALDH6A1','ATP1B1','ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5MC1','ATP5MC2','ATP5MC3','ATP5ME','ATP5MF','ATP5MG','ATP5PB','ATP5PD','ATP5PF','ATP5PO','ATP6AP1','ATP6V0B','ATP6V0C','ATP6V0E1','ATP6V1C1','ATP6V1D','ATP6V1E1','ATP6V1F','ATP6V1G1','ATP6V1H','BAX','BCKDHA','BDH2','CASP7','COX10','COX11','COX15','COX17','COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A2','COX7A2L','COX7B','COX7C','COX8A','CPT1A','CS','CYB5A','CYB5R3','CYC1','CYCS','DECR1','DLAT','DLD','DLST','ECH1','ECHS1','ECI1','ETFA','ETFB','ETFDH','FDX1','FH','FXN','GLUD1','GOT2','GPI','GPX4','GRPEL1','HADHA','HADHB','HCCS','HSD17B10','HSPA9','HTRA2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','IMMT','ISCA1','ISCU','LDHA','LDHB','LRPPRC','MAOB','MDH1','MDH2','MFN2','MGST3','MPC1','MRPL11','MRPL15','MRPL34','MRPL35','MRPS11','MRPS12','MRPS15','MRPS22','MRPS30','MTRF1','MTRR','MTX2','NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFAB1','NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFC1','NDUFC2','NDUFS1','NDUFS2','NDUFS3','NDUFS4','NDUFS6','NDUFS7','NDUFS8','NDUFV1','NDUFV2','NNT','NQO2','OAT','OGDH','OPA1','OXA1L','PDHA1','PDHB','PDHX','PDK4','PDP1','PHB2','PHYH','PMPCA','POLR2F','POR','PRDX3','RETSAT','RHOT1','RHOT2','SDHA','SDHB','SDHC','SDHD','SLC25A11','SLC25A12','SLC25A20','SLC25A3','SLC25A4','SLC25A5','SLC25A6','SUCLA2','SUCLG1','SUPV3L1','SURF1','TCIRG1','TIMM10','TIMM13','TIMM17A','TIMM50','TIMM8B','TIMM9','TOMM22','TOMM70','UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ','VDAC1','VDAC2','VDAC3'))
GLYCOLYSIS<-list(c('ABCB6','ADORA2B','AGL','AGRN','AK3','AK4','AKR1A1','ALDH7A1','ALDH9A1','ALDOA','ALDOB','ALG1','ANG','ANGPTL4','ANKZF1','ARPP19','ARTN','AURKA','B3GALT6','B3GAT1','B3GAT3','B3GNT3','B4GALT1','B4GALT2','B4GALT4','B4GALT7','BIK','BPNT1','CACNA1H','CAPN5','CASP6','CD44','CDK1','CENPA','CHPF','CHPF2','CHST1','CHST12','CHST2','CHST4','CHST6','CITED2','CLDN3','CLDN9','CLN6','COG2','COL5A1','COPB2','CTH','CXCR4','CYB5A','DCN','DDIT4','DEPDC1','DLD','DPYSL4','DSC2','ECD','EFNA3','EGFR','EGLN3','ELF3','ENO1','ENO2','ERO1A','EXT1','EXT2','FAM162A','FBP2','FKBP4','FUT8','G6PD','GAL3ST1','GALE','GALK1','GALK2','GAPDHS','GCLC','GFPT1','GLCE','GLRX','GMPPA','GMPPB','GNE','GNPDA1','GOT1','GOT2','GPC1','GPC3','GPC4','GPR87','GUSB','GYS1','GYS2','HAX1','HDLBP','HK2','HMMR','HOMER1','HS2ST1','HS6ST2','HSPA5','IDH1','IDUA','IER3','IGFBP3','IL13RA1','IRS2','ISG20','KDELR3','KIF20A','KIF2A','LCT','LDHA','LDHC','LHPP','LHX9','MDH1','MDH2','ME1','ME2','MED24','MERTK','MET','MIF','MIOX','MPI','MXI1','NANP','NASP','NDST3','NDUFV3','NOL3','NSDHL','NT5E','P4HA1','P4HA2','PAM','PAXIP1','PC','PDK3','PFKFB1','PFKP','PGAM1','PGAM2','PGK1','PGLS','PGM2','PHKA2','PKM','PKP2','PLOD1','PLOD2','PMM2','POLR3K','PPFIA4','PPIA','PPP2CB','PRPS1','PSMC4','PYGB','PYGL','QSOX1','RARS1','RBCK1','RPE','RRAGD','SAP30','SDC1','SDC2','SDC3','SDHC','SLC16A3','SLC25A10','SLC25A13','SLC35A3','SLC37A4','SOD1','SOX9','SPAG4','SRD5A3','STC1','STC2','STMN1','TALDO1','TFF3','TGFA','TGFBI','TKTL1','TPBG','TPI1','TPST1','TSTA3','TXN','UGP2','VCAN','VEGFA','VLDLR','XYLT2','ZNF292'))
MTOR<-list(c('ABCF2','ACACA','ACLY','ACSL3','ACTR2','ACTR3','ADD3','ADIPOR2','AK4','ALDOA','ARPC5L','ASNS','ATP2A2','ATP5MC1','ATP6V1D','AURKA','BCAT1','BHLHE40','BTG2','BUB1','CACYBP','CALR','CANX','CCNF','CCNG1','CCT6A','CD9','CDC25A','CDKN1A','CFP','COPS5','CORO1A','CTH','CTSC','CXCR4','CYB5B','CYP51A1','DAPP1','DDIT3','DDIT4','DDX39A','DHCR24','DHCR7','DHFR','EBP','EDEM1','EEF1E1','EGLN3','EIF2S2','ELOVL5','ELOVL6','ENO1','EPRS1','ERO1A','ETF1','FADS1','FADS2','FDXR','FGL2','FKBP2','G6PD','GAPDH','GBE1','GCLC','GGA2','GLA','GLRX','GMPS','GOT1','GPI','GSK3B','GSR','GTF2H1','HK2','HMBS','HMGCR','HMGCS1','HPRT1','HSP90B1','HSPA4','HSPA5','HSPA9','HSPD1','HSPE1','IDH1','IDI1','IFI30','IFRD1','IGFBP5','IMMT','INSIG1','ITGB2','LDHA','LDLR','LGMN','LTA4H','M6PR','MAP2K3','MCM2','MCM4','ME1','MLLT11','MTHFD2','MTHFD2L','NAMPT','NFIL3','NFKBIB','NFYC','NIBAN1','NMT1','NUFIP1','NUP205','NUPR1','P4HA1','PDAP1','PDK1','PFKL','PGK1','PGM1','PHGDH','PIK3R3','PITPNB','PLK1','PLOD2','PNO1','PNP','POLR3G','PPA1','PPIA','PPP1R15A','PRDX1','PSAT1','PSMA3','PSMA4','PSMB5','PSMC2','PSMC4','PSMC6','PSMD12','PSMD13','PSMD14','PSME3','PSMG1','PSPH','QDPR','RAB1A','RDH11','RIT1','RPA1','RPN1','RRM2','RRP9','SC5D','SCD','SDF2L1','SEC11A','SERP1','SERPINH1','SHMT2','SKAP2','SLA','SLC1A4','SLC1A5','SLC2A1','SLC2A3','SLC37A4','SLC6A6','SLC7A11','SLC7A5','SLC9A3R1','SORD','SQLE','SQSTM1','SRD5A1','SSR1','STARD4','STC1','STIP1','SYTL2','TBK1','TCEA1','TES','TFRC','TM7SF2','TMEM97','TOMM40','TPI1','TRIB3','TUBA4A','TUBG1','TXNRD1','UBE2D3','UCHL5','UFM1','UNG','USO1','VLDLR','WARS1','XBP1','YKT6'))
PI3K_AKT_MTOR <-list(c('ACACA','ACTR2','ACTR3','ADCY2','AKT1','AKT1S1','AP2M1','ARF1','ARHGDIA','ARPC3','ATF1','CAB39','CAB39L','CALR','CAMK4','CDK1','CDK2','CDK4','CDKN1A','CDKN1B','CFL1','CLTC','CSNK2B','CXCR4','DAPP1','DDIT3','DUSP3','E2F1','ECSIT','EGFR','EIF4E','FASLG','FGF17','FGF22','FGF6','GNA14','GNGT1','GRB2','GRK2','GSK3B','HRAS','HSP90B1','IL2RG','IL4','IRAK4','ITPR2','LCK','MAP2K3','MAP2K6','MAP3K7','MAPK1','MAPK10','MAPK8','MAPK9','MAPKAP1','MKNK1','MKNK2','MYD88','NCK1','NFKBIB','NGF','NOD1','PAK4','PDK1','PFN1','PIK3R3','PIKFYVE','PIN1','PITX2','PLA2G12A','PLCB1','PLCG1','PPP1CA','PPP2R1B','PRKAA2','PRKAG1','PRKAR2A','PRKCB','PTEN','PTPN11','RAC1','RAF1','RALB','RIPK1','RIT1','RPS6KA1','RPS6KA3','RPTOR','SFN','SLA','SLC2A1','SMAD2','SQSTM1','STAT2','TBK1','THEM4','TIAM1','TNFRSF1A','TRAF2','TRIB3','TSC2','UBE2D3','UBE2N','VAV3','YWHAB'))
FATTY_ACID_OX<-list(c('ACADM','ACADS','ACADVL','ADIPOR1','ADIPOR2','ALOX12','BDH2','CPT1A','CPT1B','ECH1','ECHS1','HACL1','HADHB','HAO1','HAO2','PPARA','PPARD','PPARGC1A'))





cd8.subset <- AddModuleScore(object = cd8.subset, features = proliferation.set, name = "proliferation.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = resting.set, name = "resting.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = ifn.response.set, name = "ifn.response.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = cd8.cytotoxic.set, name = "cd8.cytotoxic.set")
cd8.subset <- AddModuleScore(object = cd8.subset, features = cd8.cytokine.set, name = "cd8.cytokine.set")

cd8.subset.cc <- AddModuleScore(object = cd8.subset.cc, features = OXIDATIVE_PHOSPHORYLATION, name = "OXIDATIVE_PHOSPHORYLATION")
cd8.subset.cc <- AddModuleScore(object = cd8.subset.cc, features = GLYCOLYSIS, name = "GLYCOLYSIS")
cd8.subset.cc <- AddModuleScore(object = cd8.subset.cc, features = MTOR, name = "MTOR")
cd8.subset.cc <- AddModuleScore(object = cd8.subset.cc, features = PI3K_AKT_MTOR, name = "PI3K_AKT_MTOR")
cd8.subset.cc <- AddModuleScore(object = cd8.subset.cc, features = FATTY_ACID_OX, name = "FATTY_ACID_OX")



FeaturePlot(object = cd8.subset, features = "proliferation.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "resting.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "ifn.response.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "cd8.cytotoxic.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")
FeaturePlot(object = cd8.subset, features = "cd8.cytokine.set1",min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")

FeaturePlot(object = cd8.subset, features = "proliferation.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "resting.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)
FeaturePlot(object = cd8.subset, features = "ifn.response.set1",min.cutoff = "q10", max.cutoff = "q90", label=F)

FeaturePlot(object = cd8.subset.scores, features = c("GLYCOLYSIS1","OXIDATIVE_PHOSPHORYLATION1"), label=F,split.by = "Timepoint")
FeaturePlot(object = cd8.subset.scores, features = "OXIDATIVE_PHOSPHORYLATION1",min.cutoff = "q10", max.cutoff = "q90", label=F,split.by = "Timepoint")


FeaturePlot(object = cd8.subset, features = c("MKI67"),min.cutoff = "q10", max.cutoff = "q90", label=F, split.by = "Timepoint")


FeaturePlot(object = subset(cd8.subset,idents = c("8","9")), features = c("IFNG"),min.cutoff = "q10", max.cutoff = "q90", label=F)





filter(cd8.subset@meta.data, grepl(paste("CASSPAEMNTEAFF", collapse="|"), cd8.subset$t_cdr3s_aa))$t_barcode


DimPlot(cd8.subset.cc, label=F,  cells.highlight=filter(cd8.subset@meta.data, 
                                                      grepl(paste("CASSVAGEHEQFF", collapse="|"), 
                                                      cd8.subset$t_cdr3s_aa))$t_barcode,split.by="Timepoint")


DimPlot(cd8.subset.cc, label=F,  cells.highlight= filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode,split.by="Timepoint")


cd8.teff.1v4 <- FindMarkers(cd8.subset.cc, ident.1 = "1", ident.2="6",group.by="Timepoint",
                                        only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05, subset.ident = "CD8 TEM-1")

cd8.teff.2 <- FindMarkers(cd8.subset.cc, ident.1 = "1", ident.2="6",group.by="Timepoint",
                            only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05, subset.ident = "CD8 TEM-1")

cd8.teff.1v4 <- cbind(gene = rownames(cd8.teff.1v4), cd8.teff.1v4)
cd8.teff.1v4<- cd8.teff.1v4[, c(2,3,4,5,6,1)]


cd8.tex2 <- FindMarkers(subset(cd8.subset.cc,Timepoint ==6), ident.1 = "9",group.by="seurat_clusters",
                                        only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0.05)

cd8.tex2 <- cbind(gene = rownames(cd8.tex2), cd8.tex2)
cd8.tex2<- cd8.tex2[, c(2,3,4,5,6,1)]

EnhancedVolcano(cd8.tex2,
                lab = rownames(cd8.tex2),
                x = 'avg_log2FC',
                 y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEX-2 Timepoint 6',
                drawConnectors = TRUE,
                widthConnectors = 0.75)



EnhancedVolcano(cd8.tem1.1v4,
                lab = rownames(cd8.tem1.1v4),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEM-1 Timepoint 1 vs Timepoint 4',
                drawConnectors = TRUE,
                widthConnectors = 0.75)

EnhancedVolcano(cd8.teff.1v6,
                lab = rownames(cd8.teff.1v6),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEFF Timepoint 1 vs Timepoint 6',
                drawConnectors = TRUE,
                widthConnectors = 0.75)


EnhancedVolcano(cd8.tem1.1v6,
                lab = rownames(cd8.tem1.1v6),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEM-1 Timepoint 1 vs Timepoint 6',
                drawConnectors = TRUE,
                widthConnectors = 0.75)

EnhancedVolcano(cd8.teff.2v5,
                lab = rownames(cd8.teff.2v5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEFF Timepoint 2 vs Timepoint 5',
                drawConnectors = TRUE,
                widthConnectors = 0.75)

EnhancedVolcano(cd8.tem1.2v5,
                lab = rownames(cd8.tem1.2v5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'CD8 TEM-1 Timepoint 2 vs Timepoint 5',
                drawConnectors = TRUE,
                widthConnectors = 0.75)



### cell cycle data

DimPlot(cd8.subset.cc, label= TRUE)

FeaturePlot(cd8.subset.cc, 
            features = c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1","IL7R","IFNG","SELL","CD28","CD27","CD69","FAS","CD44","LAMP1","HLA-DRA","ICOS"), 
            cols = c("grey", "blue"),min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(cd8.subset.cc,  c("SELL", "CCR7"))
plot_density(cd8.subset.cc, c("SELL", "CCR7"), joint = TRUE)

plot_density(cd8.subset.cc, c("FOS", "JUN"))
plot_density(cd8.subset.cc, c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1","KLRD1","IL7R","IFNG","SELL","CCR7","CD28","CD27","CD69","FOS","JUN","CD44","LAMP1","HLA-DRA","ICOS"))
plot_density(cd8.subset.cc, c("PDCD1"))


plot_density(cd8.subset.cc, features="OXIDATIVE_PHOSPHORYLATION1")
plot_density(cd8.subset.cc, features="GLYCOLYSIS1")
plot_density(subset(cd8.subset.cc,Timepoint==1), features="MTOR1")
plot_density(cd8.subset.cc, features="PI3K_AKT_MTOR1")
plot_density(cd8.subset.cc, features="FATTY_ACID_OX1")

for (x in 1:1) {
plot_density(subset(cd8.subset.cc,Timepoint==x), features=c("MTOR1","PI3K_AKT_MTOR1","OXIDATIVE_PHOSPHORYLATION1","GLYCOLYSIS1"))
plot_density(subset(cd8.subset.cc,Timepoint==x), features=c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1","KLRD1","IFNG"))
plot_density(subset(cd8.subset.cc,Timepoint==x), features=c("IL7R","SELL","CCR7","CD28","CD27","CD44"))
plot_density(subset(cd8.subset.cc,Timepoint==x), features=c("CD69","FOS","JUN","LAMP1","HLA-DRA","ICOS"))

}


p<-plot_density(subset(cd8.subset.cc,idents=c("CD8 TEX-1","CD8 TEX-2"),Timepoint==2), features=c("MTOR1","PI3K_AKT_MTOR1","OXIDATIVE_PHOSPHORYLATION1","GLYCOLYSIS1"))
p1<-plot_density(subset(cd8.subset.cc,idents=c("CD8 TEX-1","CD8 TEX-2"),Timepoint==2), features=c("GZMA","GZMH","GZMB","GZMK","PRF1","GNLY","KLRG1","KLRF1","KLRC1","KLRD1","IFNG",'MKI67'))
p2<-plot_density(subset(cd8.subset.cc,idents=c("CD8 TEX-1","CD8 TEX-2"),Timepoint==2), features=c("IL7R","SELL","CCR7","CD28","CD27","CD44"))
p3<-plot_density(subset(cd8.subset.cc,idents=c("CD8 TEX-1","CD8 TEX-2"),Timepoint==2), features=c("CD69","FOS","JUN","LAMP1","HLA-DRA","ICOS")) 
# p4<-plot_density(subset(cd8.subset.cc,Timepoint==2), features=c("PDCD1","EOMES","TOX","CD244","LAG3","CTLA4",'HAVCR2','TBX21','TCF7')) 

p + plot_annotation(
  title = 'Timepoint 2')
p1 + plot_annotation(
  title = 'Timepoint 2')
p2 + plot_annotation(
  title = 'Timepoint 2')
p3 + plot_annotation(
  title = 'Timepoint 2')
p4 + plot_annotation(
  title = 'Timepoint 5')

VlnPlot(subset(cd8.subset.cc,Timepoint==1), features=c("PDCD1","EOMES","TOX","CD244","LAG3","CTLA4",'HAVCR2','TBX21','TCF7')) 


TCR_4 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSQGTGYTNTEAFF")$t_barcode
TCR_20 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CASLNQAGTALIF")$t_barcode
TCR_9 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CASLNQAGTALIF;TRB:CASSQGTGYTNTEAFF")$t_barcode
TCR_36 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode
TCR_22 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CSARDSAASTDTQYF")$t_barcode
TCR_32 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSPKDYNNEQFF")$t_barcode
TCR_11 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSVAGEHEQFF")$t_barcode
TCR_43 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLGLRESEQFF")$t_barcode
TCR_210 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLTTGNEQFF")$t_barcode
TCR_75 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLLQGKINEQFF")$t_barcode
TCR_1475 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASGPGTYGYTF")$t_barcode
TCR_374 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSHSATGESYEQYF")$t_barcode
TCR_59 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSLGGATDTQYF")$t_barcode
TCR_3801 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRA:CAVSPMNTGFQKLVF;TRB:CASSPAEMNTEAFF")$t_barcode
TCR_187 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CSARDPGLAGKWDTQYF")$t_barcode
TCR_3873 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSPAEMNTEAFF")$t_barcode
TCR_881 <- filter(cd8.subset.cc@meta.data, t_cdr3s_aa == "TRB:CASSHGGGNYEQYF")$t_barcode


DimPlot(cd8.subset.cc, label=F,  cells.highlight= list(TCR_4,TCR_20,TCR_9,TCR_36,TCR_22,TCR_32,TCR_11,TCR_43,TCR_210,TCR_75,TCR_1475,TCR_374,TCR_59,TCR_3801,TCR_187,TCR_3873,TCR_881),split.by="Timepoint")  + 
  scale_color_manual(labels = c("unselected","4_TCR","20_TCR","9_TCR","36_TCR","22_TCR","32_TCR","11_TCR","43_TCR","210_TCR","75_TCR","1475_TCR","374_TCR","59_TCR","3801_TCR","187_TCR","3873_TCR","881_TCR"), 
                     values = c("darkgrey", "blue","darkolivegreen3","brown4","black","deeppink2","blueviolet","red","orange","darkgoldenrod3","darkseagreen","deepskyblue","cyan","coral","cadetblue4","burlywood","chartreuse1","darkslateblue")) +
  labs(color = "Top TCR's Across Timepoint") 


DimPlot(cd8.subset.cc, label=F,  cells.highlight= list(TCR_4,TCR_20,TCR_9,TCR_36,TCR_22,TCR_32,TCR_11,TCR_43,TCR_210,TCR_75,TCR_1475,TCR_374,TCR_59,TCR_3801,TCR_187,TCR_3873,TCR_881),split.by="Timepoint")  
  
DimPlot(cd8.subset.cc, label=F,  cells.highlight= list(TCR_4,TCR_20,TCR_9,TCR_36,TCR_22,TCR_32,TCR_11,TCR_43,TCR_210,TCR_75,TCR_1475,TCR_374,TCR_59,TCR_3801,TCR_187,TCR_3873,TCR_881),split.by="Timepoint")  + 
  scale_color_manual(labels = c("unselected","210_TCR","43_TCR","11_TCR","32_TCR","22_TCR","36_TCR","9_TCR","20_TCR","881_TCR","3873_TCR","187_TCR","3801_TCR","59_TCR","374_TCR","1475_TCR","75_TCR","4_TCR"), 
                     values = c("darkgrey", "blue","darkolivegreen3","brown4","black","deeppink2","blueviolet","red","orange","darkgoldenrod3","darkseagreen","deepskyblue","cyan","coral","cadetblue4","burlywood","chartreuse1","darkslateblue")) +
  labs(color = "Top TCR's Across Timepoint") 

# convert to scanpy
# SaveH5Seurat(cd8.subset.cc, filename = "cd8.h5Seurat")
# Convert("cd8.h5Seurat", dest = "h5ad")

