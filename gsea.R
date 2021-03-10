#GSEA
library(msigdbr)
library(magrittr)



#set up msigdbr
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

## if DF of DEGS is from findallmarkers use this
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
cluster0.genelist <- filter(cd8.subset.markers, cluster==7 , pct.1>=.10)[,2]

## feature 2: named vector
names(cluster0.genelist) <- as.character(filter(cd8.subset.markers, cluster==7 , pct.1>=.10 )[,7])

## feature 3: decreasing order
cluster0.genelist <- sort(cluster0.genelist, decreasing = TRUE)

## if DF of DEGS is from findmarkers use this
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
cluster0.genelist <- cd8.teff.1v6[,2]
cluster0.genelist <- filter(cd8.tem1.1v4, pct.1>=.10, p_val_adj >=0.05)[,2]

## feature 2: named vector
names(cluster0.genelist) <- as.character(cd8.teff.1v6[,6])
names(cluster0.genelist) <- as.character(filter(cd8.tem1.1v4, pct.1>=.10,p_val_adj >=0.05)[,6])

## feature 3: decreasing order
cluster0.genelist <- sort(cluster0.genelist, decreasing = TRUE)

# gsea


em2 <- GSEA(cluster0.genelist, TERM2GENE = m_t2g)

# ridgeplot(em2,showCategory = 20) + ggtitle("CD8 Exhausted 2v5")

dotplot(em2,showCategory = 20,split=".sign")+ facet_grid(.~.sign)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4")

#GO
gse <- gseGO(geneList=cluster0.genelist, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb =org.Hs.eg.db, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4 GSE-GO")

# gene.test <- names(cluster0.genelist)[abs(cluster0.genelist) >1]
gene.test.up <- names(cluster0.genelist)[cluster0.genelist >= 1]
gene.test.down <- names(cluster0.genelist)[cluster0.genelist <= (-1)]
gene.test.2 <-names(cluster0.genelist)

gene.test.up<-bitr(gene.test.up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene.test.down<-bitr(gene.test.down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene.test.2<-bitr(gene.test.2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gse.up <- enrichGO(gene          = gene.test.up$ENTREZID,
                universe      = gene.test.2$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                )
dotplot(gse.up, showCategory =30)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4 GO UP")

gse.down <- enrichGO(gene          = gene.test.down$ENTREZID,
                   universe      = gene.test.2$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
)
dotplot(gse.down, showCategory =30)  + ggtitle("CD 8 TEM-1 Timepoint 1 vs Timpoint 4 GO DOWN")

