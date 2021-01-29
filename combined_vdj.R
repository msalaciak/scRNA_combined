
   
   ##########################

table(Idents(pbmc.combined))
table(Idents(pbmc.combined),pbmc.combined$Timepoint)

pbmc.combined.cellrep <-pbmc.combined

## loading tcr data

# 151231
timepoint1.tcell <- read.csv("/home/matthew/datatransfer/contig/151231_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 160411
timepoint2.tcell <- read.csv("/home/matthew/datatransfer/contig/160411_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 161962
timepoint3.tcell <- read.csv("/home/matthew/datatransfer/contig/161962_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 171094
timepoint4.tcell <- read.csv("/home/matthew/datatransfer/contig/171094_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 171642
timepoint5.tcell <- read.csv("/home/matthew/datatransfer/contig/171642_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 180251
timepoint6.tcell <- read.csv("/home/matthew/datatransfer/contig/180251_tcell/filtered_contig_annotations.csv", stringsAsFactors = F)

## loading bcr data

# 151231
timepoint1.bcell <- read.csv("/home/matthew/datatransfer/contig/151231_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 160411
timepoint2.bcell <- read.csv("/home/matthew/datatransfer/contig/160411_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 161962
timepoint3.bcell <- read.csv("/home/matthew/datatransfer/contig/161962_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 171094
timepoint4.bcell <- read.csv("/home/matthew/datatransfer/contig/171094_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 171642
timepoint5.bcell <- read.csv("/home/matthew/datatransfer/contig/171642_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)
# 180251
timepoint6.bcell <- read.csv("/home/matthew/datatransfer/contig/180251_bcell/filtered_contig_annotations.csv", stringsAsFactors = F)

#combined the lists
contig_list_tcell <- list(timepoint1.tcell,timepoint2.tcell,timepoint3.tcell,timepoint4.tcell,timepoint5.tcell,timepoint6.tcell)

contig_list_bcell <- list(timepoint1.bcell,timepoint2.bcell,timepoint3.bcell,timepoint4.bcell,timepoint5.bcell,timepoint6.bcell)

# my own script to add vdj data to seurat

#rename list names to timepoints
names(contig_list_tcell) <- c("timepoint-1", "timepoint-2","timepoint-3", "timepoint-4" ,"timepoint-5", "timepoint-6")
names(contig_list_bcell) <- c("timepoint-1", "timepoint-2","timepoint-3", "timepoint-4" ,"timepoint-5", "timepoint-6")

#add timepoint suffix
contig_list_tcell$`timepoint-1`$barcode = paste(contig_list_tcell$`timepoint-1`$barcode,"_1",sep = "")
contig_list_tcell$`timepoint-2`$barcode = paste(contig_list_tcell$`timepoint-2`$barcode,"_2",sep = "")
contig_list_tcell$`timepoint-3`$barcode = paste(contig_list_tcell$`timepoint-3`$barcode,"_3",sep = "")
contig_list_tcell$`timepoint-4`$barcode = paste(contig_list_tcell$`timepoint-4`$barcode,"_4",sep = "")
contig_list_tcell$`timepoint-5`$barcode = paste(contig_list_tcell$`timepoint-5`$barcode,"_5",sep = "")
contig_list_tcell$`timepoint-6`$barcode = paste(contig_list_tcell$`timepoint-6`$barcode,"_6",sep = "")

contig_list_bcell$`timepoint-1`$barcode = paste(contig_list_bcell$`timepoint-1`$barcode,"_1",sep = "")
contig_list_bcell$`timepoint-2`$barcode = paste(contig_list_bcell$`timepoint-2`$barcode,"_2",sep = "")
contig_list_bcell$`timepoint-3`$barcode = paste(contig_list_bcell$`timepoint-3`$barcode,"_3",sep = "")
contig_list_bcell$`timepoint-4`$barcode = paste(contig_list_bcell$`timepoint-4`$barcode,"_4",sep = "")
contig_list_bcell$`timepoint-5`$barcode = paste(contig_list_bcell$`timepoint-5`$barcode,"_5",sep = "")
contig_list_bcell$`timepoint-6`$barcode = paste(contig_list_bcell$`timepoint-6`$barcode,"_6",sep = "")

#flatten list into one df and delete extra column
df.tcr <- ldply (contig_list_tcell, data.frame)
df.bcr <- ldply (contig_list_bcell, data.frame)

#filter so we only have productive = true and consensus id != none
df.tcr <- filter(df.tcr, productive == "True" , raw_consensus_id != "None")
df.bcr <- filter(df.bcr, productive == "True",  raw_consensus_id != "None")

df.tcr$type <- 'tcell'
df.bcr$type <- 'bcell'

df.tcr <- select(df.tcr, -is_cell)
df.bcr <- select(df.bcr, -is_cell)




# merge dataframes
df.merge <-  rbind(df.tcr, df.bcr)


#check for duplicates

duplicates <- data.frame(table(df.merge$barcode))

#rownames(tcr) <- tcr[,1]




# uc david attempt

df.tcr.1<-add_clonotype("/home/matthew/datatransfer/contig/151231_tcell/",'t')
df.tcr.2<-add_clonotype("/home/matthew/datatransfer/contig/160411_tcell/",'t')
df.tcr.3<-add_clonotype("/home/matthew/datatransfer/contig/161962_tcell/",'t')
df.tcr.4<-add_clonotype("/home/matthew/datatransfer/contig/171094_tcell/",'t')
df.tcr.5<-add_clonotype("/home/matthew/datatransfer/contig/171642_tcell/",'t')
df.tcr.6<-add_clonotype("/home/matthew/datatransfer/contig/180251_tcell/",'t')

df.tcr.1$t_barcode = paste(df.tcr.1$t_barcode,"_1",sep = "")
df.tcr.2$t_barcode = paste(df.tcr.2$t_barcode,"_2",sep = "")
df.tcr.3$t_barcode = paste(df.tcr.3$t_barcode,"_3",sep = "")
df.tcr.4$t_barcode = paste(df.tcr.4$t_barcode,"_4",sep = "")
df.tcr.5$t_barcode = paste(df.tcr.5$t_barcode,"_5",sep = "")
df.tcr.6$t_barcode = paste(df.tcr.6$t_barcode,"_6",sep = "")


df.merge.tcr <-  rbind(df.tcr.1, df.tcr.2,df.tcr.3,df.tcr.4,df.tcr.5,df.tcr.6)
rownames(df.merge.tcr) <- df.merge.tcr[,1]

df.bcr.1<-add_clonotype("/home/matthew/datatransfer/contig/151231_bcell/",'b')
df.bcr.2<-add_clonotype("/home/matthew/datatransfer/contig/160411_bcell/",'b')
df.bcr.3<-add_clonotype("/home/matthew/datatransfer/contig/161962_bcell/",'b')
df.bcr.4<-add_clonotype("/home/matthew/datatransfer/contig/171094_bcell/",'b')
df.bcr.5<-add_clonotype("/home/matthew/datatransfer/contig/171642_bcell/",'b')
df.bcr.6<-add_clonotype("/home/matthew/datatransfer/contig/180251_bcell/",'b')

df.bcr.1$b_barcode = paste(df.bcr.1$b_barcode,"_1",sep = "")
df.bcr.2$b_barcode = paste(df.bcr.2$b_barcode,"_2",sep = "")
df.bcr.3$b_barcode = paste(df.bcr.3$b_barcode,"_3",sep = "")
df.bcr.4$b_barcode = paste(df.bcr.4$b_barcode,"_4",sep = "")
df.bcr.5$b_barcode = paste(df.bcr.5$b_barcode,"_5",sep = "")
df.bcr.6$b_barcode = paste(df.bcr.6$b_barcode,"_6",sep = "")


df.merge.bcr <-  rbind(df.bcr.1, df.bcr.2,df.bcr.3,df.bcr.4,df.bcr.5,df.bcr.6)
rownames(df.merge.bcr) <- df.merge.bcr[,1]

pbmc.combined.cellrep <- AddMetaData(object=pbmc.combined.cellrep, metadata=df.merge.tcr)
pbmc.combined.cellrep <- AddMetaData(object=pbmc.combined.cellrep, metadata=df.merge.bcr)


tcr.freq <- data.frame(table(pbmc.combined.cellrep@meta.data$t_clonotype_id))


#frequency t cell
tcr.freq<-table(pbmc.combined.cellrep@meta.data$t_clonotype_id,pbmc.combined.cellrep$Timepoint)

tcr.freq<-as.data.frame.matrix(tcr.freq)
tcr.freq <- cbind(tcr.freq = rownames(tcr.freq), tcr.freq)
rownames(tcr.freq) <- 1:nrow(tcr.freq)

tcr.freq <- tcr.freq %>%
  rename(
    Identity = tcr.freq)

tcr.freq1<-table(pbmc.combined.cellrep@meta.data$t_clonotype_id,Idents(pbmc.combined.cellrep))

tcr.freq1<-as.data.frame.matrix(tcr.freq1)
tcr.freq1 <- cbind(tcr.freq1 = rownames(tcr.freq1), tcr.freq1)
rownames(tcr.freq1) <- 1:nrow(tcr.freq1)

tcr.freq1 <- tcr.freq1 %>%
  rename(
    Identity = tcr.freq1)

tcr.freq1 <- tcr.freq1[, -c(2,3,6,10,11,12,15,16,17,21)] 

#### tcr per timpeoint per identity
tp1_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "1")
tp2_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "2")
tp3_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "3")
tp4_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "4")
tp5_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "5")
tp6_tcr <- subset(x = pbmc.combined.cellrep, subset = Timepoint == "6")

tcrTimePoint <- function(seuratobj) {
  tp1_tcr.freq<-table(seuratobj@meta.data$t_clonotype_id,Idents(seuratobj))
  
  tp1_tcr.freq<-as.data.frame.matrix(tp1_tcr.freq)
  tp1_tcr.freq <- cbind(tp1_tcr.freq = rownames(tp1_tcr.freq), tp1_tcr.freq)
  rownames(tp1_tcr.freq) <- 1:nrow(tp1_tcr.freq)
  
  tp1_tcr.freq <- tp1_tcr.freq %>%
    rename(
      Identity = tp1_tcr.freq)
  
  tp1_tcr.freq <- tp1_tcr.freq[, -c(2,3,6,10,11,12,15,16,17,21)] 
  
  return(tp1_tcr.freq)
  
}

tcr.freq.tp1 <-tcrTimePoint(tp1_tcr)
tcr.freq.tp2 <-tcrTimePoint(tp2_tcr)
tcr.freq.tp3 <-tcrTimePoint(tp3_tcr)
tcr.freq.tp4 <-tcrTimePoint(tp4_tcr)
tcr.freq.tp5 <-tcrTimePoint(tp5_tcr)
tcr.freq.tp6 <-tcrTimePoint(tp6_tcr)

write.csv(tcr.freq.tp1,"tcrtp1.csv")
write.csv(tcr.freq.tp2,"tcrtp2.csv")
write.csv(tcr.freq.tp3,"tcrtp3.csv")
write.csv(tcr.freq.tp4,"tcrtp4.csv")
write.csv(tcr.freq.tp5,"tcrtp5.csv")
write.csv(tcr.freq.tp6,"tcrtp6.csv")









####

# frequency b cell

bcr.freq <- data.frame(table(pbmc.combined.cellrep@meta.data$b_clonotype_id))
bcr.freq<-table(pbmc.combined.cellrep@meta.data$b_clonotype_id,pbmc.combined.cellrep$Timepoint)

bcr.freq<-as.data.frame.matrix(bcr.freq)
bcr.freq <- cbind(bcr.freq = rownames(bcr.freq), bcr.freq)
rownames(bcr.freq) <- 1:nrow(bcr.freq)

bcr.freq <- bcr.freq %>%
  rename(
    Identity = bcr.freq)

#clonotype test 1 plots
DimPlot(pbmc.combined.cellrep, group.by = "Timepoint")

clonotype1.subset <- subset(x = pbmc.combined.cellrep, subset = t_clonotype_id == "clonotype1")

timepoint_1.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "1")
timepoint_2.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "2")
timepoint_3.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "3")
timepoint_4.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "4")
timepoint_5.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "5")
timepoint_6.subset_tcr <- subset(x = clonotype1.subset, subset = Timepoint == "6")


DimPlot(timepoint_1.subset_tcr,reduction = "umap", label = TRUE ,repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 1") 
DimPlot(timepoint_2.subset_tcr,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 2") 
DimPlot(timepoint_3.subset_tcr,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 3") 
DimPlot(timepoint_4.subset_tcr,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 4") 
DimPlot(timepoint_5.subset_tcr,reduction = "umap" , label = TRUE, repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 5") 
DimPlot(timepoint_6.subset_tcr,reduction = "umap" , label = TRUE,  repel = TRUE, label.size = 3.2) +ggtitle("Timepoint 6") 








# this is scREP trial
#data
head(contig_list_tcell[[1]])
head(contig_list_bcell[[1]])

#prepare combinedTCR
combined_tcr <- combineTCR(contig_list_tcell,samples = c("timepoint-1", "timepoint-2", "timepoint-3", "timepoint-4", "timepoint-5","timepoint-6"),
                           ID = c("1", "2", "3", "4", "5", "6"),cells = c("T-AB"))

#remove prefix to match seurat obj
combined_tcr$`timepoint-1_1`$barcode <- gsub("timepoint-1_1_","",combined_tcr$`timepoint-1_1`$barcode)
combined_tcr$`timepoint-2_2`$barcode <- gsub("timepoint-2_2_","",combined_tcr$`timepoint-2_2`$barcode)
combined_tcr$`timepoint-3_3`$barcode <- gsub("timepoint-3_3_","",combined_tcr$`timepoint-3_3`$barcode)
combined_tcr$`timepoint-4_4`$barcode <- gsub("timepoint-4_4_","",combined_tcr$`timepoint-4_4`$barcode)
combined_tcr$`timepoint-5_5`$barcode <- gsub("timepoint-5_5_","",combined_tcr$`timepoint-5_5`$barcode)
combined_tcr$`timepoint-6_6`$barcode <- gsub("timepoint-6_6_","",combined_tcr$`timepoint-6_6`$barcode)

#add suffix to match seurat obj
combined_tcr$`timepoint-1_1`$barcode = paste(combined_tcr$`timepoint-1_1`$barcode,"_1",sep = "")
combined_tcr$`timepoint-2_2`$barcode = paste(combined_tcr$`timepoint-2_2`$barcode,"_2",sep = "")
combined_tcr$`timepoint-3_3`$barcode = paste(combined_tcr$`timepoint-3_3`$barcode,"_3",sep = "")
combined_tcr$`timepoint-4_4`$barcode = paste(combined_tcr$`timepoint-4_4`$barcode,"_4",sep = "")
combined_tcr$`timepoint-5_5`$barcode = paste(combined_tcr$`timepoint-5_5`$barcode,"_5",sep = "")
combined_tcr$`timepoint-6_6`$barcode = paste(combined_tcr$`timepoint-6_6`$barcode,"_6",sep = "")

pbmc.combined.cellrep <- combineExpression(combined_tcr, pbmc.combined.cellrep, groupBy = "sample")


# BCR

combined_bcr <- combineBCR(contig_list_bcell,samples = c("timepoint-1", "timepoint-2", "timepoint-3", "timepoint-4", "timepoint-5","timepoint-6"),
            ID = c("1", "2", "3", "4", "5", "6"))
#remove prefix to match seurat obj
combined_bcr$`timepoint-1_1`$barcode <- gsub("timepoint-1_1_","",combined_bcr$`timepoint-1_1`$barcode)
combined_bcr$`timepoint-2_2`$barcode <- gsub("timepoint-2_2_","",combined_bcr$`timepoint-2_2`$barcode)
combined_bcr$`timepoint-3_3`$barcode <- gsub("timepoint-3_3_","",combined_bcr$`timepoint-3_3`$barcode)
combined_bcr$`timepoint-4_4`$barcode <- gsub("timepoint-4_4_","",combined_bcr$`timepoint-4_4`$barcode)
combined_bcr$`timepoint-5_5`$barcode <- gsub("timepoint-5_5_","",combined_bcr$`timepoint-5_5`$barcode)
combined_bcr$`timepoint-6_6`$barcode <- gsub("timepoint-6_6_","",combined_bcr$`timepoint-6_6`$barcode)

#add suffix to match seurat obj
combined_bcr$`timepoint-1_1`$barcode = paste(combined_bcr$`timepoint-1_1`$barcode,"_1",sep = "")
combined_bcr$`timepoint-2_2`$barcode = paste(combined_bcr$`timepoint-2_2`$barcode,"_2",sep = "")
combined_bcr$`timepoint-3_3`$barcode = paste(combined_bcr$`timepoint-3_3`$barcode,"_3",sep = "")
combined_bcr$`timepoint-4_4`$barcode = paste(combined_bcr$`timepoint-4_4`$barcode,"_4",sep = "")
combined_bcr$`timepoint-5_5`$barcode = paste(combined_bcr$`timepoint-5_5`$barcode,"_5",sep = "")
combined_bcr$`timepoint-6_6`$barcode = paste(combined_bcr$`timepoint-6_6`$barcode,"_6",sep = "")

#combine both bcr and tcr into one list
combined_bcr_tcr <- append(combined_tcr, combined_bcr) 

#flatten list into one df and delete extra column
df <- ldply (combined_bcr_tcr, data.frame)
df <- select(df, -.id)


#find duplicates (this should not happen since tcr/bcr are exclusive BUT rarely one cell could have both )
n_occur <- data.frame(table(df$barcode))

write.csv(df,"df.csv", row.names = FALSE)


dups <- c('ACTTTCATCAAACCAC-1_5', 'CAAGTTGAGCTAAGAT-1_2', 'CACCAGGTCTCATTCA-1_2', 'CCTCAGTTCAACACCA-1_6', 'CCTTACGTCCTCTAGC-1_5', 'CGTCACTGTCACCTAA-1_4', 'CTCGTCAGTGTAACGG-1_6', 'GATGAGGCAGGGAGAG-1_2', 'GGGCACTTCCTTTCGG-1_3', 'GTAACTGCACGCATCG-1_6', 'TAGGCATCAGTACACT-1_5', 'TCACGAAGTCCCTTGT-1_6', 'TCTCATATCCACGTTC-1_5' , 'TTGGCAATCTTTCCTC-1_3')

duplicate.bcr.tcr <-filter(pbmc.combined.cellrep@meta.data, barcode %in% dups)

pbmc.combined.cellrep <- combineExpression(df, pbmc.combined.cellrep, groupBy = "sample")




## prepare tcr data for scanpy integration
timepoint1.tcell$barcode = paste(timepoint1.tcell$barcode,"_1",sep = "")
timepoint2.tcell$barcode = paste(timepoint2.tcell$barcode,"_2",sep = "")
timepoint3.tcell$barcode = paste(timepoint3.tcell$barcode,"_3",sep = "")
timepoint4.tcell$barcode = paste(timepoint4.tcell$barcode,"_4",sep = "")
timepoint5.tcell$barcode = paste(timepoint5.tcell$barcode,"_5",sep = "")
timepoint6.tcell$barcode = paste(timepoint6.tcell$barcode,"_6",sep = "")

write.csv(timepoint1.tcell,'timepoint1-tcell.csv')
write.csv(timepoint2.tcell,'timepoint2-tcell.csv')
write.csv(timepoint3.tcell,'timepoint3-tcell.csv')
write.csv(timepoint4.tcell,'timepoint4-tcell.csv')
write.csv(timepoint5.tcell,'timepoint5-tcell.csv')
write.csv(timepoint6.tcell,'timepoint6-tcell.csv')


t1<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "1")
t2<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "2")
t3<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "3")
t4<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "4")
t5<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "5")
t6<- subset(x = pbmc.combined.timepoint, subset = Timepoint == "6")



sceasy::convertFormat(t1, from="seurat", to="anndata",
                      outFile='t1.h5ad')
sceasy::convertFormat(t2, from="seurat", to="anndata",
                      outFile='t2.h5ad')
sceasy::convertFormat(t3, from="seurat", to="anndata",
                      outFile='t3.h5ad')
sceasy::convertFormat(t4, from="seurat", to="anndata",
                      outFile='t4.h5ad')
sceasy::convertFormat(t5, from="seurat", to="anndata",
                      outFile='t5.h5ad')
sceasy::convertFormat(t6, from="seurat", to="anndata",
                      outFile='t6.h5ad')


rm(t1)
rm(t2)
rm(t3)
rm(t4)
rm(t5)
rm(t6)






