


View(filter(tp1_tcr@meta.data, old.ident == "CD8 Exhausted" & !is.na(t_cdr3s_aa))$t_cdr3s_aa)


View(subset(tp1_tcr@meta.data, grepl(paste(filter(tp1_tcr@meta.data, old.ident == "CD8 Exhausted" & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse= "|"), 
                                     tp1_tcr@meta.data$t_cdr3s_aa)))

View(subset(subset(tp1_tcr@meta.data, grepl(paste(filter(tp1_tcr@meta.data, old.ident == "CD8 Exhausted" & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse= "|"), 
                                            tp1_tcr@meta.data$t_cdr3s_aa)), !(t_clonotype_id %in%  c("clonotype1","clonotype153","clonotype29","clonotype377","clonotype38","clonotype377","clonotype4","clonotype46",
                                                                                                     "clonotype480","clonotype483","clonotype257"))))

CloneTrack <- function(seurat.obj.1, seurat.obj.2,cell.1,fullList,plot,title) {

  if(plot == TRUE) {
   print( DimPlot(seurat.obj.2, label=F,  cells.highlight= filter(seurat.obj.2@meta.data,
           grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
           seurat.obj.2@meta.data$t_cdr3s_aa))$t_barcode, cols.highlight = c("darkblue", "darkred"), cols= "grey") +ggtitle(title))
  }

  if(fullList == TRUE) {

  return(filter(seurat.obj.2@meta.data,
          grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
          seurat.obj.2@meta.data$t_cdr3s_aa)))

  } else {

  return(filter(seurat.obj.2@meta.data,
         grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
        seurat.obj.2@meta.data$t_cdr3s_aa))$t_cdr3s_aa)
  }
}


# CloneTrack <- function(seurat.obj.1, seurat.obj.2,cell.1,fullList,plot,title) {
# 
#   if(plot == TRUE) {
#     print( DimPlot(seurat.obj.2, label=F,
#                    cells.highlight= filter(seurat.obj.2@meta.data, t_cdr3s_aa %in% filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa)$t_barcode
#                    , cols.highlight = c("darkblue", "darkred"), cols= "grey") +ggtitle(title))
#   }
# 
# 
#   if(fullList == TRUE) {
# 
#     return(filter(seurat.obj.2@meta.data, t_cdr3s_aa %in% filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa)
#     )
# 
#   } else {
# 
#     return(filter(seurat.obj.2@meta.data, t_cdr3s_aa %in% filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa)
#            $t_cdr3s_aa)
#   }
# }

CloneTrack <- function(seurat.obj.1, seurat.obj.2,cell.1,fullList,plot,title) {
  
  if(plot == TRUE) {
    
    filter(seurat.obj.2@meta.data,
           grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
                 seurat.obj.2@meta.data$t_cdr3s_aa))$t_barcode
    
    
    
    print( DimPlot(pbmc.combined.cellrep, label=F,  cells.highlight= , cols.highlight = c("darkblue", "darkred"), cols= "grey") +ggtitle(title))
  }
  
  if(fullList == TRUE) {
    
    return(filter(seurat.obj.2@meta.data,
                  grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
                        seurat.obj.2@meta.data$t_cdr3s_aa)))
    
  } else {
    
    return(filter(seurat.obj.2@meta.data,
                  grepl(paste(filter(seurat.obj.1@meta.data, old.ident == cell.1 & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
                        seurat.obj.2@meta.data$t_cdr3s_aa))$t_cdr3s_aa)
  }
}


matchedTCR<-filter(pbmc.combined.cellrep@meta.data,
       grepl(paste(filter(pbmc.combined.cellrep@meta.data, old.ident == "CD8 Exhausted" & !is.na(t_cdr3s_aa))$t_cdr3s_aa, collapse="|"),
             pbmc.combined.cellrep@meta.data$t_cdr3s_aa))

time1<- filter(matchedTCR,Timepoint=="1")
time2<-filter(matchedTCR,Timepoint=="2")
time3<-filter(matchedTCR,Timepoint=="3")
time4<-filter(matchedTCR,Timepoint=="4")
time5<-filter(matchedTCR,Timepoint=="5")
time6<-filter(matchedTCR,Timepoint=="6")


time1.join <- merge(time1,time2,by="t_cdr3s_aa")

DimPlot(pbmc.combined.cellrep, label=F,  cells.highlight= list(time1$t_barcode,time2$t_barcode), cols.highlight = c("darkblue", "darkred"), cols= "grey",split.by="Timepoint")






View(CloneTrack(tp1_tcr,tp1_tcr," CD8 t-cell 1", TRUE,TRUE,"test"))

View(CloneTrack(tp1_tcr,tp1_tcr,"CD4 Naive", TRUE,TRUE,"test"))

CloneTrack(tp1_tcr,tp2_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 1")
CloneTrack(tp2_tcr,tp2_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 2")
CloneTrack(tp3_tcr,tp3_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 3")
CloneTrack(tp4_tcr,tp4_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 4")
CloneTrack(tp5_tcr,tp5_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 5")
CloneTrack(tp6_tcr,tp6_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 6")


CloneTrack(tp1_tcr,tp2_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 1")




FindTimePoint2.5 <-function(logTimePoint, markers){
  test.df<-NULL
  
  for (i in 1:nrow(logTimePoint)) {
    
    if(
      ((logTimePoint[i,1] < logTimePoint[i,2])
       &&
       (logTimePoint[i,2] > logTimePoint[i,3])
       &&
       (logTimePoint[i,2] > logTimePoint[i,4])
       &&
       (logTimePoint[i,2] > logTimePoint[i,6])
      )&& 
      ((logTimePoint[i,1] < logTimePoint[i,5])
       &&
       (logTimePoint[i,5] > logTimePoint[i,3])
       &&
       (logTimePoint[i,5] > logTimePoint[i,4])
       &&
       (logTimePoint[i,5] > logTimePoint[i,6])
      )
    ) {
      if(length(filter(markers, gene == logTimePoint[i,7])$pct.1) == 0) {
        
      } else if (as.numeric(filter(markers, gene == logTimePoint[i,7])$pct.1) >.1){
        test.df<-rbind(test.df,logTimePoint[i,])
      }
    }
  }
  
  return (test.df)
}




View(filter(tp1_tcr@meta.data, t_cdr3s_aa %in% filter(tp1_tcr@meta.data, old.ident == "CD8 Exhausted" & !is.na(t_cdr3s_aa))$t_cdr3s_aa))

filter(tp1_tcr@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode


DimPlot(pbmc.combined.cellrep, label=F,  cells.highlight= filter(pbmc.combined.cellrep@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode,split.by="Timepoint")

DimPlot(pbmc.combined.cellrep, label=F,  cells.highlight= filter(pbmc.combined.cellrep@meta.data, t_cdr3s_aa == "TRB:CASSHSATGESYEQYF")$t_barcode,split.by="Timepoint")

DimPlot(pbmc.combined.cellrep, label=F,  cells.highlight= filter(pbmc.combined.cellrep@meta.data, t_cdr3s_aa == "TRB:CASSQGTGYTNTEAFF")$t_barcode,split.by="Timepoint")

DefaultAssay(tp2_tcr) <- "RNA"

FeaturePlot(tp1_tcr, features=c("TOX"),cells=filter(tp1_tcr@meta.data, t_cdr3s_aa == "TRB:CASSQGTGYTNTEAFF")$t_barcode,label=T)
FeaturePlot(tp6_tcr, features=c("TCF7"),cells=filter(tp6_tcr@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode,label=T)
FeaturePlot(tp2_tcr, features=c("GZMK"),cells=WhichCells(object = tp2_tcr, idents ="CD8 t-cell 1"),label=T)


FeaturePlot(pbmc.combined.cellrep, features=c("GZMK"),cells=filter(pbmc.combined.cellrep@meta.data, t_cdr3s_aa == "TRB:CASSLTGTTYNEQFF")$t_barcode,label=T,split.by="Timepoint")
FeaturePlot(subset(pbmc.combined.cellrep,idents ="CD8 t-cell 1"), features=c("GZMK"),label=T,split.by="Timepoint")


table(Idents(pbmc.combined.cellrep), pbmc.combined.cellrep$Timepoint)
table(pbmc.combined.cellrep$Timepoint)




