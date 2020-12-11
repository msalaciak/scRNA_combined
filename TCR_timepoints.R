


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







View(CloneTrack(tp1_tcr,tp1_tcr," CD8 t-cell 1", TRUE,TRUE,"test"))

View(CloneTrack(tp2_tcr,tp2_tcr,"CD8 TEM", TRUE,TRUE,"test"))

CloneTrack(tp1_tcr,tp1_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 1")
CloneTrack(tp2_tcr,tp2_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 2")
CloneTrack(tp3_tcr,tp3_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 3")
CloneTrack(tp4_tcr,tp4_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 4")
CloneTrack(tp5_tcr,tp5_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 5")
CloneTrack(tp6_tcr,tp6_tcr,"CD8 Exhausted", TRUE,TRUE,"Timepoint 6")







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




