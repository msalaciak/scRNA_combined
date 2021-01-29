library(SeuratWrappers)
library(monocle3)


cds <- as.cell_data_set(cd8.subset)
plot_cells(cds)

cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)


cds <- learn_graph(cds,use_partition = F)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)

colnames(subset(x = cd8.subset, subset = Timepoint == "1"))

cds <- order_cells(cds, root_cells = colnames(subset(x = cd8.subset, subset = Timepoint == "1")))

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)




integrated.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(integrated.sub, "monocle3_pseudotime")

cds <- order_cells(cds, root_cells = colnames(subset(x = cd8.subset, ident=c("3"))))

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

integrated.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(integrated.sub, "monocle3_pseudotime")


