pub<-DimPlot(timepoint_1.subset,reduction = "umap",label = TRUE ,repel = TRUE, label.size = 2.9) + theme_classic() + NoLegend() 


ggsave("pub.png",width = 5, height = 5,dpi = 300)




combo<-(t1+t2)
ggsave("combo.png",width = 9, height = 5,dpi = 300)


# cell_proportion_per <- read.csv(file = 'pbmc_combined_percent_format2.csv')
cell_proportion_per <- read.csv(file = 'cd8.csv')

# table1<- gt(data = cell_proportion_per) %>% tab_header(
#   title = "Amount of Cells Per Timepoint"
#   
# )%>%
#   tab_stubhead(label = "Cell Identity") %>%
#   cols_label(
#     N.Cells = html("N Cells"),
#     X..Cells = html("% Cells"),
#     N.Cells.1 = html("N Cells"),
#     X..Cells.1 = html("% Cells"),
#     N.Cells.2 = html("N Cells"),
#     X..Cells.2 = html("% Cells"),
#     N.Cells.3 = html("N Cells"),
#     X..Cells.3 = html("% Cells"),
#     N.Cells.4 = html("N Cells"),
#     X..Cells.4 = html("% Cells"),
#     N.Cells.5 = html("N Cells"),
#     X..Cells.5 = html("% Cells"),
#   ) %>%
#   tab_spanner(
#     label = "1",  
#     columns = vars(N.Cells, X..Cells),
#     
#   )  %>%
#   tab_spanner(
#     label = "3" ,
#     columns = vars(N.Cells.2, X..Cells.2)
#     
#   ) %>%
#   tab_spanner(
#     label = "4", 
#     columns = vars(N.Cells.3, X..Cells.3)
#     
#   )  %>%
#   tab_spanner(
#     label = "5",
#     columns = vars(N.Cells.4, X..Cells.4)
#     
#   ) %>%
#   tab_spanner(
#     label = "6",
#     columns = vars(N.Cells.5, X..Cells.5)
#     
#   )  %>%
#   tab_spanner(
#     label = "2",
#     columns = vars(N.Cells.1, X..Cells.1)
#     
#   )  %>%
#   tab_row_group(
#     group = "PBMC with T Cell Subsets",
#     rows = 3:15) %>%
#   tab_row_group(
#     group = "T Cell Summary",
#     rows = 1:2)%>%
#   tab_style(
#     style = list(
#       cell_fill(color = "#ACEACE"),
#       cell_text(weight = "bold")
#     ),
#     locations = cells_body(
# 
#       rows = c(2)
#       )
#   )%>%
#   tab_style(
#     style = list(
#       
#       cell_text(weight = "bold")
#     ),
#     locations = cells_body(
#       
#       rows = 16)
#   )
#     
# 
# 
# table1 %>% gt_theme_538()
# 
# table1 %>% gt_theme_538() %>%
#   gtsave(
#     "tab_1_new.png",expand=10,vwidth = 1628,
#     vheight = 882
#     
#   )

table1<- gt(data = cell_proportion_per) %>% tab_header(
  title = "Amount of Cells Per Timepoint"
  
)%>%
  tab_stubhead(label = "Cell Identity") %>%
  cols_label(
    N.Cells = html("N Cells"),
    X..Cells = html("% Cells"),
    N.Cells.1 = html("N Cells"),
    X..Cells.1 = html("% Cells"),
    N.Cells.2 = html("N Cells"),
    X..Cells.2 = html("% Cells"),
    N.Cells.3 = html("N Cells"),
    X..Cells.3 = html("% Cells"),
    N.Cells.4 = html("N Cells"),
    X..Cells.4 = html("% Cells"),
    N.Cells.5 = html("N Cells"),
    X..Cells.5 = html("% Cells"),
  ) %>%
  tab_spanner(
    label = "1",  
    columns = vars(N.Cells, X..Cells),
    
  )  %>%
  tab_spanner(
    label = "3" ,
    columns = vars(N.Cells.2, X..Cells.2)
    
  ) %>%
  tab_spanner(
    label = "4", 
    columns = vars(N.Cells.3, X..Cells.3)
    
  )  %>%
  tab_spanner(
    label = "5",
    columns = vars(N.Cells.4, X..Cells.4)
    
  ) %>%
  tab_spanner(
    label = "6",
    columns = vars(N.Cells.5, X..Cells.5)
    
  )  %>%
  tab_spanner(
    label = "2",
    columns = vars(N.Cells.1, X..Cells.1)
    
  )  %>%
  tab_row_group(
    group = "T Cell Subsets",
    rows = 1:12) %>%

  tab_style(
    style = list(
      
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      
      rows = 13)
  )



table1 %>% gt_theme_538()

table1 %>% gt_theme_538() %>%
  gtsave(
    "tab_1_new_cd8.png",expand=10,vwidth = 1628,
    vheight = 882
    
  )





gt_theme_538 <- function(data,...) {
  data %>%
    opt_all_caps()  %>%
    opt_table_font(
      font = list(
        google_font("Myriad Pro")
      )
    ) %>%
    tab_style(
      style = cell_borders(
        sides = "bottom", color = "transparent", weight = px(2)
      ),
      locations = cells_body(
        columns = TRUE,
        # This is a relatively sneaky way of changing the bottom border
        # Regardless of data size
        rows = nrow(data$`_data`)
      )
    )  %>% 
    tab_options(
      column_labels.background.color = "white",
      table.border.top.width = px(3),
      table.border.top.color = "transparent",
      table.border.bottom.color = "transparent",
      table.border.bottom.width = px(3),
      column_labels.border.top.width = px(3),
      column_labels.border.top.color = "transparent",
      column_labels.border.bottom.width = px(3),
      column_labels.border.bottom.color = "black",
      data_row.padding = px(3),
      source_notes.font.size = 12,
      table.font.size = 16,
      heading.align = "left",
      grand_summary_row.background.color = "#990000",
      
      ...
    ) 
}



tcrtime <- read.csv(file = 'tptcr.csv')

tcrtime <- melt(tcrtime, id.vars=c("Timepoint..clonotype"))

tcrtime[,2] <-gsub("[A-z].","",tcrtime[,2])





ggplot(tcrtime, aes(x = variable, y = value, color = Timepoint..clonotype, group = Timepoint..clonotype)) + 
  geom_point()+geom_line() +facet_wrap(vars(Timepoint..clonotype),scales = "free_y") +theme(legend.position="none",plot.title = element_text(face = "plain")) + ggtitle("Frequency of T Cell Clones Across 6 Timepoints") +
xlab("Timepoints") +
ylab("Frequency of Clones")

ggsave("tcrtime.png",width = 16, height = 5,dpi = 300)
