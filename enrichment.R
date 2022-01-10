library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)


setwd("dir/to/Gorilla_outputs")


# gorilla = read.table(file = "clusters/full_GOrilla.txt.csv", sep=",", header=T)
# gorilla$count = as.numeric(lapply(gorilla[,5], 
#                        function(x) gsub('.{1}$','',unlist(strsplit(x, split=","))[4])))
# gorilla$enrichment = as.numeric(lapply(gorilla[,5],
#                                        function(x) unlist(strsplit(x, split=" "))[1]))

# colnames(gorilla)[1] = "GO_ID"
# #for every cluster, take top 10 genes with highest enrichment score

# save_by_cluster <- function(table,c) {
#   print(c)
#   write.table(table[table$cluster == c, ], 
#               file=paste0("enrichment_heatmap/by_clusters/gorilla_",c), 
#               quote=F, sep="\t", row.names=F)
#   clust = table[(table$cluster == c & table$enrichment > 10), ]
  
#   # clust = clust[order(-clust$enrichment), ]
#   clust = clust[order(clust$qvalue), ]
#   write.table(clust, file=paste0("enrichment_heatmap/by_clusters/gorilla_",c,"_filtered"), 
#               quote=F, sep="\t", row.names=F)
#   # print(clust[1:3,])
#   print(dim(clust))
#   return(clust[1:30,])
# }

# gorilla.top.N = lapply(1:6, function(x) save_by_cluster(gorilla,x))
# gorilla.filtered = do.call(rbind, gorilla.top.N)
# gorilla.filtered = na.omit(gorilla.filtered)

# write.table(gorilla.filtered, file="enrichment_heatmap/by_clusters/gorilla_qval",
#             quote=F, sep="\t", row.names=F)



gorilla.selected = read.table(file = "enrichment_heatmap/GOrilla_select",
                              header = T, sep = "\t")
selected.list = unique(gorilla.selected$GO_ID)
gorilla.expanded = subset(gorilla, gorilla$GO_ID %in% selected.list)


# enrich_result = gorilla.filtered
enrich_result = gorilla.expanded


## tidy up the enricher result

## appropriately edit num in cluster
num.in.cluster <- (gene_cluster %>%
                     group_by(cluster) %>% 
                     summarise(n=n()))
num.in.cluster$cluster = as.integer(num.in.cluster$cluster)
save_cluster = as.integer(levels(factor(enrich_result$cluster)))
num.in.cluster = num.in.cluster[match(save_cluster, num.in.cluster$cluster),]
num.clusters = dim(num.in.cluster)[1]


## conver to matrix form
enrichment.results_wide = pivot_wider(enrich_result, 
                                      id_cols = c(Description, ontology), 
                                      names_from = cluster, 
                                      values_from=c(pvalue, qvalue, count, 
                                                    GO_ID, enrichment))


enrichment.results_wide <- arrange(enrichment.results_wide, ontology)
em.q.values <- select(enrichment.results_wide, grep("qvalue", colnames(enrichment.results_wide), value=T))
colnames(em.q.values) <- seq(1:num_cluster)
ratio.values <- select(enrichment.results_wide, grep("enrichment", colnames(enrichment.results_wide), value=T))
em.q.values <- -log10(as.matrix(em.q.values))
colnames(em.q.values) = save_cluster
max.q = max(em.q.values, na.rm=T)
min.q = min(em.q.values, na.rm=T)
col_fun = colorRamp2(c(min.q, max.q-25, max.q), c("#300101", "#e41a1c", "#e41a1c"))


## plotting

col.annot <- HeatmapAnnotation(cluster = anno_simple(colnames(em.q.values),
                                                     pt_gp = gpar(fontsize = 6),
                                                     height = unit(1, "mm"),
                                                     col = structure(brewer.pal(num.clusters, "Set3"), 
                                                                     names = colnames(em.q.values))),
                               show_annotation_name = F,
                               show_legend = F,
                               annotation_name_gp = gpar(fontsize = 6))


h <- Heatmap(em.q.values, 
             # labels
             column_title = "Enriched gene sets",
             column_title_gp = gpar(fontsize=9),
             row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 6),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize=7),
             row_labels = enrichment.results_wide$Description,
             # legends
             bottom_annotation = col.annot,
             show_heatmap_legend = F,
             
             col = col_fun,
             heatmap_legend_param = list(col_fun = col_fun,
                                         title = "-log10(q-value)",
                                         color_bar = "continuous",
                                         title_gp = gpar(fontsize = 6),
                                         labels_gp = gpar(fontsize = 6),
                                         grid_height = unit(2,"mm"),
                                         direction = "horizontal"),
             
             cluster_rows = F,
             cluster_columns = F,
             show_row_names = T,
             show_column_names = T,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(fill = "white", col = "#EEEEEE"))
               grid.circle(x = x, y = y, r = sqrt(ratio.values[i, j]) * 0.012,
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }
)

## Q-value legend ##
q_legend <- Legend(col_fun = col_fun,
                   title = "-log10(q-value)",
                   title_gp = gpar(fontsize = 7),
                   labels_gp = gpar(fontsize = 6),
                   grid_height = unit(2, "mm"),
                   direction = "horizontal")



## Ratio value legend ##
# ratio_breaks <- c(0.5, 0.4, 0.3, 0.2, 0.1)
ratio_breaks <- c(0.1, 0.08, 0.06, 0.04 ,0.02)

ratio_legend <- Legend(labels = ratio_breaks,
                       grid_height = unit(5, "mm"),
                       grid_width = unit(5, "mm"),
                       title = "Percent\nof Cluster",
                       title_gp = gpar(fontsize = 7),
                       labels_gp = gpar(fontsize = 6),
                       type = "points",
                       pch = 1, # circle
                       size = unit(sqrt(ratio_breaks)*1.5, "npc"),# <<< not the same as in cell_fun, change for consistency
                       background = "white"
)



png("enrichment_heatmap/gorilla_expanded.png", height = 1500, width = 1200, res = 300)
draw(h, padding = unit(c(8, 8, 8, 12), "mm")) # padding: empty space around heatmap, for legends, row names...
draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.75, "npc"))
dev.off() 
