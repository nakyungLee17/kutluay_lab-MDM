library(dplyr)
library(ggplot2)
library(stringr)
library(edgeR)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)



######  arguments  ######

ex_name = "MDM-ribo"
filelist = list.files(path="ribo_counts", pattern='rp')
num_clust = 12
num_rep = 100
times = c("0h", "30m", "1h", "4h")

gtf = "/path/to/annotation"

col_order = c("file","names","in","desired","orer")

outdir = paste0(ex_name,"_outputs")
dir.create(path=outdir)


######  functions  ######

conv_to_rep <- function(explist) {
  exps <- unique(explist)
  for (i in 1:length(exps)) {
    explist <- gsub(exps[i], paste('rep',toString(i),sep=''), explist)
  }
  return(explist)
}

# match each experiment title to its donor
exp_to_donor <- function(exp) {
  if (grepl("rp4|rp9", exp)) return("d1")
  else if (grepl("rp3|rp6", exp)) return("d3")
  else if (grepl("rp2|rp5", exp)) return("d4")
  else {
    print(paste0("unidentified donor for: ", exp))
    return(exp)
  }
}

# convert timepoints from string to number
conv_to_num <- function(s) {
  if (s == "1h") return(1)
  else if (s == "30m") return(0.5)
  else if (s == "4h") return(4)
  else if (s == "0h") return(0)
  else {
    print(paste0("time point undefined for: ", s))
    return(NA)
  }
}

# fix mock's corresponding element in timepoint
getMock <- function(cond,time) {
  if (cond=='mock') return('0h')
  else return(time)
}

# convert an array to Gaussian distribution
gauss <- function(x) {
  (x - mean(x)) / sd(x) 
}


# make sample matrix - filename, condition, timepoint (in num), donor, experiment

# file names in format of: experiment_condition_timepoint e.g. rp2_IFN_1h

experiments = sapply(strsplit(filelist,'_'), function(x) x[1])
condition = sapply(strsplit(filelist,'_'), function(x) x[2])
timepoints = sapply(strsplit(filelist,'_'), function(x) x[3])
timepoints = unname(mapply(getMock, condition, timepoints))
timepoints = unname(mapply(conv_to_num, timepoints))
donors = unname(mapply(exp_to_donor, experiments))
  
design = data.frame(experiments = experiments,
                    condition = condition,
                    timepoints = timepoints,
                    donors = donors,
                    row.names = filelist)


# load counts data

counts = NULL

for (file in filelist) {
  table <- read.table(paste0("ribo_counts/",file), 
                      head=FALSE, sep='\t', col.names=c('id',file))
  counts <- cbind(counts, table[,2])
}

rownames(counts) = table[ ,1]
colnames(counts) = filelist


# make DGE, convert ID, filter & normalize [B]

dge = DGEList(counts = counts)
rownames(dge$counts) = gsub('.{1}$', '', read.table(gtf, sep = " ")[ ,6])

keep = rowSums(cpm(dge) > 1) >= 2
dge = dge[keep, ,keep.lib.sizes = F]

dups = rownames(dge)[duplicated(rownames(dge))]
dge = dge[!(row.names(dge) %in% dups), ,keep.lib.sizes = F]

dge = calcNormFactors(dge)
cpm = cpm(dge, log = T, prior.count = 1)

if (save_tables) write.table(cpm, file = paste0(outdir, "/cpm"), quote = F, sep = "\t")

# load diffex genes data & E. extract union of sig genes from [D]

de_table = NULL
de_union = NULL

for (file in dir(path = "ribo_diffs", pattern = "*full")) {
  # read in table
  fname = paste0(unlist(strsplit(file,"_"))[2], "_", unlist(strsplit(file,"_"))[3])
  table = read.table(paste0("ribo_diffs/",file))
  # save the sig genes
  sig = rownames(table[abs(table[, "logFC"]) > 1 & table[, "FDR"] < 0.05, ]) 
  de_union = union(de_union, sig)
  # merge tables
  colnames(table) = paste(fname,colnames(table), sep=".")
  table[ ,"genes"] = rownames(table)
  if (is.null(de_table)) de_table = table
  else de_table = dplyr::full_join(de_table, table, by="genes")
}

# tidy up the combined de genes table
rownames(de_table) = de_table[ ,"genes"]
de_table = de_table[ ,!grepl("genes",colnames(de_table))]

# take intersect and find corresponding cpm values

union_cpm = cpm[intersect(de_union, rownames(cpm)), ]

# reorder and transform

union_cpm = union_cpm[ ,col_order]
g_table = t(apply(union_cpm,1,gauss))


# create cluster: pre-select num of rep & clusters

cluster <- ConsensusClusterPlus(t(g_table),
                                maxK = num_clust, reps = num_rep,
                                pItem = 0.85, pFeature = 1,
                                clusterAlg = "km",
                                # clusterAlg = "hc",
                                # innerLinkage = "complete",
                                # finalLinkage = "ward.D2",
                                distance = "euclidean",
                                plot = "pdf",
                                title = outdir,
                                seed = 3.1415926535
)

clust_info = cluster[[num_clust]][["consensusClass"]]

# update cluster info to g_table
g_table = data.frame(g_table)
rm(cluster)


# make fc table 

fc_table = NULL
temp = NULL

for (d in c("donor1","donor3","donor4")) {
  for (i in 2:length(times)) {
    # print(paste0(d,"_",times[i],".logFC"))
    temp = cbind(names(clust_info), rep(conv_to_num(times[i]), length(clust_info)),
                 de_table[names(clust_info), 
                          grepl(paste0(d,"_",times[i],".logFC"), colnames(de_table))],
                 rep(d, length(clust_info)), clust_info)
    # print(paste0("new temp table: ", dim(temp)))
    temp = na.omit(temp)
    # print(paste0("table after na filter: ", dim(temp)))
    fc_table = rbind(fc_table, temp)
    # print(paste0("fc table after rbind: ", dim(fc_table)))
  }
}

colnames(fc_table) = c("gene","time","logFC","donor","cluster")
fc_table = as.data.frame(fc_table)
fc_table[, c(2,3,5)] <- sapply(fc_table[, c(2,3,5)], as.numeric)

summary(fc_table)


# plot cluster 

colPalette = c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F", "#D627287F",
                "#FF7F0E7F", "#9467BD7F", "#8C564B7F", "#E377C27F",
                "#F7B6D27F")


# clusters to remove > should be marked 99
dc = c(2, 5, 6, 9, 11, 12)

pos.1 <- which(clust_info == 1)
pos.2 <- which(clust_info == 2)
pos.3 <- which(clust_info == 3)
pos.4 <- which(clust_info == 4)
pos.5 <- which(clust_info == 5)
pos.6 <- which(clust_info == 6)
pos.7 <- which(clust_info == 7)
pos.8 <- which(clust_info == 8)
pos.9 <- which(clust_info == 9)
pos.10 <- which(clust_info == 10)
pos.11 <- which(clust_info == 11)
pos.12 <- which(clust_info == 12)

clust_info[pos.1] <- 1
clust_info[pos.2] <- 99
clust_info[pos.3] <- 2
clust_info[pos.4] <- 3
clust_info[pos.5] <- 99
clust_info[pos.6] <- 99
clust_info[pos.7] <- 3
clust_info[pos.8] <- 2
clust_info[pos.9] <- 99
clust_info[pos.10] <- 3
clust_info[pos.11] <- 99
clust_info[pos.12] <- 99

clust_info = clust_info[clust_info!=99]
g_table_cut = g_table[(intersect(rownames(g_table),
                                 names(clust_info[clust_info!=99]))),]
num.clusters <- 3


if (save_tables) write.table(g_table_cut, file = paste0(outdir, "/heatmap_table"), quote = F, sep = "\t")


# plot heatmap 

cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = colPalette, col = NA)),
                                  width = unit(1, "mm"))

h2 <- Heatmap(as.matrix(g_table_cut),
                   column_title = paste0("Differentially Expressed Genes\nn = ",
                                         nrow(g_table_cut)),
                   column_title_gp = gpar(fontsize = 7),
                   show_row_names = F,
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5), # 6 for HBEC, 7 for Vero
                   name = "z-score",
                   col = colorRamp2(breaks = seq(-3, 3, length = 256), 
                                    colors = rev(colorRampPalette(brewer.pal(10,
                                                                             "RdBu"))(256))),
                   # legends
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(color_bar = "continuous",
                                               title_gp = gpar(fontsize = 6),
                                               labels_gp = gpar(fontsize = 6),
                                               grid_width = unit(2, units = "mm")),
                   # clustering
                   cluster_columns = F,
                   clustering_distance_rows = "euclidean",
                   cluster_row_slices = F,
                   show_row_dend = F,
                   split = clust_info,
                   left_annotation = cluster.rowAnnot,
                   # splitting
                   column_split = factor(split_order),
                   cluster_column_slices = F,
                   column_gap = unit(2, units = "mm"),
                   # labels
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 7),
                   # row_title = NULL,
                   # size
                   width = ncol(g_table_cut) * 0.2
)


png(paste0(outdir, "/heatmap.png"), height = 2000, width = 1000, res = 300)
h2
dev.off()


change_cluster <- function(x) {
  if (x == 1) return(1)
  else if (x %in% c(3,8)) return(2)
  else if (x %in% c(4,7,10)) return(3)
  else return(0)
}

fc_table_cut = fc_table[!(fc_table$cluster %in% dc), ]
fc_table_cut$cluster = mapply(change_cluster, fc_table_cut$cluster)

sumLabel = data.frame(table(clust_info))
clustLabel = paste0(sumLabel$clust_info, " (n=", sumLabel$Freq, ") ")
names(clustLabel) = 1:dim(sumLabel)[1]

muProfile_cut <- fc_table_cut %>%
  group_by(time, cluster) %>%
  summarise(logFC = mean(logFC))

p2 = (ggplot(as.data.frame(fc_table_cut), aes(x = as.numeric(time),
                           y = as.numeric(logFC),
                           group = interaction(gene, donor, cluster),
                           color = as.factor(cluster))) 
        + geom_line(size = 1, alpha=0.2)
        + facet_wrap(. ~ cluster,
                     ncol = 3,
                     scales = "free",
                     labeller = labeller(cluster = clustLabel))
        + geom_line(size = 1, data = muProfile_cut,
                    aes(x = time, y = logFC, 
                        group = as.factor(cluster)), color = "black")
        + geom_point(size = 1.5, data = muProfile_cut,
                     aes(x = time, y = logFC, group = NULL), color = "black")
        + scale_y_continuous(name = "log2(fold-change over mock)")
        + scale_x_continuous(name = "hours after treatment",
                             breaks = c(0, 1, 2, 3, 4),
                             expand = c(0.05, 0))
      # + scale_y_continuous(limits = c(-5, 5))
        + theme(legend.position = "none",
                axis.title.x = element_text(size=8),
                axis.title.y = element_text(size=8),
                strip.background = element_blank())
        )

png(paste0(outdir, "/profile.png"), width = 1500, height = 600, res = 300)
p2
dev.off()
