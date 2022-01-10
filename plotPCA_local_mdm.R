# ---
# title: "PCA"
# author: "NaKyung Lee"
# date: "11/16/2021"
# output: PNG file
# ---

library(RColorBrewer)
library(ggplot2)
#library(ggrepel)
library(matrixStats)
library(edgeR)


# # Set directory to the folder with counts files
# # filename format should be: replicate_condition_timepoint
# # replicate = rep1, rep2, rep3...
# # condition = IFN or mock
# # timepoint = 30m, 1h, or 4h

# # each count file contains two columns: gene id and count (with no header)

setwd('/path/to/count_files/') # with slash in the end
expname <- 'experiment_name' # prefix for naming the output png
output_dir <- "/path/to/save/outputs/" # with slash, output directory

### LOADING DATA
filelist <- list.files()
print(filelist)
data <- read.table(filelist[1], header = FALSE, sep = '\t')
row.names(data) <- data[,1]
data <- data[, -c(1, 2)]

for (file in filelist) {
  table <- read.table(file, header = FALSE, sep = '\t', col.names = c('id', file))
  data <- cbind(data, table[, 2])
}

colnames(data) <- filelist


### SETTING UP DESIGN MATRIX
conv_to_rep <- function(explist) { # converts the experiment names to "repn"
  exps <- unique(explist)
  for (i in seq_along(exps)) {
    explist <- gsub(exps[i], paste0('rep', toString(i)), explist)
  }
  return(explist)
}

saveMock <- function(cond,time) { # saves the time points of the mock samples
  if (cond == 'mock') return(cond) # saves as "mock"
  else return(time)
}

experiments <- sapply(strsplit(filelist,'_'), function(x) x[1])
experiments <- conv_to_rep(experiments) # use this line if need to convert experiment names to repN
print(experiments)
condition <- sapply(strsplit(filelist,'_'), function(x) x[2])
print(condition)
timepoints <- sapply(strsplit(filelist,'_'), function(x) x[3])
timepoints <- unname(mapply(saveMock, condition, timepoints))
print(timepoints)

design <- data.frame(experiments, condition, timepoints)
rownames(design) <- filelist



### PREPARING FOR PCA PLOT

## Data Normalization
data <- DGEList(counts = data, group = timepoints)
print(paste0("before filtering: ", nrow(data$counts)))

keep <- filterByExpr(data)
data <- data[keep, , keep.lib.sizes = F]
print(paste0("after filtering: ", nrow(data$counts)))

data <- calcNormFactors(data)

data <- cpm(data)


## selecting genes with cpm > 1 across all samples

topNdata <- data[rowSums(data > 1) == ncol(data), ]
print(nrow(topNdata))


## transforming the normalized counts into log10 scale

topNdata <- log10(topNdata)


## finding var_perc for PC1 and PC2

data_t <- as.data.frame(t(topNdata))

pca <- prcomp(data_t, center = T, scale = F) # main PCA calculation step
vars <- pca$sdev ^ 2
var_perc <- sprintf("%.2f", round((vars / sum(vars)) * 100, 2))

df <- cbind(data.frame(pca$x), design)



### PLOTTING PCA

# Custom color palette
colPalette <- c("#AEC7E8", "#98DF8A", "#FF9896", "#C49C94", "#C5B0D5", "#FFBB78", "#D62728",
                "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#1F77B4", "#2CA02C", "#F7B6D2")

## setting the x-axis and y-axis range
x_max <- max(df$PC1)
x_min <- min(df$PC1)
y_max <- max(df$PC2)
y_min <- min(df$PC2)
y.range <- y_max - y_min
x.range <- x_max - x_min

df$timepoints <- factor(df$timepoints, levels = c("mock", "30m", "1h", "4h")) # to specify the ordering of timepoint legends for a particular plot

p <- (ggplot(df, aes(x = PC1, y = PC2, shape = as.factor(experiments), color = timepoints)) # PC1 vs PC2
        + geom_point(size = 4)

        + scale_x_continuous(limits = c(x_min - x.range / 10, x_max + x.range / 10))
        + scale_y_continuous(limits = c(y_min - y.range / 10, y_max + y.range / 10))
        + scale_shape_discrete(name = "Experiment")
        # + scale_shape_manual(name = "Experiment",
        #                      values = c(19, 17, 15, 3, 7, 8, 5, 11, 10)) # custom shapes
        + scale_color_manual(values = colPalette, name = "Timepoint")

        + xlab(paste0("PC1: ", var_perc[1], "% variance"))
        + ylab(paste0("PC2: ", var_perc[2], "% variance"))

        + ggtitle(label = "Principal Component Analysis")
        + theme_bw()
        + theme(legend.position = 'right',
                panel.border = element_rect(colour = "black", fill = NA, size = 1),
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5),
                plot.caption = element_text(size = 10),
                aspect.ratio = 1)
)



### SAVING THE PLOT
png(paste0(output_dir, paste(expname, "PCA_PC1_v_PC2.png", sep = '_')),2000,1400, res = 300)
print(p)
dev.off()
