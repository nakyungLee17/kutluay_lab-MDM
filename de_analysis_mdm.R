library(edgeR)
library(RColorBrewer)
library(gplots)
library(ComplexHeatmap)
library(ggplot2)


############################################################################

# # Note that the order of count files specified in files, donor_list, and conditions has to be in the order of
# # "control, treatment, control, treatment, control, treatment..."

# path to count files folder
directory <- "/path/to/count_files" # no slash!

# names of count files being used for edgeR analysis
files <- "rep1_mock,rep1_IFN,rep2_mock,rep2_IFN,rep3_mock,rep3_IFN"

# replicate information of count files used
donor_list <- "rep1,rep1,rep2,rep2,rep3,rep3"

# experimental condition info of count files used
conditions <- "mock,IFN,mock,IFN,mock,IFN"

# output prefix
prefix <- "prefix_of_choice"

# gtf file for gene name conversion from ensembl id to gene symbols #
# this gtf file needs to be gene-only and have the same number of genes as the count files;
# this can be achieved by only selecting rows which have "gene" as the third field (instead of "exon" or anything else)
# the genes in this gtf also need to be in the same order as the genes in count files.
gtf_path <- "path/to/gene_only_gtf_file"


############################################################################

setwd(directory)
labels <- paste(donor_list, conditions, sep = "_")
print(labels)
DG <- readDGE(files, header = FALSE, labels = labels)


# Gene ID Conversion #

# Gene ID conversion using gtf used in mapping. The total number and order of genes in the gtf should be exactly the same as in the individual count files.
# Thus, can directly replace gene ids in DG with list of gene names in the gtf.
if (grepl(pattern = "ENSG", rownames(DG$counts)[1])) {
	gtf <- read.table(gtf_path, sep = " ")
	stopifnot(nrow(gtf) == nrow(DG))
	gene_names <- as.character(gtf[, 6]) # <<< gene name info; may change column index according to gtf format
	gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
	rownames(DG) <- gene_names
}



# Filtering #

keep <- rowSums(cpm(DG) > 1) >= length(unique(donor_list)) # filtering by expression level
print(paste0("before filtering: ", nrow(DG)), quote = F)
DG <- DG[keep, ,keep.lib.sizes = FALSE]
print(paste0("after filtering: ", nrow(DG)), quote = F)


# Deduplicating # added 7/22/2021

duplicate_genes <- row.names(DG)[duplicated(row.names(DG))]
DG <- DG[!(row.names(DG) %in% duplicate_genes), ,keep.lib.sizes = F]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
			 " duplicated genes removed."), quote = F)



## Normalization ##

DG <- calcNormFactors(DG)
DGgroups <- conditions #conditions is a comma separated list of conditions
DG <- DGEList(counts = DG, group = factor(DGgroups))
DG <- calcNormFactors(DG)
Donor <- factor(donor_list) #donor_list is list of the donor of each file;
Condition <- factor(conditions) #conditions is same as conditions for DGgroups;



## Dispersion Estimation, BCV Plot, and Exact Test ##

data.frame(Sample=labels, Donor, Condition) # what is this line for?
design <- model.matrix(~ Donor + Condition)
DG <- estimateDisp(DG, design, robust=TRUE)


bcv_file_name <- paste0(prefix, "_BCVplot.pdf")
pdf(bcv_file_name) #name of file to write bcv graph to;
plotBCV(DG)
dev.off()


de <- exactTest(DG, pair = unique(conditions))
de.genes <- rownames(topTags(de)$table)
diffex_file <- paste0(prefix, "_diffexgenes")
diffex_file_full <- paste0(prefix, "_diffexgenes_full")
write.table(de.genes, diffex_file, quote=FALSE, row.names=FALSE, col.names=TRUE) #file to write de.genes to
write.table(topTags(de)$table, diffex_file_full,
			quote=FALSE, row.names=TRUE, col.names=TRUE, sep='\t')

