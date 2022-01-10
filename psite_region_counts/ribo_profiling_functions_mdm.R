# Functions for ribosome profiling analysis

make_output_directories <- function(outputDir) {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "objs"))) {
    dir.create(file.path(outputDir, "objs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "figs"))) {
    dir.create(file.path(outputDir, "reports", "figs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "de_genes"))) {
    dir.create(file.path(outputDir, "reports", "de_genes"), recursive = T)
  }
  
}

make_dge <- function(counts_path, countFilter, region="cds") {
  # function to make dge object for Ribo-seq or RNA-seq and output the counts table and CPM table
  # Args: 
  #   counts_path: path for the output count tables from htseq-count
  #   countFilter: common pattern of the htseq file names in the data type
  #   region: cds, utr3, or utr5
  # Returns:
  #   null. But will save the dge object 
  
  files <- dir(path = counts_path, pattern = "*.count_reduced$")
  
  ## sort the files in asc order
  files <- sort(files)
  counts <- readDGE(files, path = counts_path, header = F)$counts

  print(noquote("counts:"))
  print(dim(counts))

  # get gene information from ensembl
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 75", dataset = "hsapiens_gene_ensembl", host = "http://feb2014.archive.ensembl.org") # human

  stopifnot(length(unique(rownames(counts))) == nrow(counts))
  
  gtf.ens <- getBM(attributes=c('ensembl_gene_id','external_gene_id',"gene_biotype"), filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl) # human

  saveRDS(gtf.ens, file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  
  genes <- readRDS(file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  print(noquote("genes:"))
  print(dim(genes))
  
  genes <- genes[match(rownames(counts), genes$ensembl_gene_id),]

  # added lines ##############################
  genes <- na.omit(genes)
  counts <- counts[genes$ensembl_gene_id, ] # human
  print(noquote("genes in count table also in gene db:"))
  print(dim(counts))
  ############################################
  
  stopifnot(genes$ensembl_gene_id == rownames(counts))

  ## get subset of columns from count table
  counts_sub <- counts[, grep(countFilter, colnames(counts))]

  ## get samples information
  samples <- read_csv("samples.csv")
  print(samples)
  # make sure order of the colnames in counts is the same as samples
  counts_sub <- counts_sub[, paste0(samples$sample, "")]
  stopifnot(colnames(counts_sub) == samples$sample)
  write.csv(counts_sub, file.path(file.path(OUTPUT, "reports"), paste0("gene_counts_", region, ".csv")))
  
  colnames(counts_sub) <- paste(samples$experiment, samples$independent_var, samples$experiment_type, sep = "_")
  dge <- DGEList(counts = counts_sub, 
                 group = paste("time", samples$independent_var, samples$experiment_type, sep = "_"),
                 genes = genes)
  dge$samples$experiment <- samples$experiment
  dge$samples$independent_var <- samples$independent_var
  dge$samples$experiment_type <- samples$experiment_type
  
  # store only protein coding genes for later analysis
  dge <- dge[dge$genes$gene_biotype == "protein_coding", ,keep.lib.sizes = F]
  dge <- na.omit(dge)
  print(noquote("protein coding: "))
  print(dim(dge$genes))
  
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(file.path(OUTPUT, "objs"), paste0("dge_", region, "_protein_coding.rds")))
  
  cpm <- cpm(dge, log = T, prior.count = 1)
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  
  write.csv(as.data.frame(cpm), file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, ".csv")), row.names = F)
  write.csv(colSums(dge$counts), file.path(file.path(OUTPUT, "reports"), paste0("colsums", region, ".csv")))
}

filter_dge <- function(dge, region, ribothreshold = 1, rseqthreshold = 1) {
  # function to filter out genes with very low counts across all libraries
  # Genes must have at least NUMREP number of ribo-seq samples with count more than
  # RIBOTHRESHOLD, and at least NUMREP number of rna-seq samples with count more than
  # RSEQTHRESHOLD, to pass the filter
  # Args
  #   dge: DGE obejct created by make_dge
  #   region: cds, utr3 or utr5
  #   ribothreshold: cut-off cpm value for riboseq samples for filtering
  #   rseqthreshold: cut-off cpm value for rnaseq samples for filtering
  print("dimension before filtering:", quote = F)
  print(dim(dge))

  numrep <- ncol(dge) / 4 # divided by 2 data types and 2 conditions
  cpm_temp <- cpm(dge)
  keep <- (rowSums(cpm_temp[, 1:(ncol(cpm_temp) / 2)] > ribothreshold) >= numrep) &
    (rowSums(cpm_temp[, (ncol(cpm_temp) / 2 + 1):(ncol(cpm_temp))] > rseqthreshold) >= numrep)
  print(table(keep))
  
  dge <- dge[keep, ,keep.lib.sizes = FALSE]
  print("dimension after filtering:", quote = F)
  print(dim(dge))

  # normalization
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(OUTPUT, "objs", paste0("dge_", region, "_filt_norm.rds")))
  cpm <- cpm(dge, log = T, prior.count = 1) # default prior.count is 2
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  write.csv(cpm, file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, "_filt.csv")), row.names = F)
}

