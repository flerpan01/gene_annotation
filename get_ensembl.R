sort_chrom <- function(d) {
  d$chr <- sub("chr", "", d$chr)

  # Sort table after chromosomes > start pos
  d$chr <- sub("MT", 100, d$chr)
  d$chr <- sub("X", 101, d$chr)
  d$chr <- sub("Y", 102, d$chr)
  d <- d[order(as.numeric(d$chr), d$start), ]
  d$chr <- sub(100, "MT", d$chr)
  d$chr <- sub(101, "X", d$chr)
  d$chr <- sub(102, "Y", d$chr)
  rownames(d) <- NULL
  d$chr <- sub("^", "chr", d$chr)
  
  return(d)
}

# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
library(biomaRt)

get_ensembl <- function(chrom = c(1:19, "MT", "X", "Y")){
  # listEnsembl() # list available datasets
  # mart <- useEnsembl(biomart="genes") # download all gene lists
  # searchDatasets(mart=mart, pattern="mus") # identify house mouse dataset
  # mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  # listAttributes(mart) # list of available attributes

  # get mart
  mart <- useEnsembl(
    biomart = "genes",
    dataset = "mmusculus_gene_ensembl") # GRCm39

  # columns to import
  #attr <- listAttributes(mart)$name # holds all columns
  cols <- c(
    "chromosome_name", 
    "start_position", 
    "end_position", 
    "strand", 
    "ensembl_gene_id", 
    "external_gene_name",
    "entrezgene_id",
    "description", 
    "gene_biotype")

  # import columns
  ens <- getBM(mart = mart, attributes = cols)
  ens <- ens[, cols]

  print(paste("Generating ensembl table..."))

  colnames(ens) <- c(
    "chr", 
    "start", 
    "end", 
    "strand", 
    "gene_id",     
    "gene_name",
    "entrez_id",  
    "gene_info", 
    "gene_type"
  )

  ens <- subset(ens, chr %in% chrom)
  ens <- sort_chrom(ens)

  # reduce gene types, collapse pseudogenes
  ens$gene_type2 <- ens$gene_type
  rows <- grep("pseudo", ens$gene_type2)
  ens$gene_type2[rows] <- "pseudogene"
  rows <- ens$gene_type %in% c(
    "snoRNA", "misc_RNA", "sRNA", "scaRNA", "snoRNA",
    "snRNA", "scRNA", "miRNA")
  ens$gene_type2[rows] <- "ncRNA"

  # 400+ IG / TR genes hid in protein coding
  rows <- grep("_gene", ens$gene_type2)
  ens$gene_type2[rows] <- "protein_coding"

  ens$size <- ens$end - ens$start

  print(paste("Total number of genes imported =", nrow(ens))) # 57287
  return(ens)
}

ens <- get_ensembl()