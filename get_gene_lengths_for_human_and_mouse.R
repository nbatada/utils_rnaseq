# Get gene lengths and save to file
# Human and Mouse
# Ensemble ID and Genesymbol

# Install necessary packages if not already installed
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))

##----------------------------------
rm(list=ls())
library(GenomicFeatures)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # #BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)  # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(org.Mm.eg.db)

# Function to get gene lengths and map IDs
get_gene_info <- function(txdb, orgdb, species_name) {
  message(paste("Processing", species_name, "..."))
  
  # Get transcripts grouped by gene
  tx_by_gene <- transcriptsBy(txdb, by="gene")
  
  # Calculate max transcript width per gene (gene length)
  gene_lens <- max(width(tx_by_gene))
  
  entrez_ids <- names(gene_lens)
  
  # Map Entrez IDs to Ensembl and Symbol
  ensembl_ids <- mapIds(orgdb, keys=entrez_ids, column="ENSEMBL", keytype="ENTREZID", multiVals="first")
  symbols <- mapIds(orgdb, keys=entrez_ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
  
  gene_info_df <- data.frame(
    Entrez_ID = entrez_ids,
    Ensembl_ID = ensembl_ids,
    Gene_Symbol = symbols,
    Gene_Length = gene_lens
  )
  
  # Remove rows where mapping failed for Ensembl or Symbol (optional, but cleans up)
  gene_info_df <- na.omit(gene_info_df)
  
  return(gene_info_df)
}

# Process Human (hg38)
human_gene_info <- get_gene_info(TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, "Human")
write.table(human_gene_info, file="human_gene_lengths.tsv", sep="\t", quote=FALSE, row.names=FALSE)


# Process Mouse (mm10)
mouse_gene_info <- get_gene_info(TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, "Mouse")
write.table(mouse_gene_info, file="mouse_gene_lengths.tsv", sep="\t", quote=FALSE, row.names=FALSE)

