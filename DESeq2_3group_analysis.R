# ============================================================================
# DESeq2 Analysis: 3-Group Adult-Only Model
# ============================================================================
# 
# Performs differential expression analysis using DESeq2 on 
# adult Drosophila samples only (excluding newly eclosed samples).
#
# Groups analyzed:
#   - A_inoA: Embryonically colonized, collected as adults (n=10)
#   - B1_reintro: Adult recolonization (n=11)
#   - B2_conA: Germ-free control, collected as adults (n=11)
#
# ============================================================================

#Install DESeq2
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")

# install.packages("tidyverse")

# Load required libraries
library(DESeq2)
library(tidyverse)



# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

# Load count data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
counts <- read.csv("salmon.merged.gene_counts_angelatian.tsv", 
                   sep="\t", row.names=1)

gene_names <- counts$gene_name
names(gene_names) <- rownames(counts)

counts <- counts[, !colnames(counts) %in% c("gene_name")]

# Remove outlier sample A_inoA_11
counts <- counts[, !colnames(counts) %in% c("A_inoA_11")]

# ============================================================================
# 2. FILTER TO ADULT SAMPLES ONLY
# ============================================================================

sample_names <- colnames(counts)

# Identify adult samples (A_inoA, B1_reintro, B2_conA)
adult_samples <- sample_names[
  grepl("^A_inoA", sample_names) | 
  grepl("^B1_reintro", sample_names) | 
  grepl("^B2_conA", sample_names)
]

counts_adult <- counts[, adult_samples]

cat("Adult samples selected:", length(adult_samples), "\n")
cat("  A_inoA:", sum(grepl("^A_inoA", adult_samples)), "\n")
cat("  B1_reintro:", sum(grepl("^B1_reintro", adult_samples)), "\n")
cat("  B2_conA:", sum(grepl("^B2_conA", adult_samples)), "\n")

# ============================================================================
# 3. CREATE SAMPLE METADATA
# ============================================================================

conditions <- sapply(adult_samples, function(x) {
  if (grepl("^A_inoA", x)) return("A_inoA")
  if (grepl("^B1_reintro", x)) return("B1_reintro")
  if (grepl("^B2_conA", x)) return("B2_conA")
})

coldata <- data.frame(
  row.names = adult_samples,
  condition = factor(conditions, levels = c("B2_conA", "B1_reintro", "A_inoA"))
)

coldata$condition <- relevel(coldata$condition, ref = "B2_conA")

# ============================================================================
# 4. CREATE DESeq2 OBJECT AND RUN ANALYSIS
# ============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_adult),
  colData = coldata,
  design = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("\nGenes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

# ============================================================================
# 5. EXTRACT PAIRWISE COMPARISONS (without lfcShrink for simplicity)
# ============================================================================

padj_threshold <- 0.05 #p-value ≤ 0.05
lfc_threshold <- 0.5 #log2 ≥ 0.5

# Comparison 1: A_inoA vs B2_conA
res_AinoA_vs_B2conA <- results(dds, 
                                contrast = c("condition", "A_inoA", "B2_conA"),
                                alpha = padj_threshold)

# Comparison 2: B1_reintro vs B2_conA
res_B1reintro_vs_B2conA <- results(dds, 
                                    contrast = c("condition", "B1_reintro", "B2_conA"),
                                    alpha = padj_threshold)

# Comparison 3: A_inoA vs B1_reintro
res_AinoA_vs_B1reintro <- results(dds, 
                                   contrast = c("condition", "A_inoA", "B1_reintro"),
                                   alpha = padj_threshold)

# ============================================================================
# 6. FILTER SIGNIFICANT DEGs
# ============================================================================

get_significant_degs <- function(res, lfc_thresh, padj_thresh) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[res_df$gene_id]

  sig <- res_df[!is.na(res_df$padj) & 
                res_df$padj <= padj_thresh & 
                abs(res_df$log2FoldChange) >= lfc_thresh, ]

  return(sig)
}

sig_AinoA_vs_B2conA <- get_significant_degs(res_AinoA_vs_B2conA, lfc_threshold, padj_threshold)
sig_B1reintro_vs_B2conA <- get_significant_degs(res_B1reintro_vs_B2conA, lfc_threshold, padj_threshold)
sig_AinoA_vs_B1reintro <- get_significant_degs(res_AinoA_vs_B1reintro, lfc_threshold, padj_threshold)

cat("A_inoA vs B2_conA:", nrow(sig_AinoA_vs_B2conA), "\n")
cat("B1_reintro vs B2_conA:", nrow(sig_B1reintro_vs_B2conA), "\n")
cat("A_inoA vs B1_reintro:", nrow(sig_AinoA_vs_B1reintro), "\n")

# ============================================================================
# 7. GENERATE GENE LISTS
# ============================================================================

genes_AinoA_vs_B2conA <- rownames(sig_AinoA_vs_B2conA)
genes_B1reintro_vs_B2conA <- rownames(sig_B1reintro_vs_B2conA)
genes_AinoA_vs_B1reintro <- rownames(sig_AinoA_vs_B1reintro)

# RECOVERED: A_inoA_vs_B2 ∩ B1_vs_B2
RECOVERED <- intersect(genes_AinoA_vs_B2conA, genes_B1reintro_vs_B2conA)

# NOT_RECOVERED: A_inoA_vs_B2 ∩ A_inoA_vs_B1
NOT_RECOVERED <- intersect(genes_AinoA_vs_B2conA, genes_AinoA_vs_B1reintro)

# UNIQUE: A_inoA_vs_B2 only
UNIQUE <- setdiff(genes_AinoA_vs_B2conA, 
                  union(genes_B1reintro_vs_B2conA, genes_AinoA_vs_B1reintro))

cat("\nGene Lists:\n")
cat("  RECOVERED:", length(RECOVERED), "\n")
cat("  NOT_RECOVERED:", length(NOT_RECOVERED), "\n")
cat("  UNIQUE:", length(UNIQUE), "\n")

# ============================================================================
# 8. SAVE RESULTS
# ============================================================================

dir.create("DESeq2_3group_results", showWarnings = FALSE)

save_results <- function(res, filename) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[res_df$gene_id]
  res_df <- res_df[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", 
                        "lfcSE", "pvalue", "padj")]
  write.csv(res_df, filename, row.names = FALSE)
}

# Save full results
save_results(res_AinoA_vs_B2conA, "DESeq2_3group_results/A_inoA_vs_B2_conA_full.csv")
save_results(res_B1reintro_vs_B2conA, "DESeq2_3group_results/B1_reintro_vs_B2_conA_full.csv")
save_results(res_AinoA_vs_B1reintro, "DESeq2_3group_results/A_inoA_vs_B1_reintro_full.csv")

# Save significant DEGs
write.csv(sig_AinoA_vs_B2conA, "DESeq2_3group_results/A_inoA_vs_B2_conA_significant.csv", row.names = FALSE)
write.csv(sig_B1reintro_vs_B2conA, "DESeq2_3group_results/B1_reintro_vs_B2_conA_significant.csv", row.names = FALSE)
write.csv(sig_AinoA_vs_B1reintro, "DESeq2_3group_results/A_inoA_vs_B1_reintro_significant.csv", row.names = FALSE)

# Save gene lists with direction info
save_gene_list_with_direction <- function(gene_ids, res_primary, filename) {
  res_df <- as.data.frame(res_primary)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[res_df$gene_id]

  df <- res_df[res_df$gene_id %in% gene_ids, 
               c("gene_id", "gene_name", "log2FoldChange", "padj")]
  df$direction <- ifelse(df$log2FoldChange > 0, "up", "down")

  write.csv(df, filename, row.names = FALSE)

  # Also save gene names only for enrichment tools
  writeLines(df$gene_name, gsub(".csv", "_for_enrichment.txt", filename))
}

save_gene_list_with_direction(RECOVERED, res_AinoA_vs_B2conA, 
                               "RECOVERED.csv")
save_gene_list_with_direction(NOT_RECOVERED, res_AinoA_vs_B2conA, 
                               "NOT_RECOVERED.csv")
save_gene_list_with_direction(UNIQUE, res_AinoA_vs_B2conA, 
                               "UNIQUE.csv")

# ============================================================================
# 9. SAVE VST DATA FOR VISUALIZATION
# ============================================================================

vsd <- vst(dds, blind = TRUE)
vst_counts <- assay(vsd)
vst_df <- as.data.frame(vst_counts)
vst_df$gene_id <- rownames(vst_df)
vst_df$gene_name <- gene_names[vst_df$gene_id]
write.csv(vst_df, "vst_counts.csv", row.names = FALSE)

# Save PCA data
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE, ntop=500)
write.csv(pca_data, "pca_data.csv", row.names = FALSE)

# Save sample metadata
write.csv(coldata, "sample_metadata.csv")


# ============================================================================
# 10. SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("I need this\n")
cat("============================================================\n")
cat("Total samples:", ncol(counts_adult), "\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("\nDEGs by comparison:\n")
cat("  A_inoA vs B2_conA:", nrow(sig_AinoA_vs_B2conA), 
    "(", sum(sig_AinoA_vs_B2conA$log2FoldChange > 0), "up,",
    sum(sig_AinoA_vs_B2conA$log2FoldChange < 0), "down )\n")
cat("  B1_reintro vs B2_conA:", nrow(sig_B1reintro_vs_B2conA),
    "(", sum(sig_B1reintro_vs_B2conA$log2FoldChange > 0), "up,",
    sum(sig_B1reintro_vs_B2conA$log2FoldChange < 0), "down )\n")
cat("  A_inoA vs B1_reintro:", nrow(sig_AinoA_vs_B1reintro),
    "(", sum(sig_AinoA_vs_B1reintro$log2FoldChange > 0), "up,",
    sum(sig_AinoA_vs_B1reintro$log2FoldChange < 0), "down )\n")
cat("\nGene lists:\n")
cat("  RECOVERED:", length(RECOVERED), "\n")
cat("  NOT_RECOVERED:", length(NOT_RECOVERED), "\n")
cat("  UNIQUE:", length(UNIQUE), "\n")
