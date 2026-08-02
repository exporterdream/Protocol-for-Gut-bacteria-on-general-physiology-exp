# ============================================================================
# DESeq2 Analysis: 2-Group Newly-Eclosed Model
# ============================================================================
#
# Performs differential expression analysis using DESeq2 on
# newly eclosed Drosophila samples only (excluding adult samples).
#
# Similiar to DESeq2_3group_analysis.R (same filtering,
# significance thresholds, and output structure), restricted to:
#
# Groups analyzed:
#   - A_inoR: Embryonically colonized, collected at newly eclosed stage (n=14)
#   - B_conR: Germ-free control, collected at newly eclosed stage (n=11)
#
# ============================================================================

# Load required libraries
library(DESeq2)
library(tidyverse)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

# Load count data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #no idea why we need this line, but without it, the script doesn't work
counts <- read.csv("salmon.merged.gene_counts.tsv", 
                   sep="\t", row.names=1)

gene_names <- counts$gene_name
names(gene_names) <- rownames(counts)

counts <- counts[, !colnames(counts) %in% c("gene_name")]

# Remove outlier sample A_inoA_11 (high duplication rate: 72.3%)
counts <- counts[, !colnames(counts) %in% c("A_inoA_11")]

# ============================================================================
# 2. FILTER TO NEWLY ECLOSED SAMPLES ONLY
# ============================================================================

sample_names <- colnames(counts)

# Identify newly eclosed samples (A_inoR, B_conR)
newly_samples <- sample_names[
  grepl("^A_inoR", sample_names) |
  grepl("^B_conR", sample_names)
]

counts_newly <- counts[, newly_samples]

cat("Newly eclosed samples selected:", length(newly_samples), "\n")
cat("  A_inoR:", sum(grepl("^A_inoR", newly_samples)), "\n")
cat("  B_conR:", sum(grepl("^B_conR", newly_samples)), "\n")

# ============================================================================
# 3. CREATE SAMPLE METADATA
# ============================================================================

conditions <- sapply(newly_samples, function(x) {
  if (grepl("^A_inoR", x)) return("A_inoR")
  if (grepl("^B_conR", x)) return("B_conR")
})

coldata <- data.frame(
  row.names = newly_samples,
  condition = factor(conditions, levels = c("B_conR", "A_inoR"))
)

coldata$condition <- relevel(coldata$condition, ref = "B_conR")

# ============================================================================
# 4. CREATE DESeq2 OBJECT AND RUN ANALYSIS
# ============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_newly),
  colData = coldata,
  design = ~ condition
)

keep <- rowSums(counts(dds)) >= 50 # same filter as 3-group script
dds <- dds[keep, ]
cat("\nGenes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

# ============================================================================
# 5. EXTRACT PAIRWISE COMPARISON (without lfcShrink, as in 3-group script)
# ============================================================================

padj_threshold <- 0.05 # p-value <= 0.05
lfc_threshold <- 0.5   # |log2FC| >= 0.5

# Comparison: A_inoR vs B_conR
res_AinoR_vs_BconR <- results(dds,
                              contrast = c("condition", "A_inoR", "B_conR"),
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

sig_AinoR_vs_BconR <- get_significant_degs(res_AinoR_vs_BconR, lfc_threshold, padj_threshold)

cat("A_inoR vs B_conR:", nrow(sig_AinoR_vs_BconR), "\n")

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

dir.create("DESeq2_2group_results", showWarnings = FALSE)

save_results <- function(res, filename) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_names[res_df$gene_id]
  res_df <- res_df[, c("gene_id", "gene_name", "baseMean", "log2FoldChange",
                        "lfcSE", "pvalue", "padj")]
  write.csv(res_df, filename, row.names = FALSE)
}

# Save full results
save_results(res_AinoR_vs_BconR, "DESeq2_2group_results/A_inoR_vs_B_conR_full.csv")

# Save significant DEGs
write.csv(sig_AinoR_vs_BconR, "DESeq2_2group_results/A_inoR_vs_B_conR_significant.csv", row.names = FALSE)

# Save significant gene names only for enrichment tools
writeLines(sig_AinoR_vs_BconR$gene_name,
           "DESeq2_2group_results/A_inoR_vs_B_conR_significant_for_enrichment.txt")

# ============================================================================
# 8. SAVE VST DATA FOR VISUALIZATION
# ============================================================================

vsd <- vst(dds, blind = TRUE)
vst_counts <- assay(vsd)
vst_df <- as.data.frame(vst_counts)
vst_df$gene_id <- rownames(vst_df)
vst_df$gene_name <- gene_names[vst_df$gene_id]
write.csv(vst_df, "DESeq2_2group_results/vst_counts_newly.csv", row.names = FALSE)

# Save PCA data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
write.csv(pca_data, "DESeq2_2group_results/pca_data_newly.csv", row.names = FALSE)

# Save sample metadata
write.csv(coldata, "DESeq2_2group_results/sample_metadata_newly.csv")

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("2-Group Newly Eclosed Analysis Summary\n")
cat("============================================================\n")
cat("Total samples:", ncol(counts_newly), "\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("\nDEGs:\n")
cat("  A_inoR vs B_conR:", nrow(sig_AinoR_vs_BconR),
    "(", sum(sig_AinoR_vs_BconR$log2FoldChange > 0), "up,",
    sum(sig_AinoR_vs_BconR$log2FoldChange < 0), "down )\n")
