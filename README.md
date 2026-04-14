# Protocol-for-Gut-bacteria-on-general-physiology-exp
Use DESeq2_3group_analysis.R to perform differential expression analysis using DESeq2 on adult Drosophila samples only (excluding newly eclosed samples). The file is salmon_merged.gene_counts.

The output significant DEGs file can be used as input for sep_up_downregulated. This script takes a CSV file containing significantly differentially expressed genes (DEGs) from a pairwise comparison and splits the genes into two text files:
Upregulated genes (log2FoldChange > 0)
Downregulated genes (log2FoldChange < 0)
Each output file contains one gene name per line (no header).

The input CSV must contain at least the following columns (standard DESeq2 output):
gene_name – gene identifier (can also be gene_id; adjust script if needed)
log2FoldChange – numeric value, positive = upregulated, negative = downregulated

Open the script and change the comparison variable to match your file prefix.

Output files:
Two files will be created in the same folder:
{comparison}_upregulated.txt
{comparison}_downregulated.txt


visualization_3group.ipynb

This script generates:
1. PCA plot (using VST-transformed data from DESeq2)
2. Volcano plots for all 3 pairwise comparisons (dot same size, no special meaning)
3. Gene list comparison (Venn diagram / UpSet plot)

Prerequisites: Run DESeq2_3group_analysis.R first to generate input files. All input for this script is output of DESeq2_3group_analysis.R.

volcano.ipynb loads DESeq2 results for all three comparisons and gene categories to create a volcano plot with:
    - Upregulated genes in red
    - Downregulated genes in blue
    - Dot size reflects expression level: size = log10(baseMean + 1) × 30

mass_length_analysis_complete_code.ipynb
ANALYSIS CODE FOR FLY MASS AND LENGTH DATA
This script performs:
1. Male vs Female comparison within each group
2. Male vs Female total
3. All pairwise t-tests between groups and U Tests
