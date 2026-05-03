# ============================================================================
# R Workflow Master Script
# Consolidated workflows from .github/workflows/*.R
# ============================================================================
# This master script organizes all R-based analysis workflows into functions
# for easy execution and management.
#
# Available Workflows:
# 1. differential_expression_analysis() - RNA-seq DESeq2 analysis
# 2. multi_genome_analysis() - Comparative genomics & synteny
# 3. virulence_factor_analysis() - Fungal virulence gene identification
# 4. visualization_plots() - Publication-quality visualizations
# 5. qc_pipeline_rnaseq() - QC pipeline for RNA-seq data
# 6. sporothrix_genomic_analysis() - Sporothrix species analysis
# ============================================================================

# ============================================================================
# 1. DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ============================================================================
differential_expression_analysis <- function() {
  cat("\n=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
  
  # Load libraries
  library(DESeq2)
  library(airway)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  
  # Load data
  data(airway)
  counts <- assay(airway)
  metadata <- as.data.frame(colData(airway))
  
  # ============================================
  # Step 1: Create DESeq2 object and run analysis
  # ============================================
  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("dex", "trt", "untrt"))
  res <- res[order(res$padj), ]
  
  # ============================================
  # Step 2: Explore the results
  # ============================================
  summary(res)
  sig_genes <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  n_sig <- nrow(sig_genes)
  cat("Number of significant genes (|log2FC| > 1, padj < 0.05):", n_sig, "\n")
  
  top10 <- head(res, 10)
  cat("Top 10 most significant genes:\n")
  print(top10)
  
  # ============================================
  # Step 3: Visualization
  # ============================================
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$significant <- ifelse(!is.na(res_df$padj) & 
                                 res_df$padj < 0.05 & 
                                 abs(res_df$log2FoldChange) > 1, 
                               "Significant", "Not Significant")
  
  # Volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, 
                                     y = -log10(padj), 
                                     color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Not Significant" = "gray", 
                                  "Significant" = "red")) +
    theme_minimal() +
    labs(title = "Differential Expression: Dexamethasone Treatment",
         x = "log2 Fold Change (Treatment vs Control)",
         y = "-log10(Adjusted P-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    theme(legend.position = "bottom")
  
  print(volcano_plot)
  ggsave("DE_volcano_plot.png", volcano_plot, width = 10, height = 8)
  
  # MA plot
  ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c("Not Significant" = "gray", 
                                  "Significant" = "red")) +
    theme_minimal() +
    labs(title = "MA Plot",
         x = "Mean of Normalized Counts (log scale)",
         y = "log2 Fold Change") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    theme(legend.position = "bottom")
  
  print(ma_plot)
  ggsave("DE_ma_plot.png", ma_plot, width = 10, height = 8)
  
  # Heatmap
  top30_genes <- rownames(head(sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ], 30))
  vst_counts <- assay(vst(dds))
  heatmap_data <- vst_counts[top30_genes, ]
  annotation_col <- data.frame(
    Treatment = metadata$dex,
    CellLine = metadata$cell
  )
  rownames(annotation_col) <- colnames(vst_counts)
  
  png("DE_heatmap.png", width = 1200, height = 1000)
  pheatmap(heatmap_data,
           scale = "row",
           annotation_col = annotation_col,
           show_rownames = TRUE,
           main = "Top 30 Differentially Expressed Genes",
           fontsize_row = 8,
           clustering_distance_rows = "correlation",
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  
  # ============================================
  # Step 4: Export results
  # ============================================
  write.csv(as.data.frame(res), "DE_results_all_genes.csv")
  write.csv(as.data.frame(sig_genes), "DE_results_significant.csv")
  
  cat("\n✓ Differential expression analysis complete!\n")
  cat("Files saved:\n")
  cat("  - DE_volcano_plot.png\n")
  cat("  - DE_ma_plot.png\n")
  cat("  - DE_heatmap.png\n")
  cat("  - DE_results_all_genes.csv\n")
  cat("  - DE_results_significant.csv\n")
}

# ============================================================================
# 2. MULTI-GENOME ANALYSIS (Comparative Genomics & Synteny)
# ============================================================================
multi_genome_analysis <- function(genome_dir = "") {
  cat("\n=== STARTING MULTI-GENOME ANALYSIS ===\n")
  
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(reshape2)
  library(tidyr)
  
  if (genome_dir != "") setwd(genome_dir)
  
  # Find genome files
  genome_files <- list.files(pattern = "\\.(gff|gtf|gff3)$", full.names = TRUE)
  
  if (length(genome_files) == 0) {
    cat("No GFF/GTF files found!\n")
    return(NULL)
  }
  
  cat("Found", length(genome_files), "genome files\n")
  
  # Read genome annotation function
  read_genome_data <- function(gff_file) {
    cat("Processing:", basename(gff_file), "\n")
    tryCatch({
      gff <- rtracklayer::import(gff_file)
      genes <- gff[gff$type == "gene", ]
      
      if (length(genes) == 0) {
        genes <- gff[grep("gene|mRNA|CDS", gff$type), ]
      }
      
      gene_names <- if (!is.null(genes$Name)) genes$Name 
                    else if (!is.null(genes$gene_id)) genes$gene_id 
                    else if (!is.null(genes$ID)) genes$ID 
                    else paste0("gene_", seq_along(genes))
      
      gene_df <- data.frame(
        chrom = as.character(seqnames(genes)),
        start = start(genes),
        end = end(genes),
        strand = as.character(strand(genes)),
        name = gene_names,
        stringsAsFactors = FALSE
      )
      
      gene_df <- na.omit(gene_df)
      cat("  Loaded", nrow(gene_df), "genes\n")
      return(gene_df)
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Load all genomes
  genome_data <- list()
  genome_names <- gsub("\\.(gff|gtf|gff3)$", "", basename(genome_files))
  
  for (i in seq_along(genome_files)) {
    data <- read_genome_data(genome_files[i])
    if (!is.null(data)) {
      genome_data[[genome_names[i]]] <- data
    }
  }
  
  if (length(genome_data) < 2) {
    cat("Need at least 2 genomes for analysis\n")
    return(NULL)
  }
  
  # Calculate gene statistics
  stats <- list()
  stats$gene_counts <- sapply(genome_data, nrow)
  stats$total_genes <- sum(stats$gene_counts)
  stats$avg_genes_per_genome <- mean(stats$gene_counts)
  
  # Gene content similarity
  genomes <- names(genome_data)
  similarity_matrix <- matrix(0, nrow = length(genomes), ncol = length(genomes))
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- genomes
  
  for (i in 1:length(genomes)) {
    for (j in 1:length(genomes)) {
      if (i != j) {
        genes_i <- genome_data[[genomes[i]]]$name
        genes_j <- genome_data[[genomes[j]]]$name
        shared_genes <- length(intersect(genes_i, genes_j))
        similarity_matrix[i, j] <- shared_genes / min(length(genes_i), length(genes_j))
      } else {
        similarity_matrix[i, j] <- 1
      }
    }
  }
  
  # Visualizations
  cat("\nGenerating visualizations...\n")
  
  # Gene counts plot
  p1 <- ggplot(data.frame(genome = names(stats$gene_counts), count = stats$gene_counts),
               aes(x = reorder(genome, -count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Gene Counts per Genome",
         x = "Genome", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  ggsave("multi_genome_gene_counts.png", p1, width = 12, height = 8)
  
  # Similarity heatmap
  similarity_melt <- melt(similarity_matrix)
  p2 <- ggplot(similarity_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "darkred", midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Gene Content Similarity Matrix",
         x = "Genome", y = "Genome", fill = "Similarity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  ggsave("multi_genome_similarity_heatmap.png", p2, width = 10, height = 8)
  
  # Save results
  save(genome_data, stats, similarity_matrix, file = "multi_genome_analysis_results.RData")
  write.csv(as.data.frame(stats), "multi_genome_statistics.csv")
  
  cat("\n✓ Multi-genome analysis complete!\n")
  cat("Files saved:\n")
  cat("  - multi_genome_gene_counts.png\n")
  cat("  - multi_genome_similarity_heatmap.png\n")
  cat("  - multi_genome_analysis_results.RData\n")
  cat("  - multi_genome_statistics.csv\n")
}

# ============================================================================
# 3. VIRULENCE FACTOR ANALYSIS
# ============================================================================
virulence_factor_analysis <- function(genome_dir = "") {
  cat("\n=== STARTING VIRULENCE FACTOR ANALYSIS ===\n")
  
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(rtracklayer)
  library(GenomicRanges)
  
  if (genome_dir != "") setwd(genome_dir)
  
  # Find genome files
  gtf_files <- list.files(pattern = "\\.gtf$", full.names = TRUE)
  gff_files <- c(list.files(pattern = "\\.gff$", full.names = TRUE),
                 list.files(pattern = "\\.gff3$", full.names = TRUE))
  
  all_files <- c(gtf_files, gff_files)
  cat("Found", length(all_files), "genome files\n")
  
  if (length(all_files) == 0) {
    cat("No GTF/GFF files found!\n")
    return(NULL)
  }
  
  # Parse genome files
  parse_genome_file <- function(file_path) {
    cat("Processing:", basename(file_path), "\n")
    tryCatch({
      file_ext <- tools::file_ext(file_path)
      genome_name <- tools::file_path_sans_ext(basename(file_path))
      
      data <- read.table(file_path, sep = "\t", header = FALSE, 
                        stringsAsFactors = FALSE, comment.char = "#",
                        quote = "", fill = TRUE, na.strings = ".",
                        col.names = c("seqid", "source", "type", "start", "end", 
                                      "score", "strand", "phase", "attributes"))
      
      data$genome <- genome_name
      data$file_type <- file_ext
      
      # Parse attributes
      data <- data %>%
        mutate(
          gene_id = str_extract(attributes, "gene_id[= ]\"?([^;\"]+)\"?"),
          gene_name = str_extract(attributes, "gene_name[= ]\"?([^;\"]+)\"?"),
          product = str_extract(attributes, "product[= ]\"?([^;\"]+)\"?"),
          description = str_extract(attributes, "description[= ]\"?([^;\"]+)\"?"),
          name = str_extract(attributes, "Name[= ]\"?([^;\"]+)\"?"),
          id = str_extract(attributes, "ID[= ]\"?([^;\"]+)\"?"),
          gene_id = coalesce(
            str_remove_all(gene_id, 'gene_id[= ]"?|"|;$'),
            str_remove_all(id, 'ID[= ]"?|"|;$')
          ),
          product = coalesce(
            str_remove_all(product, 'product[= ]"?|"|;$'),
            str_remove_all(description, 'description[= ]"?|"|;$')
          ),
          gene_name = coalesce(
            str_remove_all(gene_name, 'gene_name[= ]"?|"|;$'),
            str_remove_all(name, 'Name[= ]"?|"|;$'),
            gene_id
          )
        )
      
      cat("Successfully processed", nrow(data), "features\n")
      return(data)
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      return(NULL)
    })
  }
  
  genome_data <- map_dfr(all_files, parse_genome_file)
  
  if (is.null(genome_data) || nrow(genome_data) == 0) {
    cat("No data processed\n")
    return(NULL)
  }
  
  gene_data <- genome_data %>% 
    filter(type %in% c("gene", "mRNA", "CDS", "exon"))
  
  # Identify virulence factors
  virulence_keywords <- c(
    "protease", "peptidase", "lipase", "hydrolase", "chitinase",
    "toxin", "effector", "adhesin", "biofilm", "catalase",
    "secreted", "invasin", "siderophore", "kinase", "virulence"
  )
  
  search_data <- gene_data %>%
    mutate(
      search_text = paste(coalesce(product, ""),
                          coalesce(description, ""),
                          coalesce(gene_name, ""),
                          sep = " | "),
      search_text = tolower(search_text)
    )
  
  virulence_pattern <- paste(virulence_keywords, collapse = "|")
  
  virulence_results <- search_data %>%
    mutate(
      is_virulence = str_detect(search_text, regex(virulence_pattern, ignore_case = TRUE))
    ) %>%
    filter(is_virulence) %>%
    mutate(
      virulence_category = case_when(
        str_detect(search_text, "protease|peptidase|lipase|hydrolase") ~ "Hydrolases",
        str_detect(search_text, "toxin|effector") ~ "Toxins/Effectors",
        str_detect(search_text, "adhesin|biofilm") ~ "Adhesion/Biofilm",
        str_detect(search_text, "catalase|oxidative|stress") ~ "Stress Response",
        str_detect(search_text, "secreted|invasin") ~ "Secreted/Invasion",
        TRUE ~ "Other Virulence Factors"
      )
    )
  
  cat("Found", nrow(virulence_results), "potential virulence factors\n")
  
  # Summary statistics
  summary_stats <- virulence_results %>%
    group_by(genome, virulence_category) %>%
    summarise(count = n(), .groups = 'drop')
  
  genome_summary <- virulence_results %>%
    group_by(genome) %>%
    summarise(
      total_virulence_genes = n(),
      total_genes = nrow(filter(gene_data, genome == first(genome))),
      virulence_percentage = (total_virulence_genes / total_genes) * 100
    )
  
  # Visualizations
  p1 <- ggplot(summary_stats, aes(x = genome, y = count, fill = virulence_category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Virulence Factor Distribution",
         x = "Genome", y = "Number of Virulence Factors",
         fill = "Category") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  ggsave("virulence_distribution.png", p1, width = 12, height = 8)
  
  p2 <- ggplot(genome_summary, aes(x = reorder(genome, -virulence_percentage), 
                                   y = virulence_percentage)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Virulence Factor Percentage by Genome",
         x = "Genome", y = "Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  ggsave("virulence_percentage.png", p2, width = 12, height = 8)
  
  # Save results
  write.csv(virulence_results, "virulence_factors_detailed.csv", row.names = FALSE)
  write.csv(summary_stats, "virulence_summary_by_category.csv", row.names = FALSE)
  write.csv(genome_summary, "virulence_summary_by_genome.csv", row.names = FALSE)
  
  cat("\n✓ Virulence factor analysis complete!\n")
  cat("Files saved:\n")
  cat("  - virulence_distribution.png\n")
  cat("  - virulence_percentage.png\n")
  cat("  - virulence_factors_detailed.csv\n")
}

# ============================================================================
# 4. QC PIPELINE FOR RNA-SEQ
# ============================================================================
qc_pipeline_rnaseq <- function() {
  cat("\n=== STARTING QC PIPELINE FOR RNA-SEQ ===\n")
  
  library(DESeq2)
  library(airway)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  
  # Load data
  data(airway)
  counts <- assay(airway)
  metadata <- as.data.frame(colData(airway))
  
  # Log transformation
  log_counts <- log2(counts + 1)
  cat("✓ Log transformation complete\n")
  
  # Top 10 genes by expression
  gene_means <- rowMeans(counts)
  top_genes <- names(sort(gene_means, decreasing = TRUE)[1:10])
  cat("✓ Top 10 genes identified\n")
  
  # TP53 expression boxplot
  tp53_counts <- data.frame(
    SampleID = colnames(counts),
    Count = counts["ENSG00000141510", ]
  ) %>%
    left_join(metadata %>% rownames_to_column("SampleID"), by = "SampleID")
  
  p1 <- ggplot(tp53_counts, aes(x = dex, y = log2(Count + 1), fill = dex)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    theme_minimal() +
    labs(title = "TP53 Expression: Treated vs Untreated",
         x = "Treatment", y = "log2(TP53 Expression + 1)") +
    scale_fill_manual(values = c("untrt" = "lightblue", "trt" = "lightcoral"))
  
  print(p1)
  ggsave("QC_TP53_boxplot.png", p1, width = 8, height = 6)
  
  # Replicate correlation
  rep1 <- counts[, "SRR1039508"]
  rep2 <- counts[, "SRR1039509"]
  correlation <- cor(rep1, rep2, method = "pearson")
  cat("Replicate correlation:", round(correlation, 4), "\n")
  
  replicate_corr_df <- data.frame(
    Replicate1 = log2(rep1 + 1),
    Replicate2 = log2(rep2 + 1)
  )
  
  p2 <- ggplot(replicate_corr_df, aes(x = Replicate1, y = Replicate2)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    theme_minimal() +
    annotate("text", x = 5, y = 15, 
             label = paste("r =", round(correlation, 3)), 
             size = 5, color = "darkred") +
    labs(title = "Technical Replicate Correlation",
         x = "log2(SRR1039508 Expression + 1)",
         y = "log2(SRR1039509 Expression + 1)")
  
  print(p2)
  ggsave("QC_replicate_correlation.png", p2, width = 8, height = 6)
  
  # PCA plot
  gene_vars <- apply(log_counts, 1, var)
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])
  log_counts_subset <- log_counts[top_var_genes, ]
  
  pca_result <- prcomp(t(log_counts_subset), scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$SampleID <- rownames(pca_df)
  pca_df$Treatment <- metadata$dex
  pca_df$CellLine <- metadata$cell
  
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = CellLine)) +
    geom_point(size = 5, alpha = 0.8) +
    theme_minimal() +
    scale_color_manual(values = c("untrt" = "blue", "trt" = "red")) +
    labs(title = "PCA of RNA-seq Samples",
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme(legend.position = "bottom")
  
  print(p3)
  ggsave("QC_PCA_plot.png", p3, width = 10, height = 8)
  
  cat("\n✓ QC pipeline complete!\n")
  cat("Files saved:\n")
  cat("  - QC_TP53_boxplot.png\n")
  cat("  - QC_replicate_correlation.png\n")
  cat("  - QC_PCA_plot.png\n")
}

# ============================================================================
# 5. SPOROTHRIX GENOMIC ANALYSIS
# ============================================================================
sporothrix_genomic_analysis <- function(genome_dir = "") {
  cat("\n=== STARTING SPOROTHRIX GENOMIC ANALYSIS ===\n")
  
  library(Biostrings)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(ape)
  
  if (genome_dir != "") setwd(genome_dir)
  
  # Create project directories
  if (!dir.exists("results")) dir.create("results")
  if (!dir.exists("plots")) dir.create("plots")
  
  # Read genome sequences
  sporothrix_genomes <- list()
  genome_files <- list.files(pattern = "\\.fna$|\\.fasta$")
  
  for(file in genome_files) {
    species_name <- gsub("\\.fna$|\\.fasta$", "", file)
    sporothrix_genomes[[species_name]] <- readDNAStringSet(file)
    cat("Loaded", species_name, "-", length(sporothrix_genomes[[species_name]]), "sequences\n")
  }
  
  # Read annotation files
  sporothrix_annotations <- list()
  gff_files <- list.files(pattern = "\\.gff$|\\.gff3$")
  
  for(file in gff_files) {
    species_name <- gsub("\\.gff$|\\.gff3$", "", file)
    sporothrix_annotations[[species_name]] <- rtracklayer::import(file)
    cat("Loaded", species_name, "annotations -", 
        length(sporothrix_annotations[[species_name]]), "features\n")
  }
  
  # Calculate basic genome statistics
  genome_stats <- data.frame()
  
  for(species in names(sporothrix_genomes)) {
    genome <- sporothrix_genomes[[species]]
    annotation <- sporothrix_annotations[[species]]
    
    total_length <- sum(width(genome))
    n_contigs <- length(genome)
    gc_content <- sum(Biostrings::letterFrequency(genome, "GC")) / total_length
    n_genes <- sum(annotation$type == "gene", na.rm = TRUE)
    n_cds <- sum(annotation$type == "CDS", na.rm = TRUE)
    
    stats <- data.frame(
      Species = species,
      Total_length = total_length,
      N_contigs = n_contigs,
      GC_content = round(gc_content, 4),
      N_genes = n_genes,
      N_CDS = n_cds
    )
    
    genome_stats <- rbind(genome_stats, stats)
  }
  
  cat("\nGenome Statistics:\n")
  print(genome_stats)
  
  # Plot genome statistics
  p1 <- ggplot(genome_stats, aes(x = Species, y = Total_length/1e6, fill = Species)) +
    geom_bar(stat = "identity") +
    labs(title = "Genome Size Comparison", y = "Genome Size (Mb)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("plots/genome_sizes.png", p1, width = 10, height = 8)
  
  p2 <- ggplot(genome_stats, aes(x = Species, y = GC_content, fill = Species)) +
    geom_bar(stat = "identity") +
    labs(title = "GC Content Comparison", y = "GC Content") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("plots/gc_content.png", p2, width = 10, height = 8)
  
  # Save results
  write.csv(genome_stats, "results/genome_statistics.csv", row.names = FALSE)
  
  cat("\n✓ Sporothrix genomic analysis complete!\n")
  cat("Files saved in results/ and plots/ directories\n")
}

# ============================================================================
# MASTER WORKFLOW EXECUTION FUNCTION
# ============================================================================
run_all_workflows <- function(data_dir = "") {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║           R WORKFLOW MASTER EXECUTION SCRIPT                  ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  
  # Run each workflow
  tryCatch(differential_expression_analysis(), 
           error = function(e) cat("DESeq2 analysis skipped:", e$message, "\n"))
  
  tryCatch(multi_genome_analysis(data_dir),
           error = function(e) cat("Multi-genome analysis skipped:", e$message, "\n"))
  
  tryCatch(virulence_factor_analysis(data_dir),
           error = function(e) cat("Virulence factor analysis skipped:", e$message, "\n"))
  
  tryCatch(qc_pipeline_rnaseq(),
           error = function(e) cat("QC pipeline skipped:", e$message, "\n"))
  
  tryCatch(sporothrix_genomic_analysis(data_dir),
           error = function(e) cat("Sporothrix analysis skipped:", e$message, "\n"))
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                    ALL WORKFLOWS COMPLETE                     ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
}

# ============================================================================
# WORKFLOW GUIDE
# ============================================================================
workflow_guide <- function() {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                    R WORKFLOW GUIDE                           ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  cat("Available workflows:\n\n")
  
  cat("1. differential_expression_analysis()\n")
  cat("   - DESeq2-based RNA-seq differential expression analysis\n")
  cat("   - Generates volcano plots, MA plots, heatmaps\n")
  cat("   - Usage: differential_expression_analysis()\n\n")
  
  cat("2. multi_genome_analysis(genome_dir = '')\n")
  cat("   - Comparative genomics analysis\n")
  cat("   - Analyzes GFF/GTF files for synteny and gene content\n")
  cat("   - Usage: multi_genome_analysis('/path/to/genomes')\n\n")
  
  cat("3. virulence_factor_analysis(genome_dir = '')\n")
  cat("   - Identifies virulence factors in fungal genomes\n")
  cat("   - Categorizes by function (hydrolases, toxins, etc.)\n")
  cat("   - Usage: virulence_factor_analysis('/path/to/genomes')\n\n")
  
  cat("4. qc_pipeline_rnaseq()\n")
  cat("   - Quality control pipeline for RNA-seq data\n")
  cat("   - Includes PCA, correlation analysis, expression plots\n")
  cat("   - Usage: qc_pipeline_rnaseq()\n\n")
  
  cat("5. sporothrix_genomic_analysis(genome_dir = '')\n")
  cat("   - Specialized analysis for Sporothrix species\n")
  cat("   - Calculates genome statistics and visualizations\n")
  cat("   - Usage: sporothrix_genomic_analysis('/path/to/genomes')\n\n")
  
  cat("6. run_all_workflows(data_dir = '')\n")
  cat("   - Executes all workflows in sequence\n")
  cat("   - Usage: run_all_workflows('/path/to/data')\n\n")
  
  cat("7. workflow_guide()\n")
  cat("   - Displays this help message\n\n")
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Print welcome message and guide
cat("\n✓ R Workflows loaded successfully!\n")
cat("Type workflow_guide() to see available functions\n")
