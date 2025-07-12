# Load Libraries

library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ashr", force = TRUE)

# Set your working directory; need to have a raw_counts folder in it 
setwd("/Users/orzes/OneDrive/Desktop/Humber/BINF5503DAS/poster draft")

# Function to read in files; need to have an unzipped raw_counts file (or replace with your file name) in the working directory
merge_counts <- function(path = ".", pattern = "*.multiple.count") {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  merged <- NULL
  
  for (f in files) {
    sample_id <- sub(".multiple.count", "", basename(f)) #strip file suffix
    df <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c("gene_id", sample_id)
    
    if(is.null(merged)) {
      merged <- df
    } else {
      merged <- merge(merged, df, by = "gene_id", all = TRUE)
    }
  }
  
  merged[is.na(merged)] <- 
    return(merged)
  
}

# Call the function to make our zebrafish & kilifish count dfs
zebrafish_counts <- merge_counts("raw_counts/",  pattern = "^DR.*\\.multiple\\.count$")
kilifish_counts <- merge_counts("raw_counts/", pattern = "^JM.*\\.multiple\\.count$")

# Map nfu and ensembl id
ortholog_raw <- read.delim('GSE66712_nfu_counts_185_samples.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(ortholog_raw, 3)

ortholog_map <- ortholog_raw[, c(1, 2, 5)]
colnames(ortholog_map) <- c("nfu_gene_id", "dre_gene_id", "gene_symbol")

# Merge
dre_joined <- merge(ortholog_map, zebrafish_counts, by.x = "dre_gene_id", by.y = "gene_id")
nfu_joined <- merge(ortholog_map, kilifish_counts, by.x = "nfu_gene_id", by.y = "gene_id")

dre_clean <- dre_joined[, c("gene_symbol", grep("^DR", names(dre_joined), value = TRUE))]
nfu_clean <- nfu_joined[, c("gene_symbol", grep("^JM", names(nfu_joined), value = TRUE))]

# Get rid of duplicates
dre_clean <- aggregate(. ~ gene_symbol, data = dre_clean, FUN = sum)
nfu_clean <- aggregate(. ~gene_symbol, data = nfu_clean, FUN = sum)

# Join them into one df

total_counts <- full_join(dre_clean, nfu_clean, by = "gene_symbol")
total_counts[is.na(total_counts)] <- 0

# Gene symbols set to rownames
rownames(total_counts) <- total_counts$gene_symbol
count_matrix <- total_counts[, -1]

# Create metadata table
sample_names <- colnames(count_matrix)[colnames(count_matrix) != "gene_id"]
sample_names_cleaned <- gsub("_mzm.*|\\.bam|\\.multiple.*", "", sample_names)

metadata <- data.frame(file_column_name = sample_names_cleaned, stringsAsFactors = FALSE)
pattern <- "^(.*?)_(dre|nfu)_(liver|brain|skin)_([0-9]+)_(weeks|months)"

metadata <- metadata %>% 
  mutate(
    sample_id = sub(pattern, "\\1", file_column_name), 
    species = sub(pattern, "\\2", file_column_name), 
    tissue = sub(pattern, "\\3", file_column_name),
    age = as.integer(sub(pattern, "\\4", file_column_name)),
    age_units = sub(pattern, "\\5", file_column_name)
  )


# Normalize ages; 

metadata <- metadata %>% 
  mutate(age_in_days = ifelse(age_units == "months", age * 30, age * 7))

lifespan_days <- c("dre" = 1260, "nfu" =  150)
metadata <- metadata %>% 
  mutate(age_norm = age_in_days / lifespan_days[species])

# Define age groups
metadata <- metadata %>% 
  mutate(
    age_group = cut(
      age_norm,
      breaks = c(-Inf, 0.56 , 0.75, Inf),
      labels = c("young", "middle", "old")
    )) %>% 
  mutate(age_group = factor(age_group, levels = c("young", "middle", "old")))

# Run this if the above still doesn't work
metadata <- metadata %>% 
  mutate(age_group = case_when(age_in_days == 189 ~ "middle", 
                               TRUE ~ as.character(age_group)),
         age_group = factor(age_group, levels = c("young", "middle", "old"))
  )

# Finalize and Align with Count Matrix
colnames(count_matrix) <- gsub("(_mzm.*|\\.bam|\\.multiple.*)$", "", colnames(count_matrix))
rownames(metadata) <- metadata$file_column_name
count_matrix <- count_matrix[,metadata$file_column_name]
all(colnames(count_matrix) %in% metadata$file_column_name)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = metadata,
  design = ~ tissue + age_group + tissue:age_group
)

keep <- rowSums(counts(dds)) >= 10
# can count how many in keep to see how many genes were filtered out
dds <- dds[keep, ]
# dds <- DESeq(dds)

# Loop through age and tissue combinations
tissues <- c("liver", "brain", "skin")
species <- c("dre", "nfu")

# Circadian Clock Gene
clock_genes <- c("bmal1a", "bmal1b", "clock", "cry1", "cry2", "per1a", "per1b", "per2", "per3", "roraa", "rorab", "rorb", "rorca", "rorcb" )


plot_volcano_ggplot <- function(res, clock_genes, title = "") {
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  res_df <- res_df %>%
    mutate(
      significance = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 1 & gene %in% clock_genes ~ "Significant Clock Gene",
        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
        gene %in% clock_genes ~ "Not Significant Clock Gene",
        TRUE ~ "Not Significant"
      )
    )
  
  # Set factor levels
  res_df$significance <- factor(res_df$significance, levels = c(
    "Significant Clock Gene", "Significant", "Not Significant"
  ))
  
  # Separate layers
  base_layer <- res_df %>% filter(significance != "Significant Clock Gene")
  clock_layer <- res_df %>% filter(significance == "Significant Clock Gene")
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    # Non-clock genes base layer
    geom_point(data = base_layer,
               aes(color = significance, shape = significance),
               size = 0.8, alpha = 0.7) +
    
    # Significant clock genes on top
    geom_point(data = clock_layer,
               aes(color = significance, shape = significance),
               size = 2, alpha = 1) +
    
    # Label only significant clock genes
    geom_text(data = subset(clock_layer, padj < 0.05 & abs(log2FoldChange) > 1),
              aes(label = gene),
              vjust = 1.5, size = 3, color = "black") +
    
    scale_color_manual(values = c(
      "Significant Clock Gene" = "purple",      
      "Significant" = "lightpink",              
      "Not Significant" = "grey"
    )) +
    scale_shape_manual(values = c(
      "Significant Clock Gene" = 18, # diamond
      "Significant" = 19, # circle
      "Not Significant" = 15 # square
    )) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  return(p)
}



# Use in loop

all_results <- list()

dir.create("volcano_plots", showWarnings = FALSE)
for (sp in species) {
  
  dds_sub <- dds[, dds$species == sp]
  dds_sub$age_group <- droplevels(dds_sub$age_group)
  dds_sub$age_group <- relevel(dds_sub$age_group, ref = "young")
  
  design(dds_sub) <- ~ tissue + age_group + tissue:age_group
  keep <- rowSums(counts(dds_sub)) >= 10
  dds_sub <- dds_sub[keep, ]
  dds_sub <- DESeq(dds_sub)
  
  res_names <- resultsNames(dds_sub)
  
  for (tissue in tissues) {
    if (tissue == "brain") {
      contrast <- c("age_group_old_vs_young")
    } else {
      interaction_term <- paste0("tissue", tissue, ".age_groupold")
      if (!interaction_term %in% res_names) {
        message("Interaction term not found: ", interaction_term)
        next
      }
      contrast <- c("age_group_old_vs_young", interaction_term)
    }
    
    res <- results(dds_sub, contrast = list(contrast))
    res <- lfcShrink(dds_sub, contrast = contrast, res = res, type = "ashr")
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$species <- sp
    res_df$tissue <- tissue
    
    key <- paste(sp, tissue, sep = "_")
    all_results[[key]] <- res_df
    
    p <- plot_volcano_ggplot(res, clock_genes, title = key)
    
    ggsave(filename = paste0("volcano_plots/", key, ".png"),
           plot = p,
           width = 8,
           height = 6, 
           dpi = 300
    )
  }
}



ggsave(paste0("volcano_", sp, "_", tissue, ".pdf"), plot = p, width = 6, height = 5)
