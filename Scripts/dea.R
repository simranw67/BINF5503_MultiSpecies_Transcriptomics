# Load Libraries

library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)


# Set your working directory; need to have a raw_counts folder in it 
setwd("~/Humber/BINF 5503 Dara Analytics & Storytelling")

# Function to read in files
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
      breaks = c(-Inf, 0.4 , 0.75, Inf),
      labels = c("young", "middle", "old")
    )
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
dds <- dds[keep, ]
dds <- DESeq(dds)

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
      significiance = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 1 & gene %in% clock_genes ~ "Significant Clock Gene",
        padj < 0.05 & abs(log2FoldChange) >= 1 ~"Significant",
        TRUE ~ "Not Significant"
      ),
      is_clock = gene %in% clock_genes
      
    )
  #Plot
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom+point(aes(color = significance), alpha = 0.6) + 
    scale_color_manual(values = c("
        Clock Gene (Significant)" = "#D55E00",
        "Significant" = "#0072B2",
        "Not Significant" = "gray70")) +
    geom_text(
      data = subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1 & is_clock),
      aes(label = gene),
      vjust = 1.5, size = 3.5, color = "black") +
    labs(
      title = title, 
      x = "Log2 Fold Change",
      y =  "-Log10 Adjusted P-Value"
    ) +
    theme_minima()+
    theme(legend.title = eleemnt_blank())
  
}

# Use in loop

all_results <- list()

for (sp in species) {
  dds_sub <- dds[, dds$species == sp]
  dds_sub$age_group <- relevel(dds_sub$age_group, ref = "young")
  design(dds_sub) <- ~ tissue + age_group + tissue:age_group
  
  keep <- rowSums(counts(dds_sub)) >= 10
  dds_sub <- dds_sub[keep, ]
  dds_sub <- DESeq(dds_sub)
  
  for (tissue in tissue_list) {
    contrast <- list(c("age_group_old", paste0("tissue", tissue, ".age_group_old")))
    res <- results(dds_sub, contrast = contrast)
    res <- lfcShrink(dds_sub, contrast = contrast[[1]], res = res, type = "normal")
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$species <- sp
    res_df$tissue <- tissue
    
    all_results[[paste(sp, tissue, sep = "_")]] <- res_df
    
    plot_volcano_ggplot(res, clock_genes, title =  paste(sp, tissue, sep = " - "))
  }
}




# ggsave(paste0("volcano_", sp, "_", t, ".pdf"), width = 6, height = 5) -> if we need to save it
