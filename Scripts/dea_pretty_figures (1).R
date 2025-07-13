
#---- Setup ----

# Load Libraries
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ashr)

# install ashr
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ashr", force = TRUE)

# Set your working directory; need to have a raw_counts folder in it 
setwd("/Users/orzes/OneDrive/Desktop/Humber/BINF5503DAS/poster draft")

# Function to read in files; need to have an unzipped raw_counts file (or insert your file name) in the working directory
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


# Call the function to make our zebrafish & kilifish count dfs; pattern using regex
zebrafish_counts <- merge_counts("raw_counts/",  pattern = "^DR.*\\.multiple\\.count$")
kilifish_counts <- merge_counts("raw_counts/", pattern = "^JM.*\\.multiple\\.count$")


#---- Data Preprocessing ----

### ---- Map nfu and ensemblid (killifish & zebrafish respectively) to gene codes ---###

#Read in raw file that will be used for mapping
ortholog_raw <- read.delim('GSE66712_nfu_counts_185_samples.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(ortholog_raw, 3)

# Extract nfu, esnembl and gene symbol columns to make orthology mapping table
ortholog_map <- ortholog_raw[, c(1, 2, 5)]
colnames(ortholog_map) <- c("nfu_gene_id", "dre_gene_id", "gene_symbol")

# Merge each species dfs with orthology mapping table
dre_joined <- merge(ortholog_map, zebrafish_counts, by.x = "dre_gene_id", by.y = "gene_id")
nfu_joined <- merge(ortholog_map, kilifish_counts, by.x = "nfu_gene_id", by.y = "gene_id")

# Remove nfu and ensembl id columns
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

# Convert to a count matrix (ready to be put into deseq2); with gene symbols set to rownames
count_matrix <- total_counts[, -1]


###--- Create metadata table ---###

# Extract column names of the count matrix (filenames from raw counts)
sample_names <- colnames(count_matrix)[colnames(count_matrix) != "gene_id"]
# Use regex to remove file extensions
sample_names_cleaned <- gsub("_mzm.*|\\.bam|\\.multiple.*", "", sample_names)

#Use sample names cleaned as the file column name (identifiers) in metadata
metadata <- data.frame(file_column_name = sample_names_cleaned, stringsAsFactors = FALSE)

# Define regex pattern; will be used to extract metadata information from file names of raw counts
pattern <- "^(.*?)_(dre|nfu)_(liver|brain|skin)_([0-9]+)_(weeks|months)"

# Create columns in metadata using information from file names - sample ID, species, tissue and age
metadata <- metadata %>% 
            mutate(
              sample_id = sub(pattern, "\\1", file_column_name), 
              species = sub(pattern, "\\2", file_column_name), 
              tissue = sub(pattern, "\\3", file_column_name),
              age = as.integer(sub(pattern, "\\4", file_column_name)),
              age_units = sub(pattern, "\\5", file_column_name)
            )


# Normalize the ages 
metadata <- metadata %>% 
            mutate(age_in_days = ifelse(age_units == "months", age * 30, age * 7)) # create an age in days column

# Define lifespan for zebrafish and killifish in days
lifespan_days <- c("dre" = 1260, "nfu" =  150)

# Create age normalization column
metadata <- metadata %>% 
            mutate(age_norm = age_in_days / lifespan_days[species])

# Define age groups according to age normalization column
metadata <- metadata %>% 
            mutate(
              age_group = cut(
              age_norm,
              breaks = c(-Inf, 0.56 , 0.75, Inf),
              labels = c("young", "middle", "old")
            )) %>% 
            mutate(age_group = factor(age_group, levels = c("young", "middle", "old"))) 

# Manually label 189 days to middle for killifish
metadata <- metadata %>% 
            mutate(age_group = case_when(age_in_days == 189 ~ "middle", 
                                         TRUE ~ as.character(age_group)),
                   age_group = factor(age_group, levels = c("young", "middle", "old"))
                  )

# Finalize and check that colnames of count matrix and rownames of metdata align
colnames(count_matrix) <- gsub("(_mzm.*|\\.bam|\\.multiple.*)$", "", colnames(count_matrix)) # remove the file extensions from count matrix
rownames(metadata) <- metadata$file_column_name
count_matrix <- count_matrix[,metadata$file_column_name]
all(colnames(count_matrix) %in% metadata$file_column_name) # returns true



#---- Differential Expression Analysis -----

# Check that design matrix works and can be implemented
dds <- DESeqDataSetFromMatrix(
       countData = round(count_matrix),
       colData = metadata,
       design = ~ tissue + age_group + tissue:age_group
       )

# Check how many genes will be filtered (gene counts < 10)
keep <- rowSums(counts(dds)) >= 10
table(keep)

# Create species and tissues lists (to be able to run contrasts)
tissues <- c("liver", "brain", "skin")
species <- c("dre", "nfu")

# Create list of Circadian Clock Genes of interest
clock_genes <- c("bmal1a", "bmal1b", "clock", "cry1", "cry2", "per1a", "per1b", "per2", "per3", "roraa", "rorab", "rorb", "rorca", "rorcb" )


###--- Function for volcano plot using ggplot2 ----###
plot_volcano_ggplot <- function(res, clock_genes, title = "") {
                       res_df <- as.data.frame(res) #input of results from deseq; turned to dataframe
                       res_df$gene <- rownames(res_df)
  
                       res_df <- res_df %>%
                       mutate( # define significance levels as significant and clock gene, just significant and not significant
                        significance = case_when(
                        padj < 0.05 & abs(log2FoldChange) >= 1 & gene %in% clock_genes ~ "Significant Clock Gene",
                        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
                        gene %in% clock_genes ~ "Not Significant Clock Gene",
                        TRUE ~ "Not Significant"
                        )
                       )
  
                       # Set significane levels to factor
                       res_df$significance <- factor(res_df$significance, levels = c("Significant Clock Gene", "Significant", "Not Significant"))
  
                       # Separate layers
                       base_layer <- res_df %>% filter(significance != "Significant Clock Gene")
                       clock_layer <- res_df %>% filter(significance == "Significant Clock Gene")
  
                       # Save ggplot to p
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
    
                            scale_color_manual(values = c("Significant Clock Gene" = "purple",      
                                                          "Significant" = "lightpink",              
                                                          "Not Significant" = "grey")) +
                            scale_shape_manual(values = c("Significant Clock Gene" = 18, # diamond
                                                          "Significant" = 19, # circle
                                                          "Not Significant" = 15 # square
                                                          )) +
                            labs(title = title,
                                 x = "Log2 Fold Change",
                                 y = "-Log10 Adjusted P-Value") +
                            theme_minimal() +
                            theme(legend.title = element_blank())
  
                       return(p) # let function return the ggplot  
}



# Use plot_volcano_ggplot in a loop, running the various contrasts to do DEA (for kilifish and zebrafish separately)

all_results <- list() #empty list to store results 

dir.create("volcano_plots", showWarnings = FALSE) # Create a folder in your working directory to save volcano plots

for (sp in species) {
  
                      dds_sub <- dds[, dds$species == sp] # Will run through each species at a time
                      dds_sub$age_group <- droplevels(dds_sub$age_group)
                      dds_sub$age_group <- relevel(dds_sub$age_group, ref = "young") # Put young as reference; just to make everything the same
  
                      design(dds_sub) <- ~ tissue + age_group + tissue:age_group
                      
                      keep <- rowSums(counts(dds_sub)) >= 10
                      dds_sub <- dds_sub[keep, ]
                      dds_sub <- DESeq(dds_sub)
  
                      res_names <- resultsNames(dds_sub)
  
                          for (tissue in tissues) {
                            if (tissue == "brain") {
                              contrast <- c("age_group_old_vs_young")
                              } 
                            else {
                              interaction_term <- paste0("tissue", tissue, ".age_groupold")
                          
                              if (!interaction_term %in% res_names) {
                                message("Interaction term not found: ", interaction_term)
                                next
                                }
                            
                              contrast <- c("age_group_old_vs_young", interaction_term)
                              }
    
                            res <- results(dds_sub, contrast = list(contrast))
                            res <- lfcShrink(dds_sub, contrast = contrast, res = res, type = "ashr") #need to use ashr instead of normal because we're using contrasts
                            res_df <- as.data.frame(res)
                            res_df$gene <- rownames(res_df)
                            res_df$species <- sp
                            res_df$tissue <- tissue
    
                            key <- paste(sp, tissue, sep = "_")
                            all_results[[key]] <- res_df # save the results from dea into the empty list we made earlier
    
                            p <- plot_volcano_ggplot(res, clock_genes, title = key) # run volcano plot function for each res generated
    
                            ggsave(filename = paste0("volcano_plots/", key, ".png"), #save plots into volcano plots folder (we can change filetype)
                                   plot = p,
                                   width = 8,
                                   height = 6, 
                                   dpi = 300)
                        }
                    }

#---- Conservation of Clock Genes across the species ----

# Create a df to store the correlations for each tissue
clock_correlations_df <- data.frame()

# Loop through each tissue
for (tissue in tissues) {
                          key_dre <- paste("dre", tissue, sep = "_")
                          key_nfu <- paste("nfu", tissue, sep = "_")
  
                          if(!(key_dre %in% names(all_results)) || !(key_nfu %in% names(all_results))) {
                            next
                            }
  
                          res_dre <- all_results[[key_dre]] %>% filter(gene %in% clock_genes)
                          res_nfu <- all_results[[key_nfu]] %>% filter(gene %in% clock_genes)

                          # Merge l2fc by gene
                          merged_lf2c <- inner_join(res_dre, res_nfu, by = "gene", suffix = c("_dre", "_nfu"))
  
                          # If there's not enough clock genes found, skip normality check & correlation tests
                          if (nrow(merged_lf2c) < 3) {
                            warning("Too few circadian genes for tissue ", tissue)
                            next
                            }
  
                          # Shapiro wilk-test
                          shapiro_dre <- shapiro.test(merged_lf2c$log2FoldChange_dre)
                          shapiro_nfu <- shapiro.test(merged_lf2c$log2FoldChange_nfu)
  
                          # Correlations: 
                          pearson_corr <- cor(merged_lf2c$log2FoldChange_dre, merged_lf2c$log2FoldChange_nfu, method = "pearson")
                          spearman_corr <- cor(merged_lf2c$log2FoldChange_dre, merged_lf2c$log2FoldChange_nfu, method = "spearman")
  
                          # Create a table with shapiro wilk test and correlation results
                          clock_correlations_df <- rbind(clock_correlations_df, data.frame(
                                                         tissue = tissue,
                                                         n_genes = nrow(merged_lf2c), # number of clock genes that are being compared
                                                         pearson = round(pearson_corr, 3),
                                                         spearman = round(spearman_corr, 3),
                                                         shapiro_p_dre = round(shapiro_dre$p.value, 4),
                                                         shapiro_p_nfu = round(shapiro_nfu$p.value, 4)
                                                         ))
  
                          }

#Print out the correaltions df
print(clock_correlations_df)


#---- PCA of core clock gene expressions ----

# variance stabilizing transformation
vsd <- vst(dds, blind = TRUE) 
vsd_mat <- assay(vsd)

# subsetting clock genes
clock_genes_present <- clock_genes[clock_genes %in% rownames(vsd_mat)]
vsd_clock <- vsd_mat[clock_genes_present, ]

# transpose vsd_clock
vsd_clock_t <- t(vsd_clock)

# actual PCA and combining with metadata
pca_res <- prcomp(vsd_clock_t, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- rownames(pca_df)
pca_df <- merge(pca_df, as.data.frame(colData(vsd)), by.x = "sample", by.y = "row.names")

# GG that plot
pca = ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue, shape = species)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = c(
                         "brain" = "#FFC0CB",    
                         "skin"  = "#CD6090",  
                         "liver" = "#FDB863"   
                         )) +
      labs(title = "PCA of Clock Gene Expression") +
      theme_minimal() +
      theme(legend.title = element_blank())

ggsave(filename = "volcano_plots/PCA_clock_genes.png",
       plot = pca,
       width = 8,
       height = 6,
       dpi = 300
      )
