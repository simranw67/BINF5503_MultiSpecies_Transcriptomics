#Load libraries

library(tidyverse)
library(dplyr)

#Set your own working directory
setwd("~/Humber/BINF 5503 Dara Analytics & Storytelling")

#Read in the files
nfu_ensembl_df <- read.delim('GSE66712_nfu_counts_185_samples.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_ontology_df <- read.csv("go_terms_for_protein_coding_genes.full.150922.csv", sep = "\t", header = TRUE)


kilifish_brain <- read.csv("~/Humber/BINF 5503 Dara Analytics & Storytelling/raw_counts/JM11_nfu_brain_39_weeks_mzm_tophat2_mapped.bam.multiple.count", sep = "\t", header = FALSE)

zebrafish_brain <- read.csv("~/Humber/BINF 5503 Dara Analytics & Storytelling/raw_counts/DR106_dre_liver_42_months.multiple.count", sep = "\t", header = FALSE)

#Name the columns for kilifish_brain 
colnames(kilifish_brain) <- c("Gene", "Counts")
colnames(gene_ontology_df)
colnames(nfu_ensembl_df)

#Homogenize the column names for all the dfs.
names(nfu_ensembl_df)[names(nfu_ensembl_df) == "Nfu_gene_id"] <- "Gene"

#Because the GO database maps the same nfu gene id to multiple go ids, we need to make some sort of selection criteria for easier data processing
# This first way -> take the first occurence of the Nfu_gene_id and filter out the dataframe for that

gene_ontology_dedup <- gene_ontology_df %>%
                       group_by(Gene) %>%
                       slice(1) %>%
                       ungroup() #first occurence

dim(gene_ontology_dedup) #dimensions are a lot more manageable now

# Filter the dataframe, so it only keeps those values that are in the kilifish_brain nfu id

gene_ontology_dedup <- gene_ontology_dedup %>% 
                       filter(Gene %in% kilifish_brain$Gene)

dim(gene_ontology_dedup)

kilifish_go_mapped <- kilifish_brain %>% 
                      left_join(gene_ontology_dedup, by = "Gene")
dim(kilifish_go_mapped)

# Second way -> can sort it by the lowest GO id for example
gene_ontology_dedup2 <- gene_ontology %>%
                        arrange(Gene, GO.Id) %>%
                        group_by(Gene) %>%
                        slice(1) %>%
                        ungroup() #deterministic method; lowest ensemble id for ex

# And then repeat the same steps as above.

# Can clean up df for only the columns we need (counts; nfu id, go id)


# Map to ensmbl id -> we're using the data frame that the prof supplied; as the dimensions match quite well

kilifish_ensmbl_mapped <- kilifish_brain %>% 
                          left_join(nfu_ensembl_df, by = "Gene")
dim(kilifish_ensmbl_mapped) #we lose some observations -> about 3000 not in df (when you loook at the values that are empty)

# If we're using GO; would also need to map zebrafish to GO. I have to look for the dataframe we would use




 
