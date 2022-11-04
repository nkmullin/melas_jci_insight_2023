#object_path: path to RData object
#clust_num_column: final_cluster_labels or equivalent (column of meta.data with cluster numbers)
#subset_logical: TRUE or FALSE (should a subset of the object be performed before dge?)
#subset_set_ident: "celltype" or equivalent (column name of meta.data for the subset). Unique values in this column will be called with the array index
#dge_ident: "final_cluster_labels" , "region", or equivalent: identity for the object before differential expression.
#dge_ident_1: NA or "fovea" (or equivalent). If performing dge on a biological variable, include the identity of the first group. 
#dge_ident_2: NA or "peripheral" (or equivalent). If performing dge on a biological variable, include the identity of the second group. 
#save_prefix: "clust_vs_all" or equivalent - prefix of each dge file to be saved
#save_suffix_ident: "final_cluster_labels" or "celltype" - column on meta.data that categorizes the cells into the groupings in which dge was performed
#save_directory: character string of the filepath to the directory for results saving

#############################
##### Load in Libraries #####
#############################
print("loading libraries")

library(crayon, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(backports, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(withr, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(rstudioapi, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(cli, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(Seurat, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(Signac, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(labeling, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(tidyverse, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

## I had troubles installing scaR in a custom location.
## The default libPath should be accessible to everyone
## If problems come up, the scaR binary is available in /Shared/IVR/apvoigt/programs/scaR_0.1.0.tar.gz  
library(scaR)

args <- commandArgs(trailingOnly = TRUE)

csv_args <- read_csv(as.character(args[2])) ## The first argument supplied to the R script is a csv file of all arguments
array_job_index <- as.numeric(args[3]) ## the array job index is the second argument supplied to the R script

########################
##### Read in data #####
########################
print("loading data")
file <- csv_args %>% filter(variable == "object_path") %>% dplyr::select("value") %>% as.character()
write(paste0("Reading file... ", file), stderr())
assign('seurat_object', get(load(file = file)))
DefaultAssay(seurat_object) <- "RNA"

#Identify suffix with which to save dge (before subsetting); then subset object (if requested)
##################
##### SUBSET #####
##################

save_suffix_ident <- csv_args %>% filter(variable == "save_suffix_ident") %>% dplyr::select("value") %>% as.character()

Idents(seurat_object) <- save_suffix_ident
subset_levels <- seurat_object@active.ident %>% unique() %>% as.character()
subset_levels <- str_sort(subset_levels, numeric = TRUE)
save_suffix <- str_c("cluster", "_", subset_levels[array_job_index])

### Subset (if requested)
subset_logical <- csv_args %>% filter(variable == "subset_logical") %>% dplyr::select("value") %>% as.logical()

if(subset_logical == TRUE){
  print("subsetting data")
  subset_set_ident <- csv_args %>% filter(variable == "subset_set_ident") %>% dplyr::select("value") %>% as.character()

  Idents(seurat_object) <- subset_set_ident
  
  if(subset_set_ident == "final_cluster_labels") {
    subset_levels <- seurat_object@active.ident %>% unique() %>% as.character()
    subset_levels <- str_sort(subset_levels, numeric = TRUE) # in case you are sorting stupidly-named
                                                             # final_cluster_labels like 8A and 8B
  } else {
    subset_levels <- seurat_object@active.ident %>% unique() %>% as.character() %>% sort()
  }
  

  save_suffix <- subset_levels[array_job_index]
  seurat_object <- subset(seurat_object, idents = subset_levels[array_job_index])
} 

#######################
##### Perform dge #####
#######################
print("performing DGE")
dge_ident <- csv_args %>% filter(variable == "dge_ident") %>% dplyr::select("value") %>% as.character()
dge_ident_1 <- csv_args %>% filter(variable == "dge_ident_1") %>% dplyr::select("value") %>% as.character()
dge_ident_2 <- csv_args %>% filter(variable == "dge_ident_2") %>% dplyr::select("value") %>% as.character()

if(is.na(dge_ident_2)) {
  if(is.na(dge_ident_1)) {
    Idents(seurat_object) <- dge_ident
    subset_levels <- seurat_object@active.ident %>% unique() %>% as.character()
    subset_levels <- str_sort(subset_levels, numeric = TRUE) 
    dge <- sc_dge(seurat_object, dge_ident, ident.1 = subset_levels[array_job_index], min.pct = 0, logfc.threshold = 0)
  } else {
    dge <- sc_dge(seurat_object, dge_ident, ident.1 = dge_ident_1, min.pct = 0, logfc.threshold = 0)
  }
} else {
  dge <- sc_dge(seurat_object, dge_ident, ident.1 = dge_ident_1, dge_ident_2, min.pct = 0, logfc.threshold = 0)
}

#####################################
##### Combine with GTF and save #####
#####################################

save_prefix <- csv_args %>% filter(variable == "save_prefix") %>% dplyr::select("value") %>% as.character()
save_directory <- csv_args %>% filter(variable == "save_directory") %>% dplyr::select("value") %>% as.character()

if(!file.exists(save_directory)) {dir.create(save_directory)}

dge %>% 
  as_tibble() %>%
  left_join(biomart_querry, by = c("gene" = "cellranger_gene_name")) %>%
  left_join(hg38gtf, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::select(-gene.y) %>%
  mutate(gene = gene.x) %>%
  dplyr::select(-gene.x) %>%
  dplyr::select(gene, dplyr::everything()) %>%
  write_csv(file.path(str_c(save_directory, "/", save_prefix, "_", save_suffix, ".csv")))

