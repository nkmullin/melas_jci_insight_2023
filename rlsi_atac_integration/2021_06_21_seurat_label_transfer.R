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

args <- commandArgs(trailingOnly = TRUE)

rna_path <- as.character(args[2]) ## The first argument supplied to the R script is a the path to the Seurat RNA object RData file
                                  ## THIS OBJECT MUST BE ANNOTATED WITH A COLUMN IN METADATA CALLED 'celltype' (lowercase)
atac_path <- as.character(args[3]) ## The second argument supplied to the R script is a the path to the signac (integrated) object RData file
                                  ## THIS OBJECT MUST HAVE AN RNA SLOT. LSI must be computed on the integrated data. 

save_path <- as.character(args[4]) ## Full save path (including .RData extension) where integrated file should be saved

########################
##### LOAD OBJECTS #####
########################
assign('rna_object', get(load(file = rna_path)))
DefaultAssay(rna_object) <- "RNA"
rna_object <- FindVariableFeatures(
  object = rna_object,
  nfeatures = 10000
)

assign('atac_object', get(load(file = atac_path)))
DefaultAssay(atac_object) <- "RNA"

###########################
##### TRANFER ANCHORS #####
###########################

rna_object$celltype <- Idents(object = rna_object)

transfer.anchors <- FindTransferAnchors(
  reference = rna_object,
  query = atac_object,
  reduction = 'cca',
  dims = 1:40,
  verbose = TRUE
)

##########################
##### PREDICT LABELS #####
##########################

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_object$celltype,
  weight.reduction = atac_object[['integrated_lsi']],
  dims = 2:30
)

integrated <- AddMetaData(object = atac_object, metadata = predicted.labels)
Idents(integrated) <- "predicted.id"


################
##### SAVE #####
################
save(integrated, file = save_path)
