## Block 1 
library(Seurat)
setwd("/Users/kurtisstefan/Documents/Byrd/NK_scRNA")
# 
# matrix_file <- "merge_5001/matrix.mtx.gz"
# features_file <- "merge_5001/features.tsv.gz"
# barcodes_file <- "merge_5001/barcodes.tsv.gz"
# 
# #scp -r stefankj@arcc2.uc.edu:/N/lustre/project/proj-481/ANUSHA_1/CSTX_4-results_human/raw_matrix ./CSTX_4-results_human  
# #scp -r stefankj@arcc2.uc.edu:/N/lustre/project/proj-481/ANUSHA_1/CSTX_7-results_human/raw_matrix ./CSTX_7-results_human   
# # scp -r stefankj@arcc2.uc.edu:/N/lustre/project/proj-481/ANUSHA_1/HMY2_4-results_human/raw_matrix ./HMY2_4-results_human  
# # scp -r stefankj@arcc2.uc.edu:/N/lustre/project/proj-481/ANUSHA_1/HMY2_7-results_human/raw_matrix ./HMY2_7-results_human  
# 
# loop = c("CSTX_4-results_human", "CSTX_7-results_human", "HMY2_4-results_human", "HMY2_7-results_human")
# for (i in loop){
#   print(i)
#   matrix_file = paste0(i, "/matrix.mtx.gz")
#   features_file = paste0(i, "/features.tsv.gz")
#   barcodes_file = paste0(i, "/barcodes.tsv.gz")
#   # Read the data into Seurat
#   seurat_object_name <- paste0("seurat_", i)
#   # Create the Seurat object
#   assign(seurat_object_name, 
#          CreateSeuratObject(
#            counts = ReadMtx(
#              mtx = matrix_file,
#              features = features_file,
#              cells = barcodes_file
#            )
#          ))
# }
# 
# pbmc.big <- merge(`seurat_CSTX_4-results_human`, 
#                   y = c(`seurat_CSTX_7-results_human`, 
#                         `seurat_HMY2_4-results_human`,
#                         `seurat_HMY2_7-results_human`), add.cell.ids = c("C4", "C7", "H4", "H7"), project = "PBMC15K")
# 
# pbmc.big[["percent.mt"]] <- PercentageFeatureSet(pbmc.big, pattern = "^MT-")
# seurat_object <- subset(
#   pbmc.big,
#   subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
# )
# 
# saveRDS(seurat_object, file="seuratobject.RDS")
seurat_object <- readRDS("seuratobject.RDS")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

rm(`seurat_HMY2_4-results_human`)
rm(`seurat_HMY2_7-results_human`)
rm(`seurat_CSTX_4-results_human`)
rm(`seurat_CSTX_7-results_human`)
rm(pbmc)
rm(pbmc.big)

#seurat_object
head(seurat_object@meta.data)
seurat_object$Method <- sub("_.*", "", row.names(seurat_object@meta.data))
seurat_object$Feeder <- sub("^(.).*", "\\1", row.names(seurat_object@meta.data))
table(seurat_object$Feeder)
seurat_object
# obj <- seurat_object
# obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj)
# obj <- ScaleData(obj)
# obj <- RunPCA(obj)
# obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
# obj <- FindClusters(obj, resolution = .5, cluster.name = "unintegrated_clusters")

options(future.globals.maxSize = 8000 * 1024^2)  # Set to 8 GB
# cluster without integration first
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = "pca")
seurat_object <- FindClusters(seurat_object, resolution = 0.8)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(seurat_object, reduction = "umap.unintegrated", split.by="Method", group.by="seurat_clusters") 

seurat_object <- JoinLayers(seurat_object)



markers <- FindMarkers(  object = seurat_object,  ident.1 = "C",  ident.2 = "H",  group.by = "Feeder",  test.use = "wilcox" )
saveRDS(markers, "CvsH.RDS")
markers <- readRDS("CvsH.RDS")
markers
?VlnPlot
saveRDS(seurat_object, "so.RDS")
seurat_object <- readRDS("so.RDS")
VlnPlot(object = seurat_object, features=c("GZMB", "GZMA"), group.by="Feeder")
# H have more Granzyme A
# C ha
library(tidyverse)


tw %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(desc(abs(avg_log2FC)))
markers






## BLOCK 2 
options(future.globals.maxSize = 8000 * 1024^2)  # Set to 8 GB
# cluster without integration first
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = "pca")
seurat_object <- FindClusters(seurat_object, resolution = 0.8)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(seurat_object, reduction = "pca", split.by="Method", group.by="seurat_clusters") 



#### Block 2
obj <- IntegrateLayers(object = seurat_object, method = RPCAIntegration,  #  CCAIntegration
                                 orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = .8, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(obj, reduction = "umap.cca", split.by = c("Method"), group.by="cca_clusters", combine = FALSE, label.size = 2)

# Block 3
# Summarize counts of each Method within each cca_clusters
library(tidyverse)
method_counts <- obj@meta.data %>%
  group_by(cca_clusters, Method) %>%
  tally(name = "count") %>%
  ungroup() %>% 
  group_by(Method) %>% 
  mutate(total = (count / sum(count))) %>%
  ungroup()

library(ggplot2)
library(RColorBrewer)
# Define the color palette
num_clusters <- length(unique(method_counts$cca_clusters))
num_colors <- 14

# Generate a custom color palette
custom_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_colors)
# Create the stacked bar plot
ggplot(method_counts, aes(x = Method, y = total, fill = cca_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Method", y = "Percentage (%)", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  


# Block 4 
install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
library(metap)
obj <- JoinLayers(obj)
obj$Method
obj@meta.data$cca_clusters
Idents(obj) <- "cca_clusters"
Idents(obj)

## We can find conserved unchanged markers that don't change based on condition here if desired 
#one.markers <- FindConservedMarkers(obj, ident.1 = 3, grouping.var="Method", verbose = FALSE)

tw <- FindMarkers(obj, ident.1 = "12", ident.2 =NULL , verbose = FALSE)
library(tidyverse)
tw
tw %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC > 1) %>% 
  arrange(p_val_adj)

