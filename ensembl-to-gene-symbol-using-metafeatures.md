# This code is a short hand way to achieve the conversion by using the in-built meta.features. If you don't have this info, use the other code which extracts the gene symbols from clusterProfiler:bitr() function and then constructs the new seurat object with updated rownames.
```
ref <- readRDS("/path/to/your/dataset.rds")

new_names <- ref[["RNA"]]@meta.features$feature_name

current_ids <- rownames(ref)
missing_mask <- is.na(new_names) | new_names == ""
new_names[missing_mask] <- current_ids[missing_mask]

new_names <- make.unique(as.character(new_names))

counts_matrix <- LayerData(ref, assay = "RNA", layer = "counts")
data_matrix   <- LayerData(ref, assay = "RNA", layer = "data")

rownames(counts_matrix) <- new_names
rownames(data_matrix)   <- new_names

# Rebuild the Assay with Updated Row Names
new_assay <- CreateAssay5Object(counts = counts_matrix, data = data_matrix)

# Replace the existing 'RNA' assay
ref[["RNA"]] <- new_assay
```
