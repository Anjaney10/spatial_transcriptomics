```
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Seurat v5 helper: safely rename features by rebuilding the assay
RenameFeaturesSeurat5 <- function(object, assay = "RNA", mapping_df,
                  from_col = "ENSEMBL", to_col = "SYMBOL") {
  stopifnot(inherits(object, "Seurat"))
  stopifnot(assay %in% Assays(object))
  stopifnot(all(c(from_col, to_col) %in% colnames(mapping_df)))

  # Clean mapping: valid target symbols; 1:1 on source keys
  map <- mapping_df %>%
    dplyr::select(dplyr::all_of(c(from_col, to_col))) %>%
    dplyr::filter(!is.na(.data[[to_col]]), .data[[to_col]] != "") %>%
    dplyr::distinct(.data[[from_col]], .keep_all = TRUE)

  # Helper: build new names from old names using the mapping
  # old_names may still have version suffixes, so strip for matching
  build_new_names <- function(old_names) {
    stripped <- sub("\\.\\d+$", "", old_names)
    new <- map[[to_col]][match(stripped, map[[from_col]])]
    # Keep the original name when no mapping found
    new[is.na(new) | new == ""] <- old_names[is.na(new) | new == ""]
    make.unique(as.character(new))
  }

  # Extract existing layer matrices
  counts_mat <- tryCatch(LayerData(object, assay = assay, layer = "counts"), error = function(e) NULL)
  data_mat   <- tryCatch(LayerData(object, assay = assay, layer = "data"),   error = function(e) NULL)

  # Rename rownames on the extracted matrices
  if (!is.null(counts_mat)) {
    rownames(counts_mat) <- build_new_names(rownames(counts_mat))
  }
  if (!is.null(data_mat)) {
    rownames(data_mat) <- build_new_names(rownames(data_mat))
  }

  # Rebuild the Assay5 object from scratch with renamed matrices
  if (!is.null(counts_mat) && !is.null(data_mat)) {
    new_assay <- CreateAssay5Object(counts = counts_mat, data = data_mat)
  } else if (!is.null(counts_mat)) {
    new_assay <- CreateAssay5Object(counts = counts_mat)
  } else if (!is.null(data_mat)) {
    new_assay <- CreateAssay5Object(data = data_mat)
  } else {
    warning("No counts or data layer found in assay '", assay, "'")
    return(object)
  }

  # Replace the assay
  object[[assay]] <- new_assay

  # Update VariableFeatures if any
  vf <- tryCatch(VariableFeatures(object[[assay]]), error = function(e) character(0))
  if (length(vf) > 0) {
    vf_new <- build_new_names(vf)
    VariableFeatures(object[[assay]]) <- vf_new
  }

  object
}

# ---- Build ENSEMBL -> SYMBOL mapping from current RNA features (Seurat v5) ----
# Prefer counts features; fallback to data
ref_features <- tryCatch(
  rownames(LayerData(ref, assay = "RNA", layer = "counts")),
  error = function(e) rownames(LayerData(ref, assay = "RNA", layer = "data"))
)

# If Ensembl IDs contain version suffix (e.g., ENSG... .12), strip it for mapping
ref_features_stripped <- sub("\\.\\d+$", "", ref_features)

# Remove duplicates before bitr call (bitr errors on duplicate input)
ref_features_unique <- unique(ref_features_stripped)

mapping <- clusterProfiler::bitr(
  ref_features_unique,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

# bitr can return 1:many — keep first symbol per ENSEMBL
mapping <- mapping %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

# ---- Apply rename by rebuilding the assay ----
ref <- RenameFeaturesSeurat5(
  object = ref,
  assay = "RNA",
  mapping_df = mapping,
  from_col = "ENSEMBL",
  to_col = "SYMBOL"
)

# ---- Minimal validation ----
counts_rn <- tryCatch(rownames(LayerData(ref, "RNA", "counts")), error = function(e) NULL)
data_rn   <- tryCatch(rownames(LayerData(ref, "RNA", "data")),   error = function(e) NULL)

if (!is.null(counts_rn)) {
  message("Duplicates in counts: ", sum(duplicated(counts_rn)))
  message("NA in counts: ", sum(is.na(counts_rn)))
}
if (!is.null(data_rn)) {
  message("Duplicates in data: ", sum(duplicated(data_rn)))
  message("NA in data: ", sum(is.na(data_rn)))
}
if (!is.null(counts_rn) && !is.null(data_rn)) {
  message("counts/data rownames identical: ", isTRUE(all.equal(counts_rn, data_rn)))
}
```
