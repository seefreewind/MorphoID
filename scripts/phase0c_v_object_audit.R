#!/usr/bin/env Rscript

## Phase 0C-V author-corrected Seurat object audit.
## Descriptive only: no clustering, prediction, residualization, or model fitting.

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

script_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
raw <- file.path(root, "data/raw/phase0c_v")
out <- file.path(root, "results/phase0/phase0c_v")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rds <- file.path(raw, "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")
remove_ids_file <- gzfile(file.path(raw, "GSE250346_remove_nuclei_VUILD105MA1.csv.gz"), open = "rt")
removed_ids <- read.csv(remove_ids_file, stringsAsFactors = FALSE, check.names = FALSE)
close(remove_ids_file)
remove_xy_file <- gzfile(file.path(raw, "GSE250346_remove_coordinates_VUILD105MA1.csv.gz"), open = "rt")
removed_xy <- read.csv(remove_xy_file, stringsAsFactors = FALSE, check.names = FALSE)
close(remove_xy_file)

cat("Loading corrected Seurat object...\n")
obj <- readRDS(rds)
if (!inherits(obj, "Seurat")) stop("Corrected RDS is not a Seurat object")
## Access the stored metadata directly. The user-facing x[[]] accessor expands
## the full 1.63-million-row object and exceeded this R runtime's vector limit.
meta0 <- slot(obj, "meta.data")
cells0 <- rownames(meta0)
if (is.null(cells0) || anyDuplicated(cells0)) stop("Cell metadata rownames are missing or duplicated")
needed <- c("sample", "patient", "tma", "run", "disease_status", "sample_affect",
            "final_CT", "final_lineage", "cell_id", "full_cell_id", "x_centroid", "y_centroid",
            "CNiche", "TNiche", "nCount_RNA", "nFeature_RNA")
missing <- setdiff(needed, colnames(meta0))
if (length(missing)) stop(paste("Missing required metadata fields:", paste(missing, collapse = ",")))

rna <- obj@assays[["RNA"]]
assay_slots <- slotNames(rna)
if ("layers" %in% assay_slots) {
  layer_list <- slot(rna, "layers")
  layers <- names(layer_list)
  if (!"counts" %in% layers) stop(paste("RNA counts layer missing; found", paste(layers, collapse = ",")))
  counts <- layer_list[["counts"]]
} else if ("counts" %in% assay_slots) {
  layers <- "counts"
  counts <- slot(rna, "counts")
} else stop("RNA assay has neither a counts layer nor counts slot")
if (is.null(colnames(counts))) stop("RNA counts cell IDs are missing")
if (nrow(counts) != 343L) stop(paste("Expected 343 corrected panel genes; found", nrow(counts)))
if (anyDuplicated(colnames(counts))) stop("RNA counts cell IDs are duplicated")

cross <- read.delim(file.path(out, "vannan_sample_crosswalk.tsv"), check.names = FALSE,
                    stringsAsFactors = FALSE)
primary_samples <- cross$sample_id[cross$published_primary == "YES"]
primary_donors <- unique(cross$donor_id[cross$published_primary == "YES"])
if (length(primary_samples) != 28L || length(primary_donors) != 19L) {
  stop("Frozen TMA1-4 cohort does not resolve to 28 sections / 19 donors")
}

# Apply the author-supplied VUILD105MA1 removal lists before creating any joined cell table.
target_sample <- "VUILD105MA1"
id_col <- intersect(c("cell_id", "cell", "nucleus_id"), colnames(removed_ids))[1]
if (is.na(id_col)) stop("Cannot identify cell ID column in removed-nuclei list")
id_values <- as.character(removed_ids[[id_col]])
hit_ids <- meta0$sample == target_sample & meta0$cell_id %in% id_values
xy_x <- intersect(c("X", "x", "x_centroid"), colnames(removed_xy))[1]
xy_y <- intersect(c("Y", "y", "y_centroid"), colnames(removed_xy))[1]
if (is.na(xy_x) || is.na(xy_y) || nrow(removed_xy) < 3L) stop("Cannot read author coordinate exclusion polygon")
poly <- as.matrix(removed_xy[, c(xy_x, xy_y)])
point_in_polygon <- function(x, y, polygon) {
  inside <- rep(FALSE, length(x))
  j <- nrow(polygon)
  for (i in seq_len(nrow(polygon))) {
    yi <- polygon[i, 2]; yj <- polygon[j, 2]
    xi <- polygon[i, 1]; xj <- polygon[j, 1]
    cross <- ((yi > y) != (yj > y)) &
      (x < (xj - xi) * (y - yi) / ifelse(yj == yi, .Machine$double.eps, (yj - yi)) + xi)
    inside <- xor(inside, cross)
    j <- i
  }
  inside
}
target <- meta0$sample == target_sample
hit_poly <- rep(FALSE, nrow(meta0))
hit_poly[target] <- point_in_polygon(as.numeric(meta0$x_centroid[target]),
                                     as.numeric(meta0$y_centroid[target]), poly)
hit_any <- hit_ids | hit_poly
keep <- !hit_any
meta <- meta0[keep, needed, drop = FALSE]

if (!all(meta$sample %in% cross$sample_id)) {
  stop("One or more Seurat sample IDs do not match the corrected GEO sample crosswalk")
}
meta$sample <- as.character(meta$sample)
meta$patient <- as.character(meta$patient)
meta$tma <- as.character(meta$tma)
meta$run <- as.character(meta$run)
meta$disease_status <- as.character(meta$disease_status)
meta$sample_affect <- as.character(meta$sample_affect)
meta$final_CT <- as.character(meta$final_CT)
meta$final_lineage <- as.character(meta$final_lineage)
meta$x_centroid <- as.numeric(meta$x_centroid)
meta$y_centroid <- as.numeric(meta$y_centroid)

join <- match(meta$sample, cross$sample_id)
if (anyNA(join)) stop("Unexpected unmatched sample during donor join")
if (any(as.character(meta$patient) != cross$donor_id[join])) {
  mismatch <- unique(data.frame(sample = meta$sample[meta$patient != cross$donor_id[join]],
                                patient = meta$patient[meta$patient != cross$donor_id[join]],
                                GEO_donor = cross$donor_id[join][meta$patient != cross$donor_id[join]]))
  write.table(mismatch, file.path(out, "object_crosswalk_mismatches.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  stop("Seurat patient field does not match the GEO donor crosswalk")
}

primary <- meta$tma %in% paste0("TMA", 1:4)
if (length(unique(meta$sample[primary])) != 28L || length(unique(meta$patient[primary])) != 19L) {
  stop("Corrected Seurat object does not cover the frozen TMA1-4 sample/donor set")
}
meta_out <- data.frame(
  cell_id = rownames(meta), full_cell_id = as.character(meta$full_cell_id),
  cell_id_local = as.character(meta$cell_id), sample_id = meta$sample,
  donor_id = meta$patient, TMA = meta$tma, run = meta$run,
  disease_status = meta$disease_status, affected_status = meta$sample_affect,
  original_cell_type = meta$final_CT, lineage = meta$final_lineage,
  x_um = meta$x_centroid, y_um = meta$y_centroid,
  CNiche = as.character(meta$CNiche), TNiche = as.character(meta$TNiche),
  nCount_RNA = as.numeric(meta$nCount_RNA), nFeature_RNA = as.numeric(meta$nFeature_RNA),
  stringsAsFactors = FALSE
)
meta_out <- meta_out[meta_out$sample_id %in% primary_samples, , drop = FALSE]
write.table(meta_out, gzfile(file.path(out, "vannan_primary_cell_metadata.tsv.gz"), "wt"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# Stable object summary and sample-level counts (all current GEO samples, with primary flag).
sample_rows <- lapply(sort(unique(meta$sample)), function(s) {
  ii <- which(meta$sample == s)
  ct <- meta$final_CT[ii]
  lin <- meta$final_lineage[ii]
  cr <- cross[match(s, cross$sample_id), ]
  data.frame(sample_id = s, donor_id = unique(meta$patient[ii])[1], TMA = unique(meta$tma[ii])[1],
             run = unique(meta$run[ii])[1], disease_status = unique(meta$disease_status[ii])[1],
             affected_status = unique(meta$sample_affect[ii])[1], published_primary = cr$published_primary,
             n_cells_before_exclusions = sum(meta0$sample == s), n_cells_removed_nucleus = sum(hit_ids & meta0$sample == s),
             n_cells_removed_coordinate = sum(hit_poly & meta0$sample == s), n_cells_after_exclusions = length(ii),
             n_cell_types = length(unique(ct)), n_lineages = length(unique(lin)),
             n_cells_with_finite_xy = sum(is.finite(meta$x_centroid[ii]) & is.finite(meta$y_centroid[ii])),
             n_count_cells = length(ii), n_count_sum = sum(meta$nCount_RNA[ii]),
             n_panel_genes = nrow(counts), stringsAsFactors = FALSE)
})
sample_df <- do.call(rbind, sample_rows)
write.table(sample_df, file.path(out, "sample_expression_annotation_qc.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

ct_rows <- as.data.frame(table(sample_id = meta$sample, original_cell_type = meta$final_CT,
                               lineage = meta$final_lineage), stringsAsFactors = FALSE)
ct_rows <- ct_rows[ct_rows$Freq > 0, ]
write.table(ct_rows, file.path(out, "celltype_counts_by_sample.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
ontology <- as.data.frame(table(original_cell_type = meta$final_CT, lineage = meta$final_lineage),
                          stringsAsFactors = FALSE)
ontology <- ontology[ontology$Freq > 0, ]
names(ontology)[names(ontology) == "Freq"] <- "n_cells"
write.table(ontology, file.path(out, "vannan_celltype_mapping.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

gene_registry <- data.frame(
  gene_symbol = rownames(counts), panel_version = "Xenium human lung base panel PD_277 (246 genes) + custom CVEVZD (97 genes); 343 total",
  n_genes_in_RNA_assay = nrow(counts), shared_panel_status = "COMMON_343_GENE_ASSAY",
  correction_check = ifelse(rownames(counts) %in% c("ACTA2", "MRC1"), "known corrected symbol present", "current author-corrected RDS"),
  stringsAsFactors = FALSE
)
write.table(gene_registry, file.path(out, "gene_panel_registry.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# Frozen deterministic spot-check: 1000 cells proportionally spread over every primary section.
primary_meta <- meta_out
sections <- sort(unique(primary_meta$sample_id))
quota <- rep(1000L %/% length(sections), length(sections))
quota[seq_len(1000L %% length(sections))] <- quota[seq_len(1000L %% length(sections))] + 1L
set.seed(20260928)
selected <- unlist(lapply(seq_along(sections), function(k) {
  ids <- which(primary_meta$sample_id == sections[k])
  if (length(ids) < quota[k]) stop(paste("Not enough cells for deterministic QC in", sections[k]))
  sample(ids, quota[k], replace = FALSE)
}), use.names = FALSE)
qc_meta <- primary_meta[selected, , drop = FALSE]
qc_idx <- match(qc_meta$cell_id, colnames(counts))
if (anyNA(qc_idx)) stop("Sampled metadata cell IDs do not join to RNA count columns")
qc_counts <- counts[, qc_idx, drop = FALSE]
qc_df <- data.frame(qc_meta,
                    expression_count_sum = as.numeric(Matrix::colSums(qc_counts)),
                    expression_detected_genes = as.integer(Matrix::colSums(qc_counts > 0)),
                    expression_column_present = !is.na(qc_idx),
                    coordinate_finite = is.finite(qc_meta$x_um) & is.finite(qc_meta$y_um),
                    sample_donor_match = qc_meta$sample_id %in% cross$sample_id &
                      qc_meta$donor_id == cross$donor_id[match(qc_meta$sample_id, cross$sample_id)],
                    stringsAsFactors = FALSE)
write.table(qc_df, file.path(out, "expression_coordinate_qc.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

audit <- data.frame(
  key = c("RDS_path", "RDS_bytes", "RDS_class", "R_version", "SeuratObject_version",
          "assays", "RNA_layers", "n_cells_before_exclusions", "n_cells_after_exclusions",
          "n_unique_cell_ids", "n_count_columns", "n_primary_sections", "n_primary_donors",
          "n_panel_genes", "n_exclusion_nuclei_ids", "n_exclusion_coordinate_vertices",
          "matched_removed_nuclei", "matched_removed_polygon_cells", "deterministic_qc_cells",
          "deterministic_qc_sections", "deterministic_qc_donors", "corrected_gene_ACTA2_present",
          "corrected_gene_MRC1_present"),
  value = c(rds, file.info(rds)$size, paste(class(obj), collapse = ","), as.character(getRversion()),
            as.character(packageVersion("SeuratObject")), paste(names(obj@assays), collapse = ","),
            paste(layers, collapse = ","), nrow(meta0), nrow(meta), length(unique(rownames(meta0))),
            ncol(counts), length(unique(meta$sample[primary])), length(unique(meta$patient[primary])),
            nrow(counts), length(id_values), nrow(poly), sum(hit_ids), sum(hit_poly), nrow(qc_df),
            length(unique(qc_df$sample_id)), length(unique(qc_df$donor_id)),
            "ACTA2" %in% rownames(counts), "MRC1" %in% rownames(counts)),
  stringsAsFactors = FALSE
)
write.table(audit, file.path(out, "seurat_object_audit.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Wrote corrected Seurat audit for", nrow(meta), "cells; primary:",
    length(unique(meta$sample[primary])), "sections /", length(unique(meta$patient[primary])),
    "donors; sampled QC cells:", nrow(qc_df), "\n")
