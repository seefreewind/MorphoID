#!/usr/bin/env Rscript

# Build Phase 1A's frozen 150-um region matrix and outcome-blind signature
# scores. No cell/region exclusions are made here except the frozen >=20
# lineage-cell eligibility rule for within-lineage targets.

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

root <- normalizePath(".", mustWork = TRUE)
cfg <- file.path(root, "configs")
v <- file.path(root, "results/phase0/phase0c_v")
phase0d <- file.path(root, "results/phase0/phase0d")
phase1 <- file.path(root, "results/phase1")
dir.create(phase1, recursive = TRUE, showWarnings = FALSE)

cohort <- fread(file.path(cfg, "frozen_vannan_primary_cohort.tsv"))
cohort <- cohort[primary_inclusion == "YES"]
if (nrow(cohort) != 26L || uniqueN(cohort$donor_id) != 19L) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: admitted cohort is not 26 sections / 19 donors")
}

micro <- fread(file.path(v, "microregion_feasibility.tsv"))
micro <- micro[grid_um == 150 & sample_id %chin% cohort$sample_id]
if (nrow(micro) != 12260L || uniqueN(micro$sample_id) != 26L) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: expected 12,260 frozen 150-um regions")
}
micro[, region_id := paste(sample_id, 150L, region_x, region_y, sep = "_")]
setkey(micro, sample_id, region_x, region_y)

cell_path <- file.path(v, "vannan_primary_cell_metadata.tsv.gz")
cells <- fread(cell_path, select = c(
  "cell_id", "sample_id", "donor_id", "disease_status", "affected_status",
  "TMA", "run", "x_um", "y_um", "nCount_RNA", "morphoid_coarse_ontology"
))
cells <- cells[sample_id %chin% cohort$sample_id]
if (anyDuplicated(cells$cell_id) || anyNA(cells[, .(sample_id, donor_id, x_um, y_um,
                                                       nCount_RNA, morphoid_coarse_ontology)])) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: cell metadata has duplicate IDs or missing required fields")
}
cells[, region_x := as.integer(floor(x_um / 150))]
cells[, region_y := as.integer(floor(y_um / 150))]
cells[, region_id := paste(sample_id, 150L, region_x, region_y, sep = "_")]
cells[, matrix_index := NA_integer_]

targets <- list()
within <- fread(file.path(cfg, "frozen_within_lineage_targets.tsv"))
within <- within[coverage_class %chin% c("PRIMARY_USABLE", "SECONDARY_USABLE")]
for (i in seq_len(nrow(within))) {
  r <- within[i]
  targets[[paste(r$lineage, gsub("[^A-Za-z0-9]+", "_", r$target), sep = "__")]] <- list(
    target = paste(r$lineage, gsub("[^A-Za-z0-9]+", "_", r$target), sep = "__"),
    gene_text = r$available_panel_genes,
    lineage = r$lineage,
    target_class = if (r$coverage_class == "PRIMARY_USABLE") "PRIMARY" else "SECONDARY_WITHIN_LINEAGE",
    min_cells = 20L
  )
}
programs <- fread(file.path(cfg, "frozen_phase1_candidate_programs.tsv"))
programs <- programs[phase1_target_class == "SECONDARY_USABLE"]
for (i in seq_len(nrow(programs))) {
  r <- programs[i]
  nm <- paste0("broad__", gsub("[^A-Za-z0-9]+", "_", r$candidate_program))
  targets[[nm]] <- list(target = nm, gene_text = r$available_panel_genes,
                        lineage = "all_cells", target_class = "SECONDARY_BROAD",
                        min_cells = 1L)
}
if (length(targets) < 9L) stop("Frozen primary/secondary target definitions are incomplete")

target_genes <- lapply(targets, function(z) unique(trimws(strsplit(z$gene_text, ";", fixed = TRUE)[[1]])))
target_names <- names(targets)
all_genes <- unique(unlist(target_genes, use.names = FALSE))

rds <- file.path(root, "data/raw/phase0c_v/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")
obj <- readRDS(rds)
assay <- obj[["RNA"]]
counts <- assay@counts
if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: expression matrix lacks gene/cell identifiers")
}
missing_genes <- setdiff(all_genes, rownames(counts))
coverage <- rbindlist(lapply(target_names, function(nm) {
  g <- target_genes[[nm]]
  data.table(target = nm, target_class = targets[[nm]]$target_class,
             lineage = targets[[nm]]$lineage, n_frozen_panel_genes = length(g),
             n_genes_found_in_RNA_assay = sum(g %chin% rownames(counts)),
             genes_found = paste(g[g %chin% rownames(counts)], collapse = ";"),
             genes_missing = paste(setdiff(g, rownames(counts)), collapse = ";"),
             status = if (all(g %chin% rownames(counts))) "PASS" else "HOLD_TARGET_COVERAGE")
}))
fwrite(coverage, file.path(phase1, "target_gene_coverage.tsv"), sep = "\t", na = "NA")
if (length(missing_genes)) stop("HOLD_TARGET_COVERAGE: frozen signature genes absent from corrected RNA assay")

matrix_cells <- colnames(counts)
cells[, matrix_index := match(cell_id, matrix_cells)]
if (anyNA(cells$matrix_index)) {
  stop(sprintf("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: %d frozen primary cells missing from expression object",
               sum(is.na(cells$matrix_index))))
}
if (any(cells$nCount_RNA <= 0)) stop("HOLD_TARGET_COVERAGE: nonpositive RNA library size in frozen cells")

# Validate IDs and library-size metadata on a deterministic sample of cells.
check_idx <- unique(round(seq(1, nrow(cells), length.out = min(200L, nrow(cells)))))
check_counts <- counts[, cells$matrix_index[check_idx], drop = FALSE]
count_sums <- colSums(check_counts)
if (any(abs(count_sums - cells$nCount_RNA[check_idx]) > 0)) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: expression count sums disagree with frozen metadata")
}
rm(check_counts, count_sums)

sample_map <- cohort[, .(sample_id, donor_id, disease_status, affected_status)]
region_matrix <- merge(micro, sample_map, by = c("sample_id", "donor_id", "disease_status", "affected_status"),
                       all.x = TRUE, sort = FALSE)
if (anyNA(region_matrix$donor_id)) stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: missing frozen donor metadata")
region_matrix[, disease := fifelse(tolower(disease_status) %chin% c("control", "unaffected"),
                                   "control", "pulmonary_fibrosis")]
region_matrix[, x_center_um := (region_x + 0.5) * 150]
region_matrix[, y_center_um := (region_y + 0.5) * 150]

grid_meta <- fread(file.path("results/phase0/phase0c_vr/registered_he_150um_grid_audit.tsv"))
grid_meta <- grid_meta[sample_id %chin% cohort$sample_id]
if (nrow(grid_meta) != 26L) stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: missing registered-image grid metadata")
region_matrix <- merge(region_matrix,
                       grid_meta[, .(sample_id, width_px, height_px)],
                       by = "sample_id", all.x = TRUE, sort = FALSE)
if (anyNA(region_matrix$width_px) || anyNA(region_matrix$height_px)) {
  stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: missing frozen image dimensions")
}
region_matrix[, `:=`(width_um = width_px * 0.2125, height_um = height_px * 0.2125,
                     normalized_x = x_center_um / (width_px * 0.2125),
                     normalized_y = y_center_um / (height_px * 0.2125))]

# C is the frozen Phase 0C coarse ontology. N is the mean composition of
# nonempty cardinal neighbors exactly one 150-um tile away; diagonals are
# outside the frozen 150-um Euclidean radius. No neighborhood tuning occurs.
comp <- grep("^p_", names(region_matrix), value = TRUE)
if (length(comp) != 15L) stop("Frozen composition ontology does not contain 15 coarse categories")
region_matrix[, neighborhood_id := paste(sample_id, region_x, region_y, sep = "|")]
lookup <- setNames(seq_len(nrow(region_matrix)), region_matrix$neighborhood_id)
N <- matrix(0, nrow(region_matrix), length(comp), dimnames = list(NULL, paste0("N_", comp)))
N_n <- integer(nrow(region_matrix))
for (j in seq_len(nrow(region_matrix))) {
  gx <- region_matrix$region_x[j]; gy <- region_matrix$region_y[j]
  keys <- paste(region_matrix$sample_id[j], c(gx - 1L, gx + 1L, gx, gx),
                c(gy, gy, gy - 1L, gy + 1L), sep = "|")
  ix <- unname(lookup[keys]); ix <- ix[!is.na(ix)]
  if (length(ix)) {
    N[j, ] <- colMeans(as.matrix(region_matrix[ix, ..comp]))
    N_n[j] <- length(ix)
  }
}
region_matrix <- cbind(region_matrix, as.data.table(N))
region_matrix[, n_primary_neighborhood_regions := N_n]

scores <- data.table(region_id = region_matrix$region_id)
for (nm in target_names) {
  spec <- targets[[nm]]
  keep <- if (spec$lineage == "all_cells") rep(TRUE, nrow(cells)) else if (spec$lineage == "macrophage") {
    cells$morphoid_coarse_ontology %chin% c("macrophage", "monocyte")
  } else {
    cells$morphoid_coarse_ontology == spec$lineage
  }
  ci <- which(keep)
  gi <- match(target_genes[[nm]], rownames(counts))
  raw <- t(counts[gi, cells$matrix_index[ci], drop = FALSE])
  norm <- log1p(raw * (10000 / cells$nCount_RNA[ci]))
  keys <- cells$region_id[ci]
  # Aggregate mean log-normalized expression over the fixed genes and the
  # frozen target lineage; no outcome-dependent genes or cells are selected.
  sums <- rowsum(norm, group = keys, reorder = FALSE)
  cell_n <- as.numeric(table(factor(keys, levels = rownames(sums))))
  target_score <- rowMeans(sums) / cell_n
  names(target_score) <- rownames(sums)
  score_col <- paste0("Y_", nm)
  scores[[score_col]] <- unname(target_score[match(region_matrix$region_id, names(target_score))])
  count_col <- paste0("n_target_cells__", nm)
  n_by_region <- table(keys)
  region_matrix[[count_col]] <- as.integer(n_by_region[match(region_matrix$region_id, names(n_by_region))])
  region_matrix[is.na(get(count_col)), (count_col) := 0L]
  region_matrix[, (paste0("target_class__", nm)) := spec$target_class]
  rm(raw, norm, sums, cell_n, target_score)
  gc(verbose = FALSE)
}

region_matrix <- cbind(region_matrix, scores[, -1, with = FALSE])
image_manifest <- fread(file.path(phase0d, "d5_color_feature_extraction_manifest.tsv"))
image_manifest <- image_manifest[sample_id %chin% cohort$sample_id]
region_matrix <- merge(region_matrix,
                       image_manifest[, .(sample_id, image_file, registered_HE_sha256,
                                          fixed_scale_um_per_px)],
                       by = "sample_id", all.x = TRUE, sort = FALSE)
if (anyNA(region_matrix$image_file)) stop("TECHNICAL_HOLD_FROZEN_INPUT_MISMATCH: missing image reference")
region_matrix[, tissue_fraction := NA_real_] # populated from fixed 150-um image crops in feature extraction
region_matrix[, image_crop_reference := paste(image_file, region_x, region_y, sep = "#")]

fwrite(region_matrix, file.path(phase1, "region_analysis_matrix.tsv"),
       sep = "\t", na = "NA")

exclusions <- rbindlist(lapply(target_names, function(nm) {
  if (!identical(targets[[nm]]$target_class, "PRIMARY")) return(NULL)
  ncol <- paste0("n_target_cells__", nm)
  x <- region_matrix[get(ncol) < 20L,
                     .(region_id, sample_id, donor_id, target = nm,
                       observed_target_lineage_cells = get(ncol), exclusion_stage = "1B",
                       exclusion_reason = "BELOW_FROZEN_20_CELL_PRIMARY_THRESHOLD",
                       sensitivity_only_10_cell_threshold = ifelse(get(ncol) >= 10L, "YES", "NO"))]
  x
}), fill = TRUE)
if (is.null(exclusions) || nrow(exclusions) == 0L) {
  exclusions <- data.table(region_id = character(), sample_id = character(), donor_id = character(),
                           target = character(), observed_target_lineage_cells = integer(),
                           exclusion_stage = character(), exclusion_reason = character(),
                           sensitivity_only_10_cell_threshold = character())
}
fwrite(exclusions, file.path(phase1, "region_exclusion_log.tsv"), sep = "\t", na = "NA")

region_summary <- rbindlist(lapply(target_names, function(nm) {
  spec <- targets[[nm]]
  count_col <- paste0("n_target_cells__", nm)
  score_col <- paste0("Y_", nm)
  eligible <- if (spec$min_cells > 1L) region_matrix[get(count_col) >= spec$min_cells & is.finite(get(score_col))] else
    region_matrix[is.finite(get(score_col))]
  eligible[, .(target = nm, target_class = spec$target_class,
               n_donors = uniqueN(donor_id), n_sections = uniqueN(sample_id),
               n_regions = .N, median_target_lineage_cells = median(get(count_col)),
               minimum_cells = spec$min_cells)]
}), fill = TRUE)
fwrite(region_summary, file.path(phase1, "phase1a_target_feasibility.tsv"), sep = "\t", na = "NA")

cat(sprintf("PHASE1A_MATRIX_PASS: regions=%d; donors=%d; sections=%d; targets=%d\n",
            nrow(region_matrix), uniqueN(region_matrix$donor_id),
            uniqueN(region_matrix$sample_id), length(targets)))
if (length(warnings())) {
  cat("R_WARNINGS_BEGIN\n")
  print(warnings())
  cat("R_WARNINGS_END\n")
}
