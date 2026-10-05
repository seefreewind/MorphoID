#!/usr/bin/env python3
"""Coverage-only audit of predeclared programs against the corrected 343-gene panel."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/phase0/phase0c_v"
SRC = BASE / "source_metadata"

def read_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                result[fields[0]] = {gene.upper() for gene in fields[2:] if gene}
    return result

hallmark = read_gmt(SRC / "h.all.v2025.1.Hs.symbols.gmt")
c7 = read_gmt(SRC / "c7.all.v2025.1.Hs.symbols.gmt")
reactome = read_gmt(SRC / "c2.cp.reactome.v2025.1.Hs.symbols.gmt")
gobp = read_gmt(SRC / "c5.go.bp.v2025.1.Hs.symbols.gmt")
genes = set()
with (BASE / "gene_panel_registry.tsv").open(encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        genes.add(row["gene_symbol"].upper())
if len(genes) != 343:
    raise RuntimeError(f"Expected a frozen 343-gene panel, found {len(genes)} unique symbols")

# The eight MorphoID programs are mapped to versioned public gene sets before
# any expression summarization. These are counts only; no score is calculated.
core = [
    ("proliferation", "HALLMARK_G2M_CHECKPOINT", hallmark, "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark v2025.1.Hs"),
    ("hypoxia", "HALLMARK_HYPOXIA", hallmark, "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark v2025.1.Hs"),
    ("EMT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", hallmark, "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark v2025.1.Hs"),
    ("ECM remodeling", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", reactome, "c2.cp.reactome.v2025.1.Hs.symbols.gmt", "MSigDB Reactome v2025.1.Hs"),
    ("IFN-gamma response", "HALLMARK_INTERFERON_GAMMA_RESPONSE", hallmark, "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark v2025.1.Hs"),
    ("T-cell exhaustion", "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP", c7, "c7.all.v2025.1.Hs.symbols.gmt", "MSigDB C7 v2025.1.Hs; Wherry et al. PMID 17950003; source organism mouse"),
    ("cytotoxicity", "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", gobp, "c5.go.bp.v2025.1.Hs.symbols.gmt", "MSigDB GO BP v2025.1.Hs"),
    ("inflammatory response", "HALLMARK_INFLAMMATORY_RESPONSE", hallmark, "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark v2025.1.Hs"),
]

# These short lists are explicitly literature-marker availability panels, not
# validated directional modules and not expression scores.
lung = [
    ("alveolar epithelial injury", ["KRT8","KRT18","KRT19","SFN","CLDN4","LGALS3","ITGB6","MMP7","SOX4","NUPR1"], "Strunz et al., Nat Commun 2020; https://doi.org/10.1038/s41467-020-17358-3"),
    ("transitional epithelial state", ["KRT8","KRT18","KRT19","SFN","CLDN4","LGALS3","TACSTD2","ITGB6","SOX4","NUPR1","CDKN1A"], "Strunz et al., Nat Commun 2020; https://doi.org/10.1038/s41467-020-17358-3"),
    ("aberrant basaloid state", ["TP63","KRT17","KRT5","KRT15","LAMB3","LAMC2","VIM","CDH2","FN1","COL1A1","TNC","HMGA2","CDKN1A","CDKN2A","CCND1","CCND2","MDM2","MMP7","ITGB6","EPHB2"], "Adams et al., Sci Adv 2020; https://pmc.ncbi.nlm.nih.gov/articles/PMC7439502/"),
    ("fibroblast activation", ["CTHRC1","POSTN","COL1A1","COL1A2","COL3A1","FN1","TNC","ACTA2","TAGLN","LOX","MMP11","COL5A2","COL12A1"], "Coverage-marker panel informed by Vannan et al., Nat Genet 2025; https://doi.org/10.1038/s41588-025-02080-x"),
    ("myofibroblast", ["ACTA2","TAGLN","MYL9","CNN1","TPM1","TPM2","CALD1","MYH11","CTHRC1","POSTN"], "Coverage-marker panel informed by Vannan et al., Nat Genet 2025; https://doi.org/10.1038/s41588-025-02080-x"),
    ("ECM deposition", sorted(reactome["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]), "MSigDB Reactome v2025.1.Hs; coverage only"),
    ("macrophage inflammatory", ["IL1B","TNF","IL6","CCL2","CXCL8","NFKBIA","S100A8","S100A9","IL1RN","SPP1"], "Coverage-marker panel informed by Vannan et al., Nat Genet 2025; https://doi.org/10.1038/s41588-025-02080-x"),
    ("profibrotic macrophage", ["SPP1","GPNMB","TREM2","CCL18","APOE","LGALS3","MMP9","CHI3L1","MRC1","IL1RN"], "Coverage-marker panel informed by Vannan et al., Nat Genet 2025; https://doi.org/10.1038/s41588-025-02080-x"),
    ("interferon response", sorted(hallmark["HALLMARK_INTERFERON_GAMMA_RESPONSE"]), "MSigDB Hallmark v2025.1.Hs; coverage only"),
    ("proliferation", sorted(hallmark["HALLMARK_G2M_CHECKPOINT"]), "MSigDB Hallmark v2025.1.Hs; coverage only"),
]

definitions = []
for name, set_name, db, file_name, source in core:
    if set_name not in db:
        raise KeyError(f"Missing predeclared gene set {set_name}")
    definitions.append((name, set_name, db[set_name], "predefined_public_signature", file_name, source))
for name, members, source in lung:
    definitions.append((name, "LITERATURE_MARKER_PANEL", {g.upper() for g in members}, "coverage_only_literature_marker_panel", "manual versioned list in this script", source))

coverage_rows = []
cross_rows = []
for domain, subset in (("core", definitions[:8]), ("lung_specific", definitions[8:])):
    for name, set_name, members, kind, origin, source in subset:
        retained = sorted(members & genes)
        fraction = len(retained) / len(members) if members else 0.0
        status = "GOOD" if fraction >= 0.50 else "PARTIAL" if fraction >= 0.25 else "POOR"
        coverage_rows.append({
            "domain": domain, "program": name, "signature_id": set_name,
            "signature_type": kind, "source_file_or_definition": origin,
            "source": source, "n_total_unique_genes": len(members),
            "n_panel_genes": len(retained), "coverage_fraction": f"{fraction:.6f}",
            "coverage_status": status, "retained_genes": ";".join(retained),
            "interpretation": "availability only; no score, ranking, or expression-derived gene selection",
        })
        for symbol in sorted(members):
            cross_rows.append({"domain":domain,"program":name,"signature_id":set_name,
                               "gene_symbol":symbol,"in_panel":"YES" if symbol in genes else "NO"})

for name, rows in (("signature_coverage.tsv", coverage_rows), ("signature_gene_crosswalk.tsv", cross_rows)):
    with (BASE / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
print(f"Coverage-only gene audit complete: {len(coverage_rows)} programs; panel genes={len(genes)}")
