# MorphoID reproducibility code and frozen workflows

This version 2.0.0 clean reproducibility package accompanies the MorphoID conditional identifiability manuscript. It contains original code, frozen configurations, target/fold definitions, source-derived regional inputs, audit metadata, environment records and regeneration instructions. The public code repository is https://github.com/seefreewind/MorphoID. The complete release archive contains the regional data at their original relative paths. Scientific results have not been rerun or altered.

# Phikon-v2 dependency and licensing

Phikon-v2 is developed by Owkin. Obtain the official model and processor independently and accept the Phikon-v2 Non-Commercial License before local use. Weights and embeddings are not bundled in this default clean package. The extraction wrapper regenerates local embeddings using the unchanged frozen extractor, supplied model files and verified hashes. This package does not grant access to raw data or model weights.

Owkin has provided express written permission for aggregate manuscript results, tables and figures to be published under the journal's applicable publication licence. Embeddings, if separately distributed, remain Results under the same Phikon-v2 terms for non-commercial use, including section 3(c) ShareAlike; MIT does not apply to them. Excluding embeddings here is a packaging decision. Private correspondence is not included.

Cite Filiot A, Jacob P, Mac Kain A, Saillard C (2024), Phikon-v2, A large and public feature extractor for biomarker prediction, arXiv:2409.09173, https://doi.org/10.48550/arXiv.2409.09173.

## Reproduction flow

1. Obtain GEO GSE250346 source data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346 and retain the source terms and access receipts.
2. Obtain official Owkin Phikon-v2 and its processor from https://huggingface.co/owkin/phikon-v2/tree/2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd.
3. Accept the applicable Owkin licence directly with the provider.
4. Verify frozen revision 2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd and checkpoint SHA256 261ae680fa699b3b951597fd57aa19c02ef735805acb104b93af69b36d928569; verify processor/config hashes in docs/PHIKON_REGENERATION.md.
5. Map the historical raw-data paths to a separate writable study copy; reconstruct prerequisite matrices, crop audits and hash-baseline inputs. Regenerate embeddings locally using scripts/regenerate_phikon_embeddings.py with --model-dir, --revision and --accepted-provider-terms.
6. Run the frozen common-universe analysis and Phase1E–1J sequence documented in docs/PHIKON_REGENERATION.md. Preserve all donor folds, independent PF-only fits, negative R² and training-side residual cross-fitting.
7. Run the original figure and report scripts to regenerate tables and figures locally. The final finite-donor inference limitation remains part of the scientific conclusion.

A portable end-to-end reproduction has not been executed for this rights update. The adapter only handles the duplicated full_crop_in_bounds merge input in memory after equality validation; scientific extractor bytes, coordinates, row order and estimands remain unchanged. Raw H&E images, expression RDS objects and original access receipts must be obtained/reconstructed by the user.

## License and file inventory

Read LICENSE_SCOPE.md, LICENSE_CODE, NOTICE_PHIKON_V2.md and PACKAGE_MANIFEST.tsv. Original software is MIT; source-derived regional inputs remain under their source terms. No blanket licence covers third-party data. Aggregate result files and figure assets are not included in this code package. Their absence is a packaging decision after Owkin's permission, not an unresolved article-rights restriction.

## Historical release

The historical release https://github.com/seefreewind/MorphoID/releases/tag/v1.0.0 and DOI https://doi.org/10.5281/zenodo.23150000 contain 26 section embedding caches. The 2026-10-26 audit found explicit Owkin non-commercial/ShareAlike scope notices, with original code separately MIT: NO MATERIAL LICENSE DEFECT for the embedding licence scope. They remain public and are not deprecated or restricted. This separate clean package is a different, smaller package design; a separate standalone Zenodo record archives this clean package and does not replace the historical record.

Project attribution and archive creator: MorphoID. Funding: none. Competing interests: none.

## Release contents

The GitHub working tree contains code, configurations, environment records and documentation. The accompanying `MorphoID-v2.0.0-clean.zip` release asset contains the complete package, including regional analysis inputs and audit tables under `results/`. Extract it before running the workflows. `PACKAGE_MANIFEST.tsv` records SHA256 checksums for all package files except itself. Historical inventory tables retain their original audit scope. No model weights, embeddings, raw tissue images or private correspondence are included.
