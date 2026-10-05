# Phase 1D-R2 methods amendment

Before any Phase 1 predictive modeling, regions with zero valid H&E pixels after application of the pre-specified black-pixel filter were classified as technical image-input failures and excluded from the common analysis universe. This rule depended only on frozen image pixels and crop coordinates and did not use molecular outcomes, disease labels, donor identity, model predictions, or residuals. No imputation or modification of the frozen color-feature definition was performed.

This eligibility rule was introduced as a protocol amendment after the input-completeness audit and before any Phase 1 predictive results were obtained; it was not part of the original Phase 0D exclusion policy. Valid sampled pixels were counted using the original whole-slide sampling origin and stride, physical grid assignment, RGB normalization and black-pixel threshold. The common universe was applied identically to all directly compared model arms; target eligibility was evaluated afterward using the unchanged lineage-cell threshold and target-availability criteria.

Of 12,260 regions, 116 unique regions were technically excluded (77 coordinate-support failures and 39 zero-valid-image regions; overlap 0). The common universe contained 12144 regions from 26 sections and 19 donors. Original color rows and morphology caches were preserved.
