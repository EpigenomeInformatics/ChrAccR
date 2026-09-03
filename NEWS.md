# ChrAccR 0.9.26

## New features

* Allele-specific chromatin accessibility. The new `DsASC` class inherits from
  `DsAcc` and holds reference and alternative counts at heterozygous SNPs
  alongside the ordinary peak accessibility matrix. `buildMasterSNPs()` builds a
  shared het-site list across donors and `DsASC.gatk()` counts alleles against it
  from WASP-corrected BAM files.
* `calcASCStatistics()` tests each site for allelic imbalance with a two-sided
  binomial test against 0.5 and per-sample BH correction.
* Helpers for the steps between counting and testing: `ascDropHomozygous()`,
  `ascFilterRegions()`, `ascPoolReplicates()` and `filterForRecurrence()`.
* Cross-cell-type sharing with `ascAggregate()`, `ascSitePosterior()`,
  `ascSharing()` and `estimateSharedImbalance()`.
* Plotting: `plotASCBalance()`, `plotASCVolcano()` and `plotASCManhattan()`.
* `prepareMotifmatchr()` now supports JASPAR2020 (`motifs = "jaspar2020"`).
* New vignette, "Allele-specific chromatin accessibility", with a small example
  dataset in `inst/extdata` built by `inst/scripts/make_vignette_data.R`.

## Bug fixes

* Generics are declared once in `R/AllGenerics.R` instead of behind
  `if (!isGeneric(...))` guards in each class file. The guards resolved the name
  through the search path, so a second `devtools::load_all()` in the same session
  skipped `setGeneric()` and the following `setMethod()` failed with "no existing
  definition". `devtools::check()` is repeatable now.
* `prepareMotifmatchr()` read the species with `provider()`, which returns the
  BSgenome provider ("UCSC") rather than the organism, and `getMatrixSet()` then
  failed. It uses `BiocGenerics::organism()` instead.
* `filterForRecurrence()` averaged `log2FC_norm`, a column `calcASCStatistics()`
  does not produce. It now uses `log2FC_norm` when present and `log2FC`
  otherwise.
* `filterForRecurrence()` no longer calls `stringr::str_extract()` when parsing a
  donor out of `sampleId`; `stringr` is not a declared dependency.

## Documentation

* Added `NEWS.md` and `inst/CITATION`, so `citation("ChrAccR")` returns a proper
  entry.
* `DESCRIPTION` gains `URL` and `BugReports`, which pkgdown uses for the site and
  for cross-package links.
* Removed the duplicate `pkgdown/_pkgdown.yml`; only the root `_pkgdown.yml` was
  ever read.
* The `overview` and `singlecell` vignettes call `getChromVarDev()` with
  `motifs = "jaspar2018"`. The old `"jaspar"` no longer matches a motif set now
  that the JASPAR versions are selected by name.
