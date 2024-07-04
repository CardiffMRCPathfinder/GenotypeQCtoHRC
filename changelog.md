# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## 2024-07-04

## Added

- Beta release in GitHub (public).

## 2024-07-03

## Added

- QC3: Flagged the possibility of erroneous Y chromosome calls in female samples from old Illumina arrays.
- Custom function for SNP renaming after LiftOver brought back and improved (see **2024-04-25** and **2024-04-18**).

## Changed

- QC1: Reordered internal scripts to process samples and variants in the order preferred by PLINK (samples first).
- QC2: Wahlund effect correction (**--qc2-popmix**) is now turned off by default.
- QC5: Default reference panel (**--qc5-ref**) changed to 1000 Genomes while AADR performance in legacy arrays is investigated.

## Fixed

- QC0: Now accepts and correctly processes legacy PLINK files (pre-0.99 bed/bim/fam format).
- QC1: Corrected bug caused by phenotype data being present.
- QC3: Now correctly skipped in the absence of sex chromosome data.
- Imputation-ready files now use the correct VCF format (v4.2).
- Imputation-ready files for the X chromosome now handle heterozygous haploid hardcalls via PLINK2 **--set-invalid-haploid-missing**.
- Various typo and code corrections to HTML report.

## Removed

- Genotype Harmonizer removed from the SNP renaming target due to its handling of ambiguous SNPs.

## 2024-04-25

### Added

- Standalone flag to request SNP renaming (**--rename**). 
- Standalone flag to request harmonisation to an imputation reference (**--gh**).

### Changed

- dbSNP reference files changed to "sites-only" VCFs for compatibility with Genotype Harmonizer.
- Reordered content of *_targets.R* and *GenotypeQCtoHRC_targets.Rmd* for better readability.

### Fixed

- Integrity checks were stopping the pipeline unnecessarily, changed to a full problem report on the RMarkdown output.
- QC3: Removed monomorphic variants when computing the count of Y-linked markers to avoid a bug seen on legacy data.

### Removed

- Custom SNP rename function as it wasn't performing as well as Genotype Harmonizer. We hardly knew you!

## 2024-04-18

### Added

- Initial QC step to convert old PLINK1 files to PLINK2 format, fixing potential errors ("QC0").
- Commandline flag (**--clean**) to force the pipeline to run from scratch ignoring remaining files from previous runs.

### Changed

- Chipendium and LiftOver turned off by default.
- Short dataset name (**--shortname**) is automatically filled out from long name if not provided.

### Fixed

- RMarkdown cleanup to remove the remnants of an old unused .ped/.map conversion function.
- QC1: Reinforced dependencies between steps to avoid a potential error.
- PLINK files with "split chromosomes" (i.e. unsorted variants) can now be loaded and processed.
- PLINK files with unnamed chromosomes (i.e. chromosome code zero) can now be loaded and processed.

## 2024-04-04

### Added

- LiftOver function, based on the implementation in the rtracklayer package.
- Commandline flags (**--lo**) to enable LiftOver.
- SNP renaming target based on PLINK2 **--recover-var-ids** command, will probably see wider use in future releases.
- Some more assertions to make sure input flags are correct and consistent with each other.

### Changed

- Refactored custom functions into several scripts to make sure they don't get unwieldy.
- Refactored Chipendium closest-match code outside of .Rmd file and into its own function.
- QC5: Functions automatically work with post-LiftOver data if input data is not in build GRCh37.

### Fixed

- QC5: Code to exclude variants with missing IDs was error-prone, changed to use PLINK2 **--set-missing-var-ids** command.

## 2024-03-26

### Added

- DRAGON-Data ASCII logo (as much dragon as I could fit in 80-character lines).
- Welcome and progress messages, there *might* be an Easter Egg somewhere.
- Citation information at the end of each pipeline run.

## 2024-03-22

### Added

- Support for UNIX and Mac systems in all R functions calling PLINK.
- Dynamic naming (using **--name** flag) of pipeline HTML reports.
- QC5: Support for 1000 Genomes reference panel.

### Changed

- QC5: Dynamic legends (variable numbers of ancestries) in PCA plots.
- QC5: Dynamic LDA model results table.

### Fixed

- QC5: Better detection of whether genotype IDs were changed by GDS conversion.

## 2024-03-17

### Added

- Alpha release (internal).
