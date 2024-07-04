# DRAGON-Data "GenotypeQCtoHRC" pipeline v.2.0

## Description
This is an updated version of the "GenotypeQCtoHRC" pipeline introduced in the DRAGON-Data [publication](https://doi.org/10.1192%2Fbjo.2022.636) and supported by an MRC "Mental Health Data Pathfinder" [grant](https://gtr.ukri.org/projects?ref=MC_PC_17212) to Cardiff University. The pipeline was originally designed to perform stringent quality control of genotype data, and prepare it for imputation using the Michigan Imputation Server [platform](https://imputationserver.sph.umich.edu/).
 
The purpose of this update is to simplify the use of the pipeline by streamlining most of its functions via literate and reproducible programming techniques. While the previous version relied on a single RMarkdown file with a lot of dependencies and hardcodings, this version makes heavy use of the **{targets}** [framework](https://books.ropensci.org/targets/) for creating R workflows.

## Roadmap

Currently, this is an *BETA* version. All functions enabling genotype quality control described in the paper are working, with some additions.

- Genotype QC :white_check_mark:
- Population structure and relatedness checks :white_check_mark:
- Biogeographical ancestry inference :white_check_mark:
- Use of [1000 Genomes Phase 3](https://doi.org/10.1038/nature15393) as ancestral reference panel :white_check_mark:
- Liftover (GRCh36/GRCh37/GRCh38) :white_check_mark:
- Data preparation for Michigan Imputation Server :white_check_mark:
- Data preparation for TopMED Imputation Server :white_check_mark:
- Update of the public GitHub website :white_check_mark:
- Containerisation (Docker/Apptainer) :white_large_square:

## Making it work

The expected input is a genotype dataset in binary [PLINK1](https://www.cog-genomics.org/plink/1.9/data#make_bed) *.bed/.bim/.fam* format. This has to be in the same folder as the **[GenotypeQCtoHRC.R](GenotypeQCtoHRC.R)** file.
The expected output is a QC report in HTML format and a folder with original, curated and intermediate files. Download **[GenotypeQCtoHRC_targets_example.html](GenotypeQCtoHRC_targets_example.html)** and open it with your internet browser to see how this looks like.

The main command to run the pipeline is: *Rscript GenotypeQCtoHRC.R*

You can add the *-h* or *--help* flag to display the command line arguments and their potential options.
Required command line arguments for proper execution are:
- *--file* to indicate the stem (i.e. name without extension) of the PLINK file set to analyse.
- *--name* to indicate the name of the dataset that will be used in folder names and in the report.
- *--shortname* to indicate a shorter (i.e. 10 characters) name for the dataset to be used in plots.

You can also request some additional procedures by setting the following flags to TRUE:
- *--chip* to run the "Chipendium" method to identify genotyping platforms (only valid with raw array data).
- *--lo* to run [LiftOver](https://bioconductor.org/packages/release/workflows/html/liftOver.html) from build *--lo-in* (set to 0 if unknown) to build *--lo-out*. See help for more details.
- *--gh* to run [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) and align the markers to the *--gh-ref* imputation reference. Options for the latter are [HRC](https://imputationserver.sph.umich.edu) (GRCh37) or [TopMED](https://imputation.biodatacatalyst.nhlbi.nih.gov/) (GRCh38).
- *--rename* to run [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) and assign dbSNP RS IDs to the markers in the input dataset.
- *--clean* to remove all data from previous pipeline runs and start from scratch (useful for troubleshooting).

Everything else can be left at a default state to run the QC script with the same parameters used for the DRAGON-Data [publication](https://doi.org/10.1192%2Fbjo.2022.636). 

On a potato-grade computer, a dataset with ~2K individuals and ~500K markers should be analysed in about an hour. The lion' share of that runtime is at the moment taken by Chipendium. This uses multiple processors when available but won't run unless requested. Harmonising the dataset to the TopMED reference (also optional) can take a long time too.

For a quick test, the **[Behar et al. 2010](https://doi.org/10.1038/nature09103)** dataset used in the [example report](GenotypeQCtoHRC_targets_example.html) can be downloaded from the Estonian Biocentre [website](https://evolbio.ut.ee/jew/).

## Dependencies

The pipeline was originally developed in a Windows system using **[R](https://cran.r-project.org/)** v4.3 and [RStudio](https://posit.co/products/open-source/rstudio/) v2023.12. It has since been tested in Windows, Mac and Unix using **R v4.4**. The following packages are required to run it from source:
- CRAN: rmarkdown[^1], optparse, targets, tarchetypes, here, pbapply, tidyverse, parallelly, VGAM, AssocTests, caret, probably, tint, glue, colorspace, scales, scattermore, ggpubr.
- GitHub: [slugify](https://github.com/cannin/slugify/).
- Bioconductor: BiocParallel, GWASTools, GENESIS, SNPRelate, rtracklayer.

As some pipeline requirements (software[^2] and reference files) cannot be shared in this repository, the *dragondata_extra* folder is supplied partially filled. Please read the **README** files within for full instructions on how to populate it.
This is a provisional solution and a better method for file and depencency sharing will be implemented soon.

## Credits

- Original *GenotypeQCtoHRC* code by Leon Hubbard.
- Computational support by Lynsey Hall.
- Testing by Jack Underwood, Djenifer Kappel and Jessica Yang.
- Targets [package](https://docs.ropensci.org/targets/) by Will Landau.
- "Tint is not Tufte" [HTML theme](https://github.com/eddelbuettel/tint) by Dirk Eddelbuettel.
- Chipendium [procedure](https://www.chg.ox.ac.uk/~wrayner/strand/) by Will Rayner. 

## Notes

[^1]: The RMarkdown [setup](https://bookdown.org/yihui/rmarkdown/installation.html) can be problematic in some HPC servers. Please ensure the linked [Pandoc](https://pandoc.org/) installation is v2.5+.
[^2]: This includes PLINK (v1.9 and v2) and Genotype Harmonizer v1.4.20+, which requires Java. Please check that these files have the appropriate permissions to be executed once you copy them to your local or HPC system.
