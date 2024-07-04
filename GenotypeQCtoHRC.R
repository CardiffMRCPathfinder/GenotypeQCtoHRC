# LICENCE
#Copyright 2024 Cardiff University
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

# Setup
suppressWarnings(library(optparse))
suppressWarnings(library(targets))
suppressWarnings(library(tarchetypes))

# Parse arguments
option_list = list(
  make_option("--file", type="character", default="jew_paper_data_dbSNP-b131_pos-b37_1KG_strand", 
              help="PLINK (.bed/.bim/.fam) file stem.", metavar="character"),
  make_option("--name", type="character", default="Behar et al. 2013", 
              help="Dataset name to be used in file/folder names and the QC report. Use quotes if spaces present.", metavar="character"),
  make_option("--shortname", type="character", default=NA_character_, 
              help="Short dataset name (10 characters maximum) to be used in plots. Use quotes if spaces present.", metavar="character"),
  make_option("--tinydata", type="logical", default = F, action = "store_true", 
              help = "Disable file size checks for tiny (<1 Mb) datasets.", metavar="TRUE/FALSE"),
  make_option("--parallel", type="logical", default=T, 
              help="Parallelise pipeline whenever possible via the 'parallelly' package?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--clean", type="logical", default=F, 
              help="Remove any leftovers from previous pipeline runs? (useful for troubleshooting)\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--check", type="logical", default=T, 
              help="Perform dataset integrity checks?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--chip", type="logical", default=F, 
              help="Infer genotyping chip using the Chipendium procedure?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--rename", type="logical", default=F, 
              help="Rename markers to dbSNP RS identifiers?\n\t\tNote this will remove variants with alleles that do not match the dbSNP definition, as well as variants with only one named allele.\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--lo", type="logical", default=F, 
              help="LiftOver genotype data?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--lo-in", type="integer", default=37, 
              help="Genomic build of the input data in GRC assembly number. Set to 0 to automatically infer via Chipendium.\n\t\t[default= %default]", metavar="36-38"),
  make_option("--lo-out", type="integer", default=38, 
              help="Genomic build of the output data in GRC assembly number.\n\t\t[default= %default]", metavar="36-38"),
  make_option("--gh", type="logical", default=F, 
              help="Harmonise PLINK fileset to imputation reference?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--gh-ref", type="character", default="HRC", 
              help="Imputation reference, either HRC (GRCh37) or TopMED (GRCh38).\n\t\t[default= %default]", metavar="HRC/TopMED"),
  make_option("--qc1", type="logical", default=T, 
              help="Perform marker- and sample-level missingness checks?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc1-call-rate", type="numeric", default=0.95, 
              help="Minimum marker call rate.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc1-coverage", type="numeric", default=0.95, 
              help="Minimum sample genotyping coverage.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc2", type="logical", default=T, 
              help="Perform Hardy-Weinberg Equilibrium (HWE) checks?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc2-popmix", type="logical", default=F, 
              help="Account for population mixture (can reduce heterozygote count) in HWE test p-values?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc2-hwepval", type="numeric", default=0.000001, 
              help="Genome-wide significant threshold for deviations from HWE.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc3", type="logical", default=T, 
              help="Perform sex chromosome checks?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc3-max-female-F", type="numeric", default=0.2, 
              help="Maximum female X-linked homozygosity.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc3-min-male-F", type="numeric", default=0.8, 
              help="Minimum male X-linked homozygosity.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc4", type="logical", default=F, 
              help="Perform population structure and relatedness checks?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc4-maf", type="numeric", default=0.05, 
              help="Minor Allele Frequency used to distinguish common from rare variants.\n\t\tValues of at least MAF=1/sqrt(2*N) are reccomended, where N= number of genotyped individuals.\n\t\t[default= %default]", metavar="0-1"),
  make_option("--qc4-pcair", type="integer", default=10, 
              help="Number of PCs to extract from the PC-AiR computation.\n\t\t[default= %default]", metavar="1-100"),
  make_option("--qc4-pcrelate", type="integer", default=5, 
              help="Number of PCs used to correct PC-Relate calculations.\n\t\t[default= %default]", metavar="1-100"),
  make_option("--qc4-kin", type="numeric", default=0.09375, 
              help="Minimum value of the coefficient of relationship (kinship x 2) to flag pairs of samples as related.\n\t\t[default= %default; lower bound for third degree relatives]", metavar="0-1"),
  make_option("--qc5", type="logical", default=T, 
              help="Infer the biogeographical genetic ancestry of the samples?\n\t\t[default= %default]", metavar="TRUE/FALSE"),
  make_option("--qc5-ref", type="character", default="1KGP", 
              help="Reference sample for ancestry inference.\n\t\t[default= %default]", metavar="1KGP"),
  make_option("--extrafolder", type="character", default="dragondata_extra", 
              help="Name of the folder containing auxiliary software.", metavar="character"),
  make_option("--scriptsfolder", type="character", default="dragondata_scripts", 
              help="Name of the folder containing DRAGON-Data scripts.", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser,convert_hyphens_to_underscores = T)

# Detect whether file name and dataset name have been provided
# STOP if not
if (is.null(opt$file) || is.null(opt$name)){
  print_help(opt_parser)
  stop("At least one mandatory argument (file name or dataset name) is missing!", call.=FALSE)
}

# Define the pipeline steps to be performed
#qc.actions <-  names(which(unlist(opt)==T)) |> append(c("bfile","report")) 
#qc.actions <- qc.actions[-grep("_",qc.actions)]

# Pass the arguments into the targets pipeline
saveRDS(opt,"_opt.rds")

# Print DRAGON-Data initiation message
message("Welcome to the DRAGON-Data Genotype Quality Control pipeline.")
message("This software was made at the Centre for Neuropsychiatric Genetics and Genomics of Cardiff University.")
message("We acknowledge support from UK Research and Innovation through Medical Research Council project grant MC_PC_17212.")
message("Please see https://www.cardiff.ac.uk/centre-neuropsychiatric-genetics-genomics for more information about our work.")

# Clean pipeline remnants
if(isTRUE(opt$clean)){
  tar_destroy(ask=F)
}

# Run pipeline
#tar_make(starts_with(!!qc.actions),seconds_meta_append=15)
tar_make(seconds_meta_append=15)

# End message
Sys.sleep(2)
message("Pipeline run ended! See ", slugify::slugify(opt$name), ".qc_report.html for results.") 
message("Please cite Lynham et al. 2023 (DOI:10.1192/bjo.2022.636) if you use this software in a publication.")
