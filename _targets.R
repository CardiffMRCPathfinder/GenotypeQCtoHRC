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

## Read commandline arguments
opt <- readRDS("_opt.rds")

## Load libraries
options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(here)
library(slugify) # Package by Agustin Luna (https://github.com/cannin)
library(pbapply)
library(tidyverse)

## Library options
pboptions(type = "txt")

## Auxiliary and sourced functions
`%nin%` <- Negate(`%in%`)
source(here(opt$scriptsfolder, "qc_functions.R"))
source(here(opt$scriptsfolder, "aim_functions.R"))
source(here(opt$scriptsfolder, "misc_functions.R"))
source(here(opt$scriptsfolder, "liftover_rtracklayer.R"))

## Display DRAGON-Data logo
source(here(opt$scriptsfolder, "dragondata_logo.R"))
cat(dragondata.logo)
Sys.sleep(2) # Buzz roll

## Argument parsing
if (opt$parallel) {
  library(parallelly)
  opt$ncores <- availableCores(omit = 1)
} else {opt$ncores <- 1}
if (is.na(opt$shortname)) {
  opt$shortname <- str_trunc(slugify(opt$name, space_replace = "", tolower = F), 10, ellipsis = "")
} else {
  opt$shortname <- slugify(opt$shortname, space_replace = "", tolower = F)
}
opt$qc5_ref <- tolower(opt$qc5_ref)
opt$gh_ref <- tolower(opt$gh_ref)

## Argument overrides
if (isTRUE(opt$lo) && opt$lo_in == opt$lo_out) {
  opt$lo <- F # Abort LiftOver if both genomic builds identical
  message("LiftOver requested but input/output genomic builds identical.")
  message("LiftOver will not be carried out.")
}
if (opt$lo_in == 0 && isTRUE(opt$lo) && opt$chip == F) {
  opt$chip <- T # Chipendium needed to carry out LiftOver if input build not provided
  message("Genomic build of input data not provided but LiftOver requested.")
  message("Chipendium will be executed to detect the genomic build.")
}
if (opt$lo_in == 0 && isTRUE(opt$rename) && opt$chip == F) {
  opt$chip <- T # Chipendium needed to carry out LiftOver if input build not provided
  message("Genomic build of input data not provided but dbSNP rename requested.")
  message("Chipendium will be executed to detect the genomic build.")
}
if (opt$lo_in %in% c(36, 38) && isTRUE(opt$qc5)) {
  opt$lo <- T
  opt$lo_out <- 37 # LiftOver needed to carry out ancestry inference
  message("Input data is not GRCh37 but biogeographical ancestry inference requested.")
  message("LiftOver will be executed to map SNP data to GRCh37 for compatibility with ancestry reference panels.")
}
if (opt$lo_in %in% c(36, 38) && isTRUE(opt$gh) && opt$gh_ref == "hrc") {
  opt$lo <- T
  opt$lo_out <- 37 # LiftOver needed to carry out harmonisation
  message("Input data is not GRCh37 but harmonisation to HRC reference requested.")
  message("LiftOver will be executed to map SNP data to GRCh37 for compatibility with HRC reference panel.")
}
if (opt$lo_in %in% c(36, 37) && isTRUE(opt$gh) && opt$gh_ref == "topmed") {
  opt$lo <- T
  opt$lo_out <- 38 # LiftOver needed to carry out harmonisation
  message("Input data is not GRCh38 but harmonisation to TopMED reference requested.")
  message("LiftOver will be executed to map SNP data to GRCh38 for compatibility with TopMED reference panel.")
}

# Paths

## Software (OS-dependent)
if (get_os() == "windows"){
  PLINK1 <- here(opt$extrafolder, "plink.exe") # PLINK v1
  PLINK2 <- here(opt$extrafolder, "plink2.exe") # PLINK v2
  GH <- here(opt$extrafolder, "GH1.4.20", "GenotypeHarmonizer.bat") # Genotype Harmonizer v1.4.20
} else if (get_os() == "osx"){
  PLINK1 <- here(opt$extrafolder, "plink.osx") # PLINK v1 (renamed to avoid conflict with Unix version)
  PLINK2 <- here(opt$extrafolder, "plink2.osx") # PLINK v2 (renamed to avoid conflict with Unix version)
  GH <- here(opt$extrafolder, "GH1.4.20", "GenotypeHarmonizer.sh") # Genotype Harmonizer v1.4.20
} else {
  PLINK1 <- here(opt$extrafolder, "plink") # PLINK v1
  PLINK2 <- here(opt$extrafolder, "plink2") # PLINK v2
  GH <- here(opt$extrafolder, "GH1.4.20", "GenotypeHarmonizer.sh") # Genotype Harmonizer v1.4.20
}

## Reference and auxiliary files
CHIPENDIUM <- here(opt$extrafolder, "chipendium2019") # Chipendium SNP coordinates
LOCHAINS <- here(opt$extrafolder, "liftover-chains") # Liftover .chain files
HRCREF <- here(opt$extrafolder, "imputation-sites", "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz") # HRC r1.1 sites file (renamed from .tab.gz)
TOPMEDREF <- here(opt$extrafolder, "imputation-sites", "ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz") # TopMED freeze 5 sites file
GRC36REF <- here(opt$extrafolder, "dbSNP", "hapmap_r23a.sites.bim") # HapMap Phase 2 r23 (Gibbs 2005; PMID:26432245)
GRC37REF <- here(opt$extrafolder, "dbSNP", "all_phase3_ns.snp.sites.bim") # 1000 Phase 3 (Auton 2015; PMID:26432245)
GRC38REF <- here(opt$extrafolder, "dbSNP", "all_hg38.snp.sites.bim") # 1000 Phase 3 (Byrska-Bishop 2022; PMID:36055201)
lrld.plink <- here(opt$extrafolder, "lrld4plink.txt") # Coordinates of Long-Range LD regions (PLINK BED format)
qc5.filestem <- if(opt$qc5_ref == "aadr") {"aadr_v50_present"} else if(opt$qc5_ref == "1kgp") {"kgp3.array_snps"} # Possible ancestry references
qc5.refpath <- here(opt$extrafolder, "ancestry_reference", opt$qc5_ref, qc5.filestem) # Chosen reference for ancestry inference

## Outputs
outpath <- paste0("./", gsub("\\_+", "_", slugify(opt$name)))
if (!dir.exists(outpath)) {
  dir.create(outpath)
}
if (!dir.exists(paste0(outpath, "/genotypes"))) {
  dir.create(paste0(outpath, "/genotypes"))
}
if (!dir.exists(paste0(outpath, "/original"))) {
  dir.create(paste0(outpath, "/original"))
}
if (opt$qc5 == T && !dir.exists(paste0(outpath, "/ancestry"))) {
  dir.create(paste0(outpath, "/ancestry"))
}
if (opt$lo == T && !dir.exists(paste0(outpath, "/liftOver"))) {
  dir.create(paste0(outpath, "/liftOver"))
}
if (opt$gh == T && !dir.exists(paste0(outpath, "/imputeMe"))) {
  dir.create(paste0(outpath, "/imputeMe"))
}

# Checks

## Options
tar_assert_le(nchar(opt$shortname), 10, "Short dataset name cannot have more than 10 characters.")
if (isTRUE(opt$lo)) {
  tar_assert_in(opt$lo_in, c(36, 37, 38, 0), "Input genomic build unknown, only GRC assemblies 36, 37 and 38 currently supported.")
  tar_assert_in(opt$lo_out, c(36, 37, 38), "Output genomic build unknown, only GRC assemblies 36, 37 and 38 currently supported.")
  if(opt$lo_out == 36){tar_assert_file(GRC36REF)}
  if(opt$lo_out == 37){tar_assert_file(GRC37REF)}
  if(opt$lo_out == 38){tar_assert_file(GRC38REF)}
}

## Software
tar_assert_file(PLINK1)
tar_assert_file(PLINK2)
tar_assert_file(GH)

## Files
tar_assert_file(CHIPENDIUM)
tar_assert_file(LOCHAINS)
tar_assert_file(lrld.plink)
tar_assert_file(paste0(qc5.refpath, ".pgen"))
tar_assert_file(paste0(opt$file, ".bed"))
if (!exists("tinydata", where = opt)) {tar_assert_ge(file.size(opt$file, ".bed")/1024^2, 1, "Dataset files seem too small (< 1Mb), please check contents of .bed file. You can enable the flag --tinydata if you wish to proceed anyway.")}
tar_assert_file(paste0(opt$file, ".bim"))
opt$n_markers <- length(readLines(paste0(opt$file, ".bim")))
tar_assert_ge(opt$nmarkers, 1, "No markers found in dataset, please check contents of .bim file.")
tar_assert_file(paste0(opt$file, ".fam"))
opt$n_samples <- length(readLines(paste0(opt$file, ".fam")))
tar_assert_ge(opt$n_samples, 1, "No individuals found in dataset, please check contents of .fam file.")

# Future developments
## Not working: Implement PED/MAP compatibility for legacy data
## isPED=file.exists(paste0(opt$file, ".ped"))
## FIXME: tar_skip(pfile.input, skip=isPED==F, format="file", setNames(c(paste0(opt$file, ".ped"), paste0(opt$file, ".map")), c("ped", "map"))),
## FIXME: tar_skip(pfile.conversion, skip=isPED==F, system2(PLINK2, args=paste0("--ped ", pfile.input["ped"], "--map ", pfile.input["map"], "--make-bed --out ", opt$file))),

# Pipeline actions
Sys.sleep(2) # Buzz roll
message("Executing pipeline actions...")
list(

  # Read input files and load SNP/sample data
  tar_file(bfile.input, c(paste0(opt$file, ".bed"), paste0(opt$file, ".bim"), paste0(opt$file, ".fam"))),
  tar_target(bfile, list(bim=read_table(bfile.input[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                         fam=read_table(bfile.input[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii", na=c("na", "N/A")))),
  tar_target(bfile.nonauto, list(x=subset(bfile$bim, CHR %in% c(23, "X")),
                                 y=subset(bfile$bim, CHR %in% c(24, "Y")),
                                 par=subset(bfile$bim, CHR %in% c(25, "XY", "PAR1", "PAR2")),
                                 mito=subset(bfile$bim, CHR %in% c(26, "M", "MT")),
                                 nonchr=subset(bfile$bim, CHR %nin% c(1:26, "X", "Y", "XY", "PAR1", "PAR2", "M", "MT")))),
  
  # Extract mitochondrial DNA content
  tar_skip(bfile.mito.extract, skip=nrow(bfile.nonauto$mito)==0,
                               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", qc0.output[1], " --pvar ", qc0.output[2], " --psam ", qc0.output[3], " --chr 26 --make-bed --out ", outpath, "/genotypes/", opt$file, ".MTDNA"))),
  
  # Chipendium procedure
  tar_target(chip.manifests, if (opt$chip==F) {character(0)} else {grep("gz", list.files(CHIPENDIUM, full.names = T), value = T)}),
  tar_target(chip.matches, if (opt$chip==F) {character(0)
                         } else if (opt$ncores > 1) {message("Finding most likely genotyping array using ", length(chip.manifests), " manifest files...")
													 chip.matches.parallel(chip.manifests, chipendium_match, bimfile=bfile$bim)
                         } else {message("Finding most likely genotyping array using ", length(chip.manifests), " manifest files...")
								 pblapply(chip.manifests, chipendium_match, bimfile=bfile$bim)}),
  tar_target(chip.overlap, if (opt$chip==F) {data.frame(build=opt$lo_in)} else{chip_find_best(chip.matches)}),
  
  # QC0: Convert to PGEN to ensure proper formatting
  ## Define a reference for marker renaming (if needed)
  tar_skip(rename.reference, skip=opt$rename==F, format="file",
           if(opt$lo_in==37){GRC37REF
           }else if(opt$lo_in==38){GRC38REF
           }else if(opt$lo_in==36){GRC36REF
             # Account for missing reference
           }else if(opt$lo_in==0){if(chip.overlap$build==37){GRC37REF
           }else if(chip.overlap$build==38){GRC38REF
           }else if(chip.overlap$build==36){GRC36REF
           }}),
  ## This steps merges PAR regions (will be split later) and sorts variants according to chr:pos
  tar_target(qc0.makefile, if (isTRUE(opt$rename)) {
    system2(PLINK1, args=paste0(" --threads ", opt$ncores, " --bed ", bfile.input[1], " --bim ", bfile.input[2], " --fam ", bfile.input[3], " --make-bed --out ", outpath, "/", opt$file, ".new_format.tmp"))
    system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bfile ", outpath, "/", opt$file, ".new_format.tmp --sort-vars --make-pgen --output-chr 26 --out ", outpath, "/", opt$file, ".old_names.tmp"))
    system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pfile ", outpath, "/", opt$file, ".old_names.tmp --no-categorical --set-all-var-ids ", r"(@:#:$r:$a)", " --new-id-max-allele-len 23 missing --chr 1-26 --make-bed --output-chr 26 --out ", outpath, "/", opt$file, ".chrpos"))
    invisible(file.remove(list.files(path=paste0(outpath, "/"), pattern="*.new_format.tmp.*", full.names = T)))
    invisible(file.remove(list.files(path=paste0(outpath, "/"), pattern="*.old_names.tmp.*", full.names = T)))
    #system2(GH, args=paste0(" --update-id --keep -I PLINK_BED -R VCF --input ", outpath, "/", opt$file, ".chrpos --ref ", shQuote(rename.reference), " --output ", outpath, "/", opt$file, ".rename"))
    plink_rename_snp(bfile=paste0(outpath, "/", opt$file, ".chrpos"),
                     reffile=shQuote(rename.reference), outdir=paste0(outpath, "/"),
                     ncores=opt$ncores, PLINK2=PLINK2, PLINK1=PLINK1,
                     outfile=paste0(outpath, "/", opt$file, ".rename"))
    system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bfile ", outpath, "/", opt$file, ".rename --merge-par --make-pgen --out ", outpath, "/", opt$file, ".qc0"))
    } else {
    system2(PLINK1, args=paste0(" --threads ", opt$ncores, " --bed ", bfile.input[1], " --bim ", bfile.input[2], " --fam ", bfile.input[3], " --merge-x no-fail --make-bed --out ", outpath, "/", opt$file, ".new_format.tmp"))
    system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bfile ", outpath, "/", opt$file, ".new_format.tmp --sort-vars --make-pgen --output-chr 26 --out ", outpath, "/", opt$file, ".old_names.tmp"))
    system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pfile ", outpath, "/", opt$file, ".old_names.tmp --no-categorical --chr 1-26 --make-pgen --out ", outpath, "/", opt$file, ".qc0"))
    invisible(file.remove(list.files(path=paste0(outpath, "/"), pattern="*.new_format.tmp.*", full.names = T)))
    invisible(file.remove(list.files(path=paste0(outpath, "/"), pattern="*.old_names.tmp.*", full.names = T)))
    }),
  tar_file(qc0.output,{qc0.makefile
    c(paste0(outpath, "/", opt$file, ".qc0.pgen"), paste0(outpath, "/", opt$file, ".qc0.pvar"), paste0(outpath, "/", opt$file, ".qc0.psam"))}),
  tar_target(pvar.qc0, if(isTRUE(opt$rename)) {read_table(qc0.output[2], col_names=c("CHR", "BP", "SNP", "A2", "A1"), col_types="ciccc", skip=1)
                       } else {data.frame(SNP=character(0))}),

  # QC1: Missingness checks
  tar_target(qc1, list(plinkrun.sample=system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", qc0.output[1], " --pvar ", qc0.output[2], " --psam ", qc0.output[3], " --missing sample-only scols=-missphenos --out ", outpath, "/", opt$file, ".qc1a")),
                       plinkrun.ind=system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", qc0.output[1], " --pvar ", qc0.output[2], " --psam ", qc0.output[3], " --missing variant-only --mind ",1-opt$qc1_coverage, " --out ", outpath, "/", opt$file, ".qc1b")),
                       snp=read_table(paste0(outpath, "/", opt$file, ".qc1b.vmiss"), col_types="_c__d"),
                       ind=read_table(paste0(outpath, "/", opt$file, ".qc1a.smiss"), col_types="cc__d"))),
  tar_target(qc1.problem, list(snp=subset(qc1$snp, F_MISS > 1-opt$qc1_call_rate),
                               ind=subset(qc1$ind, F_MISS > 1-opt$qc1_coverage))),
  tar_target(qc1.makefile, {qc1.problem
                            if_else(opt$qc1==F,
                                   system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", qc0.output[1], " --pvar ", qc0.output[2], " --psam ", qc0.output[3], " --split-par b", chip.overlap$build, " --output-chr MT --min-alleles 1 --make-bed --out ", outpath, "/", opt$file, ".qc1")),
                                   system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", qc0.output[1], " --pvar ", qc0.output[2], " --psam ", qc0.output[3], " --split-par b", chip.overlap$build, " --geno ",1-opt$qc1_call_rate, " --mind ",1-opt$qc1_coverage, " --output-chr MT --make-bed --out ", outpath, "/", opt$file, ".qc1")))}),
  tar_file(qc1.output,{qc1.makefile
                       c(paste0(outpath, "/", opt$file, ".qc1.bed"), paste0(outpath, "/", opt$file, ".qc1.bim"), paste0(outpath, "/", opt$file, ".qc1.fam"))}),
  tar_target(bfile.qc1, list(bim=read_table(qc1.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                            fam=read_table(qc1.output[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii"))),
  
  # QC2: Hardy-Weinberg Equilibrium
  tar_target(qc2, list(plinkrun=system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc1.output[1], " --bim ", qc1.output[2], " --fam ", qc1.output[3], " --hardy midp --out ", outpath, "/", opt$file, ".qc2")),
                  res=if (nrow(bfile.nonauto$x)==0) {read_table(paste0(outpath, "/", opt$file, ".qc2.hardy"), col_types="_c__iii__d", col_names=c("SNP", "A1A1", "A1A2", "A2A2", "MIDP"), skip=1)
                  } else {bind_rows(read_table(paste0(outpath, "/", opt$file, ".qc2.hardy"), col_types="_c__iii__d", col_names=c("SNP", "A1A1", "A1A2", "A2A2", "MIDP"), skip=1),
                                    read_table(paste0(outpath, "/", opt$file, ".qc2.hardy.x"), col_types="_c__iii______d", col_names=c("SNP", "A1A1", "A1A2", "A2A2", "MIDP"), skip=1))})),
  tar_target(qc2.makefile, if_else(opt$qc2==F, system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc1.output[1], " --bim ", qc1.output[2], " --fam ", qc1.output[3], " --make-bed --out ", outpath, "/", opt$file, ".qc2")),
                                              ifelse(opt$qc2_popmix==F, system2(PLINK2, args=paste0(" --bed ", qc1.output[1], " --bim ", qc1.output[2], " --fam ", qc1.output[3], " --hwe ", opt$qc2_hwepval, " midp --make-bed --out ", outpath, "/", opt$file, ".qc2")),
                                                                        system2(PLINK2, args=paste0(" --bed ", qc1.output[1], " --bim ", qc1.output[2], " --fam ", qc1.output[3], " --hwe ", opt$qc2_hwepval, " midp keep-fewhet --make-bed --out ", outpath, "/", opt$file, ".qc2"))))),
  tar_file(qc2.output,{qc2.makefile
                       c(paste0(outpath, "/", opt$file, ".qc2.bed"), paste0(outpath, "/", opt$file, ".qc2.bim"), paste0(outpath, "/", opt$file, ".qc2.fam"))}),
  tar_target(bfile.qc2, list(bim=read_table(qc2.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                             fam=read_table(qc2.output[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii"))),
  tar_target(qc2.problem, bfile.qc1$bim$SNP[bfile.qc1$bim$SNP %nin% bfile.qc2$bim$SNP]),
  
  # QC2b: Remove SNPs with duplicated names for further QC steps
  tar_target(qc2b.makefile, system2(PLINK2, args=paste0(" --bed ", liftOver.output[1], " --bim ", liftOver.output[2], " --fam ", liftOver.output[3], " --rm-dup exclude-all --make-bed --out ", outpath, "/", opt$file, ".infile.tmp"))),
  tar_file(qc2b.output,{qc2b.makefile
                       c(paste0(outpath, "/", opt$file, ".infile.tmp.bed"), paste0(outpath, "/", opt$file, ".infile.tmp.bim"), paste0(outpath, "/", opt$file, ".infile.tmp.fam"))}),
  tar_target(bfile.qc2b, list(bim=read_table(qc2b.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                              fam=read_table(qc2b.output[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii"))),

  # QC3: Sex chromosome checks
  tar_file(bfile.lrld, lrld.plink),
  tar_target(bfile.makeindep, system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc2b.output[1], " --bim ", qc2b.output[2], " --fam ", qc2b.output[3], " --not-chr y mt --exclude bed1 ", shQuote(bfile.lrld), " --indep-pairwise 500kb 0.2 --maf ", opt$qc4_maf, " --out ", outpath, "/", opt$file, ".r02"))),
  tar_file(bfile.indep, {bfile.makeindep
                         c(paste0(outpath, "/", opt$file, ".r02.prune.in"), paste0(outpath, "/", opt$file, ".r02.prune.out"))}),
  tar_target(qc3.makefile, if_else(opt$qc3==F | nrow(bfile.nonauto$x)==0,
                                   qc2.makefile,
                                   system2(PLINK2, args=paste0(" --bed ", qc2.output[1], " --bim ", qc2.output[2], " --fam ", qc2.output[3], " --chr x y --min-alleles 2 --make-bed --out ", outpath, "/", opt$file, ".qc3.xy")))),
  tar_file(qc3.output,{qc3.makefile
                       if(opt$qc3==F | nrow(bfile.nonauto$x)==0){qc2.output
                       } else {c(paste0(outpath, "/", opt$file, ".qc3.xy.bed"), paste0(outpath, "/", opt$file, ".qc3.xy.bim"), paste0(outpath, "/", opt$file, ".qc3.xy.fam"))}}),
  tar_target(bfile.qc3, list(bim.x=read_table(qc3.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc")%>%subset(CHR %in% c(23, "X")),
                            bim.y=read_table(qc3.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc")%>%subset(CHR %in% c(24, "Y")))),
  tar_target(qc3, if (opt$qc3==F | nrow(bfile.nonauto$x)==0) {numeric(0)
                                                } else if (nrow(bfile.nonauto$y)==0) {list(plinkrun=system2(PLINK1, args=paste0("--allow-extra-chr --threads ", opt$ncores, " --bed ", qc3.output[1], " --bim ", qc3.output[2], " --fam ", qc3.output[3], " --exclude ", bfile.indep[2], " --check-sex ", opt$qc3_max_female_F, " ", opt$qc3_min_male_F, " --out ", outpath, "/", opt$file, ".qc3")),
                                                                                           ind=read_table(paste0(outpath, "/", opt$file, ".qc3.sexcheck"), col_types="cciicd"))
                                                } else {list(plinkrun=system2(PLINK1, args=paste0("--allow-extra-chr --threads ", opt$ncores, " --bed ", qc3.output[1], " --bim ", qc3.output[2], " --fam ", qc3.output[3], " --exclude ", bfile.indep[2], " --check-sex ycount ", opt$qc3_max_female_F, " ", opt$qc3_min_male_F, " --out ", outpath, "/", opt$file, ".qc3")),
                                                             ind=read_table(paste0(outpath, "/", opt$file, ".qc3.sexcheck"), col_types="cciicdi"))}),
  tar_skip(qc3.problem, skip=(opt$qc3==F | nrow(bfile.nonauto$x)==0), subset(qc3$ind, F>=opt$qc3_max_female_F & F<=opt$qc3_min_male_F)),
  
  # QC4: Population structure
  tar_target(bfile.makegds, packages = "SNPRelate",
             snpgdsBED2GDS(bed.fn=qc2b.output[1], bim.fn=qc2b.output[2], fam.fn=qc2b.output[3], out.gdsfn=paste0(outpath, "/", opt$file, ".qc2.gds"), option=snpgdsOption(X=23, Y=24, PAR1=25, PAR2=25, MT=26))),
  tar_file(bfile.gds, {bfile.makegds
                       c(paste0(outpath, "/", opt$file, ".qc2.gds"))}),
  tar_target(bfile.king, packages = "SNPRelate", king_from_gds(bfile.gds, opt$ncores)),
  tar_target(qc4.pcair, packages = c("GWASTools", "GENESIS"), pcair_from_gds(bfile.gds, bfile.king, bfile.indep[1], opt$ncores, opt$qc4_pcair)),
  
  # QC4: Relatedness
  tar_target(qc4.pcrelate, packages = c("GWASTools", "GENESIS", "BiocParallel"), pcrelate_from_gds(bfile.gds, qc4.pcair, bfile.indep[1], opt$qc4_pcrelate, opt$ncores)),
  
  # QC5: Biogeographical Ancestry
  ## QC5-a: Find AIMs
  tar_file(qc5.ref.fst, Sys.glob(paste0(qc5.refpath, ".*.fst.var.gz"))),
  tar_file(qc5.ref, c(paste0(qc5.refpath, ".pgen"), paste0(qc5.refpath, ".pvar"), paste0(qc5.refpath, ".psam"), paste0(qc5.refpath, ".balanced.ids"), paste0(qc5.refpath, ".norel.ids"), paste0(qc5.refpath, ".norel.pheno"))),
  tar_target(qc5.ref.bfile, list(bim=read_table(qc5.ref[2], col_names=c("CHR", "BP", "SNP", "A1", "A2", "cm"), col_types="iicccd", skip=1),
                                 fam=read_table(qc5.ref[3], col_names=c("FID", "IID", "Sex", "Pheno"), col_types="ccii", skip=1),
                                 pheno=read_table(qc5.ref[6], col_names=T, col_types="ccc"))),
  tar_target(bfile.aims, packages = "VGAM",
             if (isTRUE(opt$lo) & opt$lo_in %in% c(36,38) | (isTRUE(opt$lo) & opt$lo_in==0 & chip.overlap$build!=37)) {
			 message("LiftOver to GRCh37 requested, using lifted files for ancestry inference..")
             find_aims(bfile.liftOver$bim, qc5.ref.bfile$bim, bfile.lrld, qc5.ref.fst, opt$ncores) 
             } else {
             find_aims(bfile.qc2b$bim, qc5.ref.bfile$bim, bfile.lrld, qc5.ref.fst, opt$ncores)}
			 ),
  tar_target(bfile.makeaims, write.table(bfile.aims, file=paste0(outpath, "/ancestry/", opt$file, ".aims"), col.names = T, row.names = F, quote = F)),
  tar_file(qc5.aims,{bfile.makeaims
                     paste0(outpath, "/ancestry/", opt$file, ".aims")}),
  tar_target(qc5.indep, system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", shQuote(qc5.ref[1]), " --pvar ", shQuote(qc5.ref[2]), " --psam ", shQuote(qc5.ref[3]), " --keep ", shQuote(qc5.ref[4]), " --clump cols=-total,-bins ", qc5.aims, " --clump-snp-field ID --clump-field HUDSON_FST_P --clump-kb 500 --clump-r2 0.2 --clump-p1 1 --clump-p2 1 --out ", outpath, "/ancestry/", qc5.filestem, ".AIM.r02"))),
  tar_file(qc5.aims.indep,{qc5.indep
                               paste0(outpath, "/ancestry/", qc5.filestem, ".AIM.r02.clumps")}),
  ## QC5-b: Extract AIMs from reference and input sample 
  tar_target(qc5.ref.extract, system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pgen ", shQuote(qc5.ref[1]), " --pvar ", shQuote(qc5.ref[2]), " --psam ", shQuote(qc5.ref[3]), " --keep ", shQuote(qc5.ref[5]), " --extract ", qc5.aims.indep, " --set-missing-var-ids ", r"(@:#:$r:$a)", " --make-bed --out ", outpath, "/ancestry/", qc5.filestem, ".AIM.r02.extract"))),
  tar_file(qc5.ref.aims,{qc5.ref.extract
                           c(paste0(outpath, "/ancestry/", qc5.filestem, ".AIM.r02.extract.bed"), paste0(outpath, "/ancestry/", qc5.filestem, ".AIM.r02.extract.bim"), paste0(outpath, "/ancestry/", qc5.filestem, ".AIM.r02.extract.fam"))}),
  tar_target(qc5.input.extract,
             if (isTRUE(opt$lo) & opt$lo_in %in% c(36,38) | (isTRUE(opt$lo) & opt$lo_in==0 & chip.overlap$build!=37)) {
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", liftOver.output[1], " --bim ", liftOver.output[2], " --fam ", liftOver.output[3], " --extract ", qc5.aims.indep, " --set-missing-var-ids ", r"(@:#:$r:$a)", " --make-bed --out ", outpath, "/ancestry/", opt$file, ".AIM.r02.extract"))
             } else {
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc2b.output[1], " --bim ", qc2b.output[2], " --fam ", qc2b.output[3], " --extract ", qc5.aims.indep, " --set-missing-var-ids ", r"(@:#:$r:$a)", " --make-bed --out ", outpath, "/ancestry/", opt$file, ".AIM.r02.extract"))}),
  tar_file(qc5.input.aims,{qc5.input.extract
    c(paste0(outpath, "/ancestry/", opt$file, ".AIM.r02.extract.bed"), paste0(outpath, "/ancestry/", opt$file, ".AIM.r02.extract.bim"), paste0(outpath, "/ancestry/", opt$file, ".AIM.r02.extract.fam"))}),
  ## QC5-c: Merge reference and input sample 
  tar_target(qc5.ref.merge, system2(PLINK1, args=paste0("--allow-extra-chr --allow-no-sex --threads ", opt$ncores, " --bed ", qc5.input.aims[1], " --bim ", qc5.input.aims[2], " --fam ", qc5.input.aims[3], " --bmerge ", qc5.ref.aims[1], " ", qc5.ref.aims[2], " ", qc5.ref.aims[3], " --indiv-sort 0 --make-bed --out ", outpath, "/ancestry/", opt$file, ".", opt$qc5_ref))),
  tar_file(qc5.input.ref,{qc5.ref.merge
                          c(paste0(outpath, "/ancestry/", opt$file, ".", opt$qc5_ref, ".bed"), paste0(outpath, "/ancestry/", opt$file, ".", opt$qc5_ref, ".bim"), paste0(outpath, "/ancestry/", opt$file, ".", opt$qc5_ref, ".fam"))}),
  ## QC5-d: Load into R and fit LDA model
  tar_target(qc5.input.ref.makegds, packages = "SNPRelate",
             snpgdsBED2GDS(bed.fn=qc5.input.ref[1], bim.fn=qc5.input.ref[2], fam.fn=qc5.input.ref[3], out.gdsfn=paste0(outpath, "/ancestry/", opt$file, ".", opt$qc5_ref, ".gds"))),
  tar_file(qc5.input.ref.gds, {qc5.input.ref.makegds
                               c(paste0(outpath, "/ancestry/", opt$file, ".", opt$qc5_ref, ".gds"))}),
  tar_target(qc5.pcair, packages = c("GWASTools", "GENESIS"), pcair_project_from_gds(qc5.input.ref.gds, bfile.qc2$fam, opt$ncores)),
  tar_target(qc5.tw, packages = c("GWASTools", "GENESIS", "AssocTests"), tw(qc5.pcair$values,100, criticalpoint = 0.9793)$SigntEigenL),
  tar_target(qc5.train, packages = "caret", lda_train(qc5.pcair, qc5.tw, qc5.ref.bfile)),
  tar_target(qc5.youden, packages = "probably", lda_best_threshold(qc5.train)),
  tar_target(qc5.test, packages = "caret", lda_test(bfile.qc2$fam, qc5.pcair, qc5.tw, qc5.train, qc5.youden)),
  
  # Liftover procedure
  tar_skip(lo.reference, skip=opt$lo==F, format="file",
           if(opt$lo_out==37){GRC37REF
           }else if(opt$lo_out==38){GRC38REF
           }else if(opt$lo_out==36){GRC36REF
           }else{(character(0))}),
  tar_target(liftOver.snp, packages = "rtracklayer",
             if (opt$lo==F) {character(0)
             } else if (opt$lo_in==0) {
               liftOver_bim(bfile.qc2$bim, build.in=chip.overlap$build, build.out=opt$lo_out, chainpath=LOCHAINS)
             } else {
               liftOver_bim(bfile.qc2$bim, build.in=opt$lo_in, build.out=opt$lo_out, chainpath=LOCHAINS)
             }),
  tar_skip(liftOver.map, skip=opt$lo==F, format="file",{
    write.table(liftOver.snp, file=paste0(outpath, "/liftOver/", opt$file, ".b", opt$lo_out, ".map"), col.names = F, row.names = F, quote = F)
    paste0(outpath, "/liftOver/", opt$file, ".b", opt$lo_out, ".map")}),
  tar_target(liftOver.output, format="file",
             if(opt$lo==F){
               qc2.output
             } else {
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc2.output[1], " --bim ", qc2.output[2], " --fam ", qc2.output[3], " --update-map ", liftOver.map, " --extract ", liftOver.map, " --make-pgen --sort-vars --out ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".old_names.tmp"))
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --pfile ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".old_names.tmp --set-all-var-ids ", r"(@:#:$r:$a)", " --new-id-max-allele-len 23 missing --make-bed --out ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".chrpos.tmp"))
               plink_rename_snp(bfile=paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".chrpos.tmp"),
                                reffile=shQuote(lo.reference), outdir=paste0(outpath, "/liftOver/"),
                                ncores=opt$ncores, PLINK2=PLINK2, PLINK1=PLINK1,
                                outfile=paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out)) 
               #system2(GH, args=paste0(" --update-id --keep -I PLINK_BED -R VCF --input ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".chrpos.tmp --ref ", shQuote(lo.reference), " --output ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out))
               invisible(file.remove(list.files(path=paste0(outpath, "/liftOver/"), pattern="*.old_names.tmp.*", full.names = T)))
               invisible(file.remove(list.files(path=paste0(outpath, "/liftOver/"), pattern="*.chrpos.tmp.*", full.names = T)))
               c(paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".bed"),
                 paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".bim"),
                 paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".fam"))}),
  tar_target(bfile.liftOver,
             if(opt$lo==F){
               bfile.qc2
             } else {	
               list(bim=read_table(liftOver.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                    fam=read_table(liftOver.output[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii"))}),
  
  # Genotype Harmonisation
  tar_skip(gh.input, skip=opt$gh==F,
           if(opt$lo_in %in% c(36,38) & opt$gh_ref=="hrc") {
             paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out)
           } else if(opt$lo_in %in% c(36,37) & opt$gh_ref=="topmed") {
             system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", liftOver.output[1], " --bim ", liftOver.output[2], " --fam ", liftOver.output[3], " --output-chr chrMT --chr 1-23 --make-bed --out ", outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".prefix"))
             paste0(outpath, "/liftOver/", opt$file, ".qc2.b", opt$lo_out, ".prefix")
           } else if(opt$gh_ref=="topmed"){
             system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc2.output[1], " --bim ", qc2.output[2], " --fam ", qc2.output[3], " --output-chr chrMT --chr 1-23 --make-bed --out ", outpath, "/", opt$file, ".qc2.prefix"))
             paste0(outpath, "/", opt$file, ".qc2.prefix")
           } else {paste0(outpath, "/", opt$file, ".qc2")}),
  tar_target(gh.input.hh,
             if(opt$gh==F | nrow(bfile.nonauto$x)==0){gh.input
             } else if(opt$gh_ref=="topmed"){system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bfile ", gh.input, " --set-invalid-haploid-missing --output-chr chrMT --chr 1-23 --make-bed --out ", outpath, "/imputeMe/", opt$file, ".nohh"))
             paste0(outpath, "/imputeMe/", opt$file, ".nohh")
             } else {system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bfile ", gh.input, " --set-invalid-haploid-missing --make-bed --out ", outpath, "/imputeMe/", opt$file, ".nohh"))
             paste0(outpath, "/imputeMe/", opt$file, ".nohh")
             }),
  tar_skip(gh.reference, skip=opt$gh==F, format="file",
           if(opt$gh_ref=="hrc"){HRCREF
           }else if(opt$gh_ref=="topmed"){TOPMEDREF
           }else{(character(0))}),
  tar_target(gh.output, format="file",
           if(opt$gh==F){qc2.output
           } else {system2(GH, args=paste0(" -asf --update-id -I PLINK_BED -R VCF --input ", gh.input.hh, " --ref ", shQuote(gh.reference), " --output ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref))
           if (nrow(bfile.nonauto$x)>=1) {invisible(file.remove(list.files(path=paste0(outpath, "/imputeMe/"), pattern="*.nohh.*", full.names = T)))}
           c(paste0(outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".bed"),
             paste0(outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".bim"),
             paste0(outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".fam"))}),
  tar_target(bfile.gh,
             if(opt$gh==F){
               bfile.qc2
             } else {
               list(bim=read_table(gh.output[2], col_names=c("CHR", "SNP", "cm", "BP", "A1", "A2"), col_types="ccdicc"),
                    fam=read_table(gh.output[3], col_names=c("FID", "IID", "PatID", "MatID", "Sex", "Pheno"), col_types="ccccii"))}),
  tar_skip(gh.split.format, skip=opt$gh==F,
           if(length(unique(bfile.gh$fam$IID))==nrow(bfile.gh$fam)) {
             message("Genotype harmonisation to imputation reference panel completed.")
             message("VCF files (v4.2 format) for Imputation Server will have IID as sample identifier.")
             "iid"
           } else if (length(grep("_", c(bfile.gh$fam$FID, bfile.gh$fam$IID)))==0) {
             message("Genotype harmonisation to imputation reference panel completed.")
             message("VCF files (v4.2 format) for Imputation Server will have FID_IID as sample identifier.")
             "_"
           } else if (length(grep("-", c(bfile.gh$fam$FID, bfile.gh$fam$IID)))==0) {
             message("Genotype harmonisation to imputation reference panel completed.")
             message("VCF files (v4.2 format) for Imputation Server will have FID-IID as sample identifier.")
             "-"
           } else if (length(grep(".", c(bfile.gh$fam$FID, bfile.gh$fam$IID)))==0) {
             message("Genotype harmonisation to imputation reference panel completed.")
             message("VCF files (v4.2 format) for Imputation Server will have FID.IID as sample identifier.")
             "."
           } else {
             message("Genotype harmonisation to imputation reference panel completed.")
             message("VCF files (v4.2 format) for Imputation Server will have FID_IID as sample identifier.")
             "standard"}),
  tar_skip(gh.split, skip=opt$gh==F,
           for (chr in 1:23) {
             if (opt$gh_ref=="topmed"){
               if (gh.split.format=="iid"){
                 system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr chrMT --export vcf-4.2 bgz id-paste=", gh.split.format, " --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                         stdout = NULL, stderr = NULL)
               } else if (gh.split.format%in%c("_", "-", ".")){
                 system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr chrMT --export vcf-4.2 bgz id-delim=", gh.split.format, " --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                         stdout = NULL, stderr = NULL)}
               else {
                 system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr chrMT --export vcf-4.2 bgz --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                         stdout = NULL, stderr = NULL)} 
             } else {
             if (gh.split.format=="iid"){
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr MT --export vcf-4.2 bgz id-paste=", gh.split.format, " --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                       stdout = NULL, stderr = NULL)
             } else if (gh.split.format%in%c("_", "-", ".")){
               system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr MT --export vcf-4.2 bgz id-delim=", gh.split.format, " --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                       stdout = NULL, stderr = NULL)}
               else {
              system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", gh.output[1], " --bim ", gh.output[2], " --fam ", gh.output[3], " --chr ", chr, " --ref-allele ", shQuote(gh.reference), " 4 3 '#' --output-chr MT --export vcf-4.2 bgz --out ", outpath, "/imputeMe/", opt$file, ".gh.", opt$gh_ref, ".chr", chr),
                      stdout = NULL, stderr = NULL)}                 
               }}),

  # Cleanup QC folder
  tar_file(qc6.final, {qc1.problem
                       qc2.problem
                       qc3.problem
                       qc4.pcrelate
                       qc5.test
                       bfile.liftOver
                       gh.split
                       invisible(file.remove(list.files(path=paste0(outpath, "/"), pattern="*.infile.tmp.*", full.names = T)))
                       system2(PLINK2, args=paste0(" --threads ", opt$ncores, " --bed ", qc2.output[1], " --bim ", qc2.output[2], " --fam ", qc2.output[3], " --make-pgen --out ", outpath, "/genotypes/", opt$file, ".qc"))
                       c(paste0(outpath, "/genotypes/", opt$file, ".qc.pgen"), paste0(outpath, "/genotypes/", opt$file, ".qc.pvar"), paste0(outpath, "/genotypes/", opt$file, ".qc.psam"))}),
  
  # RMarkdown report
  tar_render(report, "GenotypeQCtoHRC_targets.Rmd", params = list(softfolder = opt$extrafolder,
                                                                  scriptsfolder = opt$scriptsfolder,
                                                                  dataset = opt$name,
                                                                  dataset.short = opt$shortname,
                                                                  dataset.filename = opt$file,
                                                                  dataset.outpath = outpath,
                                                                  check = opt$check,
                                                                  chip = opt$chip,
                                                                  rename = opt$rename,
                                                                  lo = opt$lo,
                                                                  lo.in = chip.overlap$build,
                                                                  lo.out = opt$lo_out,
                                                                  gh = opt$gh,
                                                                  gh.ref =opt$gh_ref,
                                                                  qc1 = opt$qc1,
                                                                  qc1.call.rate = opt$qc1_call_rate,
                                                                  qc1.coverage = opt$qc1_coverage,
                                                                  qc2 = opt$qc2,
                                                                  qc2.hwepval = opt$qc2_hwepval,
                                                                  qc2.popmix = opt$qc2_popmix,
                                                                  qc3 = opt$qc3,
                                                                  qc3.max.femaleF = opt$qc3_max_female_F,
                                                                  qc3.min.maleF = opt$qc3_min_male_F,
                                                                  qc4 = opt$qc4,
                                                                  qc4.maf = opt$qc4_maf,
                                                                  qc4.pcairs = opt$qc4_pcair,
                                                                  qc4.pcrelate = opt$qc4_pcrelate,
                                                                  qc4.kin = opt$qc4_kin,
                                                                  qc5 = opt$qc5,
                                                                  qc5.ref = opt$qc5_ref),
                                                    output_file=paste0(slugify(opt$name), ".qc_report.html"))
)
