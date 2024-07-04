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

# OS detection
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
} # Function by Will Lowe (https://conjugateprior.org/2015/06/identifying-the-os-from-r/)

# Chipendium procedure
chipendium_match <- function(chip_strand_file,bimfile) {
  chippy = as.data.frame(read_table(chip_strand_file,col_names=F,col_types = "c",progress = F,show_col_types = F))
  bimfile$chrpos = paste0(bimfile$CHR,":",bimfile$BP)
  overlap = bimfile[bimfile$chrpos %in% unlist(chippy$X1), ]$chrpos
  chip_overlap =as.data.frame(cbind(length(overlap), as.character(chip_strand_file)))
  colnames(chip_overlap) = c("OverlappingSNPs", "CHIP")
  chip_overlap$OverlappingSNPs = as.numeric(as.character(chip_overlap$OverlappingSNPs))
  chip_overlap = chip_overlap[order(chip_overlap$OverlappingSNPs), ]
  return(chip_overlap)
}

chip.matches.parallel <- function(chip.manifests, chipendium_match, bimfile){
  hatchlings <- parallelly::makeClusterPSOCK(opt$ncores, autoStop = TRUE, rscript_startup = quote(local({
    options(tidyverse.quiet = TRUE)
    suppressWarnings(library("tidyverse"))})))
  pblapply(chip.manifests, chipendium_match, bimfile, cl=hatchlings)
}

chip_find_best <- function(chip.matches) {
  # Compare number of overlapping SNPs with 391 genotyping platforms and builds
  chip.overlap = as.data.frame(do.call(rbind, chip.matches))
  # Report the closest matching genotyping platform
  best.guess.chip = gsub("-strand.zip", "", chip.overlap[which(chip.overlap$OverlappingSNPs == max(chip.overlap$OverlappingSNPs)),2])
  best.guess.chip = basename(best.guess.chip)
  best.guess.chip = gsub(".processed.txt.gz", "", best.guess.chip)
  best.guess.chip.n = gsub("-strand.zip", "", chip.overlap[which(chip.overlap$OverlappingSNPs == max(chip.overlap$OverlappingSNPs)),1])
  best.guess.chip.n = unique(best.guess.chip.n)
  best.guess.chip.grc = as.integer(unique(str_extract(best.guess.chip,"3[6-8]")))
  message("Based on Chipendium results, most likely genomic build is GRCh",best.guess.chip.grc)
  return(list(chip=best.guess.chip,
              snps=best.guess.chip.n,
              build=best.guess.chip.grc))
}

# Rename procedure
plink_rename_snp <- function(bfile, reffile, outfile, outdir, ncores, PLINK2, PLINK1) {
  # First pass
  pfile <- bfile
  ## Match markers with reference
  system2(PLINK2, args=paste0(" --threads ",ncores," --bfile ",bfile," --recover-var-ids ",reffile," partial force --output-chr 26 --make-pgen --out ",pfile,".TMP.r01"))
  ## Fill out missing SNP IDs (missing in reference)
  system2(PLINK2, args=paste0(" --threads ",ncores," --pfile ",pfile,".TMP.r01 --set-missing-var-ids ",r"(@:#:$r:$a)"," --new-id-max-allele-len 23 missing --rm-dup exclude-all --make-pgen --sort-vars --output-chr 26 --out ",pfile,".TMP.r1"))
  ## Read PVAR file and create a file from all non-RS markers
  pvar <- read_table(paste0(pfile,".TMP.r1.pvar"), col_names=c("CHR","BP","SNP","A2","A1"),col_types="ciccc",skip=1)
  pvar <- mutate(pvar,chrpos=paste(CHR,BP,sep = ":"))
  pvar.match <- subset(pvar, grepl("rs",SNP), select=c("SNP","chrpos")) # Matched markers
  pvar.flip <- subset(pvar, !grepl("rs",SNP) & SNP!=".", select=c("SNP","chrpos")) # Unmatched markers
  write.table(pvar.flip$SNP,file=paste0(pfile,".TMP.r1.flip"), col.names = F, row.names = F, quote = F)
  # Second pass
  ## Extract unmatched markers
  system2(PLINK2, args=paste0(" --threads ",ncores," --pfile ",pfile,".TMP.r1 --extract ",pfile,".TMP.r1.flip --output-chr 26 --make-bed --out ",pfile,".TMP.r1b"))
  ## Flip in PLINKv1
  system2(PLINK1, args=paste0(" --threads ",ncores," --bfile ",pfile,".TMP.r1b --flip ",pfile,".TMP.r1.flip --make-bed --out ",pfile,".TMP.r1c"))
  ## Match flipped markers with reference
  system2(PLINK2, args=paste0(" --threads ",ncores," --bfile ",pfile,".TMP.r1c --recover-var-ids ",reffile," partial force --make-pgen --output-chr 26 --out ",pfile,".TMP.r02"))
  ## Fill out missing SNP IDs (missing in reference)
  system2(PLINK2, args=paste0(" --threads ",ncores," --pfile ",pfile,".TMP.r02 --set-missing-var-ids ",r"(@:#:$r:$a)"," --new-id-max-allele-len 23 missing --rm-dup exclude-all --make-bed --output-chr 26 --out ",pfile,".TMP.r2"))
  ## Read PVAR file and create a file from all RS matches that did not exist in the first pass
  bim <- read_table(paste0(pfile,".TMP.r2.bim"), c("CHR","SNP","cm","BP","A1","A2"),col_types="ccdicc")
  bim <- mutate(bim,chrpos=paste(CHR,BP,sep = ":"))
  bim.keep <- subset(bim, grepl("rs",SNP) & !(SNP %in% pvar.match$SNP), select=c("SNP","chrpos"))
  write.table(bim.keep$SNP,file=paste0(pfile,".TMP.r2.flip"), col.names = F, row.names = F, quote = F)
  ## Flip RS matches again
  system2(PLINK1, args=paste0(" --threads ",ncores," --bfile ",pfile,".TMP.r2 --extract ",pfile,".TMP.r2.flip --flip ",pfile,".TMP.r2.flip --make-bed --out ",pfile,".TMP.r2b"))
  ## Extract matched markers from pass 1 and unmatched markers from pass 2 from fist pass file
  pvar.keep <- rbind(pvar.match,subset(pvar.flip,!chrpos %in% bim.keep$chrpos))
  write.table(pvar.keep$SNP,file=paste0(pfile,".TMP.r1.noflip"), col.names = F, row.names = F, quote = F)
  system2(PLINK2, args=paste0(" --threads ",ncores," --pfile ",pfile,".TMP.r1 --extract ",pfile,".TMP.r1.noflip --output-chr 26 --make-bed --out ",pfile,".TMP.r1a"))
  # Merge and convert to PLINK2 standard format
  system2(PLINK1, args=paste0(" --threads ",ncores," --bfile ",pfile,".TMP.r1a --bmerge ",pfile,".TMP.r2b --make-bed --out ",outfile,".TMP"))
  system2(PLINK2, args=paste0(" --threads ",ncores," --bfile ",outfile,".TMP -output-chr MT --make-pgen --sort-vars --out ",outfile,".TMP.2"))
  system2(PLINK2, args=paste0(" --threads ",ncores," --pfile ",outfile,".TMP.2 --make-bed --out ",outfile))
  bim.output <- read_table(paste0(outfile,".bim"), c("CHR","SNP","cm","BP","A1","A2"),col_types="ccdicc")
  bim.output.match <- bim.output$SNP[grep("rs",bim.output$SNP)]
  message(length(bim.output.match)," markers matched to dbSNP RS identifiers.")
  # Clean up temporary files
  invisible(file.remove(list.files(path=outdir,pattern="*.TMP.*",full.names = T)))
}