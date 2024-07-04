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

liftOver_bim <- function(bimfile,check.rev=T,ensembl=F,build.in,build.out,chainpath) {

# Identify genomic builds
build.in <- tolower(build.in)
build.out <- tolower(build.out)
build.in <- case_match(tolower(build.in),
                         c(18, 36,"h36","hg18","grch36") ~ "Hg18",
                         c(19, 37,"h37","hg19","grch37") ~ "Hg19",
                         c(38,"h38","hg38","grch38") ~ "Hg38")
build.out <- case_match(tolower(build.out),
                          c(18, 36,"h36","hg18","grch36") ~ "Hg18",
                          c(19, 37,"h37","hg19","grch37") ~ "Hg19",
                          c(38,"h38","hg38","grch38") ~ "Hg38") 

# Error out if builds are not expected
stopifnot("Unrecognised input genomic build"=build.in %in% c("Hg18","Hg19","Hg38"),
          "Unrecognised output genomic build"=build.out %in% c("Hg18","Hg19","Hg38"))

# Do not liftover if same build has been set for in and out arguments
if (build.in == build.out) {
  
  message("Same genomic build specified as input and output arguments. LiftOver aborted.") 
  bim.old <- select(bimfile, c("SNP","BP"))
  return(bim.old)
  
} else {

# Override argument: Do not attempt reverse liftover for GRCh36 -> GRCh38 conversions
if (build.in == "Hg18" & build.in == "Hg38") {check.rev <- F} # No .chain file available

# Munge BIM file
bim <- mutate(bimfile,
CHR=case_match(CHR, .default=CHR,
"23" ~ "X",
"24" ~ "Y",
c("25","XY","PAR1","PAR2") ~ "X",
c("26","MT") ~ "M"))

# Use BIM file data to create a GRanges object
if (isTRUE(ensembl)) {
  # Ensembl positions are 1-based
  bim.track <- GRangesForUCSCGenome(genome=tolower(build.in),chrom=paste0("chr",bim$CHR),ranges=IRanges(start=bim$BP,width=1,names=bim$SNP))
} else {
  # But UCSC positions are 0-based
  bim.track <- GRangesForUCSCGenome(genome=tolower(build.in),chrom=paste0("chr",bim$CHR),ranges=IRanges(start=bim$BP-1,width=1,names=bim$SNP))
}

if (build.in == "Hg18") {
  bim.track.qc <- bim.track # No CUP file available for GRCh36
} else {
  ## Import CUPs (Ormond et al. 2021; DOI:10.1093/bib/bbab069)
  cup <- read_table(paste0(chainpath,"/",build.in,".Ormond2021.CUPs.bed"),col_names = c("chrom","start","end"),col_types = "cii")
  if (isTRUE(ensembl)) {
    cup <- with(cup,GRangesForUCSCGenome(genome=tolower(build.in),chrom,ranges=IRanges(start+1,end+1)))
  } else {
    cup <- with(cup,GRangesForUCSCGenome(genome=tolower(build.in),chrom,ranges=IRanges(start,end)))
  }
  # Find SNPs affected
  bim.cup <- findOverlaps(bim.track,cup)
  bim.cup <- queryHits(bim.cup)
  bim.cup <- bim.track[bim.cup]
  # Remove SNPs affected
  bim.track.qc <- unlist(subtract(bim.track,bim.cup))
}

# Liftover process
## Change genomic build name and track style for Ensembl
if (isTRUE(ensembl)) {
  seqlevelsStyle(bim.track.qc) <- "Ensembl"
  build.in <- case_match(build.in, "Hg18" ~ "NCBI36", "Hg19" ~ "GRCh37", "Hg38" ~ "GRCh38")
  build.out <- case_match(build.out, "Hg18" ~ "NCBI36", "Hg19" ~ "GRCh37", "Hg38" ~ "GRCh38")
}

## Import chain files
ch = import.chain(paste0(chainpath,"/",build.in,"To",build.out,".chain"))
if (isTRUE(check.rev)) {rev.ch = import.chain(paste0(chainpath,"/",build.out,"To",build.in,".chain"))}

# Liftover
bim.lifted = unlist(liftOver(bim.track.qc, ch))

if (isTRUE(check.rev)) {
  # Reverse Liftover (could highlight more CUPs)
  bim.lifted_rev = unlist(liftOver(bim.lifted, rev.ch))
  # Find variants that survive reverse liftover (conversion-stable)
  snp.lifted <- findOverlaps(bim.track.qc,bim.lifted_rev,type="equal")
  snp.lifted <- bim.track.qc[queryHits(snp.lifted)]
  snp.lifted <- names(snp.lifted)
  # Get lifted ranges for conversion-stable variants
  bim.lifted.qc <- bim.lifted[names(bim.lifted) %in% snp.lifted]
} else {
  bim.lifted.qc <- bim.lifted
}

# Generate output data frame
if (isTRUE(ensembl)) {
  bim.new <- as.data.frame(bim.lifted.qc) %>%
    rownames_to_column() %>%
    select(c("rowname","start"))
} else {
  bim.new <- mutate(as.data.frame(bim.lifted.qc),BP=start+1) %>% # Convert positions to 1-based
    rownames_to_column() %>%
    select(c("rowname","BP"))
}
return(bim.new)
}
}