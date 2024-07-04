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

# QC4: PC-AiR including initial relatedness estimates from KING
king_from_gds <- function(gdsfile,ncores) {
  genofile <- snpgdsOpen(gdsfile)
  message("Reticulating splines...")
  kindata <- snpgdsIBDKING(genofile, num.thread=ncores)
  colnames(kindata$kinship) <- kindata$sample.id
  rownames(kindata$kinship) <- kindata$sample.id
  snpgdsClose(genofile)
  return(kindata)}

pcair_from_gds <- function(gdsfile,kindata,snpfile,ncores,qc4_pcair) {
  genodata <- GdsGenotypeReader(gdsfile)
  genodata <- GenotypeData(genodata)
  snpindep <- unlist(read.table(snpfile))
  rpcdata <- pcair(genodata,kinobj=kindata$kinship,divobj=kindata$kinship,snp.include=snpindep, 
                   num.cores=ncores, eigen.cnt=qc4_pcair,
                   algorithm=if(nscan(genodata)>=10000){"randomized"} else {"exact"})
  close(genodata)
  return(rpcdata)
}

# QC4: PC-Relate
pcrelate_from_gds <- function(gdsfile,rpcdata,snpfile,qc4_pcrelate,ncores) {
  # Set up multithreading
  if(ncores>1) {snow_hatchlings <- BiocParallel::SnowParam(workers = ncores, type = "SOCK")}
  # Load genotype data in GDS
  genodata <- GdsGenotypeReader(gdsfile)
  genodata <- GenotypeData(genodata)
  # Select only autosomal LD-independent SNPs
  snpindep <- unlist(read.table(snpfile))
  genodata.snpid <- data.frame(CHR=getChromosome(genodata),SNP=getSnpID(genodata))
  genodata.snpid <- subset(genodata.snpid,CHR%in%c(1:22))
  snpindep <- snpindep[snpindep%in%genodata.snpid$SNP]
  # Run PC-Relate
  genoiter <- GenotypeBlockIterator(genodata, snpBlock=10000, snpInclude=snpindep)
  phidata <- pcrelate(genoiter, pcs=rpcdata$vectors[,1:qc4_pcrelate], training.set=rpcdata$unrels,
                      scale="overall", small.samp.correct=(nscan(genodata)<=100),
                      BPPARAM=if(ncores>1) {snow_hatchlings} else {BiocParallel::SerialParam()})
  # Clean up
  close(genodata)
  return(phidata)
}
