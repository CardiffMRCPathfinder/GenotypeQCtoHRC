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
king_from_gds <- function(gdsfile,ncores,mom=F) {
  genofile <- snpgdsOpen(gdsfile)
  message("Reticulating splines...")
  if (isTRUE(mom)){kindata <- snpgdsIBDMoM(genofile, num.thread=ncores,kinship=T)
  } else {kindata <- snpgdsIBDKING(genofile, num.thread=ncores)}
  colnames(kindata$kinship) <- kindata$sample.id
  rownames(kindata$kinship) <- kindata$sample.id
  snpgdsClose(genofile)
  return(kindata)}

pcair_from_gds <- function(gdsfile,kindata,snpfile,ncores,qc4_pcair,mom=F) {
  genodata <- GdsGenotypeReader(gdsfile)
  genodata <- GenotypeData(genodata)
  snpindep <- unlist(read.table(snpfile))
  rpcdata <- pcair(genodata,kinobj=kindata$kinship,divobj=kindata$kinship,snp.include=snpindep,
				   kin.thresh=if(isTRUE(mom)){0.044}else{0.022},
				   div.thresh=if(isTRUE(mom)){0}else{-0.022},
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
  genodata.n <- nscan(genodata)
  genodata.blocksize <- if(genodata.n>3000){ceiling(genodata.n/ceiling(genodata.n/2999))} else {3000}
  snpindep <- snpindep[snpindep%in%genodata.snpid$SNP]
  # Run PC-Relate
  genoiter <- GenotypeBlockIterator(genodata, snpBlock=10000, snpInclude=snpindep)
  if(genodata.n>3000){message("Using ",genodata.blocksize," individuals per block to ensure equal-sized blocks")}
  phidata <- pcrelate(genoiter, pcs=rpcdata$vectors[,1:qc4_pcrelate], training.set=rpcdata$unrels,
                      scale="overall", small.samp.correct=(genodata.n<=100), 
					  sample.block.size=genodata.blocksize,
                      BPPARAM=if(ncores>1) {snow_hatchlings} else {BiocParallel::SerialParam()})
  # Clean up
  close(genodata)
  return(phidata)
}
