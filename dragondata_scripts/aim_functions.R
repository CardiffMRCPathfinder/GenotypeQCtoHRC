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

# Ancestry Inference Model: Find ancestry-informative markers
find_aims <- function(bimfile,bimfile.ref,lrld.regions.txt,fst.filevector,ncores) {
  message("Reading pairwise FST values from reference files...")
  if(ncores>1) {
    hatchlings <- parallelly::makeClusterPSOCK(opt$ncores, autoStop = TRUE, rscript_startup = quote(local({
      options(tidyverse.quiet = TRUE)
      suppressWarnings(library("tidyverse"))})))
    fst.list <- pblapply(fst.filevector, read_table, col_names=T,col_types = "iicid",
                         progress = F,show_col_types = F,na = c("NA","nan"), cl=hatchlings)  
  } else { fst.list <- pblapply(fst.filevector, read_table, col_names=T,col_types = "iicid",
                                progress = F,show_col_types = F,na = c("NA","nan"))}
  # Remove Long-range LD regions from input marker list
  lrld.regions <- read_table(lrld.regions.txt, col_names=c("CHR","START","STOP"),col_types = "iii")
  for (i in 1:nrow(lrld.regions)) {
    bimfile <- subset(bimfile, !(CHR==lrld.regions$CHR[i] & BP>=lrld.regions$START[i] & BP<=lrld.regions$STOP[i]))
  }
  # Ensure markers from input dataset and reference are in the same strand
  snp.match.strand <- distinct(rbind(inner_join(select(bimfile,c("SNP","A1","A2")),select(bimfile.ref,c("SNP","A1","A2")),by=join_by(SNP,A1,A2)),
                                     inner_join(select(bimfile,c("SNP","A1","A2")),select(bimfile.ref,c("SNP","A1","A2")),by=join_by(SNP,A1==A2,A2==A1))),
                               SNP, .keep_all=T) %>%
    subset(SNP != ".") # Remove any SNP with no name
  message("Extracting markers with top 2.5% FST values in reference dataset...")
  fst.threshold <- pblapply(fst.list, function(x){
    # Select SNPs in common between reference and test dataset
    # Make sure selected SNPs will not be dropped because of strand issues
    x.trim <- subset(x,ID%in%snp.match.strand$SNP)
    # Select AIMs with top 2.5% FST values
    fst.top <- quantile(x.trim$HUDSON_FST,0.975,type=8,na.rm=T)
    # P-value type scaling for clumping (not a real p-value)
    x.trim$HUDSON_FST_P <- 10^(VGAM::logclink(x.trim$HUDSON_FST))
    return(subset(x.trim,HUDSON_FST>=fst.top,select=c("ID","HUDSON_FST","HUDSON_FST_P")))
  })
  fst.threshold <- as_tibble(bind_rows(fst.threshold))
  aims <- distinct(arrange(fst.threshold,HUDSON_FST_P),ID,.keep_all = T)
  message("Identified ",nrow(aims)," AIMs overlapping input dataset.")
  return(aims) }

# Ancestry Inference Model: Project input sample into reference eigenvector space
pcair_project_from_gds <- function(gdsfile,famfile.input,ncores) {
  genodata <- GdsGenotypeReader(gdsfile)
  genodata <- GenotypeData(genodata)
  genodata.ids <- getScanID(genodata)
  ref.ids <- tail(genodata.ids,-nrow(famfile.input))
  #message("First reference sample identified as ",head(ref.ids,1),". Last reference sample identified as ",tail(ref.ids,1))
  message("Projecting input genotypes into reference PCA space...")
  rpcdata <- pcair(genodata,unrel.set=ref.ids, 
                   num.cores=ncores, eigen.cnt=100,
                   algorithm=if(nscan(genodata)>=10000){"randomized"} else {"exact"})
  close(genodata)
  return(rpcdata)
}

# Ancestry Inference Model: Train linear discriminant analysis model on reference data
lda_train <- function(qc5.pcair,qc5.tw,qc5.ref.bfile) {
  # Create PCA data frame
  qc5.pca <- as.data.frame(qc5.pcair$vectors)
  colnames(qc5.pca) <- paste0("PC",1:100)
  qc5.pca$IID <- row.names(qc5.pca)
  qc5.pca <- subset(qc5.pca,select=c("IID",paste0("PC",1:100))) %>%
    arrange(IID) %>%
    mutate(ref=ifelse(IID%in%qc5.pcair$unrels,T,F))
  qc5.ref.ancestry <- qc5.ref.bfile[["pheno"]] %>%
    mutate(IID2=paste(FID,IID,sep="-")) %>%
    arrange(IID2)
  # Train and tune
  lda.train.data <- as.data.frame(subset(qc5.pca,ref==T,select=paste0("PC",1:qc5.tw)))
  lda.train.class <- as.character(qc5.ref.ancestry$biogeographic_group)
  message("Training Linear Discriminant Analysis (LDA) model for ancestry inference...")
  lda.train <- train(x=lda.train.data,y=lda.train.class, method="lda", metric="Accuracy", trControl=trainControl(method="repeatedcv", number=10, repeats=5, classProbs = T, savePredictions = "final"))
  return(lda.train)}

# Ancestry Inference Model: Select best probability threshold for ancestry inference
lda_best_threshold <- function(lda_train) {
  pops <- data.frame(biogeographic_group=unique(lda_train$pred$obs),threshold=NA)
  message("Finding optimal probability threshold for assigning discrete ancestry categories...")
  for (ancestry in pops$biogeographic_group) {
    lda.train.ancestry <- subset(lda_train$pred,select=c("pred","obs",ancestry))
    colnames(lda.train.ancestry) <- c("pred","obs","prob")
    lda.train.ancestry <- mutate(lda.train.ancestry,
                                 pred=ifelse(pred==ancestry,ancestry,"OTHER"),
                                 pred=factor(pred,levels=c(ancestry,"OTHER")),
                                 obs=ifelse(obs==ancestry,ancestry,"OTHER"),
                                 obs=factor(obs,levels=c(ancestry,"OTHER")))
    # Find best threshold from 0.5 to 1
    lda.train.performance <- threshold_perf(lda.train.ancestry, obs, prob, thresholds = seq(0.5, 1, by = 0.01))
    lda.train.performance <- subset(lda.train.performance,.metric=="j_index")
    lda.train.threshold <- unlist(subset(lda.train.performance,.estimate==max(.estimate),select=".threshold"),use.names = F)
    lda.train.threshold <- min(lda.train.threshold)
    # Repeat to refine
    lda.train.performance <- threshold_perf(lda.train.ancestry, obs, prob, thresholds = seq(max(0.501,lda.train.threshold-0.01), min(lda.train.threshold+0.01,1), by = 0.001))
    lda.train.performance <- subset(lda.train.performance,.metric=="j_index")
    lda.train.threshold <- unlist(subset(lda.train.performance,.estimate==max(.estimate),select=".threshold"),use.names = F)
    lda.train.threshold <- min(lda.train.threshold)
    # Assign threshold
    pops$threshold[pops$biogeographic_group == ancestry] <- lda.train.threshold
    rm(lda.train.ancestry,lda.train.performance,lda.train.threshold)}
  return(pops)}

# Ancestry Inference Model: Apply linear discriminant analysis model to input data
lda_test <- function(famfile,qc5.pcair,qc5.tw,qc5.train,qc5.youden) {
  # Create PCA data frame
  qc5.pca <- as.data.frame(qc5.pcair$vectors)
  colnames(qc5.pca) <- paste0("PC",1:100)
  qc5.pca$IID <- row.names(qc5.pca)
  qc5.pca <- subset(qc5.pca,select=c("IID",paste0("PC",1:100))) %>%
    arrange(IID) %>%
    mutate(ref=ifelse(IID%in%qc5.pcair$unrels,T,F))
  # Rearrange FAM file
  if (max(famfile$IID %in% qc5.pca$IID)==1) {
    fam <- mutate(famfile, IID2=IID) %>%
      select(c("FID","IID","IID2"))
  } else {
    fam <- mutate(famfile, IID2=paste(FID,IID,sep="-")) %>%
      select(c("FID","IID","IID2"))
  }
  # Predict ancestries
  lda.test.data <- as.data.frame(subset(qc5.pca,ref==F,select=paste0("PC",1:qc5.tw)))
  row.names(lda.test.data) <- unlist(subset(qc5.pca,ref==F,select="IID"),use.names = F)
  message("Applying LDA model to input genotype data...")
  lda.test <- predict(object=qc5.train,newdata=lda.test.data,type="prob")
  lda.test <- round(lda.test,3) # Restrict output to three decimal places
  lda.test$IID2 <- rownames(lda.test)
  lda.test <- right_join(fam,lda.test,by="IID2")
  lda.test$biogeographic_prediction <- "UNK"
  for (i in 1:nrow(qc5.youden)) {
    bga <- as.character(qc5.youden$biogeographic_group[i])
    lda.test$biogeographic_prediction[lda.test[bga]>=qc5.youden$threshold[i]] <- bga
  }
  lda.test <- mutate(lda.test, 
                     biogeographic_prediction=factor(biogeographic_prediction,
                                                     labels=c("Unclassified","African American / Afro-Caribbean","American","East Asian","European","Latino","Near Eastern / North African","Oceanian","Central / South Asian","Sub-Saharan African"),
                                                     levels=c("UNK","AAC","AME","EAS","EUR","LAT","NEA","OCE","SAS","SSA")))
  return(lda.test)}
