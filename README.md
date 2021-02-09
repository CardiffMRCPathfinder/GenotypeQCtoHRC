# GenotypeQCtoHRC

Overview: This automated Rmd pipeline is designed to perform stringent quality control of genotype data in ped/map or bed/bim/fam format, align against the Haplotype Reference Consortium V1.1 reference panel and output in compressedd vcf.gz format ready for imputation. Each pipeline produces an Rmd QC report, which gives you information regarding number of individuals SNPs excluded based on user parameters, sex discrepencies, relatedness adjusted PCAs and ancestry information. 




# Required R Packages 
library("tidyverse") # Package version: 1.2.1  
library("ggthemes") # Package version: 4.0.1  
library("scales") # Package version 1.0.0  
library("RColorBrewer") # Package version 1.1-2  
library("evaluate") # devtools::install_github('hadley/evaluate')  
library("data.table") # Package version 1.12.0  
library("rvest") # Package version 0.3.2   
library("kableExtra") # Package version: 1.0.0   
library("gridGraphics") # Package version: 0.3-3   
library("cowplot") # Package version 0.9.4   
library("AssocTests") # Package version:0.0-4    
library("caret") # Package version: 6.0-81   
  
library("BiocManager") # Use to install/update Bioconductor packages. Bioconductor version 3.8. Package version 1.30.4   
library("GWASTools") # Bioconductor package : Package version: 1.28.0   
library("SNPRelate") # Bioconductor package. Package version: 1.16.0   
library("GENESIS") # Bioconductor package. Package version: 2.12.2   
library("snpStats") # Bioconductor package. Package version: 1.32.0   

library("knitr")  
library("prettydoc")  

# Walkthrough

1) Download 1000 genomes genotype data
2) 




5) The generated HTML QC report is in this file: 
