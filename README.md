# GenotypeQCtoHRC

This automated Rmd pipeline is designed to perform stringent quality control of genotype data in ped/map or bed/bim/fam format, align against the Haplotype Reference Consortium V1.1 reference panel and output in compressed vcf.gz format ready for imputation using the Michigan Imputation Server. Each pipeline produces an Rmd QC report, which gives you information regarding number of individuals SNPs excluded based on user parameters, sex discrepencies, possible sex chromosome abnormalities, relatedness information, relatedness adjusted PCAs and ancestry probabilities.  
  
We demonstrate the utilty of the pipeline using publically available data: https://www.nature.com/articles/nature09103  
Genotype data available for download here: https://evolbio.ut.ee/jew/.  

The automated report can be found for the Genome Wide Survey of the Jewish People is located in this file: GenotypeQCtoHRC.V1.5.GWS_JP.html.

# Pipeline Order Of Operations  
The pipeline performs the following operations:
1) Converts ped / map into bed/bim/fam (if applicable).  
2) Checks bim file for number of chromosomes in data. If < 22 then stops. If >=22 then continues.
3) Checks SNP overlap against a range of available genotyping platforms. It returns the genotyping platform that contains the highest match based on overlapping values(chr+bp).
4) Identifies Y/Mito SNPs and saves to a separate bed/bim/fam (if applicable).  
5) Performs sex checks against biological sex and gender in the fam file. If there is a discordance, individuals are flagged and removed. If sex is set to missing, these individuals are retained. 
6) Runs genome-harmoniser using the Haplotype Reference Consortium v1.1 reference panel.
i) Runs genome-harmoniser and produces a report for why each SNP is excluded.  
ii) After harmonisation, checks for duplicate RS IDs.   
iii) Checks for total number of SNPs removed by GH. If the total number of SNPs retained is > 2/3, the pipeline continues otherwise runs liftover to update build from hg18 to hg19 and creates new bed/bim/fam.  
iv) Re-runs genome-harmoniser. If the number of SNPs removed is >2/3 the pipeline stops.    
7) Identifies SNPs low genotyping rates < specified threshold, writes to file and generates plot.  
8) Identifies individual low genotyping rates < specified threshold, writes to file and generates plot.  
9) Identifies SNPs failing Hardy-Weinberg equilibrium < specified threshold, writes to file and generates plot.
10) Sex chromosome checks and abnormalities.  
i) Checks for individuals with X (F) and Y (Call Rate) call rates < specific thresholds so identify outlyers indictive of sex chr abnormalities.  
ii) Writes these individuals to file, but does not exclude them.
11) Saves to vcf format and sorts the resulting file. This is now in the correct format for imputation of autosomal and X chrs using Michigan Imputation Server.  
12) Relatedness / PCA - PC-AiR performs a Principal Components Analysis on genome-wide SNP data for the detection of population structure in a sample that may contain known or cryptic relatedness. Unlike standard PCA, PC-AiR accounts for relatedness in the sample to provide accurate ancestry inference that is not confounded by family structure.
13) Global biogeographical ancestry inference - This analysis is based on the set of Ancestry Informative Markers (AIMs) contained in the EUROFORGEN (PMID: [24631693](https://www.sciencedirect.com/science/article/pii/S1872497314000404), [27208666](https://www.sciencedirect.com/science/article/pii/S1872497316300643)) and Kidd (PMID: [24508742](https://www.sciencedirect.com/science/article/pii/S1872497314000039)) panels. The original sets contain a grand total of 167 AIMs, but this will be lower depending on the number of SNPs in your dataset. As a reference panel, a combination of Human Genome Diversity Project (PMID: [18292342](http://science.sciencemag.org/content/319/5866/1100)) and South Asian Genome Project (PMID: [25115870](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102645)) is used. Ancestry inference is based on Fisher's Linear Discriminant Analysis algorithm (PMID: [16193326](https://link.springer.com/article/10.1007%2Fs00439-005-0012-1)), which has been applied to the first `r paste0(pca.aim1.tw)` principal components of the HGDP+SAGP dataset. This optimal number of PCs has been determined using the Tracy-Widom test for eigenvalues (PMID: [17194218](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190)).  
14) The html report is generated and data is saved in appropriate sub-directories. 

# Usage 

To run the pipeline you will need to modify the accompanying Rmd and .sh script for your own purposes. 
The following lines (87-98) can be easily modified for your own QC parameters:  
qc0.max.femaleF <- 0.2 # Maximum female X-linked homozygosity  
qc0.min.maleF <- 0.8 # Minimum male X-linked homozygosity  
qc1.call.rate <- 0.95 # Minimum marker call rate  
qc2.coverage <- 0.95 # Minimum individual genotyping coverage  
qc3.hwepval <- 1e-6 # Maximum mid p-value for the HWE test  
qc4.maf <- 0.01 # Minimum MAF for PCA  
qc4.pcairs <- 5 # Number of PCs used to correct kinship matrix  
qc4.kin <- 0.1 # Minimum kinship to include pairs of samples in report  
max.aimpcs <- 32 #   
min.prob.aim1 <- 0.8 # Minimal ancestry percentage to classify one person  
min.prob.admx <- 0.4 # Minimal ancestry percentage to send a person  
min.prob.aim2 <- 0.8 # Minimal ancestry percentage to classify one person  


# R Packages

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

# Other software 
plink1.9  
plink2  
GenotypeHarmonizer  
vcf-sort  

# Other required files
SNP coordinates / builds across platforms: Please contact hubbardl@cardiff.ac.uk initially - this will be made available in a repository.  
HRC v1.1 SNP information: ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
