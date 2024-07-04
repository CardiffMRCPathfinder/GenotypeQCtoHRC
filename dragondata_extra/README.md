# Software dependencies

## Required files
Please download here the PLINK and PLINK2 binaries appropriate for your system. 
You can find these at:

https://www.cog-genomics.org/plink/1.9/

https://www.cog-genomics.org/plink/2.0/

Our scripts were tested with **PLINK1 beta 7.2** and **PLINK2 alpha 5.11**, but later versions might work.

## Note for MAC users
Please rename your plink binaries to either *plink.osx* (PLINK1) or *plink2.osx"* (PLINK2).

## Folder organisation
With all files downloaded and tarballs decompressed, the folder should look like this:
```
dragondata_extra/
├── .gitignore
├── lrld.txt
├── lrld4plink.txt
├── plink
├── plink2
├── ancestry_reference/
│   └── 1kgp/
│       ├── kgp3.array_snps.balanced.ids
│       ├── kgp3.array_snps.id
│       ├── kgp3.array_snps.norel.ids
│       ├── kgp3.array_snps.norel.pheno
│       ├── kgp3.array_snps.pgen
│       ├── kgp3.array_snps.pheno
│       ├── kgp3.array_snps.psam
│       ├── kgp3.array_snps.pvar
│       ├── kgp3.array_snps.snplist
│       ├── kgp3.array_snps.Huddart2018.AAC.EAS.fst.var.gz
│       ├── kgp3.array_snps.Huddart2018.AAC.EUR.fst.var.gz
│       ├── ...
│       └── README.md
├── chipendium2019/
│   ├── Axiom_UKB.processed.txt.gz
│   ├── BDCHP-1X10-HUMANHAP550_11218540_C-b36.processed.txt.gz
│   ├── ...
│   └── README.md
├── dbSNP/
│   ├── all_hg38.snp.sites.bim
│   ├── all_phase3_ns.snp.sites.bim
│   └── hapmap_r23a.sites.bim
├── GH1.4.20/
│   ├── GenotypeHarmonizer.bat
│   ├── GenotypeHarmonizer.jar
│   ├── GenotypeHarmonizer.sh
│   ├── HarmonizeBinaryPlinkExample.bat
│   ├── HarmonizeBinaryPlinkExample.sh
│   ├── HarmonizeShapeit2HapsExample.bat
│   ├── HarmonizeShapeit2HapsExample.sh
│   ├── LICENSE
│   ├── README.md
│   ├── exampleData/
│   │   ├── 1000gCeuChr20Mb6.vcf.gz
│   │   ├── 1000gCeuChr20Mb6.vcf.gz.tbi
│   │   └── ...
│   └── lib/
│       ├── annotations-1.3.2.jar
│       ├── ant-1.8.4.jar
│       └── ...
├── imputation-sites/
│   ├── ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
│   ├── ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz.tbi
│   ├── HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
│   └── HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz.tbi
└── liftover-chains/
    ├── GRCh37ToGRCh38.chain
    ├── GRCh37ToNCBI36.chain
    ├── GRCh38ToGRCh37.chain
    ├── Hg18ToHg19.chain
    ├── Hg18ToHg38.chain
    ├── Hg19.Ormond2021.CUPs.bed
    ├── Hg19ToHg18.chain
    ├── Hg19ToHg38.chain
    ├── Hg38.Ormond2021.CUPs.bed
    ├── Hg38ToHg19.chain
    ├── NCBI36ToGRCh37.chain
    ├── NCBI36ToGRCh38.chain
    └── README.md
```
