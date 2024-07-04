# Biogeographic Ancestry Estimation: 1000 Genomes Project Reference Panel

## Required files
Please download a copy of the 1000 Genomes reference panel (build 37, release 2016-05-05) from:
https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg

Then use PLINK2 to [*--keep*](https://www.cog-genomics.org/plink/2.0/filter#sample) the individuals in **kgp3.array_snps.id** and [*--extract*](https://www.cog-genomics.org/plink/2.0/filter#variant) the markers in **kgp3.array_snps.snplist**.
Write the output with [*--make-pgen*](https://www.cog-genomics.org/plink/2.0/data#make_pgen) to a set of pgen/pvar/psam files with the stem **kgp3.array_snps**.

Finally, download and decompress cross-population F<sub>ST</sub> estimates from here: https://walters.psycm.cf.ac.uk/dragondata/ancestry_reference_1kgp.tar

MD5 checksum: https://walters.psycm.cf.ac.uk/dragondata/ancestry_reference_1kgp.tar.md5

## Last update
This 1000 Genomes release was last updated by the PLINK authors on **2022-04-13**.
