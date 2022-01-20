#!/bin/bash --login
###
#job name
#SBATCH --job-name=
#job stdout file
#SBATCH --output=AFFY.out
#job stderr file
#SBATCH --error=AFFY.err
#maximum job time in D-HH:MM
#SBATCH --time=0-24:00
#number of parallel processes (tasks) you are requesting - maps to MPI processes
#SBATCH --mem-per-cpu=20000
#SBATCH -n 1
#SBATCH -p c_compute_neuro1
module load R

### Path to R libraries
export R_LIBS_USER=/home/username/R/x86_64-pc-linux-gnu-library/3.5

### Path to pandoc
export RSTUDIO_PANDOC="/home/c.c1002680/Software/pandoc-2.5/bin"

### Please replace GWS_JP with the name of your bed/bim/fam
sed s/"NAME_OF_DATASET"/GWS_JP/g GenotypeQCtoHRC.V1.5.Individual.DatasetQC.Rmd > GenotypeQCtoHRC.V1.5.GWS_JP.Rmd
Rscript -e "rmarkdown::render(\"GenotypeQCtoHRC.V1.5.GWS_JP\", params = list(DATASET = \"GWS_JP\"))"
