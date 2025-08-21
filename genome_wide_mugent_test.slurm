#!/bin/bash

#SBATCH --job-name=noah_genome_wide_mugent_test
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guptaa12@ccf.org

module load R/4.3.3
cd /home/guptaa12/beegfs/guptaa12/v8de8bueja89e/mugent_test
# download AD GWAS summary statistics
afr_loc="https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_B2_ALL_afr_leave23andme_20220403.tsv.gz"
eur_loc="https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz"
[ ! -f "HOSPITALIZED_vs_POPULATION_AFR.tsv.gz" ] && wget -O HOSPITALIZED_vs_POPULATION_AFR.tsv.gz ${afr_loc}
[ ! -f "HOSPITALIZED_vs_POPULATION_EUR.tsv.gz" ] && wget -O HOSPITALIZED_vs_POPULATION_EUR.tsv.gz ${eur_loc}

# download pre-computed LD matrices
wget https://raw.githubusercontent.com/noahlorinczcomi/gent/refs/heads/main/LD.fetch
chmod +x LD.fetch
[ ! -d "ld_matrices/EUR" ] && ./LD.fetch EUR ld_matrices
[ ! -d "ld_matrices/AFR" ] && ./LD.fetch AFR ld_matrices

Rscript genome_wide_mugent_test.R > log.txt 2>&1
