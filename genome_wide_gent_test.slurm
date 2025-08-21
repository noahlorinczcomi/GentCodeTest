#!/bin/bash

#SBATCH --job-name=noah_genome_wide_gent_test
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guptaa12@ccf.org

module load R/4.3.3
cd /home/guptaa12/beegfs/guptaa12/v8de8bueja89e/gent_test
# download AD GWAS summary statistics
ad_fp="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt"
[ ! -f "Kunkle_etal_Stage1_results.txt" ] && wget ${ad_fp}

# download pre-computed LD matrices
wget https://raw.githubusercontent.com/noahlorinczcomi/gent/refs/heads/main/LD.fetch
chmod +x LD.fetch
[ ! -d "ld_matrices/EUR" ] && ./LD.fetch EUR ld_matrices

Rscript genome_wide_gent_test.R > log.txt 2>&1
