cd /home/guptaa12/beegfs/guptaa12/
mkdir v8de8bueja89e && cd v8de8bueja89e
Rscript -e "remotes::install_github('noahlorinczcomi/gent',lib=getwd())"
mkdir -p gent_test mugent_test
mv genome_wide_gent_test.* gent_test
mv genome_wide_mugent_test.* mugent_test
sbatch gent_test/genome_wide_gent_test.slurm
sbatch mugent_test/genome_wide_mugent_test.slurm
