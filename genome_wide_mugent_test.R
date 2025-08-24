setwd('/home/guptaa12/beegfs/guptaa12/v8de8bueja89e')
library(gent,lib=getwd())
#source('https://raw.githubusercontent.com/noahlorinczcomi/gent/refs/heads/main/genomewide.R')
library(dplyr)
library(ggplot2)
setwd('mugent_test')

# load full GWAS summary statistics for European (eur_gwas) and African (afr_gwas) populations
eur_gwas=data.table::fread('HOSPITALIZED_vs_POPULATION_EUR.tsv.gz') # European COVID GWAS
afr_gwas=data.table::fread('HOSPITALIZED_vs_POPULATION_AFR.tsv.gz') # African COVID GWAS
# add Z-scores
eur_gwas=eur_gwas %>% mutate(z=all_inv_var_meta_beta/all_inv_var_meta_sebeta)
afr_gwas=afr_gwas %>% mutate(z=all_inv_var_meta_beta/all_inv_var_meta_sebeta)

# Genome-wide multi-population gene-based association testing (MuGenT)
covid_result=mugent_genomewide(
    gwas_list = list(EUR=eur_gwas, AFR=afr_gwas),     # full GWAS summary statistics in each population
    ld_population_list = list(EUR='EUR', AFR='AFR'),  # populations of the LD references to use (matching GWAS order)
    ld_directory = 'ld_matrices',                     # directory in which all LD matrices were saved (second argument to LD.fetch)
    KbWindow = 50,                                    # Kb window size of SNP-gene set assignment
    snp_list = list(EUR='rsid', AFR='rsid'),          # column names of unique SNP identifiers
    chromosome_list = list(EUR='#CHR', AFR='#CHR'),   # column names of chromosomes
    position_list = list(EUR='POS', AFR='POS'),       # column names of SNP base pair positions (hg19)
    effect_allele_list = list(EUR='ALT', AFR='ALT'),  # effect allele column names
    z_statistic_list = list(EUR='z', AFR='z'),        # Z-statistic column is not present in either GWAS so set NULL
    verbose = TRUE)                                   # TRUE indicates progress should be printed to the console 

# Plot summaries
#p1=gent_manhattan(covid_result)
#p2=gent_qq(covid_result) + ggtitle('QQ plot of COVID-19 MuGenT (EUR and AFR)')

# Fine-mapping
finemapping_results_clumping=gent_finemap(
    gent_results=covid_result$MuGenT,       # output of gent_genomewide(), mugent_genomewide(), or xgent_genomewide()
    ld_population='EUR',                    # LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR'
    gwas_n=1e6+5e5,                         # GWAS sample size
    index_genes=NULL,                       # if NULL, will find loci automatically for fine-mapping
    clump_p=0.05/12727,                     # clumping P-value threshold for significance. Only loci with P<this threshold may be fine-mapped
    clump_r2=0.01,                          # genes whose squared correlation is below this threshold are considered statistically independent
    window_kb_width=2000,                   # size of locus windows (kilobases) in which gene-level fine-mapping will be performed
    verbose=TRUE)                           # if TRUE, messages will be printed to the console

finemapping_results_setgenes=gent_finemap(
    gent_results=covid_result$MuGenT,                   # output of gent_genomewide(), mugent_genomewide(), or xgent_genomewide()
    ld_population='EUR',                                # LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR
    gwas_n=1e6+5e5,                                     # GWAS sample size
    index_genes=c('MAPT','ACE'),                        # genes indexing test loci
    window_kb_width=2000,                               # size of locus windows (kilobases) in which gene-level fine-mapping will be performed
    verbose=TRUE)                                       # if TRUE, messages will be printed to the console

# save data and plots
saveRDS(covid_result,'covid_result.Rds')
saveRDS(finemapping_results_clumping,'finemapping_results_clumping.Rds')
saveRDS(finemapping_results_setgenes,'finemapping_results_setgenes.Rds')
#ggsave('mugent_manhattan.pdf',p1,width=6,height=4)
#ggsave('mugent_qq.pdf',p2,width=4,height=4)


