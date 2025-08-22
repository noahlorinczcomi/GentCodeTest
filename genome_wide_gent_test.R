setwd('/home/guptaa12/beegfs/guptaa12/v8de8bueja89e')
library(gent,lib=getwd())
source('https://raw.githubusercontent.com/noahlorinczcomi/gent/refs/heads/main/genomewide.R')
library(dplyr)
library(ggplot2)
setwd('gent_test')

# load Alzheimer's disease GWAS data
ad_gwas=data.table::fread('Kunkle_etal_Stage1_results.txt') # Alzheimer's disease GWAS
ad_gwas %>% mutate(z=Beta/SE)

# Genome-wide gene-based association testing (GenT)
ad_result=gent_genomewide(
    gwas=ad_gwas,                  # full GWAS summary statistics
    KbWindow=50,                   # Kb window size. SNPs which are within this window will be tested in gene-specific sets
    ld_population='EUR',           # population of the LD reference to use. must match your GWAS population
    ld_directory='ld_matrices',    # directory in which LD matrices were saved (second argument to LD.fetch)
    snp='MarkerName',              # column name of unique SNP identifier in ad_gwas
    chromosome='Chromosome',       # column name of chromosome in ad_gwas
    position='Position',           # column name of SNP base pair position (hg19) in ad_gwas
    effect_allele='Effect_allele', # column name of SNP effect allele in ad_gwas
    z='z',                        # Z-statistic column is not present in ad_gwas so set NULL
    verbose=TRUE)                  # TRUE indicates that progress should be printed to the console

# Plot summaries
#p1=gent_manhattan(ad_result) + coord_cartesian(ylim=c(0,100))
#p2=gent_qq(ad_result) + ggtitle('QQ plot of GenT (AD)')

# Fine-mapping
finemapping_results_clumping=gent_finemap(
    gent_results=ad_result,      # output of gent_genomewide(), mugent_genomewide(), or xgent_genomewide()
    ld_population='EUR',            # LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR'
    gwas_n=158000,                  # GWAS sample size
    index_genes=NULL,               # if NULL, will find loci automatically for fine-mapping
    clump_p=0.05/12727,             # clumping P-value threshold for significance. Only loci with P<this threshold may be fine-mapped
    clump_r2=0.01,                  # genes whose squared correlation is below this threshold are considered statistically independent
    window_kb_width=2000,           # size of locus windows (kilobases) in which gene-level fine-mapping will be performed
    verbose=TRUE)                   # if TRUE, messages will be printed to the console

finemapping_results_setgenes=gent_finemap(
    gent_results=ad_result,                             # output of gent_genomewide(), mugent_genomewide(), or xgent_genomewide()
    ld_population='EUR',                                # LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR
    gwas_n=158000,                                      # GWAS sample size
    index_genes=c('BIN1','APOE','PCSK9'),               # genes indexing test loci
    window_kb_width=2000,                               # size of locus windows (kilobases) in which gene-level fine-mapping will be performed
    verbose=TRUE)                                       # if TRUE, messages will be printed to the console

# save data and plots
saveRDS(ad_result,'ad_result.Rds')
saveRDS(finemapping_results_clumping,'finemapping_results_clumping.Rds')
saveRDS(finemapping_results_setgenes,'finemapping_results_setgenes.Rds')
#ggsave('gent_manhattan.pdf',p1,width=6,height=4)
#ggsave('gent_qq.pdf',p2,width=4,height=4)

