# 2 June 2017
# Genomic landscape of divergence in south Pacific silvereyes
# R script used to reorder vcf file based on satsuma output

# ------------------------------------------------------------

setwd("~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/chromosome")

# STEP 1: Read in and massage the satsuma and scaffold data


library(tidyverse)
#created in previous script
scaff_ord <- read_csv("GCRF_scaffold_order_from_ZFinch_Un.csv") %>%
  rename(scaffold = sca)
scaff_lengths <- read_tsv("scaffold_lengths_clean",col_names =F) %>% rename(scaffold=X1,length=X2)

missing_lengths <- anti_join(scaff_ord, scaff_lengths)

missing_ords <- anti_join(scaff_lengths, scaff_ord)


scaffs <- scaff_ord %>%
  left_join(scaff_lengths) %>%
  rename(ZFCHROM = chr, CHROM = scaffold)
scaffs



#vcf <- read_tsv("../output/PABU_clean.SNP.rm_hets.RegionFID.imuted_gt.fix_name.new_pheno_inclLA.noNA.related_bslmm.sample5mil_bi500K.param.txt", comment = "##", progress = FALSE) #%>%
#rename(CHROM = chr, POS=ps)
#vcf <- read_tsv("../output/PABU_clean.SNP.rm_hets.RegionFID.sort.fix_name.new_pheno_inclLA.noNA.related_bslmm.sample5mil_bi500K.param.txt", comment = "##", progress = FALSE) %>% rename(CHROM = chr, POS=ps)

#do this for the vcf pre gemma, can use this on the PIP file (params)
vcf <- read_tsv("~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/bslmm/GCRF_feather_PC.bslmm.param.txt", comment = "##", progress = FALSE) %>% rename(CHROM = chr, POS=ps)

#vcf<-read_tsv("../vcf_allplates/amke-removed-hets.rmPCA.resident_allotherpops.weir.fst.90outlier.bed_exintreg_genes", comment = "##", progress = FALSE,col_names = F) %>% rename(CHROM=X1,POS=X3,FST=X4,GENE=X8,REGION=X9) %>% dplyr::select(CHROM,POS,FST,GENE,REGION)

#vcf<-read_tsv("~/Dropbox/BGP/PABU/GWAS_mrMLM/Kinship_calc_in_mrMLM_100kb_chr1-n/Manhattan_FASTmrMLM.chr1-n.4ZForder.txt" ,comment = "##", progress = FALSE) %>% rename(CHROM = chrom, POS=pos) 

#this is prepped for the params file
#vcf<-read_tsv("~/Dropbox/BGP/PABU/GWAS_mrMLM/Manhattan_FASTmrMLM.chr1-n.100kb.intermed.4ZForder.txt",comment = "##", progress = FALSE) %>% rename(CHROM = chr, POS=ps)

#vcf<-read_tsv("~/Documents/CurrentProjects/UCLA_UCSC/BGP/PABU/GWAS/GEMMA_inputfile/PABU_clean.SNP.rm_hets.RegionFID.imuted_gt.fix_name.new_pheno_inclLA.noNA.chr1-n.012.4ZF.mrMLM.genotype.txt",comment = "##", progress = FALSE) %>% rename(CHROM = CHR)

#vcf<-read_tsv("~/Documents/CurrentProjects/UCLA_UCSC/BGP/PABU/GWAS/GEMMA_inputfile/PABU_clean.SNP.rm_hets.RegionFID.imuted_gt.fix_name.new_pheno_inclLA.01.map",comment = "##", progress = FALSE,col_names = F) %>% rename(CHROM=X1,POS=X4,SNP=X2) %>% dplyr::select(CHROM,POS,SNP)

#pip10 <- vcf[vcf$gamma>=0.10,]
#vcf <- pip10

combo <- scaffs %>%
  left_join(vcf)
tail(combo)
#combo <- combo[complete.cases(combo),]

#same order as ped file
#combo <- vcf %>% left_join(scaffs)
#tail(combo)
options(scipen = 999)

zf_ified <- combo %>%
  mutate(ZFPOS = {
    ML = floor(mean.loc)  # some temp variables to make it easier to express
    Lo2 = floor(length/2)
    L = length
    ifelse(sca.ori == 1,
           ML - Lo2 + POS,           # forward orientation
           ML - Lo2 + (L - POS))     # reverse orientation
  }) %>%
  dplyr::select(ZFCHROM, ZFPOS, everything()) %>%
  mutate(ZFCHROM = factor(ZFCHROM, levels = unique(ZFCHROM))) %>%   # this is to get them to sort correctly
  arrange(ZFCHROM, ZFPOS)

tail(zf_ified)


zf_ified$CHROM <- NULL
zf_ified$POS <- NULL
zf_ified$mean.loc <- NULL
zf_ified$sca.ori <- NULL
zf_ified$length <- NULL

colnames(zf_ified)[1] <- "CHROM"
colnames(zf_ified)[2] <- "POS"
dim(zf_ified)
zf_ified %>% dplyr::select(ZFCHROM,SNP,ZFPOS) %>% write.table("~/GCRF_clean.SNP.rm_hets.RegionFID.imuted_gt.fix_name.new_pheno_inclLA.ZFchr.map",row.names=F,quote=F,sep="\t")


#Filter so that negative positions and unmapped scaffolds are removed:
zf_ified <- filter(zf_ified, POS >0)
#zf_ified <- na.omit(zf_ified) ##only if there's no NA in dataset except what you introduce
#zf_ified <- zf_ified[!grepl("_Un", zf_ified$`CHROM`),]
head(zf_ified)
dim(zf_ified)
library(pgirmess)
zf_ified %>% 
  dplyr::select(`rs`,CHROM,POS,everything()) %>% 
  write_delim("~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/bslmm/GCRF_featherPC_clean.SNP.rm_hets.RegionFID.imuted_gt.fix_name.noNA.chr1-31ZF.genotype.txt", quote = "none" ,delim = "\t")

# Note: To make into a functional vcf file need to add header. This is done by adding the
# following line to the top of the file:
# "##fileformat=VCFv4.2"
# This was done in atom.
