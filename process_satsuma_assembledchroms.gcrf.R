# Script for processing Satsuma Synteny output for Assembled Chromosomes
# This is based on a script written by Benjamin Van Doren 
# 30 January 2017
# -----------------------------------------------------------
library(tidyverse)
# Set working directory
setwd("~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/chromosome/")

#cat chromosome.1.BCRF.satsuma_summary.chained.sort.out chromosome.1A.BCRF.satsuma_summary.chained.sort.out 
#chromosome.1B.BCRF.satsuma_summary.chained.sort.out chromosome.2.1.BCRF.satsuma_summary.chained.sort.out chromosome.2.2.BCRF.satsuma_summary.chained.sort.out
#chromosome.3.BCRF.satsuma_summary.chained.sort.out chromosome.4.BCRF.satsuma_summary.chained.sort.out 
#chromosome.4A.BCRF.satsuma_summary.chained.sort.out chromosome.5.BCRF.satsuma_summary.chained.sort.out 
#chromosome.6.BCRF.satsuma_summary.chained.sort.out chromosome.7.BCRF.satsuma_summary.chained.sort.out 
#chromosome.8.BCRF.satsuma_summary.chained.sort.out chromosome.9.BCRF.satsuma_summary.chained.sort.out 
#chromosome.10.BCRF.satsuma_summary.chained.sort.out chromosome.11.BCRF.satsuma_summary.chained.sort.out 
#chromosome.12.BCRF.satsuma_summary.chained.sort.out chromosome.13.BCRF.satsuma_summary.chained.sort.out 
#chromosome.14.BCRF.satsuma_summary.chained.sort.out chromosome.15.BCRF.satsuma_summary.chained.sort.out 
#chromosome.16.BCRF.satsuma_summary.chained.sort.out chromosome.17.BCRF.satsuma_summary.chained.sort.out 
#chromosome.18.BCRF.satsuma_summary.chained.sort.out chromosome.19.BCRF.satsuma_summary.chained.sort.out 
#chromosome.20.BCRF.satsuma_summary.chained.sort.out chromosome.21.BCRF.satsuma_summary.chained.sort.out 
#chromosome.22.BCRF.satsuma_summary.chained.sort.out chromosome.23.BCRF.satsuma_summary.chained.sort.out 
#chromosome.24.BCRF.satsuma_summary.chained.sort.out chromosome.25.BCRF.satsuma_summary.chained.sort.out 
#chromosome.26.BCRF.satsuma_summary.chained.sort.out chromosome.27.BCRF.satsuma_summary.chained.sort.out 
#chromosome.28.BCRF.satsuma_summary.chained.sort.out chromosome.LG2.BCRF.satsuma_summary.chained.sort.out 
#chromosome.LG5.BCRF.satsuma_summary.chained.sort.out chromosome.LGE22.BCRF.satsuma_summary.chained.sort.out 
#chromosome.MT.BCRF.satsuma_summary.chained.sort.out chromosome.Z.BCRF.satsuma_summary.chained.sort.out 
#nonchromosomal.fa.BCRF.satsuma_summary.chained.sort.out > BCRF.satsuma_summary.chained.sort.all.out
#cat BCRF.satsuma_summary.chained.sort.all.out| tr ":" "\t" | cut -f 3,4,8-|head (check tail too to make sure you have the right things) This then gets piped into a
#cat BCRF.satsuma_summary.chained.sort.all.out| tr ":" "\t" | cut -f 3,4,8- > BCRF.satsuma_summary.chained.sort.all.processed.out
#output file "processed"
 

# A necessary function
expand.indices = function(df) {
  if (nrow(df)>1) {
    df.expanded = apply(df,1,function(X) { return(X[1]:X[length(X)]) })
    if (class(df.expanded)=="matrix") {
      df.expanded = lapply(apply(df.expanded,2,list),unlist)
    }
    df.cat = do.call("c",df.expanded)
  } else {
    df.cat = df[1,1]:df[1,2]
  }
  return(df.cat)
}

# Read in original Satsuma output file. Note, that I haven't pre-checked this section. I process the sasuma_summary.chained.sort.out using a unix script in terminal rather than this. 
sat = readLines("../PABU.all_chr.satsuma_summary.chained.sort.out")

# The following shortens and breaks up each line of Satsuma Synteny output from this:
# KB442450.1_Taeniopygia_guttata_chromosome_Z_un_1040_1490_LAII01000182.1_Zosterops_lateralis_melanops_scaffold_181_437837_438271_0.793333_-
# to this:
# KB442450.1	chromosome_Z	1040	1490	LAII01000182.1	437837	438271	0.793333	-

#sat = gsub(",_whole_genome_shotgun_sequence","",sat)
#sat =gsub("_Taeniopygia_guttata","\t",sat)
#sat =gsub("_chromosome","Chr",sat)
#sat = gsub("_Zosterops_lateralis_melanops_scaffold_","\tSca_",sat)

# Now that we're done with this one-time pre-processing, save the file.
# Obviously we can now skip the above steps in the future.
fileConn = file("../PABU.all_chr.satsuma_summary.chained.sort.out_processed")

#writeLines(sat, fileConn)
#close(fileConn)


##########################

# Read in the pre-processed file, not yet as a table and split into columns by the tab character. 

library(data.table)
library(dplyr)
chr.map<-fread("BCRF.satsuma_summary.chained.sort.all.processed.out", header = F,sep="\t", fill=T)
head(chr.map)

# name the columns
colnames(chr.map) = c("chr.code","chr","chr.start","chr.end","sca.no","sca.start","sca.end","conserve.strength","sca.ori")

#I do this to sometimes check the different sections of the file
#chr.map %>% dplyr::select(sca.ori) %>% distinct()
chr.map %>% distinct(chr)
chr.map %>% distinct(sca.no)

# Make into a data.frame
chr.map = as.data.frame(chr.map)
head(chr.map)
# Load as a separate file the query scaffold names and sizes
# brown capped length of each scaffold, make this with some cut functions
# cat leucosticte_australis_final_assembly.fasta | grep ">" > scaffold_names
# cat leucosticte_australis_final_assembly.fasta | grep ">" | cut -d "_" -f 8 > lengths
# paste scaffold_names lengths > scaffold_lengths
# sed 's/^>//g' scaffold_lengths > scaffold_lengths_clean
ZLat.scafs = read.table("scaffold_lengths_clean",header=F, sep="\t")
colnames(ZLat.scafs) <- c("scaffold","length")
head(ZLat.scafs)
nrow(ZLat.scafs %>% distinct(scaffold))

#converting some of the variables
#chr.map[,1] = as.character(chr.map[,1])
#chr.map[,2] = as.character(chr.map[,2])
#chr.map[,5] = as.character(chr.map[,5])
#chr.map[,6] = as.character(chr.map[,6])
#for (i in c(3:4,6:8)) {
#  chr.map[,i] = as.numeric(as.character(chr.map[,i]))
#}

head(chr.map)

chr.map2<-chr.map %>% mutate(sca.ori2=if_else(sca.ori=="+",1,-1)) %>% dplyr::select(-sca.ori) %>% rename(sca.ori=sca.ori2)
chr.map = chr.map2

chr.map$len.sca.chunk = chr.map$sca.end-chr.map$sca.start
chr.map$sca.length = ZLat.scafs$length[match(chr.map$sca.no,ZLat.scafs$scaffold)]
chr.map$prop.of.sca = chr.map$len.sca.chunk/chr.map$sca.length

head(chr.map)
# "chr.map" should now look something like this:
#X chr.code   chr    chr.start chr.end   sca.code        sca.no sca.start  sca.end  conserve.strength sca.ori len.sca.chunk sca.length  prop.of.sca
#1 CM000515.1 Chr_1  23946267  23947763  LAII01000001.1  Sca_0  125        1631     0.818182          1       1506          15146310    9.943016e-05
#2 CM000515.1 Chr_1  23947845  23948483  LAII01000001.1  Sca_0  1719       2349     0.724138          1       630           15146310    4.159429e-05
#3 CM000515.1 Chr_1  23948622  23949952  LAII01000001.1  Sca_0  2515       3854     0.874436          1       1339          15146310    8.840437e-05
#4 CM000515.1 Chr_1  23949954  23950127  LAII01000001.1  Sca_0  3852       4025     0.612717          1       173           15146310    1.142192e-05
#5 CM000515.1 Chr_1  23950213  23951375  LAII01000001.1  Sca_0  4097       5371     0.854561          1       1274          15146310    8.411290e-05
#6 CM000515.1 Chr_1  23951376  23951400  LAII01000001.1  Sca_0  5391       5415     1.000000          1       24            15146310    1.584544e-06

# The last column, "prop.of.sca," tells you how much of the scaffold is mapped by that mapping
# It might be low, but what matters is when you add things up.

# Write to csv
write.csv(chr.map, "~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/chromosome/chr.mapAssembled.BCRF.csv")


#------------------------------
#STEP 3: CREATING SCA.CHROM.MAP
library(dplyr)
# load in necessary files:
chr.map = read.csv("chr.mapAssembled.BCRF.csv",stringsAsFactors = F)
#ZLat.scafs = read.table("scaffold_length",header=F,stringsAsFactors = F) %>% rename(scaffold=V1,length=V2)
ZLat.scafs = read.table("scaffold_length_clean",header=F, sep="\t")
colnames(ZLat.scafs) <- c("scaffold","length")
# head(chr.map)

# and necessary function:
expand.indices = function(df) {
  if (nrow(df)>1) {
    df.expanded = apply(df,1,function(X) { return(X[1]:X[length(X)]) })
    if (class(df.expanded)=="matrix") {
      df.expanded = lapply(apply(df.expanded,2,list),unlist)
    }
    df.cat = do.call("c",df.expanded)
  } else {
    df.cat = df[1,1]:df[1,2]
  }
  return(df.cat)
}


# This loop could take a while (a minute or two or three?) but will speed up as you get to the smaller scaffolds
sca.chrom.map = data.frame(ZLat.scaffold=ZLat.scafs$scaffold,best.ZFinch.chrom=NA,mean.loc.ZFinch.chrom=NA,ori=NA,prop.of.scaf=NA)

split.scaffolds = list()

for (i in 1:nrow(ZLat.scafs)) {
  sca = ZLat.scafs$scaffold[i]
  chr.sca = chr.map[chr.map$sca.no==sca,]
  chr.sca = chr.sca[order(chr.sca$sca.start),]
  # get chromosome with most coverage of scaffold
  chr.cov = tapply(chr.sca$prop.of.sca,chr.sca$chr,sum)
  chr.cov = chr.cov[chr.cov>0.2]
  # If more than one chromosome mapped to over 20% of the scaffold's length, save it here for examination
  if (length(chr.cov)>1) { 
    split.scaffolds = c(split.scaffolds,list(rbind(paste(unique(chr.sca$sca)),round(chr.cov,2))))
  } 
  if (length(chr.cov)==0) {
    sca.chrom.map$best.ZFinch.chrom[i] = sca.chrom.map$mean.loc.ZFinch.chrom[i] = sca.chrom.map$ori[i] = sca.chrom.map$prop.of.scaf[i] = NA
  } else { # all good
    chr.max = names(chr.cov[chr.cov==max(chr.cov)])
    # take weighted mean of position by coverage on chromsome
    chr.sca = chr.sca[chr.sca$chr==chr.max,]
    mean.pos = mean(expand.indices(chr.sca[,c("chr.start","chr.end")]))
    stopifnot(sca.chrom.map$ZLat.scaffold[i]==sca)
    sca.chrom.map$best.ZFinch.chrom[i] = chr.max
    sca.chrom.map$mean.loc.ZFinch.chrom[i] = mean.pos
    sca.chrom.map$prop.of.scaf[i] = max(chr.cov)
    ori.mean = weighted.mean(chr.sca$sca.ori,chr.sca$prop.of.sca)
    if (ori.mean>0) { 
      sca.chrom.map$ori[i] = 1
    } else {
      sca.chrom.map$ori[i] = -1
    }
  }
  print(nrow(ZLat.scafs)-i)
}



sca.chrom.map = sca.chrom.map[complete.cases(sca.chrom.map),]
head(sca.chrom.map$best.ZFinch.chrom)
unique(sca.chrom.map$best.ZFinch.chrom)
# The next lines hack together a way of ordering the names of the Zebra Finch chromosomes such that
# ones of "unknown" position ("un") come after the known ones (e.g., chr_1 vs chr_1_un)
# and the Z and LGE22 chromosomes are last
o = substring(regmatches(sca.chrom.map$best.ZFinch.chrom,regexpr("_.*",sca.chrom.map$best.ZFinch.chrom)),2)
o = gsub("A",".1",o)
o = gsub("B",".12",o)
o = gsub("_random","101",o)
o = gsub("Z","98",o)
o = gsub("LGE22","99",o)
o = gsub("Un","100",o)
o = as.numeric(o)

sca.chrom.map = sca.chrom.map[order(o,sca.chrom.map$mean.loc.ZFinch.chrom),]
write.csv(sca.chrom.map,"GCRF_scaffold_chrom_map_ZFinch_Un.csv",row.names=F)


#------------------------------
#STEP 4: ORDERING SCA.CHROM.MAP & SAVING OUTPUT FILE

# measuring length of successfully mapped scaffolds and comparing it to the total length of the included scaffolds
sum(ZLat.scafs$length[match(as.character(sca.chrom.map$ZLat.scaffold),ZLat.scafs$scaffold)])/sum(ZLat.scafs$length)
# For me, this value was 0.9447883 <- proportion mapped to Zebra Finch Genome
head(sca.chrom.map)
# Now, you may want to output the "sca.chrom.map" object for further use. 
# I pared it down further to just the three columns I wanted before saving it: 
setwd("~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/chromosome/")
scaffold.order = sca.chrom.map[,c("best.ZFinch.chrom","ZLat.scaffold","mean.loc.ZFinch.chrom","ori")]  
colnames(scaffold.order) = c("chr","sca","mean.loc","sca.ori")
head(scaffold.order)

unique(scaffold.order$chr)
scaffold.order$chr[scaffold.order$chr == '1A'] <- '1.1'
#scaffold.order$chr[scaffold.order$chr == '1B'] <- '1.2'
scaffold.order$chr[scaffold.order$chr == '4A'] <- '4.1'
scaffold.order$chr[scaffold.order$chr == 'Z'] <- '29'
scaffold.order$chr[scaffold.order$chr == 'LGE22'] <- '30'
#scaffold.order$chr[scaffold.order$chr == '3_random'] <- '31.1'
#scaffold.order$chr[scaffold.order$chr == '4_random'] <- '31.2'
scaffold.order$chr[scaffold.order$chr == '8_random'] <- '31.3'
scaffold.order$chr[scaffold.order$chr == '13_random'] <- '31.4'
scaffold.order$chr[scaffold.order$chr == '25_random'] <- '31.5'
scaffold.order$chr[scaffold.order$chr == '26_random'] <- '31.6'
scaffold.order$chr[scaffold.order$chr == 'Un'] <- '32'
unique(scaffold.order$chr)

write.csv(scaffold.order,"GCRF_scaffold_order_from_ZFinch_Un.csv",row.names=F)


