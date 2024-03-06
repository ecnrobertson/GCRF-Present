#Following tutorial put up online https://visoca.github.io/popgenomworkshop-gwas_gemma/

#First, we are going to play with the hyperparameter estimates, which will inform us of the genetic architecture of the trait.

setwd("GenomicPipeline/data/bslmm")

# Load hyperparameter file
# ==============================================================================
hyp.params<-read.table("GCRF_bill_PC.bslmm.hyp.txt",header=T)
# ==============================================================================

# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
hyp.params.table
# write table to file
write.table(hyp.params.table, file="hyperparameters.dsv", sep="\t", quote=F)
# ==============================================================================



# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file="hyperparameters.pdf", width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# No gamma
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$pi), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------
dev.off()
# ==============================================================================

#Okay so these graphs don't look too great, so might need to rerun the bslmm and see if I get "better mixing".

#pulling out the numbers from the scaffolds so I have a numeric to work with for sorting and whatnot.
# cat GCRF_bill_PC.bslmm.param.txt | cut -d "_" -f 2 > scaffold_nums.txt
# sed '1d' scaffold_nums.txt > scaffold_nums_clean_billPC.txt
# head scaffold_nums_clean_billPC.txt 



# library to speed up loading of big tables

setwd("GenomicPipeline/data/bslmm")

library(data.table)

#nums <- read.table("scaffold_nums_clean_billPC.txt")
#olnames(nums) <- "chr"

# Load parameters output
# ==============================================================================
params<-fread("GCRF_billPC_clean.chr1-31ZF.genotype.params.txt",header=T,sep="\t", data.table=F)
# ==============================================================================

#params <- data.frame(nums, params[,-1])

# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
#26
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
#252
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
#2516
# ------------------------------------------------------------------------------

# write tables
#write.table(top1, file="top1eff.dsv", quote=F, row.names=F, sep="\t")
#write.table(top01, file="top0.1eff.dsv", quote=F, row.names=F, sep="\t")
#write.table(top001, file="top0.01eff.dsv", quote=F, row.names=F, sep="\t")

# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01.bill<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10.bill<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25.bill<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50.bill<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
#write.table(pip01, file="pip01.dsv", quote=F, row.names=F, sep="\t")
#write.table(pip10, file="pip10.dsv", quote=F, row.names=F, sep="\t")
#write.table(pip25, file="pip25.dsv", quote=F, row.names=F, sep="\t")
#write.table(pip50, file="pip50.dsv", quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------

# ==============================================================================
# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================
# Prepare data
# ------------------------------------------------------------------------------

# add linkage group column (chr)
chr<-params$ZFCHROM

# sort by linkage group and position
params.sort<-params[order(as.numeric(params$ZFCHROM), params$ZFPOS),]

# get list of linkage groups/chromosomes
chrs<-sort(unique(as.numeric(chr)))
# ------------------------------------------------------------------------------

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file="../bslmm/new_pip_plot_bill0.01_chr.png", width=11.7,height=8.3,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,1),ylab="PIP",xlab="Chromosome", xaxt="n")

# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$ZFCHROM==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"lightgrey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}
# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP 
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/abs(log(z))
# plot
symbols(x,y,circles=z, bg="grey14",inches=1/5, fg=NULL,add=T)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h=0.01,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.01],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.01]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.01]
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="tomato",inches=1/5,fg=NULL,add=T)
# ------------------------------------------------------------------------------

# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.01], adj=c(0,0), cex=0.8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# close device
dev.off()
# ==============================================================================