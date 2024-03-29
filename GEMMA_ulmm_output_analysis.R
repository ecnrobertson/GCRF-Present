library(tidyverse)
library(ggplot2)
library(geiger)
library(ggrepel)
library(ggpubr)
library(data.table)
setwd('~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/ulmm/')

# >>>> WING <<<< 0 SNP @<6.88e-9 0 SNP @<5e-8
################################# 

width <-read_delim('wing_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_wing_lmm.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_wing_lmm.outlier_tally.txt",row.names=F,quote=F,sep="\t")


# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()



################################# 
# >>>> TARSUS <<<< 0 SNP @<6.88e-9 0 SNP @5e-8
#################################
width <-read_delim('tarsus_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
#c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_tarsus_lmm.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_tarsus_lmm.outlier_tally.txt",row.names=F,quote=F,sep="\t")


# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()


################################# 
# >>>> TAIL LENGTH <<<< 0 SNP @<6.88e-9 1 SNP @5e-8
################################
width <-read_delim('tail_length_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_tail_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_tail_lmm_5e-8.outlier_tally.txt", row.names=F,quote=F,sep="\t")

#after reorganizing at chromosome level
width <-read_delim('GCRF_tail_ulmm_snps_ZFchr_assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
out %>% write.table("GCRF_tail_lmm_5e-8.outlier_tally.txt", row.names=F,quote=F,sep="\t")

# Make the plot
y_breaks <- pretty(range(-log10(c$p_wald)), n=5)

library(qqman)
test.plot <- manhattan(c, chr="ZFCHROM", bp="ZFPOS", p="p_wald", snp="rs")

g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald), color=as.factor(ZFCHROM))) +
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("grey80","grey40", "grey20"), 200 )) +
  # custom X axis:
  scale_x_continuous(label = c$ZFCHROM) +
  scale_y_continuous(breaks=y_breaks, expand = c(0, 0) ) +     # remove space between plot area and x axis
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="tomato2", size=0.9) +
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_wald)")

png("GCRF.tail_length_lmm.assoc.png",width=11.7,height=10.3,units="in", res=200)
ggarrange(g2)
dev.off()


################################# 
# >>>> NARE LENGTH <<<< 0 SNP @<6.88e-9 0 SNP @5e-8
################################
width <-read_delim('nare_length_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_nare_length_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_nare_length_lmm_5e-8.outlier_tally.txt", row.names=F,quote=F,sep="\t")

# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()


################################# 
# >>>> BRISTLE LENGTH <<<< 0 SNP @<6.88e-9 0 SNP @5e-8
################################
width <-read_delim('bristle_length_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_bristle_length_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")


# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()



################################# 
# >>>> BEAK WIDTH <<<< 0 SNP @<6.88e-9 0 SNP @5e-8
################################
width <-read_delim('beak_width_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_beak_width_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_beak_width_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")

# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()



################################# 
# >>>> BEAK DEPTH <<<< 0 SNP @<6.88e-9 0 SNP @5e-8
################################
width <-read_delim('beak_depth_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_beak_depth_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")


# Make the plot
g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(chr))) +
  
  # Show all points
  geom_point(alpha=0.8, size=0.9,shape=1) + 
  #geom_hline(yintercept = 0.0278) +
  scale_color_manual(values = rep(c("#6600CC","#9933FF", "#CC99FF"), 200 )) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(c, outlier=="outlier"), color="orange", size=0.9) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(x="Scaffold",y="-log10(p_lrt)")

png("GCRF.wing_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()



################################# 
# >>>> CULMEN-END LENGTH <<<< 0 SNP @<6.88e-9 1 SNP @5e-8
################################
width <-read_delim('culmen_end_length_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out2<-c %>% filter(outlier=="outlier")
dim(out2)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_culmen_end_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_culmen_end_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")

###### after rearranging
width <-read_delim("GCRF_culmen_end_length_ulmm_snps_ZFchr_assoc.txt",delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width[order(as.numeric(width$ZFCHROM)),]
c <- c %>% mutate(win_num = 1:nrow(c))
#c2 <- c[sample(nrow(c), size=200000, replace = F),]
#c2 <- c2[order(as.numeric(c2$ZFCHROM)),] %>% mutate(win_num = 1:nrow(c2))
c2 <- c[sample(nrow(c), size=200000, replace = F),]

dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
out %>% write.table("GCRF_culmen_lmm_5e-8.outlier_tally.txt", row.names=F,quote=F,sep="\t")

# Make the plot

frequency <- table(c$ZFCHROM)
start <- 1
lab.pos<-vector()
for (i in seq_along(frequency)){
  size <- frequency[i]
  txtpos <- start+size/2
  lab.pos <- c(lab.pos, txtpos)
  start <- start+size
}

specific_breaks <- lab.pos
#specific_labels <- c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", 
#                       "8", "9", "10", "11", "12", "13", "14",
#                       "15", "17", "18", "19", "20", "21",
#                       "22", "23", "24", "25", "26", "27", "28",
#                       "Z", "LGE22", "MT")
specific_labels <- c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", 
                     "8", "9", "10", "11", "12", "13", "14",
                     " ", "16", " ", " ", "20", " ",
                     " ", " ", " ", " ", " ", " ", "28",
                     "Z", "LGE22", "")

g2<-ggplot(c, aes(x=win_num, y=-log10(p_wald),color=as.factor(ZFCHROM))) +
  geom_point(alpha=1, size=0.9,shape=19) + 
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", alpha=.8) +
  scale_color_manual(values = rep(c("grey80","grey40", "grey20"), 200 )) +
  scale_x_continuous(breaks = specific_breaks, labels = specific_labels, ) +
  #scale_x_continuous(breaks = waiver())+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,9) ) +     # remove space between plot area and x axis
  geom_point(data=subset(c, outlier=="outlier"), color="tomato2", size=1.2) +
  geom_text(x=1950000, y=7.9, label="CTDSPL", color="black", angle=45, size=3)+
  #geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=label), size=2,segment.color = "black",segment.size = .3,min.segment.length = unit(0, 'lines'),nudge_y = .03) +
  theme_classic() +
  theme( 
    legend.position="none")
    #panel.border = element_blank()+
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank()) +
  labs(x="Chromosome",y="-log10(p-value)")

png("../../../figures/GCRF.culmen_lmm.assoc.newxaxis.short.label.png",width=11.7,height=5,units="in", res=200)
ggarrange(g2)
dev.off()

pdf("../../../figures/GCRF.culmen_lmm.assoc.final.pdf",width=11.7,height=5)
ggarrange(g2)
dev.off()