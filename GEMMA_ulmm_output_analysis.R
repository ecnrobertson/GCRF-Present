library(tidyverse)
library(ggplot2)
library(geiger)
library(ggrepel)
library(ggpubr)
library(data.table)
setwd('~/Desktop/GCRF/GCRF-Present/GenomicPipeline/data/ulmm/')

# >>>> WING <<<< 1754 SNP @<6.25e-9
################################# 

width <-read_delim('wing_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<6.25e-9,"outlier","neutral")) 
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
# >>>> TARSUS <<<< 0 SNP @<6.25e-9 0 SNP @5e-8
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
# >>>> TAIL LENGTH <<<< 1 SNP @<6.25e-9 13 SNP @5e-8
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
  labs(x="Scaffold",y="-log10(p_wald)")

png("GCRF.tail_length_lmm.assoc.png",width=11.7,height=8.3,units="in", res=200)
ggarrange(g2)
dev.off()


################################# 
# >>>> NARE LENGTH <<<< 0 SNP @<6.25e-9 2 SNP @5e-8
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
# >>>> BRISTLE LENGTH <<<< 0 SNP @<6.25e-9 0 SNP @5e-8
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
# >>>> BEAK WIDTH <<<< 1 SNP @<6.25e-9 2 SNP @5e-8
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
# >>>> BEAK DEPTH <<<< 0 SNP @<6.25e-9 0 SNP @5e-8
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
# >>>> CULMEN-END LENGTH <<<< 0 SNP @<6.25e-9 5 SNP @5e-8
################################
width <-read_delim('culmen_end_length_lmm.assoc.txt',delim="\t") %>%  mutate(outlier=if_else(p_wald<5e-8,"outlier","neutral")) 
c <- width %>%  mutate(win_num = 1:nrow(width)) 
dim(c)
out<-c %>% filter(outlier=="outlier")
dim(out)
cg_fl<-c %>% filter(label!="neutral")
c %>% distinct(chr)
cbPalette <- c("#999999", "#56B4E9")
c %>% distinct(chr) %>% tally()
c %>% group_by(chr) %>% tally() %>% write.table("GCRF_culmen_end_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")
out %>% write.table("GCRF_culmen_end_lmm_5e-8.outlier_tally.txt",row.names=F,quote=F,sep="\t")

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
