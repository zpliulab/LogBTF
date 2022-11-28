
#################################################################################################
#################################################################################################


## 20221128 Box_GroupPlotSimple
## Visualize results of the simulated single-cell and bulk RNA-seq dataset.


rm(list = ls())

# Color Setting -----------------------------------------------------------
colors8 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors6 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#cab2d6", "#b2df8a", "#FEC071")
colors4 <- c("#29aeff", "#cab2d6", "#ff77a8", "#b2df8a")



# function def ------------------------------------------------------------
mydfdata <- function(df){
  df %>% mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
    select(-c("Group1","Group2")) %>% 
    reshape2::melt(var.ids="new_col") -> df1
  
  df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
  df1$new_col <- factor(df1$new_col,
                        levels = x_level)
  return(df1)
}


# Result 2 : on single-cell data --------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(ggprism)


# AUROC -------------------------------------------------------------------
setwd("/Users/lilingyu/E/PhD/R/Boolean/Result")  
df <- readxl::read_excel("DREAM3scResult.xlsx", sheet="AUROC")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)


pAUROC <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors8)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.1, 0.70),
                     breaks = seq(0.1, 0.70, by=0.1))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUROC")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=8,y=-0.15,yend=-0.15)+
  annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=9,xend=16,y=-0.15,yend=-0.15)+
  annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=17,xend=24,y=-0.15,yend=-0.15)+
  annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")


pAUROC
LogBTF <- list()
LogBTF[[1]] <- pAUROC
# ggsave(pAUROC, filename = 'Figuresc/AUROC.pdf', width = 8, height = 6, device = cairo_pdf)


# StAcc -------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pStAcc <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors7)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.1, 0.9),
                     breaks = seq(0.1, 0.9, by=0.1))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("StAcc")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  
  geom_segment(x=1,xend=7,y=-0.24,yend=-0.24)+
  geom_segment(x=8,xend=14,y=-0.24,yend=-0.24)+
  geom_segment(x=15,xend=21,y=-0.24,yend=-0.24)+
  coord_cartesian(clip = "off")

pStAcc
LogBTF[[2]] <- pStAcc
# ggsave(pStAcc, filename = 'Figuresc/StAcc.pdf', width = 8, height = 6, device = cairo_pdf)


# AUPRsc10 -------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="AUPRsc10")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)


pAUPR10 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors4)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.1, 0.5),
                     breaks = seq(0.1, 0.5, by=0.1))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUPR")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=4,y=-0.07,yend=-0.07)+
  geom_segment(x=5,xend=8,y=-0.07,yend=-0.07)+
  geom_segment(x=9,xend=12,y=-0.07,yend=-0.07)+
  coord_cartesian(clip = "off")

pAUPR10
LogBTF[[3]] <- pAUPR10
# ggsave(pAUPR10, filename = 'Figuresc/pAUPR10.pdf', width = 3, height = 6, device = cairo_pdf)



# AUPRsc -------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="AUPRsc50100")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUPR50 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors4)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.01, 0.09),
                     breaks = seq(0.01, 0.09, by=0.02),
                     position = "right")+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUPR")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=4,y=-0.024,yend=-0.024)+
  geom_segment(x=5,xend=8,y=-0.024,yend=-0.024)+
  coord_cartesian(clip = "off")

pAUPR50
LogBTF[[4]] <- pAUPR50
# ggsave(pAUPR50, filename = 'Figuresc/pAUPR50.pdf', width = 6, height = 6, device = cairo_pdf)


# Pre10 ---------------------------------------------------------------------

df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="Pre")[1:7,]
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)
pPre10 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors7)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.0, 0.35),
                     breaks = seq(0.0, 0.35, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("Pre")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=7,y=-0.16,yend=-0.16)+
  coord_cartesian(clip = "off")

pPre10
LogBTF[[5]] <- pPre10
# ggsave(pPre10, filename = 'Figuresc/Pre10.pdf', width = 8, height = 6, device = cairo_pdf)


# Pre 50---------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="Pre")[8:22,]
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pPre50 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors7)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.00, 0.1),
                     breaks = seq(0.00, 0.1, by=0.02),
                     position = "right")+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("Pre")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=7,y=-0.046,yend=-0.046)+
  geom_segment(x=8,xend=14,y=-0.046,yend=-0.046)+
  coord_cartesian(clip = "off")

pPre50
LogBTF[[6]] <- pPre50
# ggsave(pPre10, filename = 'Figuresc/Pre50.pdf', width = 8, height = 6, device = cairo_pdf)


# AUROC(SIGN=1) ---------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="AUROC(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUROC1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.25, 0.7),
                     breaks = seq(0.25, 0.7, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUROC (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=0.182, yend=0.182)+
  geom_segment(x=3,xend=4,y=0.182, yend=0.182)+
  geom_segment(x=5,xend=6,y=0.182, yend=0.182)+
  coord_cartesian(clip = "off")

pAUROC1
LogBTF[[7]] <- pAUROC1
# ggsave(pAUROC1, filename = 'Figuresc/AUROC(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)


# StAcc(SIGN=1) ---------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pStAcc1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.45, 0.85),
                     breaks = seq(0.45, 0.85, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("StAcc (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=0.39, yend=0.39)+
  geom_segment(x=3,xend=4,y=0.39, yend=0.39)+
  geom_segment(x=5,xend=6,y=0.39, yend=0.39)+
  coord_cartesian(clip = "off")

pStAcc1
LogBTF[[8]] <- pStAcc1
# ggsave(pStAcc1, filename = 'Figuresc/StAcc(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)


# AUROC+StAcc -------------------------------------------------------------
df2<-readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc(SIGN=1)")
df <- cbind(df1[,c(1:2)], df1[,c(3:7)] + df2[,c(3:7)])
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUROC_StAcc <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.8, 1.40),
                     breaks = seq(0.8, 1.40, by=0.1))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUROC + StAcc (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=0.7, yend=0.7)+
  geom_segment(x=3,xend=4,y=0.7, yend=0.7)+
  geom_segment(x=5,xend=6,y=0.7, yend=0.7)+
  coord_cartesian(clip = "off")

pAUROC_StAcc
LogBTF[[9]] <- pAUROC_StAcc
# ggsave(pAUROC_StAcc, filename = 'Figuresc/AUROC_StAcc(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)


# AUPR(SIGN=1) ---------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="AUPR(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUPR1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 0.35, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUPR (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=-0.06, yend=-0.06)+
  geom_segment(x=3,xend=4,y=-0.06, yend=-0.06)+
  geom_segment(x=5,xend=6,y=-0.06, yend=-0.06)+
  coord_cartesian(clip = "off")

pAUPR1
LogBTF[[10]] <- pAUPR1
# ggsave(pAUPR1, filename = 'Figuresc/pAUPR(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)



# DyAcc ---------------------------------------------------------------------
df<-readxl::read_excel("DREAM3scResult.xlsx", sheet="DyAcc")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pLogBTF <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#9F8CFF","#F266E2"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.84, 1.00),
                     breaks = seq(0.84, 1.00, by=0.02))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUC and DyAcc")) +
  
  geom_segment(x=1,xend=2,y=0.18, yend=0.18)+
  geom_segment(x=3,xend=4,y=0.18, yend=0.18)+
  geom_segment(x=5,xend=6,y=0.18, yend=0.18)+
  coord_cartesian(clip = "off")

pLogBTF
LogBTF[[11]] <- pLogBTF
# ggsave(pLogBTF, filename = 'Figuresc/DyAcc.pdf', width = 2.5, height = 4, device = cairo_pdf)




# Plot --------------------------------------------------------------------

# pdf("Figuresc/Pre.pdf",width = 8, height = 6)
cowplot::plot_grid(plotlist = LogBTF, nrow = 3)
# dev.off()


## AUC + DyAcc, AUROC, AUPR , SIGN = 1    
# pdf("Figuresc/sc2Method.pdf",width = 16, height = 5)
cowplot::plot_grid(plotlist = list(LogBTF[[11]], LogBTF[[7]], LogBTF[[10]] ), 
                   nrow = 1, rel_widths = c(1.0, 2.3, 2.3))
# dev.off()

## AUROC, AUPR, PRE
## save 18*6
# pdf("Figuresc/sc8Method.pdf",width = 18, height = 6)
cowplot::plot_grid(plotlist = list(LogBTF[[1]], LogBTF[[3]], LogBTF[[4]], LogBTF[[5]], LogBTF[[6]]), 
                   nrow = 1, rel_widths = c(2 , 0.55, 0.9, 0.8, 1.5))
# dev.off()


# StAcc just for LogBTF ---------------------------------------------------

# StAcc(SIGN=1) no use in our paper !!! ---------------------------------------------------------------------
data <-readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc(SIGN=1)")[c(1,3,5),-1]

data$median<- c(mean(as.matrix(data[1,c(2:6)])), mean(as.matrix(data[2,c(2:6)])), mean(as.matrix(data[3,c(2:6)])))
data$lower<- c(min(data[1,c(2:6)]), min(data[2,c(2:6)]),min(data[3,c(2:6)]))
data$upper<- c(max(data[1,c(2:6)]), max(data[2,c(2:6)]),max(data[3,c(2:6)]))
data$Group2 <- c("Size 10", "Size 50", "Size 600")

pStAcc <- ggplot(data,aes(factor(Group2),median))+
  geom_col(aes(fill=Group2),position = 'dodge',width = 0.8) +
  geom_errorbar(aes(factor(Group2),group=Group2,ymin=lower,ymax=upper),
                position = position_dodge(width = 0.8),
                width=0,cex=1) +
  scale_y_continuous(limits = c(0,0.9),
                     breaks = seq(0,0.9,0.1),
                     expand = c(0,0)) +
  labs(x='Network size',y='StAcc (SIGN=1)')+
  theme_test(base_size = 15)+
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'))

# pdf("Figuresc/StAccLogBTF.pdf",width = 4.5, height = 4.5)
pStAcc
# dev.off()




# StAcc for SIGN0 and SIGN1 -----------------------------------------------

data <-readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc(SIGN=1)")[c(1,3,5),-1]
data2 <- readxl::read_excel("DREAM3scResult.xlsx", sheet="StAcc")[c(1,8,15),-1]


data$median<- c(mean(as.matrix(data[1,c(2:6)])), mean(as.matrix(data[2,c(2:6)])), mean(as.matrix(data[3,c(2:6)])))
data$lower<- c(min(data[1,c(2:6)]), min(data[2,c(2:6)]),min(data[3,c(2:6)]))
data$upper<- c(max(data[1,c(2:6)]), max(data[2,c(2:6)]),max(data[3,c(2:6)]))
data$Group2 <- c("Size 10", "Size 50", "Size 600")


data2$median<- c(mean(as.matrix(data2[1,c(2:6)])), mean(as.matrix(data2[2,c(2:6)])), mean(as.matrix(data2[3,c(2:6)])))
data2$lower<- c(min(data2[1,c(2:6)]), min(data2[2,c(2:6)]),min(data2[3,c(2:6)]))
data2$upper<- c(max(data2[1,c(2:6)]), max(data2[2,c(2:6)]),max(data2[3,c(2:6)]))
data2$Group2 <- c("Size 10", "Size 50", "Size 600")


data12 <- rbind(data[1,], data2[1,],
                data[2,], data2[2,],
                data[3,], data2[3,])

data12$SIGN <- c("SIGN=1", "SIGN=0", 
                 "SIGN=1", "SIGN=0", 
                 "SIGN=1", "SIGN=0")

pStAcc01 <- ggplot(data12,aes(factor(Group2),median))+
  geom_col(aes(fill=SIGN),position = 'dodge',width = 0.8) +
  geom_errorbar(aes(factor(Group2),group=SIGN,ymin=lower,ymax=upper),
                position = position_dodge(width = 0.8),
                width=0,cex=1) +
  scale_y_continuous(limits = c(0,0.9),
                     breaks = seq(0,0.9,0.1),
                     expand = c(0,0)) +
  scale_fill_manual(values = c('#fec79e','#8ec4cb'))+
  labs(x='Network size',y='StAcc')+
  theme_test(base_size = 15) +
  theme(#legend.position = 'none',
    panel.border = element_rect(size=2,fill = 'transparent'),
    axis.text = element_text(color='black')) #+

# pdf("Figuresc/StAccLogBTF01.pdf",width = 6, height = 4.5)
pStAcc01
# dev.off()





#################################################################################################
#################################################################################################


## 202211028 Bar_BoublePlotSimple
## plot AUROC and Runtimes on Real data 2 and 3.


rm(list = ls())

colors8 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7ARA <- c("#4daf4a", "#4dbbd5", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
variable <- rep("method", 7)

library(tidyverse)  
library(reshape2) 

# Real Data1 --------------------------------------------------------------
setwd("/Users/lilingyu/E/PhD/R/Boolean/Result")  
data <- as.data.frame(read.csv("RealscResult.csv", header = T, sep = ","))
data$Method <- factor(data$Method)

p3 <- ggplot(data)+
  geom_col(aes(Method,AUROC,fill=Method),
           position = 'dodge',width = 0.6)+
  scale_fill_manual(values = colors7ARA)+
  labs(x='Method',y='AUROC') +
  theme_test(base_size = 15)+
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, 
                                   hjust = 0.5, vjust = 0.5, size = 11)) +
  scale_y_continuous(limits = c(0.0,0.6),
                     breaks = seq(0.0, 0.6, 0.1),
                     expand = c(0,0),
                     sec.axis = sec_axis(~./(0.6/160),
                                         name = 'Runtime',
                                         breaks = seq(0,160,40))) +
  geom_point(data = data,
             aes(Method, Runtime*(0.6/160)+0.01,
                 color=variable),
             size=4) +
  geom_line(data = data,
            aes(Method, Runtime*(0.6/160)+0.01,
                color=variable, group=variable),
            cex=1.3)+
  scale_color_manual(values = '#1e8b9b') 

# pdf('Figuresc/RealdataAUROC.pdf', width = 5.5, height = 4.5)
p3 
# dev.off()


# Real Data3 --------------------------------------------------------------
setwd("/Users/lilingyu/E/PhD/R/Boolean/Result")  
data <- as.data.frame(read.csv("RealscResult2.csv", header = T, sep = ","))
data$Method <- factor(data$Method)

p4 <- ggplot(data)+
  geom_col(aes(Method,AUROC,fill=Method),
           position = 'dodge',width = 0.6)+
  scale_fill_manual(values = colors7ARA)+
  labs(x='Method',y='AUROC') +
  theme_test(base_size = 15) +
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, 
                                   hjust = 0.5, vjust = 0.5, size = 11)) +
  scale_y_continuous(limits = c(0.0,0.6),
                     breaks = seq(0.0, 0.6, 0.1),
                     expand = c(0,0),
                     sec.axis = sec_axis(~./(0.6/2550),
                                         name = 'Runtime',
                                         breaks = seq(0,2550,500))) +
  geom_point(data = data,
             aes(Method, Runtime*(0.6/2550)+0.01,
                 color=variable),
             size=4) +
  geom_line(data = data,
            aes(Method, Runtime*(0.6/2550)+0.01,
                color=variable, group=variable),
            cex=1.3)+
  scale_color_manual(values = '#1e8b9b') 

# pdf('Figuresc/Realdata2AUROC.pdf', width = 5.5, height = 4.5)
p4
# dev.off()


