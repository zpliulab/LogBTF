## 20221128
## 2023.2.25 Newly add GNIPLR method
## Visualize results of the simulated single-cell and bulk RNA-seq dataset.


rm(list = ls())

# Color Setting -----------------------------------------------------------
colors8 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#c6cdd7", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#c6cdd7", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
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



# Result 2 : on Matrix data --------------------------------------------------------


library(dplyr)
library(ggplot2)
library(stringr)
library(ggprism)


# AUROC -------------------------------------------------------------------
setwd("/Users/lilingyu/E/PhD/R/Boolean/Result")  
# df <- readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="AUROC")
df <- readxl::read_excel("SERGIOSimu_scResult_Revise6583CleanSCODE.xlsx", sheet="AUROC")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)


pAUROC <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F) +
  scale_fill_manual(values = colors8) +
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.4, 0.65),
                     breaks = seq(0.4, 0.65, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.8,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("AUROC")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=9,y=0.28,yend=0.28)+
  annotate("text",x=3.5, y=0.1, label="Node 20", vjust=0.4)+
  geom_segment(x=10,xend=18,y=0.28,yend=0.28)+
  annotate("text",x=9,  y=-0.1, label="Node 100",vjust=0.4)+
  # geom_segment(x=19,xend=27,y=-0.15,yend=-0.15)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")


pAUROC
LogBTF <- list()
LogBTF[[1]] <- pAUROC
# ggsave(pAUROC, filename = 'FigurescBioRevised/AUROC_SERGIO_Matrix.pdf', width = 8, height = 6, device = cairo_pdf)


# StAcc -------------------------------------------------------------------
## SCODE 效果最好
# df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="StAcc")
df <- readxl::read_excel("SERGIOSimu_scResult_Revise6583CleanSCODE.xlsx", sheet="StAcc")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pStAcc <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors7)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.7, 1.0),
                     breaks = seq(0.7, 1.0, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.8,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("StAcc")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +


  geom_segment(x=1,xend=8,y=0.54,yend=0.54)+
  geom_segment(x=9,xend=16,y=0.54,yend=0.54)+
  # geom_segment(x=15,xend=21,y=-0.24,yend=-0.24)+
  coord_cartesian(clip = "off")

pStAcc
LogBTF[[2]] <- pStAcc
# ggsave(pStAcc, filename = 'Figuresc/StAcc.pdf', width = 8, height = 6, device = cairo_pdf)


# AUPRsc10 -------------------------------------------------------------------
# df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583CleanSCODE.xlsx", sheet="AUPRsc")
df <- readxl::read_excel("SERGIOSimu_scResult_Revise6583CleanSCODE.xlsx", sheet="AUPRsc")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)


pAUPR10 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors4)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.0, 0.3),
                     breaks = seq(0.0, 0.3, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.8,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("AUPR")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=4,y=-0.17,yend=-0.17)+
  geom_segment(x=5,xend=8,y=-0.17,yend=-0.17)+
  # geom_segment(x=9,xend=12,y=-0.07,yend=-0.07)+
  coord_cartesian(clip = "off")

pAUPR10
LogBTF[[3]] <- pAUPR10
# ggsave(pAUPR10, filename = 'Figuresc/pAUPR10.pdf', width = 3, height = 6, device = cairo_pdf)



# Pre1 ---------------------------------------------------------------------

# df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="Pre")
df <- readxl::read_excel("SERGIOSimu_scResult_Revise6583CleanSCODE.xlsx", sheet="Pre")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)


pPre10 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = colors7)+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.0, 0.3),
                     breaks = seq(0.0, 0.3, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.8,
              base_fontface = "plain",
              # base_family = "serif",
  ) +
  labs(x=NULL, y=expression("Pre")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=8,y=-0.17,yend=-0.17)+
  geom_segment(x=9,xend=16,y=-0.17,yend=-0.17)+
  # geom_segment(x=9,xend=12,y=-0.07,yend=-0.07)+
  coord_cartesian(clip = "off")

pPre10
LogBTF[[5]] <- pPre10
# ggsave(pPre10, filename = 'Figuresc/Pre10.pdf', width = 8, height = 6, device = cairo_pdf)


# AUROC(SIGN=1) ---------------------------------------------------------------------
df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="AUROC(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUROC1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.45, 0.55),
                     breaks = seq(0.45, 0.55, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.8,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
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
df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="StAcc(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pStAcc1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.70, 1.00),
                     breaks = seq(0.70, 1.00, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("StAcc (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=0.39, yend=0.39)+
  geom_segment(x=3,xend=4,y=0.39, yend=0.39)+
  geom_segment(x=5,xend=6,y=0.39, yend=0.39)+
  coord_cartesian(clip = "off")

pStAcc1
LogBTF[[8]] <- pStAcc1
# ggsave(pStAcc1, filename = 'Figuresc/StAcc(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)




# AUPR(SIGN=1) ---------------------------------------------------------------------
df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="AUPR(SIGN=1)")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pAUPR1 <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#cab2d6", "#b2df8a"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, by=0.05))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("AUPR (SIGN=1)")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=-0.06, yend=-0.06)+
  geom_segment(x=3,xend=4,y=-0.06, yend=-0.06)+
  # geom_segment(x=5,xend=6,y=-0.06, yend=-0.06)+
  coord_cartesian(clip = "off")

pAUPR1
LogBTF[[10]] <- pAUPR1
# ggsave(pAUPR1, filename = 'Figuresc/pAUPR(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)



# DyAcc ---------------------------------------------------------------------
df<-readxl::read_excel("SERGIOSimu_scResult_Revise6583Clean.xlsx", sheet="DyAcc")
x_level <- paste(df$Group1,df$Group2,sep="_")
df1 <- mydfdata(df)

pLogBTF <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#9F8CFF","#F266E2"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.45, 0.85),
                     breaks = seq(0.45, 0.85, by=0.1))+
  theme_prism(axis_text_angle =0,
              base_line_size = 0.1,
              base_fontface = "plain",
              # base_family = "serif",
              ) +
  labs(x=NULL, y=expression("AUC and DyAcc")) +
  
  geom_segment(x=1,xend=2,y=0.38, yend=0.38)+
  geom_segment(x=3,xend=4,y=0.38, yend=0.38)+
  # geom_segment(x=5,xend=6,y=0.38, yend=0.38)+
  coord_cartesian(clip = "off")

pLogBTF
LogBTF[[11]] <- pLogBTF
# ggsave(pLogBTF, filename = 'Figuresc/DyAcc.pdf', width = 2.5, height = 4, device = cairo_pdf)



# Picture  -- Matrix,  Clear -- dont seem good ----------------------------
# pdf("FigurescBioRevised/scSERGIO.pdf",width = 18, height = 6.5)
cowplot::plot_grid(plotlist = list(LogBTF[[1]],LogBTF[[3]],LogBTF[[5]]), 
                   nrow = 1, 
                   rel_widths = c(1.4, 0.7, 1.3),
                   labels = c('(a)','(b)','(c)'),
                   scale = c(0.99))
# dev.off()

