## 2022.8.31 将DREAM3的结果做成箱线图
## https://www.jianshu.com/p/834240e559ec
## https://www.nature.com/articles/s41588-021-00831-0



rm(list = ls())


setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="AUROC")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                    levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#e64b35", "#4daf4a", "#4dbbd5", 
                               "#cab2d6", "#b2df8a", "#FEC071"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.1, 0.70),
                     breaks = seq(0.1, 0.70, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("AUROC")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.3, 0.2),'cm')) +
  
  geom_segment(x=1,xend=6,y=-0.15,yend=-0.15)+
  annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=7,xend=12,y=-0.15,yend=-0.15)+
  annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=13,xend=18,y=-0.15,yend=-0.15)+
  annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\AUROC.pdf', width = 5.5, height = 5, device = cairo_pdf)



# StAcc -------------------------------------------------------------------

rm(list = ls())

setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="StAcc")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                      levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#e64b35", "#4daf4a", "#4dbbd5", 
                               "#cab2d6", "#b2df8a", "#FEC071"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(0.1, 0.9),
                     breaks = seq(0.1, 0.9, by=0.08))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("StAcc")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.3, 0.2),'cm')) +
  
  geom_segment(x=1,xend=6,y=-0.24,yend=-0.24)+
  # annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=7,xend=12,y=-0.24,yend=-0.24)+
  # annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=13,xend=18,y=-0.24,yend=-0.24)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  # annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  # annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\StAcc.pdf', width = 5.5, height = 5, device = cairo_pdf)



# Pre ---------------------------------------------------------------------

rm(list = ls())

setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="Pre")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                      levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1, aes(fill=group), show.legend = F)+
  scale_fill_manual(values = c("#e64b35", "#4daf4a", "#4dbbd5", 
                               "#cab2d6", "#b2df8a", "#FEC071"))+
  scale_x_discrete(labels = str_split_fixed(x_level, "_", 2)[,1],
                   guide = "prism_offset")+
  scale_y_continuous(limits = c(-0.05, 0.35),
                     breaks = seq(-0.05, 0.35, by=0.05))+
  theme_prism(axis_text_angle = 45,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif") +
  labs(x=NULL, y=expression("Pre")) +
  theme(plot.margin = unit(c(0.2 ,0.2, 1.3, 0.2),'cm')) +
  
  geom_segment(x=1,xend=6,y=-0.22,yend=-0.22)+
  # annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=7,xend=12,y=-0.22,yend=-0.22)+
  # annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=13,xend=18,y=-0.22,yend=-0.22)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  # annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  # annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\Pre.pdf', width = 5.5, height = 5, device = cairo_pdf)


# AUROC(SIGN=1) ---------------------------------------------------------------------

rm(list = ls())

setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="AUROC(SIGN=1)")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                      levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
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
  # annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=3,xend=4,y=0.182, yend=0.182)+
  # annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=5,xend=6,y=0.182, yend=0.182)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  # annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  # annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\AUROC(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)


# StAcc(SIGN=1) ---------------------------------------------------------------------

rm(list = ls())

setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="StAcc(SIGN=1)")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                      levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
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
  # annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=3,xend=4,y=0.39, yend=0.39)+
  # annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=5,xend=6,y=0.39, yend=0.39)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  # annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  # annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\StAcc(SIGN=1).pdf', width = 7, height = 5, device = cairo_pdf)


# DyAcc ---------------------------------------------------------------------

rm(list = ls())

setwd("D:/E/博士/R_程序/Boolean/Result")
## AUROC    StAcc
df<-readxl::read_excel("DREAM3result.xlsx", sheet="DyAcc")


## 宽格式数据，如果使用ggplot2作图需要转换成长格式
library(dplyr)
df %>% 
  mutate(new_col=paste(Group1,Group2,sep="_")) %>% 
  select(-c("Group1","Group2")) %>% 
  reshape2::melt(var.ids="new_col") -> df1

head(df1)


## ggplot2 作图

library(ggplot2)
library(stringr)
# install.packages('ggprism')
library(ggprism)

x_level <- paste(df$Group1,df$Group2,sep="_")
x_level
# df1$group <- str_sub(df1$new_col,5,7)
df1$group <- str_split_fixed(df1$new_col, "_", 2)[,1]
df1$new_col <- factor(df1$new_col,
                      levels = x_level)

# paranr = length(x_level) 
p <- ggplot(df1, aes(x=new_col, y=value))+
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
  # theme(plot.margin = unit(c(0.2 ,0.2, 4, 0.2),'cm')) +
  
  geom_segment(x=1,xend=2,y=0.18, yend=0.18)+
  # annotate("text",x=3.5, y=-0.1, label="Node 10", vjust=0.4)+
  geom_segment(x=3,xend=4,y=0.18, yend=0.18)+
  # annotate("text",x=9,  y=-0.1, label="Node 50",vjust=0.4)+
  geom_segment(x=5,xend=6,y=0.18, yend=0.18)+
  # annotate("text",x=16.5,y=-0.1, label="Node 100",vjust=0.4)+
  # annotate("text",x=1,y=-0.1,label="P1", hjust=2, vjust=0.2)+
  # annotate("text",x=1,y=-0.1,label="P2", hjust=2, vjust=0.4)+
  coord_cartesian(clip = "off")

p
ggsave(p, filename = 'Figure\\DyAcc.pdf', width = 2.5, height = 4, device = cairo_pdf)
