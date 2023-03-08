## 2022.8.31 ??DREAM3?Ľ???????????ͼ
## 2023.2.25 Newly add GNIPLR method
## plot AUROC and Runtimes on Real data 2.



rm(list = ls())


colors8 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#c6cdd7", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#c6cdd7", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7ARA <- c("#4daf4a", "#4dbbd5", "#c6cdd7", "#29aeff", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors7 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#cab2d6", "#ff77a8", "#b2df8a", "#FEC071")
colors6 <- c("#e64b35", "#4daf4a", "#4dbbd5", "#cab2d6", "#b2df8a", "#FEC071")
colors4 <- c("#29aeff", "#cab2d6", "#ff77a8", "#b2df8a")

library(tidyverse) #ggplot2包等
library(reshape2) #用到melt()函数将宽数据转换成长数据



# Real Data1 --------------------------------------------------------------


setwd("/Users/lilingyu/E/PhD/R/Boolean/Result")  
## AUROC    StAcc
data <- as.data.frame(read.csv("RealscResult_BioRevised.csv", header = T, sep = ","))


data[,c(1:3)]

data$Method <- factor(data$Method)
data$AUROC

ggplot(data,aes(Method, AUROC)) +
  geom_point(aes(fill=Method),position = 'dodge',width = 0.8) +
  # ylim(0.4, 0.6) +
  scale_y_continuous(limits = c(0.0,0.55),
                     breaks = seq(0.0, 0.55, 0.1)) +
  scale_fill_manual(values = colors7ARA) + 
  labs(x='Method',y='AUROC') +
  theme_test(base_size = 15) +
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black')) 



p1 <- ggplot(data,aes(Method, AUROC, color=variable,group=variable)) +
  geom_point(size=4) +
  geom_line(cex=1.3) +
  scale_y_continuous(limits = c(0.47,0.51),
                     breaks = seq(0.47, 0.51, 0.01)) +
  scale_fill_manual(values = '#ff8c3e') + 
  labs(x='Method',y='AUROC') +
  theme_test(base_size = 15) +
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, hjust = 0.5, vjust = 0.5))


variable <- rep("method", 8)

p2 <- ggplot(data,aes(Method, Runtime,color=variable,group=variable)) +
  geom_point(size=4) +
  geom_line(cex=1.3) +
  scale_color_manual(values = '#1e8b9b')+
  theme_test(base_size = 15) +
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, hjust = 0.5, vjust = 0.5))






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




# pdf('FigurescBioRevised/RealdataAUROC_Data2.pdf', width = 5.5, height = 5.5)
p3 
# dev.off()



