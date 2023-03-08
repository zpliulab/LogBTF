## 2022.12.2 PNAS 2018 RegNetwork -- 72 node 139 edges
## 2022.12.3 SERGIO simulated network analysis -- De-noised_100G_9T_300cPerT_4_DS1
## 2022.12.5 SERGIO simulated network -- add regulator SIGN

## 2022.12.23 PNAS 2008 TLGL_PNAS T-LGL Survival Signaling Network.
## 2022.12.23 PNAS 2021 Breast Cancer network FROM pnas.2022598118
## 2023.1.9 Delete some drugs and get BRCAdrugone.csv

## 2023.2.14 Add Simulated SERGIO - MElanoma data


################################# PNAS 2018 ###########################################
#######################################################################################

rm(list = ls())


setwd("/Users/lilingyu/E/PhD/R/Landscape/Data/PNASData")


## load data
# data <- readxl::read_excel("pnas.1722609115.sd01.xlsx", sheet="Sheet1")[c(1:72),]
# data <- readxl::read_excel("Real data.xlsx")
# data <- read.csv("Real data.csv", header=TRUE, sep = ',')


## load data
data <- readxl::read_excel("pnas.1722609115.sd02.xlsx", sheet="Sheet1")


intersect(which(is.na(data[,2]) == "TRUE"), which(is.na(data[,3]) == "TRUE"))
data[intersect(which(is.na(data[,2]) == "TRUE"), which(is.na(data[,3]) == "TRUE")), 1]



## define a dunction to seprate the list and get a list with 1 or -1
library(dplyr)
library(tidyr)
library(stringr)


############# 记住：第一列 Node 是Target；Promote 和 Inhibition 才是 TF
seplist <- function(data, l, label){
  
  ## data is three cols
  ## l is the acti (2) or inhi (3)
  ## label is the acti (1) or inhi (-1)
  
  # l <- 2
  # label <- 1
  
  exprset3 <- data
  a <- tibble(exprset3[,c(l,1)])
  
  test1 <- apply(a,1, function(x){
    str_split(x[1],',',simplify=T)       
  })
  
  test2 <- apply(a, 1, function(x){         
    paste(str_split(x[1],', ', simplify=T), x[2],sep = "---")
  })
  
  x <- tibble(unlist(test2))       
  colnames(x) <- "lala"
  x2 <- separate(x,lala,c("TF", "Target"),sep = '---') 
  
  ## delete NA and add acti or inhi
  x3 <- x2[-which(x2[,1] == "NA"),]
  x3$Sign <- rep(label, dim(x3)[1])
  
  return(x3)

}


actlist <- seplist(data, l=2, label =1)
inhilist <- seplist(data, l=3, label = -1)
reglist <- rbind(actlist, inhilist)



## original TF and Target is inverse, now correct them
## by inverstigate the Figure S1 and the Table S2
# library(dplyr)
# reglist = reglist %>% select(Target, TF)
# colnames(reglist) <- c("TF", "Target")

reglist[1,1]
reglist[1,2]

## save network and sign
# write.csv(reglist[,c(1,2)], file = "PNASnetwork3col.csv", row.names = F, quote = F)
# write.csv(reglist, file = "PNASnetworkSign.csv", row.names = F, quote = F)



## statistic information
## self roof 3
reglist <- reglist[-which(reglist[,1] == reglist[,2]),]


## Tf statistic 
## https://www.jianshu.com/p/2cf5ba339945?ivk_sa=1024320u
TFgene <- data.frame(rle(sort(reglist$TF))[2], rle(sort(reglist$TF))[1])
colnames(TFgene) <- c("TF", "Number")
str(rle(sort(reglist$TF)))
colnames(TFgene)



## target statistic
Targene <- data.frame(rle(sort(reglist$Target))[2], rle(sort(reglist$Target))[1])
colnames(Targene) <- c("Target", "Number")
str(rle(sort(reglist$Target)))
colnames(Targene)


colnames(TFgene)
colnames(Targene)


## plot 
data = TFgene
num = 8
step = 1


library(ggplot2)
pTF <- ggplot(data, aes(x=TF, y=Number)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_minimal()+
  coord_flip() +
  theme(legend.position = 'none',
        # panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', 
                                   # angle = 45, hjust = 0.5, vjust = 0.5, 
                                   size = 11)) +
  scale_y_continuous(limits = c(0.0,num),
                     breaks = seq(0.0, num, step),
                     expand = c(0,0)
  ) 

pTF



data = Targene
num = 10
step = 1


pTar <- ggplot(data, aes(x=Target, y=Number)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_minimal()+
  coord_flip() +
  theme(legend.position = 'none',
        # panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', 
                                   # angle = 45, hjust = 0.5, vjust = 0.5, 
                                   size = 11)) +
  scale_y_continuous(limits = c(0.0,num),
                     breaks = seq(0.0, num, step),
                     expand = c(0,0)
  ) 

pTar



# pdf("PNASnetworkTFnum.pdf", width = 4, height = 8)
# p
# dev.off()

# pdf("PNASnetworkTarnum.pdf", width = 4, height = 8)
# p
# dev.off()



# 双向柱形图 ------------------------------------------------------------
## https://blog.csdn.net/qq_35294674/article/details/124669111

length(unique(TFgene$TF))
length(unique(Targene$Target))
TFalone <- setdiff(TFgene$TF, Targene$Target)


## data merge
TFTar <- merge(TFgene, Targene, by.x = "TF", by.y = "Target")
colnames(TFTar) <- c("Gene", "TF num", "Target num")

which(TFgene[,1] %in% TFalone)
addlist <- TFgene[which(TFgene[,1] %in% TFalone), ]
colnames(addlist) <- c("Gene", "TF num")
addlist$`Target num` <- rep(0, length(which(TFgene[,1] %in% TFalone)))

TFTargene <- rbind(TFTar, addlist)



library(ggplot2)
library(reshape2)


## dounle direction bar plot
## https://blog.csdn.net/qq_35294674/article/details/124669111

df <- TFTargene
df <- melt(df)   

p<- ggplot(df, aes(
  # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  x = factor(Gene,levels = unique(Gene)),   
  # 判断分组情况，将两个柱子画在0的两侧
  y = ifelse(variable == "TF num", value, -value),  
  fill = variable)) +
  # 画柱形图
  geom_bar(stat = 'identity')+   
  # x轴与y轴互换位置
  coord_flip()+
  # 在图形上加上数字标签
  geom_text(                                                  
    aes(label=value, 
        # 标签的值（数据框的第三列）
        # 垂直位置。如果没有coord_flip()，则可以取消这行注释
        # vjust = ifelse(variable == "Up", -0.5, 1), 
        # 水平位置 and # 标签大小
        hjust = ifelse(variable == "TF num", -0.4, 1.1)),size=2) +
  # 调整y轴    # 刻度设置为绝对值
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1)))+                
  scale_fill_manual(values = c('#fec79e','#8ec4cb'))+
  labs(x='Gene',y='Number') +
  # 在y轴的两侧，留下一部分的空白位置，防止加标签的时候，显示不全
  theme_test(base_size = 10) 

p
# pdf("PNASnetworkTFTarnum.pdf", width = 5, height = 8)
# p
# dev.off()




################################# SERGIO ##############################################
#######################################################################################

rm(list = ls())


setwd("/Users/lilingyu/E/PhD/Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLY")


## load data
# data <- readxl::read_excel("pnas.1722609115.sd01.xlsx", sheet="Sheet1")[c(1:72),]
# data <- readxl::read_excel("Real data.xlsx")
# data <- read.csv("Real data.csv", header=TRUE, sep = ',')


## load data
data <- read.csv("gt_GRN.csv", header=TRUE, sep = ',')
## TF - 10
length(unique(data[,1]))
## Target - 93
length(unique(data[,2]))
length(setdiff(unique(data[,2]), unique(data[,1])))


## different cols
# inter <- read.table("Interaction_cID_4.txt", header = F, sep = ',')
# inter <- read.csv("Interaction_cID_4.csv", header = F, sep = ',')

## 9 cell type
reg <- read.table("Regs_cID_4.txt", header = F, sep = ',')
## 2 cell type
reg2 <- reg[,c(1:3)]
# write.table(reg2, "Regs_cID_4_2type.txt", row.names = F, col.names = F, sep = ',')

## 3 cell type
reg3 <- reg[,c(1:4)]
# write.table(reg3, "Regs_cID_4_3type.txt", row.names = F, col.names = F, sep = ',')



################################# SERGIO RegNetwork SIGN ##############################
#######################################################################################

rm(list = ls())


setwd("/Users/lilingyu/E/PhD/Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLY")


## different cols
inter <- read.csv("Interaction_cID_4.csv", header = F, sep = ',')

RegSIGN <- c()
list <- c()

for (i in 1:dim(inter)[1]) {
  # i <- 1
  for (j in 1:inter[,2][i]) {
    # j <- 1
    list <- cbind(inter[i,2+j], inter[i,1], inter[i,2 + inter[,2][i] + j]) 
    RegSIGN <- rbind(RegSIGN, list)
  }

}

dim(RegSIGN)


RegSIGN[which(RegSIGN[,3] > 0),3] <- 1
RegSIGN[which(RegSIGN[,3] < 0),3] <- 2
colnames(RegSIGN) <- c("from", "to", "type")
RegSIGNoutput <- as.data.frame(RegSIGN)
RegSIGNoutput[2,2]
RegSIGNoutput[2,3]
RegSIGNoutput[4,]
RegSIGNoutput$from <- as.character(RegSIGNoutput$from)
RegSIGNoutput$to <- as.character(RegSIGNoutput$to)
# write.csv(RegSIGNoutput, "Interaction_cID_4SIGN.csv", row.names = F)


## draft bellow ###########################################################
i <- 1
inter[,2][i]

j <- 1
inter[l,2+j]
inter[l,2 + inter[,2][i] + j]
cbind(inter[l,2+j], inter[l,1], inter[l,2 + inter[,2][i] + j]) 

inter[l,1]
###########################################################################



################################# PNAS 2021 Breast Cancer network #####################
#######################################################################################

rm(list = ls())

setwd("/Users/lilingyu/E/PhD/Python/bmodel/Data/Network/BreastCancerModel")


data <- as.data.frame(read.csv("BreastCancerModel.txt", header = T))
data[1,]

library(stringr)
str_split(data[1,],"=", 2)



## split by = to source[2] and target[1] nodes
source <- c()
target <- c()

for (i in 1:dim(data)[1]) {
  # i <- 1
  target <- rbind(target, str_split(data[i,], "=", 2)[[1]][1])
  source <- rbind(source, str_split(data[i,], "=", 2)[[1]][2]) 
}


## delete * in target
## https://www.mianshigee.com/note/detail/137235vvj/
target <- gsub('.{1}$','',target)


## split by not to active and inhibit
active <- c()
inhibt <- c()
for (i in 1:dim(data)[1]) {
  # i <- 3
  active <- rbind(active, str_split(source[i,], "not", 2)[[1]][1])
  inhibt <- rbind(inhibt, str_split(source[i,], "not", 2)[[1]][2]) 
}

# colnames(data0)
net <- cbind(target, active, inhibt)
colnames(net) <- c("Node", "Promoters", "Inhibitors")
# write.csv(net, file = "BreastCancerDrugNet.csv", row.names = F)




# Boolean rule to Node-Promote-Inhibition ---------------------------------


## 接下来，手动整理.xls文件 !!!---- 注意 空格、括号 是不是去干净 ！！！！
## not(A or B) 和 not(A and B) --  A 和 B 均为 inhibition 


## load data
rm(list = ls())

setwd("/Users/lilingyu/E/PhD/Python/bmodel/Data/Network/BreastCancerModel")
data <- readxl::read_excel("BreastCancerDrugNet.xlsx", sheet="Sheet1")


## define a dunction to seprate the list and get a list with 1 or -1
library(dplyr)
library(tidyr)
library(stringr)
############# 记住：第一列 Node 是Target；Promote 和 Inhibition 才是 TF
seplist <- function(data, l, label){
  
  ## data is three cols
  ## l is the acti (2) or inhi (3)
  ## label is the acti (1) or inhi (-1)
  
  # l <- 2
  # label <- 1
  
  exprset3 <- data
  a <- tibble(exprset3[,c(l,1)])
  
  test1 <- apply(a,1, function(x){
    str_split(x[1],',',simplify=T)       
  })
  
  test2 <- apply(a, 1, function(x){         
    paste(str_split(x[1],', ', simplify=T), x[2],sep = "---")
  })
  
  x <- tibble(unlist(test2))       
  colnames(x) <- "lala"
  x2 <- separate(x,lala,c("TF", "Target"),sep = '---') 
  
  ## delete NA and add acti or inhi
  x3 <- x2[-which(x2[,1] == "NA"),]
  x3$Sign <- rep(label, dim(x3)[1])
  
  return(x3)
  
}


actlist <- seplist(data, l=2, label =1)
inhilist <- seplist(data, l=3, label = -1)
reglist <- rbind(actlist, inhilist)


reglist[1,1]
reglist[1,2]


## statistic information
length(which(reglist[,1] == reglist[,2]))    # 22
## self roof 3
reglist <- reglist[-which(reglist[,1] == reglist[,2]),]


setdiff(union(as.matrix(reglist[,1]), as.matrix(reglist[,2])), data$Node)

## save network and sign -- no self-loof
# write.csv(reglist[,c(1,2)], file = "PNASnetwork.csv", row.names = F, quote = F)
# write.csv(reglist, file = "PNASnetworkSign.csv", row.names = F, quote = F)


## delete repeat pairs
reglist <- reglist[!duplicated(reglist),]    # 203

## 保存成 Python 可以用的 .csv 文件
reglist[which(reglist[,3] > 0),3] <- 1
reglist[which(reglist[,3] < 0),3] <- 2
colnames(reglist) <- c("from", "to", "type")
union(as.matrix(reglist[,1]), as.matrix(reglist[,2]))    # 80
# write.csv(reglist, file = "BRCAdrug.csv", row.names = F, quote = F)
## 改变存储路径  --  python 的路径
# write.csv(reglist, file = "/Users/lilingyu/E/PhD/R/Landscape/Data/BRCAdata/BRCAdrug.csv", row.names = F, quote = F)


# delete drug and keep only PI3K drug -------------------------------------------------------------

drug6 <- c("Fulvestrant", "Everolimus", "Trametinib", "Ipatasertib",
           "Palbociclib", "Neratinib")

reglist1drug <- reglist
reglist1drug[which(as.matrix(reglist1drug[,1]) %in% drug6), ]
reglist1drug <- reglist1drug[-which(as.matrix(reglist1drug[,1]) %in% drug6),]


union(as.matrix(reglist1drug[,1]), as.matrix(reglist1drug[,2]))    # 75

setdiff(as.matrix(union(as.matrix(reglist[,1]), as.matrix(reglist[,2]))),
        as.matrix(union(as.matrix(reglist1drug[,1]), as.matrix(reglist1drug[,2]))))

# write.csv(reglist1drug, file = "BRCAdrugone.csv", row.names = F, quote = F)
# # 改变存储路径  --  python 的路径
# write.csv(reglist1drug, file = "/Users/lilingyu/E/PhD/R/Landscape/Data/BRCAdata/BRCAdrugone.csv", row.names = F, quote = F)


# delete drug and cpmbina PI3K and PI3K_2 -------------------------------------------------------------

drug7 <- c("Fulvestrant", "Alpelisib", "Everolimus", "Trametinib", "Ipatasertib",
          "Palbociclib", "Palbociclib")

drug6 <- c("Fulvestrant", "Everolimus", "Trametinib", "Ipatasertib",
          "Palbociclib", "Palbociclib")

reglist1drug <- reglist
reglist1drug <- reglist1drug[-which(as.matrix(reglist1drug[,1]) %in% drug6),]

lab1 <- which(as.matrix(reglist1drug[,1])%in% "PI3K_2")
lab2 <- which(as.matrix(reglist1drug[,2])%in% "PI3K_2")

reglist1drug[lab1, 1] <- "PI3K"
reglist1drug[lab2, 2] <- "PI3K"

reglist1DRUG <- reglist1drug[!duplicated(reglist1drug),]

# write.csv(reglist1DRUG, file = "BRCAdrug.csv", row.names = F, quote = F)
## 改变存储路径  --  python 的路径
# write.csv(reglist1DRUG, file = "/Users/lilingyu/E/PhD/R/Landscape/Data/BRCAdata/BRCAdrug.csv", row.names = F, quote = F)




## Tf statistic  ######################  Plot figure #########################
## https://www.jianshu.com/p/2cf5ba339945?ivk_sa=1024320u
TFgene <- data.frame(rle(sort(reglist$TF))[2], rle(sort(reglist$TF))[1])
colnames(TFgene) <- c("TF", "Number")
str(rle(sort(reglist$TF)))
colnames(TFgene)


## target statistic
Targene <- data.frame(rle(sort(reglist$Target))[2], rle(sort(reglist$Target))[1])
colnames(Targene) <- c("Target", "Number")
str(rle(sort(reglist$Target)))
colnames(Targene)


colnames(TFgene)
colnames(Targene)


## plot 
data = TFgene
num = 10
step = 1


library(ggplot2)
pTF <- ggplot(data, aes(x=TF, y=Number)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_minimal()+
  coord_flip() +
  theme(legend.position = 'none',
        # panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', 
                                   # angle = 45, hjust = 0.5, vjust = 0.5, 
                                   size = 11)) +
  scale_y_continuous(limits = c(0.0,num),
                     breaks = seq(0.0, num, step),
                     expand = c(0,0)
  ) 

pTF



data = Targene
num = 10
step = 1


pTar <- ggplot(data, aes(x=Target, y=Number)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_minimal()+
  coord_flip() +
  theme(legend.position = 'none',
        # panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', 
                                   # angle = 45, hjust = 0.5, vjust = 0.5, 
                                   size = 11)) +
  scale_y_continuous(limits = c(0.0,num),
                     breaks = seq(0.0, num, step),
                     expand = c(0,0)
  ) 

pTar



# pdf("PNASnetworkTFnum.pdf", width = 4, height = 8)
# p
# dev.off()

# pdf("PNASnetworkTarnum.pdf", width = 4, height = 8)
# p
# dev.off()



#################  双向柱形图  #############################################

length(unique(TFgene$TF))
length(unique(Targene$Target))
TFalone <- setdiff(TFgene$TF, Targene$Target)


## data merge
TFTar <- merge(TFgene, Targene, by.x = "TF", by.y = "Target")
colnames(TFTar) <- c("Gene", "TF num", "Target num")

which(TFgene[,1] %in% TFalone)
addlist <- TFgene[which(TFgene[,1] %in% TFalone), ]
colnames(addlist) <- c("Gene", "TF num")
addlist$`Target num` <- rep(0, length(which(TFgene[,1] %in% TFalone)))

TFTargene <- rbind(TFTar, addlist)



library(ggplot2)
library(reshape2)


## dounle direction bar plot
## https://blog.csdn.net/qq_35294674/article/details/124669111

df <- TFTargene
df <- melt(df)   

p<- ggplot(df, aes(
  # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  x = factor(Gene,levels = unique(Gene)),   
  # 判断分组情况，将两个柱子画在0的两侧
  y = ifelse(variable == "TF num", value, -value),  
  fill = variable)) +
  # 画柱形图
  geom_bar(stat = 'identity')+   
  # x轴与y轴互换位置
  coord_flip()+
  # 在图形上加上数字标签
  geom_text(                                                  
    aes(label=value, 
        # 标签的值（数据框的第三列）
        # 垂直位置。如果没有coord_flip()，则可以取消这行注释
        # vjust = ifelse(variable == "Up", -0.5, 1), 
        # 水平位置 and # 标签大小
        hjust = ifelse(variable == "TF num", -0.4, 1.1)),size=2) +
  # 调整y轴    # 刻度设置为绝对值
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1)))+                
  scale_fill_manual(values = c('#fec79e','#8ec4cb'))+
  labs(x='Gene',y='Number') +
  # 在y轴的两侧，留下一部分的空白位置，防止加标签的时候，显示不全
  theme_test(base_size = 10) 

p
# pdf("PNASnetworkTFTarnum.pdf", width = 5, height = 8)
# p
# dev.off()



############################ Pathway information ##############################

pathway <- read.table("BreastCancerDrugPathway.txt", header = F, sep = "\t")

library(stringr)


Tf <- c()
Tf0 <- c()
path <- c()
path0 <- c()
for (i in 1:dim(pathway)[1]) {
  Tf <- rbind(Tf, str_split_fixed(pathway[i,], " = ", 2)[[1]][1])
  Tf0 <- rbind(Tf0, str_split_fixed(Tf[i,], "att", 2)[2])
  path <- rbind(path, str_split_fixed(pathway[i,], "type:", 2)[2])
  path0 <- rbind(path0, str_split_fixed(path[i,], " width", 2)[1])
}

## delete [  ] in target
## https://qa.1r1g.com/sf/ask/2457948741/
## 对于第一/最后一个字符,匹配string(^.)开头的字符或字符串()的结尾,.$并替换为''.
gsub('^.|.$', '', Tf0)
gsub('.$', '', path0)
Pathwayname <- cbind(gsub('^.|.$', '', Tf0), gsub('.$', '', path0))
colnames(Pathwayname) <- c("TF", "pathway")

## save
# write.csv(Pathwayname, "Pathwayname.csv", row.names = F)

## merge
reglistpathway <- merge(reglist, Pathwayname, by = "TF")
# write.csv(reglistpathway, "PNASnetworkSignPath.csv", row.names = F, quote = F)




################################# SERGIO -Melanoma RegNetwork SIGN ##############################
#######################################################################################

rm(list = ls())


setwd("/Users/lilingyu/E/PhD/Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised")


## different cols
inter <- read.csv("Interaction_cID_4_17.csv", header = F, sep = ',')

RegSIGN <- c()
list <- c()

for (i in 1:dim(inter)[1]) {
  # i <- 1
  for (j in 1:inter[,2][i]) {
    # j <- 1
    list <- cbind(inter[i,2+j], inter[i,1], inter[i,2 + inter[,2][i] + j]) 
    RegSIGN <- rbind(RegSIGN, list)
  }
  
}

dim(RegSIGN)


RegSIGN[which(RegSIGN[,3] > 0),3] <- 1
RegSIGN[which(RegSIGN[,3] < 0),3] <- 2
colnames(RegSIGN) <- c("from", "to", "type")
RegSIGNoutput <- as.data.frame(RegSIGN)
RegSIGNoutput[2,2]
RegSIGNoutput[2,3]
RegSIGNoutput[4,]
RegSIGNoutput$from <- as.character(RegSIGNoutput$from)
RegSIGNoutput$to <- as.character(RegSIGNoutput$to)
# write.csv(RegSIGNoutput, "Interaction_cID_4_17SIGN.csv", row.names = F)


################################# SERGIO 20 node RegNetwork SIGN ##############################
#######################################################################################

rm(list = ls())


setwd("/Users/lilingyu/E/PhD/Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised")


## different cols
inter <- read.csv("Interaction_cID_4_20node.csv", header = F, sep = ',')

RegSIGN <- c()
list <- c()

for (i in 1:dim(inter)[1]) {
  # i <- 1
  for (j in 1:inter[,2][i]) {
    # j <- 1
    list <- cbind(inter[i,2+j], inter[i,1], inter[i,2 + inter[,2][i] + j]) 
    RegSIGN <- rbind(RegSIGN, list)
  }
  
}

dim(RegSIGN)


RegSIGN[which(RegSIGN[,3] > 0),3] <- 1
RegSIGN[which(RegSIGN[,3] < 0),3] <- 2
colnames(RegSIGN) <- c("from", "to", "type")
RegSIGNoutput <- as.data.frame(RegSIGN)
RegSIGNoutput[2,2]
RegSIGNoutput[2,3]
RegSIGNoutput[4,]
RegSIGNoutput$from <- as.character(RegSIGNoutput$from)
RegSIGNoutput$to <- as.character(RegSIGNoutput$to)
# write.csv(RegSIGNoutput, "Interaction_cID_4SIGN_20node.csv", row.names = F)

