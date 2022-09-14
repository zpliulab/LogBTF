rm(list = ls())


setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")


# data generation  --------------------------------------------------------

## parameter setting
p <- 9
l <- 2^p 

## 十进制转二进制，并补全9位
library(R.utils)
aa <- matrix("0", p, 1)
## aa 第1列为0，即可从1开始
for (i in 1:(l-1)) {
  # i <- 2
  num <- intToBin(i)
  # 拆分字符串
  numchara <-strsplit(num,'')
  numnum <- data.frame(as.numeric(numchara[[1]]))
  sub <- data.frame(rep(0, p-length(as.numeric(numchara[[1]]))))
  colnames(numnum) <- colnames(sub) <- c("var")
  numnew <- rbind(sub, numnum)
  aa <- cbind(aa, numnew) 
}
aa
dim(aa)
View(aa[,1:10])
View(aa[,dim(aa)[2]])


## 添加列名
namenum <- rep(1:l,1)
name <- c()
library(stringr)
for (j in 1:l) {
  # j <- 2
  name <- cbind(name, str_c("var_",namenum[j]))
}
name
View(name)


## 添加行名，并对0/1进行数值化
b <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
ab <- cbind(b,aa)
colnames(ab) <- c("gene_name", name)
dim(ab)
View(ab[,1:10])
ab[2,2]
bb <- apply(ab[,2:(l+1)], 2, as.numeric)
dim(bb)
bb[1,2]
View(bb[,1:10])
rownames(bb) <- b
## 输出
# write.csv(bb, file = "dataAllList.csv", row.names = F)
# install.packages("xlsx")
library(xlsx)
# write.xlsx(bb, file = "dataAllList.xlsx")



# Data progress -----------------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\Atta")


# load Rdata use read.table, 此时0/1均为character
setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")
myfile = list.files("Atta")    
dir = paste("./Atta/", myfile, sep="")  
n = length(dir) 

att0 = as.matrix(read.table(file = dir[1], header=T)[,-1])
## for cycle
for (i in 2:n) {
  # i <- 2
  att = as.matrix(read.table(file = dir[i], header=T)[,-1]  )
  att0 <- rbind(att0, att)
}
class(att0)
att0[1,1]
# att0[att0] <- 1
att0[which(att0 == "True")] <- 1
att0[which(att0 == "False")] <- 0
att0hat <- apply(att0,2, as.numeric)


## 分别提取x和y
xyselect <- function(x){
  att1 = as.matrix(read.table(file = dir[1], header=T)[x,-1])
  for (i in 2:n) {
    att = as.matrix(read.table(file = dir[i], header=T)[x,-1])
    att1 <- rbind(att1, att)
  }
  att1[which(att1 == "True")] <- 1
  att1[which(att1 == "False")] <- 0
  att1hat <- apply(att1,2, as.numeric)
  return(att1hat)
}

## initial
att1hat <- xyselect(1)
## t+1 时刻
att2hat <- xyselect(2)


## glm model - 只取9个真实的 X and Y
x <- data.frame(att1hat) 
library(glmnet)

coefall <- c()
for (j in 1:(dim(att)[2])) {
  # j <- 1
  y <- att2hat[,j]
  glm.fit <- glm(y~., data = x, family = binomial, control = list(maxit = 100))
  coef <- round(as.vector(glm.fit$coefficients),3)
  coefall <- rbind(coefall, coef)
}

colnames(coefall) <- c("theta", colnames(att))
rownames(coefall) <- colnames(att)
# write.csv(coefall, file="coefall.csv")




# Robust analysis ------------------------------------------------------------------


rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\Atta")


# load Rdata use read.table, 此时0/1均为character
setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")
myfile = list.files("Atta")    
dir = paste("./Atta/", myfile, sep="")  
n = length(dir) 

att0 = as.matrix(read.table(file = dir[1], header=T)[,-1])
## for cycle
for (i in 2:n) {
  # i <- 2
  att = as.matrix(read.table(file = dir[i], header=T)[,-1]  )
  att0 <- rbind(att0, att)
}
class(att0)
att0[1,1]
# att0[att0] <- 1
att0[which(att0 == "True")] <- 1
att0[which(att0 == "False")] <- 0
att0hat <- apply(att0,2, as.numeric)


## 分别提取x和y
xyselect <- function(x){
  att1 = as.matrix(read.table(file = dir[1], header=T)[x,-1])
  for (i in 2:n) {
    att = as.matrix(read.table(file = dir[i], header=T)[x,-1])
    att1 <- rbind(att1, att)
  }
  att1[which(att1 == "True")] <- 1
  att1[which(att1 == "False")] <- 0
  att1hat <- apply(att1,2, as.numeric)
  return(att1hat)
}

## initial
att1hat <- xyselect(1)
## t+1 时刻
att2hat <- xyselect(2)

## random flipping  
## https://blog.csdn.net/sinat_38918949/article/details/82144344
## https://www.datasciencemadesimple.com/sample-function-in-r/
n <- dim(att1hat)[1]
p <- dim(att1hat)[2]

np <- dim(att1hat)[1]*dim(att1hat)[2]
flip <- 0.01

index <- c(1:np)    # c(1:(n*p)) different c(1:n*p)
set.seed(123)
index <- sample(index, round(flip*np), replace = FALSE, prob = NULL)  ## 向下取整
att1hat[index[which(att1hat[index] == 1)]] <- 0
att1hat[index[which(att1hat[index] == 0)]] <- 1


## glm model - 只取9个真实的 X and Y
x <- data.frame(att1hat) 
library(glmnet)

coefall <- c()
for (j in 1:(dim(att)[2])) {
  # j <- 1
  y <- att2hat[,j]
  glm.fit <- glm(y~., data = x, family = binomial, control = list(maxit = 100))
  coef <- round(as.vector(glm.fit$coefficients),3)
  coefall <- rbind(coefall, coef)
}

colnames(coefall) <- c("theta", colnames(att))
rownames(coefall) <- colnames(att)
View(coefall)
# write.csv(coefall, file="coefall.csv")
 

## updata coeff
adj_matrix <- t(coefall)[-1,]

## Final ranked list, AUROC and AUPR
adjstand <- t(apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} ))
adjstand[adjstand <= 0.5] = 0
adj_matrix <- adjstand * t(sign(adj_matrix))
## 把NAN变成0
adj_matrix[is.na(adj_matrix)] <- 0
# adj_matrix[adj_matrix > 0] = 1



## performance


SIGN = 1

## Gold standard GRN 
setwd("D:/E/博士/R_程序/Boolean/Data/AttaNetwork")
## 列为target
adj_gold <- as.matrix(read.csv("gold standard.csv",sep=",", header = F))
## 转置一下系数


# no self-loof, direct
if(SIGN == 1){
  # adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  adj_matrix <- t(coefall)[-1,]
}
# undirect
if(SIGN == 0){
  adj_gold[which(adj_gold != 0)] <- 1
  adj_matrix <- abs(t(coefall)[-1,])
}


## Final ranked list, AUROC and AUPR
adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
adjstand[adjstand <= 0.5] = 0
# adjstand[adjstand >0] = 1     ## for No noise
adj_matrix <- adjstand * sign(adj_matrix)
## 把NAN变成0
adj_matrix[is.na(adj_matrix)] <- 0
# adj_matrix[adj_matrix > 0] = 1


## load function
source('D:\\E\\博士\\R_程序/Boolean/R/PermanceFunction.R')
library(pracma)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Acc <- performance$Acc
Recall <- performance$Recall
Pre <- performance$Pre
FPR <- performance$FPR
Fmeasure <- performance$Fmeasure
Result <- cbind(Acc,Recall,Pre,FPR,Fmeasure)


ResultAll <- as.matrix(cbind(AUC, Result))
colnames(ResultAll) <- c("AUROC", "AUPR", "Acc","Recall","Pre","FPR","Fmeasure")
ResultAll



# Noise Data progress -----------------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\Atta")


# load Rdata use read.table, 此时0/1均为character
setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")
myfile = list.files("Atta")    
dir = paste("./Atta/", myfile, sep="")  
n = length(dir) 

att0 = as.matrix(read.table(file = dir[1], header=T)[,-1])
## for cycle
for (i in 2:n) {
  # i <- 2
  att = as.matrix(read.table(file = dir[i], header=T)[,-1]  )
  att0 <- rbind(att0, att)
}
class(att0)
att0[1,1]
# att0[att0] <- 1
att0[which(att0 == "True")] <- 1
att0[which(att0 == "False")] <- 0
att0hat <- apply(att0,2, as.numeric)


## 分别提取x和y
xyselect <- function(x){
  att1 = as.matrix(read.table(file = dir[1], header=T)[x,-1])
  for (i in 2:n) {
    att = as.matrix(read.table(file = dir[i], header=T)[x,-1])
    att1 <- rbind(att1, att)
  }
  att1[which(att1 == "True")] <- 1
  att1[which(att1 == "False")] <- 0
  att1hat <- apply(att1,2, as.numeric)
  return(att1hat)
}

## initial
att1hat <- xyselect(1)
## t+1 时刻
att2hat <- xyselect(2)


## glm model - 只取9个真实的 X and Y
x <- data.frame(att1hat) 


## 给X添加扰动
## noise
res <- att1hat
noiseLevel = 1e-5
res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))

## 将原来为0的变成0
for (i in 1:dim(att1hat)[1]) {
  # i <- 1
  for (j in 1:dim(att1hat)[2]) {
    # j <- 2
    if (att1hat[i,j] == 0)
      res[i,j] <- 0
  }
}
x <- data.frame(res)


library(glmnet)
coefall <- c()
for (j in 1:(dim(att)[2])) {
  # j <- 1
  y <- att2hat[,j]
  glm.fit <- glm(y~., data = x, family = binomial, control = list(maxit = 100))
  coef <- round(as.vector(glm.fit$coefficients),3)
  coefall <- rbind(coefall, coef)
}

colnames(coefall) <- c("theta", colnames(att))
rownames(coefall) <- colnames(att)
# write.csv(coefall, file="coefall.csv")


# train and test, AUC -----------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\Atta")


# load Rdata use read.table, 此时0/1均为character
setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")
myfile = list.files("Atta")    
dir = paste("./Atta/", myfile, sep="")  
n = length(dir) 

att0 = as.matrix(read.table(file = dir[1], header=T)[,-1])
## for cycle
for (i in 2:n) {
  # i <- 2
  att = as.matrix(read.table(file = dir[i], header=T)[,-1]  )
  att0 <- rbind(att0, att)
}
class(att0)
att0[1,1]
# att0[att0] <- 1
att0[which(att0 == "True")] <- 1
att0[which(att0 == "False")] <- 0
att0hat <- apply(att0,2, as.numeric)


## 分别提取x和y
xyselect <- function(x){
  att1 = as.matrix(read.table(file = dir[1], header=T)[x,-1])
  for (i in 2:n) {
    att = as.matrix(read.table(file = dir[i], header=T)[x,-1])
    att1 <- rbind(att1, att)
  }
  att1[which(att1 == "True")] <- 1
  att1[which(att1 == "False")] <- 0
  att1hat <- apply(att1,2, as.numeric)
  return(att1hat)
}

## initial
att1hat <- xyselect(1)
## t+1 时刻
att2hat <- xyselect(2)


## glm model - 只取9个真实的 X and Y
x <- data.frame(att1hat) 


## x1 for example
## 观察y的0-1
y <- att2hat[,2]
sum(y)    # 384 个 1
length(y)-sum(y)  ## 128个0

# install.packages('caret', dependencies = TRUE)
library(caret)
library(tidyverse)
library(pROC)

aucallrep <- c()


aucall <- c()
for (k in 1:10) {
  
  set.seed(123*k)
  
  training_samples <- y %>% createDataPartition(p = 0.02, list = FALSE)
  train_data  <- x[training_samples, ]
  train_label <- y[training_samples ]
  
  glm.fit <- glm(train_label~., data = train_data, family = binomial, control = list(maxit = 100))
  p <- predict(glm.fit, x, type = "response")
  A_train <- data.frame(as.matrix(p), y)
  names(A_train)<- c("p", "outcome")
  aucplot <- plot.roc(A_train$outcome, A_train$p, print.auc=T)
  auc <- aucplot$auc
  aucall <- cbind(aucall, auc)
}
aucall
mean(aucall)
sd(aucall)

aucallrep <- rbind(aucallrep, aucall)


num <- rep(1:10)

# pernumall <- c()
# for (l in 1:10) {
#   # l <- 1
#   pernum <- str_c("Percent", num[l])
#   pernumall <- cbind(pernumall, pernum)
# }
# rownames(aucallrep) <- pernumall


strnumall <- c()
for (l in 1:10) {
  # l <- 1
  strnum <- str_c("Seed", num[l])
  strnumall <- cbind(strnumall, strnum)
}
colnames(aucallrep) <- strnumall

## save data
# write.csv(aucallrep, file = "AUCall\\aucallrep9.csv")




# error bar  -----------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\AUCall")


data1 <- read.csv("aucallrep1.csv", header = T, sep=",", row.names = 1)
data2 <- read.csv("aucallrep2.csv", header = T, sep=",", row.names = 1)
data3 <- read.csv("aucallrep3.csv", header = T, sep=",", row.names = 1)
data4 <- read.csv("aucallrep4.csv", header = T, sep=",", row.names = 1)
data5 <- read.csv("aucallrep5.csv", header = T, sep=",", row.names = 1)
data6 <- read.csv("aucallrep6.csv", header = T, sep=",", row.names = 1)
data7 <- read.csv("aucallrep7.csv", header = T, sep=",", row.names = 1)
data8 <- read.csv("aucallrep8.csv", header = T, sep=",", row.names = 1)
data9 <- read.csv("aucallrep9.csv", header = T, sep=",", row.names = 1)


alllist <- list(data1,data2,data3,data4,data5,data6,data7,data8, data9)
name <- c("data1","data2","data3","data4","data5","data6","data7","data8", "data9")
dim(alllist[[2]])[1]
length(alllist)

# p <- as.numeric(c("10", "15", "20", "10"))
# p[1]

col1 <- col2 <- col3 <- col4 <- col5 <- c()
for (i in 1:length(alllist)) {
  # i <- 2
  n1 <- as.matrix(rep(name[i], dim(alllist[[i]])[1]))
  n2 <- as.matrix(rep(1: dim(alllist[[i]])[1]))
  n3 <- as.matrix(apply(alllist[[i]],1,mean))
  n4 <- as.matrix(apply(alllist[[i]],1,sd))
  n5 <- as.matrix(apply(alllist[[i]],1,sd)/sqrt(dim(alllist[[i]])[1]))
  
  col1 <- rbind(col1, n1)
  col2 <- rbind(col2, n2)
  col3 <- rbind(col3, n3)
  col4 <- rbind(col4, n4)
  col5 <- rbind(col5, n5)
}

# col1 <- rbind(as.matrix(rep("x1",p1)), as.matrix(rep("x2", p2)))
# col2 <- rbind(as.matrix(rep(1:p1)), as.matrix(rep(1:p2)))
# apply(data1,1,mean) #1代表显示行的平均值
# apply(data1,1,sd)/sqrt(10) # 标准误（Standard Error of Mean，SE） sd/sqrt(n)
# col3 <- as.numeric(rbind(as.matrix(apply(data1,1,mean)), as.matrix(apply(data2,1,mean))))
# col4 <- as.numeric(rbind(as.matrix(apply(data1,1,sd)), as.matrix(apply(data2,1,sd))))
# col5 <- as.numeric(rbind(as.matrix(apply(data1,1,sd)/sqrt(p1)), as.matrix(apply(data2,1,sd)/sqrt(p2))))


data <- cbind(col1, col2, as.numeric(col3), as.numeric(col4), as.numeric(col5))
colnames(data) <- c("Node", "Percent", "AUC", "SD", "SE")

data ## https://blog.csdn.net/qq_50522851/article/details/122051267

Data <- as.data.frame(lapply(data, as.numeric))
Data <- as.matrix(data)
setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")
# write.csv(data, file = "zhexiandata.csv", row.names = F)




# Plot  -------------------------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data")

Data <- read.csv("zhexiandata.csv", header = T, sep=",")

Data$Percent

Data[5,5]

library(ggplot2)
ggplot(Data, aes(x=Percent, y=AUC, colour=Node, shape=Node))+
  geom_line(size=0.75) +
  geom_point(size=4,shape=21,fill="white")+
  geom_errorbar(aes(ymin=AUC-SE, ymax=AUC+SE),colour="black", width=.03,size=0.75)+
  scale_colour_hue(name="Node",
                   breaks=c("data1","data2","data3","data4","data5","data6","data7","data8", "data9"),
                   labels=c("x 1- 7", "x2 - 14", "x3 - 16", "x4 - 6","x5 - 11", "x6 - 2", "x7 - 9", "x8 - 29", "x9 - 5"),
                   l=40) + 
  expand_limits(y=0.6) +                  
  scale_y_continuous(breaks=c(0.5, 0.55,0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)) +
  scale_x_continuous(breaks=c(rep(1:30)))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(color = 'black',size = 14),
        axis.text = element_text(color = 'black',size = 12), # element_blank()
        legend.text = element_text(color = 'black',size = 12),
        panel.grid.minor = element_blank(),legend.justification=c(0.85,0.1),legend.position=c(0.85,0.1)) 


# ggsave(Metabolite, filename = 'Metaboliteadd.pdf', width = 6, height = 6, device = cairo_pdf)


# error bar  -----------------------------------------------------

rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\AUCall")


data1 <- read.csv("aucallrep1.csv", header = T, sep=",", row.names = 1)
data2 <- read.csv("aucallrep2.csv", header = T, sep=",", row.names = 1)
data3 <- read.csv("aucallrep3.csv", header = T, sep=",", row.names = 1)
data4 <- read.csv("aucallrep4.csv", header = T, sep=",", row.names = 1)
data5 <- read.csv("aucallrep5.csv", header = T, sep=",", row.names = 1)
data6 <- read.csv("aucallrep6.csv", header = T, sep=",", row.names = 1)
data7 <- read.csv("aucallrep7.csv", header = T, sep=",", row.names = 1)
data8 <- read.csv("aucallrep8.csv", header = T, sep=",", row.names = 1)
data9 <- read.csv("aucallrep9.csv", header = T, sep=",", row.names = 1)




## for hiplot

hiplotdata <- function(x, p, data1){
  # x <- 1
  # p <- 7
  
  lab <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
  node1 <- rep(lab[x], 10*p)
  
  Per1 <- as.matrix(rep(1, 10))
  for (j in 2:p) {
    Per1 <- rbind(Per1, as.matrix(rep(j, 10)))
  }
  auc1 <- as.matrix(data1[1,])
  for (k in 2:p) {
    auc1 <- cbind(auc1, as.matrix(data1[k,]))
  }
  data <- cbind(t(auc1), as.matrix(node1), as.matrix(Per1))
  colnames(data) <- c("AUC", "Node", "Percent")
  
  return(data)
}

setwd("D:\\E\\博士\\R_程序\\Boolean\\Data\\AUCall_Hiplot")
# hip1 <- hiplotdata(1, 8, data1)
# write.csv(hip1, file = "hip1.csv" ,row.names = F)
# hip2 <- hiplotdata(2, 15, data2)
# write.csv(hip2, file = "hip2.csv" ,row.names = F)
# hip3 <- hiplotdata(3, 17, data3)
# write.csv(hip3, file = "hip3.csv" ,row.names = F)
# hip4 <- hiplotdata(4, 7, data4)
# write.csv(hip4, file = "hip4.csv" ,row.names = F)
# hip5 <- hiplotdata(5, 12, data5)
# write.csv(hip5, file = "hip5.csv" ,row.names = F)
# hip6 <- hiplotdata(6, 3, data6)
# write.csv(hip6, file = "hip6.csv" ,row.names = F)
# hip7 <- hiplotdata(7, 10, data7)
# write.csv(hip7, file = "hip7.csv" ,row.names = F)
# hip8 <- hiplotdata(8, 30, data8)
# write.csv(hip8, file = "hip8.csv" ,row.names = F)
# hip9 <- hiplotdata(9, 6, data9)
# write.csv(hip9, file = "hip9.csv" ,row.names = F)

