## 2022.7.23  使用GENIE3


# Glmnet + Glm noise ------------------------------------------------------------------

rm(list = ls())

# options(digits = 7)

## load function
# source('D:/E/博士/R_程序/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')

##
set.seed(123)


## 
# file <- c(10,1,1,53)
# file <- c(10,2,4,59)
# file <- c(10,3,2,40)
# file <- c(50,1,10,415)
# file <- c(50,2,20,557)
# file <- c(50,3,16,417)
# file <- c(100,1,19,1497)
# file <- c(100,2,4,1700)
file <- c(100,3,42,1718)




method = 2
genenum <- file[1]
filenum <- file[2]
k <- file[3]
geneselect <- file[4]


# GENIE3 -----------------------------------------------------------------


l <- 21

# datapath <- paste("D:/E/博士/R_程序/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Yeast",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
n <- dim(datanum)[1]
p <- dim(datanum)[2]


## Kmeans ??????
library(BiTrinA)
## features*times - A n x m matrix comprising m raw measurements of n features
# binarizeMatrix(t(datanum))
datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
# write.csv(t(datahatA), "CoefAllSize10\\datadream3p1_matrix.csv")


## Run GENIE3
library(GENIE3)
library(randomForest)
library(pROC)
set.seed(123) # For reproducibility of results
## The algorithm outputs a matrix containing the weights of the putative regulatory links, 
## with higher weights corresponding to more likely regulatory links. 
# weightMatrix <- GENIE3(exprMatrix, regulators=paste("Gene", 1:5, sep=""))


# exprMatrix <- t(datanum)
exprMatrix <- t(datahatA)


# set.seed(123) 
# weightMatrix <- GENIE3(exprMatrix)


Matrix <- myGENIE3(exprMatrix)


## results
weightMatrix <- Matrix[[1]]
predMatrix <- t(Matrix[[2]])
meanAUC <- mean(Matrix[[3]])


## DyAcc
DyAccuary <- DyAccWeight(predMatrix)
DyAccuary


## edge link
linkListNum <- min(getLinkList(weightMatrix, reportMax = geneselect)[,3])
weightMatrix0 <- weightMatrix
weightMatrix0[weightMatrix0 < linkListNum] = 0
weightMatrix0[weightMatrix0 >= linkListNum] = 1


# 计算52条边的AUC ------------------------------------------------------------

## 0-无向
SIGN = 0
# SIGN = 1

## Gold standard GRN 
datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
# datapath0 <- paste("D:/E/博士/R_程序/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
adj_gold0 <- paste(datapath0," InSilicoSize",genenum,"-Yeast",filenum,"-adj .csv",sep="")
adj_gold <- as.matrix(read.csv(file = adj_gold0))


# no self-loof, direct
if(SIGN == 1){
  adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  adj_matrix <- weightMatrix0
}
# undirect
if(SIGN == 0){
  adj_gold[which(adj_gold != 0)] <- 1
  adj_matrix <- abs(weightMatrix0)
}
# adj_matrix[abs(adj_matrix) > linkListNum] = 1


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)
# AUC
# DyAccuary
# meanAUC


meanAUCAcc <- cbind(geneselect, meanAUC, DyAccuary, AUC)
colnames(meanAUCAcc) <- c("link", "MeanAUC", "DyAcc", "AUROC", "AUPR")


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Acc0 <- performance$Acc
Recall0 <- performance$Recall
Pre0 <- performance$Pre
FPR0 <- performance$FPR
Fmeasure0 <- performance$Fmeasure
Result0 <- cbind(Acc0,Recall0,Pre0,FPR0,Fmeasure0)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")



## output
output <- cbind(meanAUCAcc, Result0)
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_GENIE3/Yeast",filenum,"Node",genenum, sep=""))
# setwd(paste("D:/E/博士/R_程序/Boolean/Data/DREAM3result_GENIE3/Yeast",filenum,"Node",genenum, sep=""))

if(method == 2){
  write.csv(output, file = "ENet_output.csv")
}else if(method == 3){
  write.csv(meanAUCAccAll, file = "Lasso_output.csv")
}else if(method == 1){
  write.csv(meanAUCAccAll, file = "Ridge_output.csv")
}

print("GENIE3 method")
output

