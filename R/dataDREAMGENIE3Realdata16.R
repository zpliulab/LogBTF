
rm(list = ls())

time_start<-Sys.time()

# options(digits = 7)

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')


set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
methodname <- "GENIE3"
## data 2 or 3
# file <- c(2,1361)
file <- c(16,36003)
dataset <- file[1]
geneselect <- file[2]


## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")


# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="") # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
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


## 0-????
SIGN = 0
# SIGN = 1


## Gold standard GRN 
fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target
# diag(adj_gold) <- 0 # remove self-regulation
# ntf <- dim(X)[1] # Number of genes (total)
# diagind <- seq(1,ntf*ntf,ntf+1) # indices of the non-diagonal elements in trueA
# # View(adj_gold[1:10,1:10])


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


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)
# AUC
# DyAccuracy
# mean(AUCall0)


# meanAUCAcc <- cbind(geneselect, AUC)
# colnames(meanAUCAcc) <- c("link", "AUROC", "AUPR")
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
output


datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
fdata <- paste(datapath00,methodname,".csv",sep="") 
write.csv(output, fdata)


## times 
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))

