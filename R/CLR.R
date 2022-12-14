## CLR method 
## single-cell data are based on bulk data generated by GNW (Ecoli and Yeast)
## single-cell data has many zero, which can be see from their distribution figure
## For bulk data, get binarized data using k-means, then all zreo are zero, all no-zeros are original value, as simulated single-cell 
## For real data, the code only need to chage the input path and file, here we show all codes


##############################################################################################################################
##############################################################################################################################

## 2022.11.28 CLR method for single-cell data Ecoli

rm(list = ls())

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/MethodFuction.R')

##
set.seed(123)


## parameters
# file <- c(10,1,4,45)
# file <- c(10,2,3,48)
# file <- c(50,1,12,576)
# file <- c(50,2,6,585)
# file <- c(100,1,29,1708)
file <- c(100,2,11,1996)

method = 2
genenum <- file[1]
filenum <- file[2]
k <- file[3]
geneselect <- file[4]


l <- 21

# datapath <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Ecoli",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
n <- dim(datanum)[1]
p <- dim(datanum)[2]


## Kmeans for single-cell data generation ##########################
library(BiTrinA)                                                   #
datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])     #
datahatA <- as.matrix(datanum) ** as.matrix(datahat)               #
datahatA[which(datahatA == 1)] <- 0                                #
datanum <- datahatA                                                #
####################################################################
## If using the simulated bulk gene expression data, 
## skip the above Box "Kmeans for single-cell data generation"


## expression matrix : sample * genes
exprMatrix <- datanum
# exprMatrix <- datahatA



# CLR -----------------------------------------------------------------
library(minet)
## input data -- sample*genes
set.seed(123)
mim <- build.mim(exprMatrix,estimator="spearman")
weightMatrix <- clr(mim)


## edge link
library(GENIE3)
linkListNum <- min(getLinkList(weightMatrix, reportMax = geneselect)[,3])
weightMatrix0 <- weightMatrix
weightMatrix0[weightMatrix0 < linkListNum] = 0
weightMatrix0[weightMatrix0 >= linkListNum] = 1


## 0-????
SIGN = 0
# SIGN = 1

## Gold standard GRN 
datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
# datapath0 <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
adj_gold0 <- paste(datapath0," InSilicoSize",genenum,"-Ecoli",filenum,"-adj .csv",sep="")
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


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUC <- cbind(AUCresult$AUROC,AUCresult$AUPR)

meanAUCAcc <- cbind(geneselect, AUC)
colnames(meanAUCAcc) <- c("link", "AUROC", "AUPR")

performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output
output <- cbind(meanAUCAcc, Result0)
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_CLR01/Ecoli",filenum,"Node",genenum, sep=""))
# setwd(paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM3result_ARACNE/Ecoli",filenum,"Node",genenum, sep=""))

if(method == 2){
  write.csv(output, file = "ENet_output.csv")
}else if(method == 3){
  write.csv(meanAUCAccAll, file = "Lasso_output.csv")
}else if(method == 1){
  write.csv(meanAUCAccAll, file = "Ridge_output.csv")
}

print("CLR method")
output


##############################################################################################################################
##############################################################################################################################

## 2022.11.28 CLR method for single-cell data -- Yeast

rm(list = ls())

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/MethodFuction.R')


set.seed(123)


## parameters
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



l <- 21
# datapath <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Yeast",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
n <- dim(datanum)[1]
p <- dim(datanum)[2]


## Kmeans for single-cell data generation ##########################
library(BiTrinA)                                                   #
datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])     #
datahatA <- as.matrix(datanum) ** as.matrix(datahat)               #
datahatA[which(datahatA == 1)] <- 0                                #
datanum <- datahatA                                                #
#####################################################################
## If using the simulated bulk gene expression data, 
## skip the above Box "Kmeans for single-cell data generation"


## expression matrix : sample * genes
exprMatrix <- datanum
# exprMatrix <- datahatA


# CLR -----------------------------------------------------------------
library(minet)
## input data -- sample*genes
set.seed(123)
mim <- build.mim(exprMatrix,estimator="spearman")
weightMatrix <- clr(mim)


## edge link
library(GENIE3)
linkListNum <- min(getLinkList(weightMatrix, reportMax = geneselect)[,3])
weightMatrix0 <- weightMatrix
weightMatrix0[weightMatrix0 < linkListNum] = 0
weightMatrix0[weightMatrix0 >= linkListNum] = 1


## 0-????
SIGN = 0
# SIGN = 1

## Gold standard GRN 
datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
# datapath0 <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
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


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUC <- cbind(AUCresult$AUROC,AUCresult$AUPR)

meanAUCAcc <- cbind(geneselect, AUC)
colnames(meanAUCAcc) <- c("link", "AUROC", "AUPR")

performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output the single-cell results
output <- cbind(meanAUCAcc, Result0)
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_CLR01/Yeast",filenum,"Node",genenum, sep=""))


####################################################################################################
## out the bulk data results                                  
# setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_CLR/Ecoli",filenum,"Node",genenum, sep=""))  
####################################################################################################


if(method == 2){
  write.csv(output, file = "ENet_output.csv")
}else if(method == 3){
  write.csv(meanAUCAccAll, file = "Lasso_output.csv")
}else if(method == 1){
  write.csv(meanAUCAccAll, file = "Ridge_output.csv")
}

print("CLR method")
output



##############################################################################################################################
##############################################################################################################################

## 2022.11.28 CLR method for real data5 or data16

rm(list = ls())

time_start<-Sys.time()

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/MethodFuction.R')

set.seed(123)

## parameters
method <- 2
noDIAG = 1
methodname <- "CLR"

## data 5 or 16
# file <- c(2,1361)
# file <- c(5,2074)
file <- c(16,36003)
dataset <- file[1]
geneselect <- file[2]


## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")

# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="")  # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
n <- dim(datanum)[1]
p <- dim(datanum)[2]


## expression matrix : sample * genes
exprMatrix <- datanum
# exprMatrix <- datahatA



# CLR -----------------------------------------------------------------
library(minet)
## input data -- sample*genes
set.seed(123)
mim <- build.mim(exprMatrix,estimator="spearman")
weightMatrix <- clr(mim)


## edge link
library(GENIE3)
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
AUC <- cbind(AUCresult$AUROC,AUCresult$AUPR)
meanAUCAcc <- cbind(geneselect, AUC)
colnames(meanAUCAcc) <- c("link", "AUROC", "AUPR")


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output
output <- cbind(meanAUCAcc, Result0)
output


# datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
# fdata <- paste(datapath00,methodname,".csv",sep="") 
# write.csv(output, fdata)

## times 
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))


