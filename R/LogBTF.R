
## 2022.7.27 method 1-Ridge; 2-ENet; 3-Lasso 
## ON single-cell or bulk RNA-seq data, and real data
## 2022.11.28 Simply the code



##############################################################################################################################
##############################################################################################################################

# For single-cell (bulk RNA-seq) gene expression data Ecoli ---------------------------------------------------

rm(list = ls())

set.seed(123)
method <- 2
noDIAG = 1
genenum <- 10
filenum <- 1


## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/MethodFuction.R')



l <- 21
datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
# datapath <- paste("/Users/lilingyu/E/PhD/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Ecoli",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()


for (k in 1:run) {
  # k <- 1
  datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  
  ## Kmeans 
  library(BiTrinA)
  datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  ## Perturbation design matrix
  datahat <- PerMatrix(datahatA)
  
  ## split train and test data
  ## for train
  xglm <- as.data.frame(datahat)
  yglm <- datahatA[-1,]
  ## for validation
  xglm01 <- as.data.frame(datahatA[1:(n-1),])
  
  
  ## glmnet
  data <- datahatA
  num_time_points <- dim(data)[1]
  numGENES <- dim(data)[2]
  
  ## no scale 
  DISTANCE_matrix <- as.matrix(data)
  ## with scale
  # DISTANCE_matrix <- scale(as.matrix(data))
  
  
  ## penalty
  X_matrix <- DISTANCE_matrix[1:(num_time_points-1),]
  n <- dim(X_matrix)[1]
  p <- dim(X_matrix)[2]
  
  #LOOCV settings
  nfold <- dim(X_matrix)[1]
  
  ## LogBTF function 
  logBTFresults <- logBTFfunction(method=2, nfold, numGENES, X_matrix, DISTANCE_matrix)
  
  
  ## AUC
  AUCall0 <- logBTFresults[[1]]
  ## coef matrix
  coefall <- logBTFresults[[2]]
  ## Dynamics Acc
  ptrainall0 <- logBTFresults[[3]]
  DyAccuracy <- DyAcc(datahatA,ptrainall0)
  
  # Direct Performance Compare -----------------------------------------------------
  
  # SIGN = 0
  SIGN = 1
  
  ## Gold standard GRN 
  datapath1 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold1 <- paste(datapath1," InSilicoSize",genenum,"-Ecoli",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold1))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
    adj_matrix <- coefall[-1,]
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
    adj_matrix <- abs(coefall[-1,])
  }
  
  
  ## Final ranked list, AUROC and AUPR
  adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
  adjstand[adjstand <= 1e-5] = 0
  adj_matrix <- adjstand * sign(adj_matrix)
  ## ??NAN????0
  adj_matrix[is.na(adj_matrix)] <- 0
  # adj_matrix[adj_matrix > 0] = 1
  
  
  ## no zero coefficients and links
  numgene <- c()
  for (k in 1:p) {
    # i <- 9
    numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
  }
  
  
  ## Final ranked list, AUROC and AUPR
  library(pracma)
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
  AUROC <- AUCresult$AUROC
  AUPR <- AUCresult$AUPR
  AUC <- cbind(AUROC,AUPR)
  AUCALL <- rbind(AUCALL, AUC)
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Acc <- performance$Acc
  Recall <- performance$Recall
  Pre <- performance$Pre
  FPR <- performance$FPR
  Fmeasure <- performance$Fmeasure
  Result <- cbind(Acc,Recall,Pre,FPR,Fmeasure)
  ResultAll <- rbind(ResultAll, Result)
  
  
  meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene), mean(AUCall0), DyAccuracy)
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("min", "max", "links", "MeanAUC", "DyAcc")
  
  
  # Undirect Performance Compare -----------------------------------------------------
  
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
    adj_matrix <- coefall[-1,]
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
    adj_matrix <- abs(coefall[-1,])
  }
  
  
  ## Final ranked list, AUROC and AUPR
  adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
  adjstand[adjstand <= 1e-5] = 0
  adj_matrix <- adjstand * sign(adj_matrix)
  ## ??NAN????0
  adj_matrix[is.na(adj_matrix)] <- 0
  
  
  library(pracma)
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
  AUC0 <- cbind(AUROC,AUPR,AUCresult$AUROC,AUCresult$AUPR)
  AUCALL0 <- rbind(AUCALL0, AUC0)
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Result0 <- cbind(Acc,Recall,Pre,FPR,Fmeasure,
                   performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
  ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
  colnames(ResultAll0) <- c("Acc","Recall","Pre","FPR","Fmeasure",
                            "Acc0","Recall0","Pre0","FPR0","Fmeasure0")
  
}


## output
output <- cbind(meanAUCAccAll, AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result/Ecoli",filenum,"Node",genenum, sep=""))
# setwd(paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM3result/Ecoli",filenum,"Node",genenum, sep=""))

if(method == 2){
  write.csv(output, file = "ENet_output.csv")
}else if(method == 3){
  write.csv(meanAUCAccAll, file = "Lasso_output.csv")
}else if(method == 1){
  write.csv(meanAUCAccAll, file = "Ridge_output.csv")
}

print("method=");method
output




##############################################################################################################################
##############################################################################################################################

# For single-cell (bulk RNA-seq) gene expression data Yeast ---------------------------------------------------

rm(list = ls())

# options(digits = 7)

##
set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
genenum <- 10
filenum <- 2


## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/LogBTFmainfunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/LogBTFmainfunction.R')

## load data
l <- 21
# datapath <- paste("/Users/lilingyu/E/PhD/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Yeast",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()


for (k in 1:run) {
  # k <- 1
  datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  
  
  ## Kmeans ??????
  library(BiTrinA)
  ## features*times - A n x m matrix comprising m raw measurements of n features
  datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  ## Perturbation design matrix
  datahat <- PerMatrix(datahatA)
  
  ## split train and test data
  ## for train
  xglm <- as.data.frame(datahat)
  yglm <- datahatA[-1,]
  ## for validation
  xglm01 <- as.data.frame(datahatA[1:(n-1),])
  
  
  ## glmnet
  data <- datahatA
  num_time_points <- dim(data)[1]
  numGENES <- dim(data)[2]
  
  
  ## no scale 
  DISTANCE_matrix <- as.matrix(data)
  ## with scale
  # DISTANCE_matrix <- scale(as.matrix(data))
  
  
  ## penalty
  X_matrix <- DISTANCE_matrix[1:(num_time_points-1),]
  n <- dim(X_matrix)[1]
  p <- dim(X_matrix)[2]
  
  
  #LOOCV settings
  nfold <- dim(X_matrix)[1]
  
  ## LogBTF function 
  logBTFresults <- logBTFfunction(method=2, nfold, numGENES, X_matrix, DISTANCE_matrix)
  
  
  ## AUC
  AUCall0 <- logBTFresults[[1]]
  ## coef matrix
  coefall <- logBTFresults[[2]]
  ## Dynamics Acc
  ptrainall0 <- logBTFresults[[3]]
  DyAccuracy <- DyAcc(datahatA,ptrainall0)
  
  
  # Direct Performance Compare -----------------------------------------------------
  
  ## 0-????
  # SIGN = 0
  SIGN = 1
  
  ## Gold standard GRN 
  datapath1 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  # datapath1 <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold1 <- paste(datapath1," InSilicoSize",genenum,"-Yeast",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold1))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
    adj_matrix <- coefall[-1,]
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
    adj_matrix <- abs(coefall[-1,])
  }
  
  
  ## Final ranked list, AUROC and AUPR
  adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
  adjstand[adjstand <= 1e-5] = 0
  adj_matrix <- adjstand * sign(adj_matrix)
  ## ??NAN????0
  adj_matrix[is.na(adj_matrix)] <- 0
  
  
  ## no zero coefficients and links
  numgene <- c()
  for (k in 1:p) {
    numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
  }
  
  
  ## Final ranked list, AUROC and AUPR
  library(pracma)
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
  AUROC <- AUCresult$AUROC
  AUPR <- AUCresult$AUPR
  AUC <- cbind(AUROC,AUPR)
  AUC
  AUCALL <- rbind(AUCALL, AUC)
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Acc <- performance$Acc
  Recall <- performance$Recall
  Pre <- performance$Pre
  FPR <- performance$FPR
  Fmeasure <- performance$Fmeasure
  Result <- cbind(Acc,Recall,Pre,FPR,Fmeasure)
  ResultAll <- rbind(ResultAll, Result)
  
  
  meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene), mean(AUCall0), DyAccuracy)
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("min", "max", "links", "MeanAUC", "DyAcc")
  
  
  # Undirect Performance Compare -----------------------------------------------------
  
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
    adj_matrix <- coefall[-1,]
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
    adj_matrix <- abs(coefall[-1,])
  }
  
  ## Final ranked list, AUROC and AUPR
  adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
  adjstand[adjstand <= 1e-5] = 0
  adj_matrix <- adjstand * sign(adj_matrix)
  ## ??NAN????0
  adj_matrix[is.na(adj_matrix)] <- 0
  
  
  library(pracma)
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
  AUC0 <- cbind(AUROC,AUPR,AUCresult$AUROC,AUCresult$AUPR)
  AUCALL0 <- rbind(AUCALL0, AUC0)
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Result0 <- cbind(Acc,Recall,Pre,FPR,Fmeasure,
                   performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
  ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
  colnames(ResultAll0) <- c("Acc","Recall","Pre","FPR","Fmeasure",
                            "Acc0","Recall0","Pre0","FPR0","Fmeasure0")
  
}


## output
output <- cbind(meanAUCAccAll, AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result/Yeast",filenum,"Node",genenum, sep=""))
# setwd(paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM3result/Yeast",filenum,"Node",genenum, sep=""))

if(method == 2){
  write.csv(output, file = "ENet_output.csv")
}else if(method == 3){
  write.csv(meanAUCAccAll, file = "Lasso_output.csv")
}else if(method == 1){
  write.csv(meanAUCAccAll, file = "Ridge_output.csv")
}

print("method=");method
output





##############################################################################################################################
##############################################################################################################################


# For real RNA-seq data2 or data5 ---------------------------------------------------

## 2022.11.28 LogBTF method, for real data
## according to dataDREAMGlmRegAdjRealdata10fold.R


rm(list = ls())
time_start<-Sys.time()


## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/LogBTFmainfunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/LogBTFmainfunction.R')

set.seed(123)
method <- 3
noDIAG = 1
# dataset <- 2
# dataset <- 5
dataset <- 16
methodname <- "LogBTFs"


## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")
# setwd("/Users/lilingyu/E/PhD/Matlab/Boolean/GRISLI")

# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="")  # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
n <- dim(datanum)[1]
p <- dim(datanum)[2]


## Boolean data
datanum[datanum > 0] <- 1
datahatA <- datanum


## Perturbation design matrix
datahat <- PerMatrix(datahatA)


## split train and test data
## for train
xglm <- as.data.frame(datahat)
yglm <- datahatA[-1,]
## for validation
xglm01 <- as.data.frame(datahatA[1:(n-1),])


## glmnet
data <- datahatA
num_time_points <- dim(data)[1]
numGENES <- dim(data)[2]


## no scale 
DISTANCE_matrix <- as.matrix(data)
## with scale
# DISTANCE_matrix <- scale(as.matrix(data))


## penalty
X_matrix <- DISTANCE_matrix[1:(num_time_points-1),]
n <- dim(X_matrix)[1]
p <- dim(X_matrix)[2]


#LOOCV settings
nfold <- dim(X_matrix)[1]


## LogBTF function 
logBTFresults <- logBTFfunction(method=2, nfold, numGENES, X_matrix, DISTANCE_matrix)


## AUC
AUCall0 <- logBTFresults[[1]]
## coef matrix
coefall <- logBTFresults[[2]]
## Dynamics Acc
ptrainall0 <- logBTFresults[[3]]
DyAccuracy <- DyAcc(datahatA,ptrainall0)


## Performance Compare 

## 0 -- undirect
SIGN = 0
# SIGN = 1

## Gold standard GRN 
fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target

# no self-loof, direct
if(SIGN == 1){
  adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  adj_matrix <- coefall[-1,]
}
# undirect
if(SIGN == 0){
  adj_gold[which(adj_gold != 0)] <- 1
  adj_matrix <- abs(coefall[-1,])
}


adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
adjstand[adjstand <= 1e-5] = 0
adj_matrix <- adjstand * sign(adj_matrix)
adj_matrix[is.na(adj_matrix)] <- 0


## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}
numgene[is.na(numgene)] <- 0


library(pracma)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)

meanAUCAcc <- cbind(min(numgene[numgene!=min(numgene)]), max(numgene), sum(numgene), mean(AUCall0), DyAccuracy, AUC)
meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
colnames(meanAUCAccAll) <- c("min", "max", "links", "MeanAUC", "DyAcc", "AUROC", "AUPR")


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output
output <- cbind(meanAUCAccAll, Result0)
output


# datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
# fdata <- paste(datapath00,methodname,".csv",sep="") 
# write.csv(output, fdata)


## times 
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))



##############################################################################################################################
##############################################################################################################################


# For real RNA-seq data16 ---------------------------------------------------

## 2022.11.28 LogBTF method, for real data
## according to dataDREAMGlmRegAdjRealdata10fold.R

rm(list = ls())
time_start<-Sys.time()


## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/LogBTFmainfunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
# source('/Users/lilingyu/E/PhD/R/Boolean/R/LogBTFmainfunction.R')

set.seed(123)
method <- 2
noDIAG = 1
# dataset <- 2
# dataset <- 5
dataset <- 16
methodname <- "LogBTFs"


## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")
# setwd("/Users/lilingyu/E/PhD/Matlab/Boolean/GRISLI")

# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="")  # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
n <- dim(datanum)[1]
p <- dim(datanum)[2]

## Kmeans 
library(BiTrinA)
## features*times - A n x m matrix comprising m raw measurements of n features
datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])


## Perturbation design matrix
datahat <- PerMatrix(datahatA)


## split train and test data
## for train
xglm <- as.data.frame(datahat)
yglm <- datahatA[-1,]
## for validation
xglm01 <- as.data.frame(datahatA[1:(n-1),])


## glmnet
data <- datahatA
num_time_points <- dim(data)[1]
numGENES <- dim(data)[2]


## no scale 
DISTANCE_matrix <- as.matrix(data)
## with scale
# DISTANCE_matrix <- scale(as.matrix(data))


## penalty
X_matrix <- DISTANCE_matrix[1:(num_time_points-1),]
n <- dim(X_matrix)[1]
p <- dim(X_matrix)[2]


#LOOCV settings
# nfold <- dim(X_matrix)[1]
nfold <- 10


## LogBTF function 
logBTFresults <- logBTFfunction(method=2, nfold, numGENES, X_matrix, DISTANCE_matrix)


## AUC
AUCall0 <- logBTFresults[[1]]
## coef matrix
coefall <- logBTFresults[[2]]
## Dynamics Acc
ptrainall0 <- logBTFresults[[3]]
DyAccuracy <- DyAcc(datahatA,ptrainall0)


## Performance Compare 

## 0 -- undirect
SIGN = 0
# SIGN = 1

## Gold standard GRN 
fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target

# no self-loof, direct
if(SIGN == 1){
  adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  adj_matrix <- coefall[-1,]
}
# undirect
if(SIGN == 0){
  adj_gold[which(adj_gold != 0)] <- 1
  adj_matrix <- abs(coefall[-1,])
}


adjstand <- apply(abs(adj_matrix), 2, function(x){(x-min(x))/(max(x)-min(x))} )
adjstand[adjstand <= 1e-5] = 0
adj_matrix <- adjstand * sign(adj_matrix)
adj_matrix[is.na(adj_matrix)] <- 0


## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}
numgene[is.na(numgene)] <- 0


library(pracma)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)

meanAUCAcc <- cbind(min(numgene[numgene!=min(numgene)]), max(numgene), sum(numgene), mean(AUCall0), DyAccuracy, AUC)
meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
colnames(meanAUCAccAll) <- c("min", "max", "links", "MeanAUC", "DyAcc", "AUROC", "AUPR")


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output
output <- cbind(meanAUCAccAll, Result0)
output


# datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
# fdata <- paste(datapath00,methodname,".csv",sep="") 
# write.csv(output, fdata)


## times 
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))






##############################################################################################################################
##############################################################################################################################
