## 2022.6.2 比较glm(noise)和penalty glm(leave one)的coef
## 2022.7.19 精简code，只保留glm和RegLog
## 2022.7.30 单细胞数据，比较稀疏，二值化--结果一般 0.506  0.504， 尝试直接非0设为1

# Glmnet + Glm noise ------------------------------------------------------------------

rm(list = ls())

time_start<-Sys.time()

# options(digits = 7)

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')


set.seed(123)
## 
method <- 2
## 0-不排除对角线
noDIAG = 1
# data 2 or 3
dataset <- 3
methodname <- "LogBTFs"


# LogBTFs -----------------------------------------------------------------

## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")


# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="") # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
n <- dim(datanum)[1]
p <- dim(datanum)[2]


# sum(datanum != 0)
# sum(datahatA != 0)
# View(datanum[,1:10])
# View(datahatA[,1:10])


datanum[datanum > 0] <- 1
# adj_matrix[adj_matrix > 0] = 1
# sum(datanum-datahatA)
datahatA <- datanum

## Kmeans ??????
# library(BiTrinA)
# ## features*times - A n x m matrix comprising m raw measurements of n features
# # binarizeMatrix(t(datanum))
# datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
# # write.csv(t(datahatA), "CoefAllSize10\\datadream3p1_matrix.csv")


## noise  ??X?????е?1??????
set.seed(2022)
res <- datahatA[1:(n-1),]
noiseLevel = 1e-5
res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))
# res[1,]

## when 1, it add noise
for (i in 1:(dim(res)[1])) {
  # i <- 1
  for (j in 1:dim(res)[2]) {
    # j <- 2
    if (datahatA[i,j] == 0)
      res[i,j] <- 0
  }
}

datahat <- res
# colnames(datahat) <- c("G1","G2","G3","G4","G5","G6","G7","G8", "G9", "G10")


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


#Generate Y and X_matrix for glmnet
if(method==1){
  alphas <- 0
}else if(method==2){
  alphas <- 0.5  
}else if(method==3){
  alphas <- 1
}else if(method==4){
  alphas <- seq(0,1,0.1)
}else{
  input <- readline(' *** Please input manually the alpha values (between 0 and 1) separated by comma: ')
  alphas <- as.numeric(unlist(strsplit(input,',')))
}

## no scale 
DISTANCE_matrix <- as.matrix(data)
## with scale
# DISTANCE_matrix <- scale(as.matrix(data))


## penalty
X_matrix <- DISTANCE_matrix[1:(num_time_points-1),]
n <- dim(X_matrix)[1]
p <- dim(X_matrix)[2]

# datascale <- scale(DISTANCE_matrix)
# datascale <- DISTANCE_matrix
# DISTANCE_matrix01 <- pmax(sign(datascale), 0)
# X_matrix01 <- DISTANCE_matrix01[1:(num_time_points-1),]

#LOOCV settings
nfold <- dim(X_matrix)[1]
foldid <- 1:nfold
keep <- TRUE
pred_lambda_min <- matrix(0, nrow = numGENES+1, ncol = numGENES)
lambda_res <- vector()
alpha_res <- vector()
AUCall0 <- c()
ptrainall0 <- c()
meanAUCAccAll <- c()

# options(digits = 3)
library(glmnet)
library(pROC)
for (gi in 1:numGENES) {
  # gi <- 59
  
  AUCall <- c()
  ptrainall <- c()
  cverrorall <- c()
  
  lambda <-  vector()
  cvERROR <-  vector()
  beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
  theta <- matrix(data=0,nrow = dim(X_matrix)[2]+1,ncol = length(alphas))
  
  
  # for (test in 1:length(alphas)) {
  test <- 1
  Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
  
  # if Y exist one 1/0, use noise 0/1 data.
  # if Y exist one 1/0, use noise 0/1 data.
  if(sum(Y_vector) == 0 | sum(Y_vector) == n){
    glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
    coef <- glm.fit$coefficients
    # # coef <- round(glm.fit$coefficients,4)
    coef[is.na(coef)] <- 0
    pred_lambda_min[,gi] <- coef
    
    # pred <- predict(glm.fit, xglm, type = "response")
    # predall <- cbind(predall, pred)
    
    ptrain <- Y_vector    ## 不进化
    ptrainall <- cbind(ptrainall, ptrain)
    ptrainall0 <- cbind(ptrainall0, ptrainall)
    # ptrain <- sgn(as.matrix(xglm01) %*% coef[-1] + coef[1])
    # ptrainall <- cbind(ptrainall, ptrain)
    # ptrainall0 <- cbind(ptrainall0, ptrainall)
    
    # aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
    # auc <- aucplot$auc
    # AUCall <- cbind(AUCall, auc)
    # AUCall0 <- cbind(AUCall0, AUCall)
  }else if(sum(Y_vector) == 1 | sum(Y_vector) == (n-1)){
    glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
    coef <- glm.fit$coefficients
    # coef <- round(glm.fit$coefficients,4)
    # coef[is.na(coef)] <- 0
    pred_lambda_min[,gi] <- coef
    
    # pred <- predict(glm.fit, xglm, type = "response")
    # predall <- cbind(predall, pred)
    
    ptrain <- sgn(as.matrix(xglm01) %*% coef[-1] + coef[1])
    ptrainall <- cbind(ptrainall, ptrain)
    ptrainall0 <- cbind(ptrainall0, ptrainall)
    
    aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
    auc <- aucplot$auc
    AUCall <- cbind(AUCall, auc)
    AUCall0 <- cbind(AUCall0, AUCall)
    
  }else{
    
    ## glnnet with positive coef
    # if(noDIAG==1){
    #   CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
    #                           keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    # }else{
    #   CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
    #                           keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    # }
    
    ## glmnet
    if(noDIAG==1){
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,
                              nfolds = nfold, foldid = foldid, keep = keep, grouped = FALSE)
    }else{
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],
                              nfolds = nfold, foldid = foldid, keep = keep, grouped = FALSE)
    }
    
    plot(CV_results)
    lambda[test] <- CV_results$lambda.min
    cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
    cverrorall <- cbind(cverrorall, cvERROR)
    # coef.CV_results <- round(coef(CV_results, s='lambda.min'),3)
    coef.CV_results <- coef(CV_results, s='lambda.min')
    ##  
    beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    theta[coef.CV_results@i+1,test] = coef.CV_results@x
    
    
    theta[1,test] <- lambda*theta[1,test]
    # theta[1,test] <- theta[1,test]/lambda
    ptrain <- sgn(as.matrix(X_matrix) %*% theta[-1,test] + theta[1,test])
    # ptrain <- sgn(as.matrix(X_matrix) %*% theta[-1,test] + theta[1,test])
    # ptrain <- as.matrix(X_matrix) %*% theta[-1,test] + theta[1,test]
    ptrainall <- cbind(ptrainall, ptrain)
    # View(cbind(as.matrix(X_matrix) %*% theta[-1,test] + theta[1,test],Y_vector))
    
    aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
    auc <- aucplot$auc
    AUCall <- cbind(AUCall, auc)
    
    
    minIdx <- max(which(cvERROR==min(cvERROR)))
    lambda_res[gi] <- lambda[minIdx]
    alpha_res[gi] <- alphas[minIdx]
    # pred_lambda_min[,gi] <- beta[,minIdx]
    pred_lambda_min[,gi] <- theta[,minIdx]
    
    ptrainall0 <- cbind(ptrainall0, as.matrix(ptrainall[,minIdx]))
    AUCall0 <- cbind(AUCall0, AUCall[,minIdx])
  }
  
  print(gi)
}

mean(AUCall0)

## coef
coefall <- pred_lambda_min

## Dynamics Acc
DyAccuracy <- DyAcc(datahatA,ptrainall0)
# DyAccuracy


# Performance Compare -----------------------------------------------------

## 0-无向
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
## 把NAN变成0
adj_matrix[is.na(adj_matrix)] <- 0
# adj_matrix[adj_matrix > 0] = 1
# View(adj_matrix[,1:10])


## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  # i <- 9
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}

numgene[is.na(numgene)] <- 0


# library(ROCR)
# library(tigress)
# library(zinbwave)
# 
# pred <- prediction(abs(adj_matrix[-diagind]), adj_gold[-diagind])
# auc <- performance(pred, measure = "auc")@y.values
# auc


library(pracma)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)
# AUC
# DyAccuracy
# mean(AUCall0)


# meanAUCAcc <- cbind(sort(numgene)[2], max(numgene), sum(numgene), mean(AUCall0), DyAccuracy, AUC)
meanAUCAcc <- cbind(min(numgene[numgene!=min(numgene)]), max(numgene), sum(numgene), mean(AUCall0), DyAccuracy, AUC)
meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
colnames(meanAUCAccAll) <- c("min", "max", "links", "MeanAUC", "DyAcc", "AUROC", "AUPR")



performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Acc0 <- performance$Acc
Recall0 <- performance$Recall
Pre0 <- performance$Pre
FPR0 <- performance$FPR
Fmeasure0 <- performance$Fmeasure
Result0 <- cbind(Acc0,Recall0,Pre0,FPR0,Fmeasure0)
colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")


## output
output <- cbind(meanAUCAccAll, Result0)
output


datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
fdata <- paste(datapath00,methodname,".csv",sep="") 
write.csv(output, fdata)


## times 
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))

