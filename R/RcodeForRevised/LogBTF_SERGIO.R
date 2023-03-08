## 2022.7.27 method 1-Ridge; 2-ENet; 3-Lasso
## 2023.2.12 LogBTF method using SERGIO simulated data
## 2023.2.11 based on code dataDREAMGlmRegAdj_Ecoli



# LogBTF ------------------------------------------------------------------

rm(list = ls())


## set pathway
# pathway <- '/home/lly/'
pathway <- '/Users/lilingyu/E/PhD/'


##
set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
epsilon <- 1e-3


dataclass <- 'count_matrix'
# dataclass <- 'exper_clean'
genenum <- 100
celltyp <- 2
cellbum <- 10


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))


# LogBTF -----------------------------------------------------------------




## load data
datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", cellbum, ".csv",sep="")
Data = t(read.csv(file = dataauto, header=F))
Data[2,1]


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()


## begin our test
datanum <- Data 
n <- dim(datanum)[1]
p <- dim(datanum)[2]


# ## Kmeans ??????
# library(BiTrinA)
# ## features*times - A n x m matrix comprising m raw measurements of n features
# datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
# # write.csv(t(datahatA), "CoefAllSize10\\datadream3p1_matrix.csv")
datahatA <- Data
datahatA[datahatA>0] <-1 


## noise  ??X????????1??????
set.seed(2022)
res <- datahatA[1:(n-1),]
noiseLevel = epsilon
res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))


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


#LOOCV settings
nfold <- dim(X_matrix)[1]
foldid <- 1:nfold
keep <- TRUE
pred_lambda_min <- matrix(0, nrow = numGENES+1, ncol = numGENES)
lambda_res <- vector()
alpha_res <- vector()
AUCall0 <- c()
ptrainall0 <- c()


library(glmnet)
library(pROC)
for (gi in 1:numGENES) {
  # gi <- 1
  
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
  if(sum(Y_vector) == 0 | sum(Y_vector) == n){
    glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
    coef <- glm.fit$coefficients
    # # coef <- round(glm.fit$coefficients,4)
    coef[is.na(coef)] <- 0
    pred_lambda_min[,gi] <- coef
    
    ptrain <- Y_vector    ## ??????
    ptrainall <- cbind(ptrainall, ptrain)
    ptrainall0 <- cbind(ptrainall0, ptrainall)

  }else if(sum(Y_vector) == 1 | sum(Y_vector) == (n-1)){
    
    glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
    coef <- glm.fit$coefficients
    # coef <- round(glm.fit$coefficients,4)
    coef[is.na(coef)] <- 0
    pred_lambda_min[,gi] <- coef

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
    ptrainall <- cbind(ptrainall, ptrain)

    aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
    auc <- aucplot$auc
    AUCall <- cbind(AUCall, auc)
    
    
    minIdx <- max(which(cvERROR==min(cvERROR)))
    lambda_res[gi] <- lambda[minIdx]
    alpha_res[gi] <- alphas[minIdx]
    pred_lambda_min[,gi] <- theta[,minIdx]
    
    ptrainall0 <- cbind(ptrainall0, as.matrix(ptrainall[,minIdx]))
    AUCall0 <- cbind(AUCall0, AUCall[,minIdx])
  }
  
  print(gi)
}

mean(AUCall0)
coefall <- pred_lambda_min


## Dynamics Acc
DyAccuracy <- DyAcc(datahatA,ptrainall0)


# Direct Performance Compare -----------------------------------------------------

## 0-????
# SIGN = 0
SIGN = 1

## load ground-truth Adj matrix
setwd(datapath)
adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix.csv"))


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
adjstand[adjstand <= epsilon] = 0
adj_matrix <- adjstand * sign(adj_matrix)
## ??NAN????0
adj_matrix[is.na(adj_matrix)] <- 0
# adj_matrix[adj_matrix > 0] = 1


## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  # k <- 13
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}

## delete NA value, using numeric class dtypr ata
numgene1 <- as.numeric(numgene)
numgene1 <- na.omit(numgene1)
numgene <- numgene1
# sum(numgene)


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

## 0-????
SIGN = 0
# SIGN = 1

## Gold standard GRN 
## load ground-truth Adj matrix
adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix.csv"))


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
adjstand[adjstand <= epsilon] = 0
adj_matrix <- adjstand * sign(adj_matrix)
## ??NAN????0
adj_matrix[is.na(adj_matrix)] <- 0
# adj_matrix[adj_matrix > 0] = 1


library(pracma)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC0 <- AUCresult$AUROC
AUPR0 <- AUCresult$AUPR
AUC0 <- cbind(AUROC,AUPR,AUROC0,AUPR0)
AUC0
AUCALL0 <- rbind(AUCALL0, AUC0)


performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
Acc0 <- performance$Acc
Recall0 <- performance$Recall
Pre0 <- performance$Pre
FPR0 <- performance$FPR
Fmeasure0 <- performance$Fmeasure
Result0 <- cbind(Acc,Recall,Pre,FPR,Fmeasure,Acc0,Recall0,Pre0,FPR0,Fmeasure0)
ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
colnames(ResultAll0) <- c("Acc","Recall","Pre","FPR","Fmeasure",
                          "Acc0","Recall0","Pre0","FPR0","Fmeasure0")

  
## output
output <- t(as.matrix(c(meanAUCAccAll[,1:5], AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])))
output


## save
# setwd(paste(pathway,"R/Boolean/RShell/Data/DREAM3result_LogBTF/",sep=""))
# 
# if(method == 2){
#   write.csv(output, file = paste("ENet_output_SERGIO_",genenum,"_",celltyp, "_", cellbum,".csv",sep=""))
# }else if(method == 3){
#   write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",genenum,"_",celltyp, "_", cellbum,".csv",sep=""))
# }else if(method == 1){
#   write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",genenum,"_",celltyp, "_", cellbum,".csv",sep=""))
# }

print("method=");method

