
## 2022.7.27 method 1-Ridge; 2-ENet; 3-Lasso 
## ON single-cell or bulk RNA-seq data, and real data
## 2022.11.28 Simply the code



##############################################################################################################################
##############################################################################################################################

# For single-cell gene expression data Ecoli ---------------------------------------------------

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
dataauto <- paste(datapath,"InSilicoSize",genenum,"-Ecoli",filenum,"-nonoise-trajectories.tsv",sep="")
Data = as.matrix(read.table(file = dataauto, header=T))
run <- dim(Data)[1]/l


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()

for (k in 1:run) {
  
  # k<-11
  datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]

  
  ## Kmeans 
  library(BiTrinA)
  datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  datahatA <- as.matrix(datanum) ** as.matrix(datahat)
  datahatA[which(datahatA == 1)] <- 0
  datanum <- datahatA
  
  
  # later scale -------------------------------------------------------------
  data <- datanum
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
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  
  lambda_res <- vector()
  alpha_res <- vector()
  
  library(glmnet)
  for (gi in 1:numGENES) {
    # gi <- 1
    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
    
    for (test in 1:length(alphas)) {
      # test <- 1
      Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
      
      # if Y exist one 1/0, use noise 0/1 data.
      if(sum(Y_vector != 0) > 1){
        
        set.seed(123)
        
        if(noDIAG==1){
          CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                                  keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
        }else{
          CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                                  keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
        }
        lambda[test] <- CV_results$lambda.min
        cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
        # coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min') 
        # beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
        # 2022.5.14 Lingyu Li corrected
        coef.CV_results <- coef(CV_results, s='lambda.min')
        beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
        # coef.CV_results <- predict(CV_results, X_matrix, s='lambda.min')
        # beta[coef.CV_results[-1],test] = coef.CV_results[-1]
      }else{
        
        set.seed(123)
        glm.fit <- lm(Y_vector~X_matrix, data.frame(X_matrix))
        coef <- glm.fit$coefficients
        coef[is.na(coef)] <- 0
        beta[,test] = coef[-1]
      }
      
      minIdx <- max(which(cvERROR==min(cvERROR)))
      lambda_res[gi] <- lambda[minIdx]
      alpha_res[gi] <- alphas[minIdx]
      pred_lambda_min[,gi] <- beta[,minIdx]
      
      print(gi)
      
    }
    
  }
  
  coef <- pred_lambda_min
  coef[is.na(coef)] <- 0
  
  # Direct Performance Compare -----------------------------------------------------

  # SIGN = 0
  SIGN = 1
  
  
  library(ppcor)
  if(SIGN==1){
    parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  coef[is.na(coef)]
  adj_matrix <- pred_lambda_min
  

  ## no zero coefficients and links
  numgene <- c()
  for (k in 1:p) {
    # i <- 9
    numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
  }
  numgene[is.na(numgene)] <- 0
  
  ## Gold standard GRN 
  SIGN = 1
  
  ## Gold standard GRN 
  datapath1 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold1 <- paste(datapath1," InSilicoSize",genenum,"-Ecoli",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold1))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
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
  
  
  # meanAUCAcc <- cbind(mean(AUCall0), DyAcc)
  # meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  # colnames(meanAUCAccAll) <- c("MeanAUC", "DyAcc")
  
  meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene))
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("min", "max", "links")
  
  
  # Undirect Performance Compare -----------------------------------------------------
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  library(ppcor)
  if(SIGN==1){
    parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  adj_matrix <- pred_lambda_min
  
  
  ## Gold standard GRN 
  datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold0 <- paste(datapath0," InSilicoSize",genenum,"-Ecoli",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold0))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
  }
  
  
  ## Final ranked list, AUROC and AUPR
  library(pracma)
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
  AUROC0 <- AUCresult$AUROC
  AUPR0 <- AUCresult$AUPR
  AUC0 <- cbind(AUROC,AUPR,AUROC0,AUPR0)
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
  
}


## output
output <- cbind(meanAUCAccAll, AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_SICNERITIES01/Ecoli",filenum,"Node",genenum, sep=""))


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

# For single-cell gene expression data Yeast ---------------------------------------------------

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


## load data
l <- 21
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
  
  # k<-42
  datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  
  ## Kmeans 
  library(BiTrinA)
  datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  datahatA <- as.matrix(datanum) ** as.matrix(datahat)
  datahatA[which(datahatA == 1)] <- 0
  datanum <- datahatA
  
  
  data <- datanum
  num_time_points <- dim(data)[1]
  numGENES <- dim(data)[2]
  
  
  # Generate Y and X_matrix for glmnet
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
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  
  lambda_res <- vector()
  alpha_res <- vector()
  
  library(glmnet)
  for (gi in 1:numGENES) {
    # gi <- 1
    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
    
    for (test in 1:length(alphas)) {
      # test <- 1
      Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
      
      # if Y exist one 1/0, use noise 0/1 data.
      if(sum(Y_vector != 0) > 1){
        
        set.seed(123)
        
        if(noDIAG==1){
          CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                                  keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
        }else{
          CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                                  keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
        }
        lambda[test] <- CV_results$lambda.min
        cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
        # coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min') 
        # beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
        # 2022.5.14 Lingyu Li corrected
        coef.CV_results <- coef(CV_results, s='lambda.min')
        beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
        # coef.CV_results <- predict(CV_results, X_matrix, s='lambda.min')
        # beta[coef.CV_results[-1],test] = coef.CV_results[-1]
      }else{
        
        set.seed(123)
        glm.fit <- lm(Y_vector~X_matrix, data.frame(X_matrix))
        coef <- glm.fit$coefficients
        coef[is.na(coef)] <- 0
        beta[,test] = coef[-1]
      }
      
      minIdx <- max(which(cvERROR==min(cvERROR)))
      lambda_res[gi] <- lambda[minIdx]
      alpha_res[gi] <- alphas[minIdx]
      pred_lambda_min[,gi] <- beta[,minIdx]
      
      print(gi)
      
      
    }
    
  }
  
  coef <- pred_lambda_min
  coef[is.na(coef)] <- 0
  
  
  # Direct Performance Compare -----------------------------------------------------
  # SIGN = 0
  SIGN = 1
  
  
  library(ppcor)
  if(SIGN==1){
    parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  adj_matrix <- pred_lambda_min
  
  
  ## no zero coefficients and links
  numgene <- c()
  for (k in 1:p) {
    # i <- 9
    numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
  }
  
  
  ## Gold standard GRN 
  SIGN = 1
  
  ## Gold standard GRN 
  datapath1 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold1 <- paste(datapath1," InSilicoSize",genenum,"-Yeast",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold1))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
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

  meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene))
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("min", "max", "links")
  
  
  # Undirect Performance Compare -----------------------------------------------------
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  library(ppcor)
  if(SIGN==1){
    parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  adj_matrix <- pred_lambda_min
  
  
  ## Gold standard GRN 
  datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold0 <- paste(datapath0," InSilicoSize",genenum,"-Yeast",filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold0))
  
  
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
  }
  
  
  ## Final ranked list, AUROC and AUPR
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
  
}


## output
output <- cbind(meanAUCAccAll, AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])
setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_SICNERITIES01/Yeast",filenum,"Node",genenum, sep=""))

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


# For real RNA-seq data5 ---------------------------------------------------

## 2022.11.28 LogBTF method, for real data
## according to dataDREAMGlmRegAdjRealdata10fold.R


rm(list = ls())
time_start<-Sys.time()


## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/LogBTFmainfunction.R')


set.seed(123)
method <- 2
noDIAG = 1
methodname <- "SINCERITIES"
## data 5 
file <- c(5,2074)
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


data <- datanum
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
pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)

lambda_res <- vector()
alpha_res <- vector()

library(glmnet)
for (gi in 1:numGENES) {
  # gi <- 1
  lambda <-  vector()
  cvERROR <-  vector()
  beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
  
  for (test in 1:length(alphas)) {
    # test <- 1
    Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
    if(noDIAG==1){
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                              keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    }else{
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                              keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    }
    lambda[test] <- CV_results$lambda.min
    cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
    # coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min') 
    # beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    # 2022.5.14 Lingyu Li corrected
    coef.CV_results <- coef(CV_results, s='lambda.min')
    beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    # coef.CV_results <- predict(CV_results, X_matrix, s='lambda.min')
    # beta[coef.CV_results[-1],test] = coef.CV_results[-1]
  } 
  
  minIdx <- max(which(cvERROR==min(cvERROR)))
  lambda_res[gi] <- lambda[minIdx]
  alpha_res[gi] <- alphas[minIdx]
  pred_lambda_min[,gi] <- beta[,minIdx]
  
  print(gi)
  
}

coef <- pred_lambda_min


# Direct Performance Compare -----------------------------------------------------

# SIGN = 0
SIGN = 0


## Gold standard GRN 
fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) 


library(ppcor)
if(SIGN==1){
  parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
  pred_lambda_min <- coef*sign(parcorr_matrix)
}
adj_matrix <- pred_lambda_min



## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  # i <- 9
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}
numgene[is.na(numgene)] <- 0


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)


# meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene),AUC)
meanAUCAcc <- cbind(sort(numgene)[2], max(numgene), sum(numgene),AUC)
meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
colnames(meanAUCAccAll) <- c("min", "max", "links", "AUROC","AUPR")



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


# da# datapath00 <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
datapath00 <- paste("/home/lly/R/Boolean/Data/DREAM3_RealData",dataset,"/",sep="")
fdata <- paste(datapath00,methodname,".csv",sep="")
write.csv(output, fdata)## times 
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


set.seed(123)

## 
method <- 2
noDIAG = 1
methodname <- "SINCERITIES"
## data 16
file <- c(16,0)
dataset <- file[1]
geneselect <- file[2]


## set path
setwd("/home/lly/MATLAB/Boolean/GRISLI")


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()


# Load data
datapath <- paste("./SCODE-master/data",dataset,"/",sep="")
fdata <- paste(datapath,"data.txt",sep="") # file of expression matrix
X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
datanum <- t(X)
n <- dim(datanum)[1]
p <- dim(datanum)[2]


data <- datanum
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
# nfold <- dim(X_matrix)[1]
nfold <- 10


foldid <- 1:nfold
keep <- TRUE
pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)

lambda_res <- vector()
alpha_res <- vector()

library(glmnet)
for (gi in 1:numGENES) {
  # gi <- 1
  lambda <-  vector()
  cvERROR <-  vector()
  beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
  
  for (test in 1:length(alphas)) {
    # test <- 1
    Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
    
    if(noDIAG==1){
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                              keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    }else{
      CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                              keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    }
    lambda[test] <- CV_results$lambda.min
    cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
    # coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min') 
    # beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    # 2022.5.14 Lingyu Li corrected
    coef.CV_results <- coef(CV_results, s='lambda.min')
    beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    # coef.CV_results <- predict(CV_results, X_matrix, s='lambda.min')
    # beta[coef.CV_results[-1],test] = coef.CV_results[-1]
  } 
  
  minIdx <- max(which(cvERROR==min(cvERROR)))
  lambda_res[gi] <- lambda[minIdx]
  alpha_res[gi] <- alphas[minIdx]
  pred_lambda_min[,gi] <- beta[,minIdx]
  
  print(gi)
  
}

coef <- pred_lambda_min


# Direct Performance Compare -----------------------------------------------------

## 0-????
# SIGN = 0
SIGN = 0


## Gold standard GRN 
fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) 


library(ppcor)
if(SIGN==1){
  parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
  pred_lambda_min <- coef*sign(parcorr_matrix)
}
adj_matrix <- pred_lambda_min



## no zero coefficients and links
numgene <- c()
for (k in 1:p) {
  # i <- 9
  numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
}
numgene[is.na(numgene)] <- 0


library(pracma)
library(Matrix)
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUC <- cbind(AUROC,AUPR)


# meanAUCAcc <- cbind(min(numgene), max(numgene), sum(numgene),AUC)
meanAUCAcc <- cbind(sort(numgene)[2], max(numgene), sum(numgene),AUC)
meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
colnames(meanAUCAccAll) <- c("min", "max", "links", "AUROC","AUPR")



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
##############################################################################################################################
##############################################################################################################################
