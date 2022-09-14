
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
methodname <- "SINCERITIES"
## data 2 or 3
file <- c(16,0)
dataset <- file[1]
geneselect <- file[2]



# LogBTFs -----------------------------------------------------------------

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

# datascale <- scale(DISTANCE_matrix)
# datascale <- DISTANCE_matrix
# DISTANCE_matrix01 <- pmax(sign(datascale), 0)
# X_matrix01 <- DISTANCE_matrix01[1:(num_time_points-1),]

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
# AUC
# DyAccuracy
# mean(AUCall0)


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