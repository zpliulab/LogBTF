
rm(list = ls())

time_start<-Sys.time()


## load function
source('D:/E/博士/R_程序/Boolean/R/SINCERITIES functions/PermanceFunction.R')
# source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')



method <- 2
## 0-不排除对角线
noDIAG = 1
SIGN = 1




AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()



setwd("D:/E/博士/山东大学/刘治平/Beginning/2022.1.12_曾健明SIGNET/2021_SIGNET/Ref17/Pseudotime-network-inference-master")
# setwd("/home/lly/E/Beginning/2022.1.12_曾健明SIGNET/2021_SIGNET/Ref17/Pseudotime-network-inference-master")


data <- read.table("binary_expression_LMPP.txt",header = T)
n <- dim(data)[1]
p <- dim(data)[2]
# k <- 4


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



# Performance Compare -----------------------------------------------------

## 0-无向
# SIGN = 0


## Gold standard GRN 
## tf- 1 -Target   ----   adjmatrix
tomaru2 <- read.table("startingNetworkParCor.txt", header = T)
type_regulation <- tomaru2[,2]
netINFO <- tomaru2[,-2]
adj_ref <- matrix(0,nrow = numGENES, ncol = numGENES)

for (i in 1:dim(netINFO)[1]) {
  # i <- 1
  idxGENEsource <- match(netINFO[i,1],colnames(data))
  idxGENEtarget <- match(as.character(t(netINFO[i,-1])),colnames(data))
  if(SIGN==1){
    adj_ref[idxGENEsource,idxGENEtarget] <- type_regulation[i]
  }else{
    adj_ref[idxGENEsource,idxGENEtarget] <- 1
  }
}
adj_gold <- as.matrix(adj_ref)


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




## Network Compare
net <- adj_matrix
colnames(net) <- rownames(net) <- colnames(data)
net[1,]
order(net[1,])

A <- net

which(A["Notch",] != 0)
which(A["Ets1",] != 0)

which(A["Gata2",] != 0)
which(A["Cbfa2t3h",] != 0)
which(A["Nfe2",] != 0)



## predict
library(tidyfst)
edgA <- mat_df(A)
edgA1 <- edgA[which(edgA[,3] != 0),]
# write.csv(edgA1, "edge_predicr_inputDegree8.csv", row.names = F, quote = F)
edgA2 <- edgA1[,-3]
colnames(edgA2) <- c("V1", "V2")

## load edges
edge <- read.table("startingNetworkParCor.txt", header = T)
edgeA <- as.data.frame(cbind(edge[,1], edge[,3], edge[,2]))
# write.csv(edgeA, "edge_given.csv", row.names = F, quote = F)
edgeA2 <- edgeA[,-3]


edg <- rbind(edgA2, edgeA2)
edgall <- edg[!duplicated(edg, fromLast=TRUE), ] 
dim(edgA2)[1]
dim(edgeA2)[1]
dim(edgall)[1]


## 43条边一样
dim(edgA2)[1] + dim(edgeA2)[1] - dim(edgall)[1]



# times
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code Time :',round(exc_time,32),'mins'))
