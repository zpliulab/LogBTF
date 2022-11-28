####################
# PACKAGE required:
# pracma
####################

auc_from_ranks_TC_sign <- function(acc_ranked,adj,n_points){
  maxWEIGHT <- max(acc_ranked)
  acc_ranked <- acc_ranked/maxWEIGHT
  adj <- adj*3
  
  fp <- matrix(0,nrow = 1,ncol = n_points+1)
  tn <- matrix(0,nrow = 1,ncol = n_points+1)
  tp <- matrix(0,nrow = 1,ncol = n_points+1)
  fn <- matrix(0,nrow = 1,ncol = n_points+1)
  fpr <- matrix(0,nrow = 1,ncol = n_points+1)
  precision <- matrix(0,nrow = 1,ncol = n_points+1)
  recall <- matrix(0,nrow = 1,ncol = n_points+1)
  
  thr <- seq(0,1,by=1/n_points)
  sign_info <- sign(acc_ranked)
  sign_info[sign_info==0] <- 1
  for (i in 1:(n_points+1)) {
    adj_m <- (abs(acc_ranked)>=thr[i])*sign_info
    compare <- abs(adj_m+adj)
    fp[1,i] <- nnzero(compare<3&compare>0)
    tp[1,i] <- nnzero(compare==4)
    fn[1,i] <- nnzero(compare==3)
    tn[1,i] <- nnzero(compare==0)
    precision[1,i]=tp[1,i]/(tp[1,i]+fp[1,i]);
    recall[1,i]=tp[1,i]/(tp[1,i]+fn[1,i]);
    fpr[1,i]=fp[1,i]/(fp[1,i]+tn[1,i]);
  }
  precision <- rev(precision)
  recall <- rev(recall)
  precision_new <- c(1,precision)
  recall_new <- c(0,recall)
  fpr <- rev(fpr)
  
  # auc <- dget("D:/E/博士/R_程序/Boolean/R/SINCERITIES functions/auc.R")
  auc <- function(x,y){
  if(length(x)!=length(y)){
    stop("The input vectors should have the same size")
  }
  n <- length(x)
  xix <- order(x)
  x_ <- x[xix]
  y_ <- y[xix]
  NAidx <- is.na(x_)|is.na(y_)
  x_ <- x_[!NAidx]
  y_ <- y_[!NAidx]
  area <- trapz(x_,y_)
  return(area)
}
  AUROC <- auc(fpr,recall)
  AUPR <- auc(recall_new,precision_new)
  
  result <- list(AUROC=AUROC,AUPR=AUPR,fpr=fpr,recall=recall,recall_new=recall_new,precision_new=precision_new)
  return(result)
}



performance_from_ranks_TC_sign <- function(acc_ranked,adj,n_points){

maxWEIGHT <- max(acc_ranked)
acc_ranked <- acc_ranked/maxWEIGHT
adj <- adj*3

fp <- matrix(0,nrow = 1,ncol = n_points+1)
tn <- matrix(0,nrow = 1,ncol = n_points+1)
tp <- matrix(0,nrow = 1,ncol = n_points+1)
fn <- matrix(0,nrow = 1,ncol = n_points+1)
fpr <- matrix(0,nrow = 1,ncol = n_points+1)
precision <- matrix(0,nrow = 1,ncol = n_points+1)
recall <- matrix(0,nrow = 1,ncol = n_points+1)
accuracy <- matrix(0,nrow = 1,ncol = n_points+1)
Fmeasure <- matrix(0,nrow = 1,ncol = n_points+1)


sign_info <- sign(acc_ranked)
sign_info[sign_info==0] <- 1

for (i in 1:(n_points+1)) {
  # i <- 1
  adj_m <- (abs(acc_ranked)>0)*sign_info
  # adj_m <- (abs(acc_ranked)>thr[i])*sign_info
  compare <- abs(adj_m+adj)
  
  ## 计算非0个数
  fp[1,i] <- nnzero(compare<3&compare>0)
  tp[1,i] <- nnzero(compare==4)   # TP  -- gold standard
  fn[1,i] <- nnzero(compare==3)   # FN  -- 实际有边，预测的没有
  tn[1,i] <- nnzero(compare==0)   # TN  -- gold standard
  
  precision[1,i] = tp[1,i]/(tp[1,i]+fp[1,i]);
  recall[1,i]    = tp[1,i]/(tp[1,i]+fn[1,i]);
  fpr[1,i]       = fp[1,i]/(fp[1,i]+tn[1,i]);
  accuracy[1,i]  = (tp[1,i] + tn[1,i])/(tp[1,i]+tn[1,i]+fp[1,i]+fn[1,i]);
  Fmeasure[1,i]  = 2*tp[1,i]/(2*tp[1,i] + fp[1,i] + fn[1,i])
}

result <- list(Acc=accuracy,Recall=recall,Pre=precision,FPR=fpr,Fmeasure=Fmeasure)
return(result)
}


## 2022.7.24
## sign function
sgn <- function(xx){
  y <- xx
  for (k in 1:dim(xx)[1]) {
    # k <- 1
    if (xx[k,1] >= 0)
      y[k,] <- 1
    else (y[k,] <- 0)
  }
  return(y) 
}


## 2022.7.24
## Dynamic Regression 
DyAcc <- function(datahatA,ptrainall){
  X <- datahatA[-1,]
  kk <- 0
  for (j in 1:dim(X)[2]) {
    for (k in 1:(dim(X)[1])) {
      if (as.numeric(ptrainall[k,j]) == as.numeric(X[k,j]))
        kk <- kk + 1
    }
  }
  acc <- kk/((dim(X)[2])*(dim(X)[1]))
  return(acc)
  
}


## 2022.7.24
## Modified Dynamic accuracy
DyAccWeight <- function(predMatrix){
  predMatrix01 <- apply(predMatrix, 2, function(x){(x-min(x))/(max(x)-min(x))})
  predMatrix01[predMatrix01 < 0.5] <- 0
  predMatrix01[predMatrix01 >= 0.5] <- 1
  DyAcc <- 1-sum((predMatrix01-exprMatrix) != 0)/(dim(exprMatrix)[1]*dim(exprMatrix)[2])
  return(DyAcc)
}


## 2022.11.28 
## Add perturbation design matrix
PerMatrix <- function(datahatA){
  ## noise matrix
  set.seed(2022)
  res <- datahatA[1:(n-1),]
  noiseLevel = 1e-5
  res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))
  
  ## when 1, it add noise
  for (i in 1:(dim(res)[1])) {
    for (j in 1:dim(res)[2]) {
      if (datahatA[i,j] == 0)
        res[i,j] <- 0
    }
  }
  return(res)
}

