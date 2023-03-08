## 2023.2.12 based on code dataDREAMGlmSINCERITIESadj_Ecoli01.R

## for count data :
## Warning messages:
# 1: In pcor(X_matrix, method = "spearman") :
#   The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.
# 2: In sqrt(1/diag(V)) : 产生了NaNs
# 3: In cov2cor(icvx) : diag(.)内有0或NA值；无限值的可靠性不高
# 4: In sqrt((n - 2 - gp)/(1 - pcor^2)) : 产生了NaNs
# 5: In pt(-abs(statistic), (n - 2 - gp)) : 产生了NaNs

## for clear data, it is available for first two
## for the third data,  Error in solve.default(cvx) : 系统计算上是奇异的: 倒条件数=1.7953e-19


# SINCERITIES method ------------------------------------------------------------


rm(list = ls())

## set pathway
pathway <- '/home/lly/'
dataLinux <- 'Data'

#pathway <- '/Users/lilingyu/E/PhD/'
#dataLinux <- "RShell/Data" 


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellbum <- 10
cellnum <- c(10,20,30,40,50)


## set seed
set.seed(123)
##
method <- 2
## 0-???ų??Խ???
noDIAG = 1
epsilon <- 1e-4


# SINCERITIES -----------------------------------------------------------------


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()



for (k in cellnum) {
  # a <- i + 1
  # print(a)
  # k <- 20
  
  ## load data
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  Data = t(read.csv(file = dataauto, header=F))
  # Data = t(read.csv(file = dataauto, header=F))[1:k,]
  Data[2,1]
  

  ## begin our test
  datanum <- Data
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  
  
  ## binarization
  if(dataclass == 'count_matrix'){
    datahatA <- Data
    datahatA[datahatA>0] <-1
    datanum <- datahatA
  }else if(dataclass == 'exper_clean_'){
    ## Kmeans ??????
    library(BiTrinA)
    ## features*times - A n x m matrix comprising m raw measurements of n features
    datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
    datahatA[datahatA>0] <-1
    datanum <- datahatA
  }
  

  
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
    # gi <- 4
    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
    
    for (test in 1:length(alphas)) {
      # test <- 1
      Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]

      # if Y exist one 1/0, use noise 0/1 data.
      # if(sum(Y_vector != 0) > 1){
      ## 2023.2.13 add. Because there are so many vector is all 0 and all 1, or only one 1 and one 0.......
      if(sum(Y_vector != 0) == n | sum(Y_vector != 0) == 0 | sum(Y_vector != 0) == 1 | sum(Y_vector != 0) == (n-1)){
        
        set.seed(123)
        glm.fit <- lm(Y_vector~X_matrix, data.frame(X_matrix))
        coef <- glm.fit$coefficients
        coef[is.na(coef)] <- 0
        beta[,test] = coef[-1]

      }else{
        
        
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
        # beta[is.na(beta)] <- 0
        # coef.CV_results <- predict(CV_results, X_matrix, s='lambda.min')
        # beta[coef.CV_results[-1],test] = coef.CV_results[-1]
        minIdx <- max(which(cvERROR==min(cvERROR)))
        lambda_res[gi] <- lambda[minIdx]
        alpha_res[gi] <- alphas[minIdx]
        pred_lambda_min[,gi] <- beta[,minIdx]
      }
      
      print(gi)
        

    }
    
  }
  
  coef <- pred_lambda_min
  coef[is.na(coef)] <- 0
  
  
  # Direct Performance Compare -----------------------------------------------------
  
  ## 0-????
  # SIGN = 0
  SIGN = 1
  
  
  library(ppcor)
  if(SIGN==1){
    parcorr_matrix <- pcor(X_matrix,method = 'spearman')$estimate
    ## 2023.2.14 add 
    parcorr_matrix[is.na(parcorr_matrix)] <- 0
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  adj_matrix <- pred_lambda_min
  
  
  
  ## no zero coefficients and links
  numgene <- c()
  for (k in 1:p) {
    # i <- 9
    numgene <- cbind(numgene, as.numeric(summary(adj_matrix[,k] != 0)[3]))
  }
  # numgene[is.na(numgene)] <- 0
  ## delete NA value, using numeric class dtypr ata
  numgene1 <- as.numeric(numgene)
  numgene1 <- na.omit(numgene1)
  numgene <- numgene1
  # sum(numgene)
  

  ## load ground-truth Adj matrix
  setwd(datapath)
  adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix.csv"))
  
  
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
    ## 2023.2.14 add 
    parcorr_matrix[is.na(parcorr_matrix)] <- 0
    pred_lambda_min <- coef*sign(parcorr_matrix)
  }
  adj_matrix <- pred_lambda_min
  
  
  ## load ground-truth Adj matrix
  setwd(datapath)
  adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix.csv"))
  
  
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
  
  print("SICNERITIES method")
}

## output
 output <- cbind(meanAUCAccAll[,1:3], AUCALL0[,1:2], ResultAll0[,1:5], AUCALL0[,3:4], ResultAll0[,6:10])
output

 # # save
# setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_SICNERITIES01/",sep=""))
# 
#  if(method == 2){
#    write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
#  }else if(method == 3){
#    write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
#  }else if(method == 1){
#    write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
#  }

  
