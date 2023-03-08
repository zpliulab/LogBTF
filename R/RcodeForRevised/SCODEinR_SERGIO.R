
## Rscript Reconstruct_dynamics.R data/init.txt out/A.txt out/dynamics.txt 100
## Rscript Reconstruct_dynamics.R data/init.txt out/A.txt out/dynamics.txt 
## 2023.2.12 copy from SCODEinR.R and test SERGIO data.



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GENIE3")


# SCODE method ------------------------------------------------------------


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
cellnum <- c(10,20,30,50,100)


## set seed
set.seed(123)
## 
method <- 2


## input LogBTF information
setwd(paste(pathway,"R/Boolean/",dataLinux,"/DREAM3result_LogBTF/",sep=""))
LogBTF <- read.csv(file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""), row.names = 1)



AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()

  
for (k in cellnum) {
  # k <- 10
  
  ## load data
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  Data = t(read.csv(file = dataauto, header=F))
  Data[2,1]
  
  
  ## begin our test
  datanum <- Data 
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  
  
  ## binarization
  if(dataclass == 'count_matrix'){
    datahatA <- Data
    datahatA[datahatA>0] <-1 
    datanum <- t(datahatA)
  }else if(dataclass == 'exper_clean_'){
    ## Kmeans ??????
    library(BiTrinA)
    ## features*times - A n x m matrix comprising m raw measurements of n features
    datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
    
    datahatA[datahatA>0] <-1 
    
    ## Ever, need to genereate simulated scData, Here we dont ues this
    # datahatA <- as.matrix(datanum) ** as.matrix(datahat)
    # datahatA[which(datahatA == 1)] <- 0
    
    datanum <- t(datahatA)
  }

  
  ## use MASS library to calculate pseudo inverse matrix.
  library(MASS)
  
  ## gene
  tfnum <- genenum
  pnum <- 4
  ## cell
  cnum <- celltyp*k
  maxite <- 100
  
  maxB <- 2.0
  minB <- -10.0
  
  mkdir <- "/home/lly/R/SingleCell/SCODE-master"
  
  
  ## gene * cell
  X <- as.matrix(datanum)
  # X <- as.matrix(read.table("data/exp_train.txt", sep="\t"))[1:tfnum,1:cnum]
  W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
  Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
  WZ <- matrix(nrow=tfnum, ncol=cnum)
  
  #read pseudo-time and normalize pseudo-time
  
  pseudotime <- seq(1,cnum,1)
  # pseudotime <- read.table("data/time_train.txt", sep="\t")[1:cnum,2]
  pseudotime <- pseudotime/max(pseudotime)
  
  new_B <- rep(0, pnum)
  old_B <- rep(0, pnum)
  
  #initialization
  RSS <- Inf
  for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
  }
  
  #function to sample Z
  sample_Z <- function(){
    for(i in 1:pnum){
      for(j in 1:cnum){
        Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
      }
    }
  }
  
  #optimize W and B iteratively
  for(ite in 1:maxite){
    #sampling B
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #for last calculation
    if(ite == maxite){
      for(i in 1:pnum){
        new_B[i] <- old_B[i]
      }
    }
    
    #sample Z from new B
    sample_Z()
    
    #regression
    for(i in 1:tfnum){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:pnum){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
      RSS <- tmp_RSS
      old_B[target] <- new_B[target]
    }
    else{
      new_B[target] <- old_B[target]
    }
  }
  
  # #output RSS
  # write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
  # 
  # #output W
  # write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #infer A
  B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
  for(i in 1:pnum){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W)
  A <- W %*% B %*% invW
  
  #output A and B
  # write.table(A, "ResultsLLY/A.txt", row.names=F, col.names=F, sep="\t")
  # write.table(B, paste(dir,"/B.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  
  ## 2023.2.12   ---  old
  weightMatrix0 <- A
  
  ## 2023.2.12   --- New
  ## edge link
  # geneselect <- 258
  geneselect <- LogBTF$links[which(cellnum == k)]
  library(GENIE3)
  linkListNum <- min(getLinkList(A, reportMax = geneselect)[,3])
  weightMatrix0 <- A
  weightMatrix0[weightMatrix0 < linkListNum] = 0
  weightMatrix0[weightMatrix0 >= linkListNum] = 1
  
  
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  ## load ground-truth Adj matrix
  setwd(datapath)
  adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix.csv"))
  
  
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
  AUCALL <- rbind(AUCALL, AUC)
  
  
  ## delect links 
  meanAUCAcc <- cbind(geneselect, AUC)
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("link", "AUROC", "AUPR")
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Acc0 <- performance$Acc
  Recall0 <- performance$Recall
  Pre0 <- performance$Pre
  FPR0 <- performance$FPR
  Fmeasure0 <- performance$Fmeasure
  Result0 <- cbind(Acc0,Recall0,Pre0,FPR0,Fmeasure0)
  ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
  colnames(ResultAll0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")
  
  
  print("SCODE method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



## save
# setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_SCODE01/",sep=""))
# if(method == 2){
#   write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }else if(method == 3){
#   write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }else if(method == 1){
#   write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }







