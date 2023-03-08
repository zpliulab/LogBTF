## 2023.2.12 copy from dataDREAMCLR01 and test SERGIO data.


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("minet")


# CLR method ------------------------------------------------------------------

rm(list = ls())


## set pathway
 pathway <- '/home/lly/'
 dataLinux <- "Data" 

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
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_LogBTF/",sep=""))
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
    datanum <- datahatA
  }else if(dataclass == 'exper_clean_'){
    ## Kmeans ??????
    library(BiTrinA)
    ## features*times - A n x m matrix comprising m raw measurements of n features
    datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
    
    datahatA[datahatA>0] <-1 
    
    ## Ever, need to genereate simulated scData, Here we dont ues this
    # datahatA <- as.matrix(datanum) ** as.matrix(datahat)
    # datahatA[which(datahatA == 1)] <- 0
    
    datanum <- datahatA
  }
  
  
  ## expression matrix : sample * genes
  exprMatrix <- datanum
  # exprMatrix <- datahatA
  
  
  # CLR -----------------------------------------------------------------
  library(minet)
  ## input data -- sample*genes
  set.seed(123)
  mim <- build.mim(exprMatrix,estimator="spearman")
  ## 2023.2.14 NAN -- 0
  mim[is.na(mim)] <- 0
  weightMatrix <- clr(mim)
  weightMatrix[is.na(weightMatrix)] <- 0
  
  ## 2023.2.12   --- New
  ## edge link
  library(GENIE3)
  geneselect <- LogBTF$links[which(cellnum == k)]
  linkListNum <- min(getLinkList(weightMatrix, reportMax = geneselect)[,3])
  weightMatrix0 <- weightMatrix
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
  
  
  print("CLR method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



 # save
 setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_CLR01/",sep=""))
 if(method == 2){
   write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }else if(method == 3){
   write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }else if(method == 1){
   write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }



