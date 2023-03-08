## 2023.2.12 copy from dataDREAMTREGRESS01 and test SERGIO data.


# install.packages('devtools')
# library(devtools)
# install_github("jpvert/tigress")


# TREGRESS method ------------------------------------------------------------------

rm(list = ls())


## set pathway
 pathway <- '/home/lly/'
 dataLinux <- "Data" 

#pathway <- '/Users/lilingyu/E/PhD/'
#dataLinux <- "RShell/Data"


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))


## Parameter setting
# # dataclass <- 'count_matrix'
# dataclass <- 'exper_clean_'
# genenum <- 100
# celltyp <- 1
# # cellbum <- 10
# cellnum <- c(10,20,30,50,100)
 
 
# ## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 20
celltyp <- 1
# cellbum <- 10
cellnum <- c(10,20,30,40,50,60)
netname <- "Interaction_cID_4SIGNmatrix_20node.csv"



## set seed
set.seed(123)
## 
method <- 2


## input LogBTF information
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_LogBTF/",sep=""))
LogBTF <- read.csv(file = paste("ENet_output_SERGIO_count_matrix",genenum,"_",celltyp,".csv",sep=""), row.names = 1)



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
  # Data = t(read.csv(file = dataauto, header=F))[1:k,]
  # Data = t(read.csv(file = dataauto, header=F))[(k+1):(2*k),]
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
  
  
  # TREGRESS --  sample*gene
  X <- t(datanum)
  rownames(X) <- seq(1,genenum,1)
  
  # TIGRESS parameters
  nstepsLARS <- 1  # 10 ?????У?????
  nsplit <- 1000
  library(tigress)
  
  
  ## Train TIGRESS on all TF
  set.seed(123)
  predTigress <- tigress(t(X), nstepsLARS=nstepsLARS, nsplit=nsplit)
  # Perf on all TF
  weightMatrix <- predTigress[[1]]
  
  
  ## 2023.2.13   --- New
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
  adj_gold <- as.matrix(read.csv(file = netname))
  
  
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
  
  
  print("TREGRESS method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



# save
 setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_TREGRESS01/",sep=""))
 if(method == 2){
   write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }else if(method == 3){
   write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }else if(method == 1){
   write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
 }
