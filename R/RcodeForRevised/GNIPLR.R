
## 2023.2.9 use Python in R
## Reference：https://blog.51cto.com/u_15264819/2886315
## https://blog.csdn.net/qq_30653631/article/details/107620137


############################################################################
# ## install package  
# py_discover_config(required_module = NULL, use_environment = NULL)
# 
# library(reticulate)
# # install_miniconda()
# py_install("scipy")
# 
# ## see the version of Python used
# py_config()
# ## install Python base
# py_install('sklearn')
# ## check whether is pandas installed?
# py_module_available("pandas")
# 
# # library(reticulate)
# # reticulate::use_virtualenv("~/pythonenvs/userenv")
# 
# 
# # Run Python envirment ----------------------------------------------------
# repl_python()
# ## Python can not work
############################################################################


# Run R code and output Data files, to input in Python  -------------------
## This code is for DREAM3 Dataset. ########################################

## dont use  classification function
## load function
# source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')


## clear workspace
rm(list = ls())
## set seed
set.seed(123)
## define function
GNIPLRsc <- function(filelist, dataname){
  # filelist <- fileEcoli
  # dataname <- "Ecoli"
  
  for (i in 1:length(filelist)) {
    # i <- 1
    
    file <- filelist[[i]]
    
    method = 2
    genenum <- file[1]
    filenum <- file[2]
    k <- file[3]
    geneselect <- file[4]
    
    
    l <- 21
    datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
    dataauto <- paste(datapath,"InSilicoSize",genenum,"-",dataname,filenum,"-nonoise-trajectories.tsv",sep="")
    Data = as.matrix(read.table(file = dataauto, header=T))
    run <- dim(Data)[1]/l
    
    
    datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
    n <- dim(datanum)[1]
    p <- dim(datanum)[2]
    
    
    ## Kmeans ??????
    library(BiTrinA)
    ## features*times - A n x m matrix comprising m raw measurements of n features
    # binarizeMatrix(t(datanum))
    datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
    # write.csv(t(datahatA), "CoefAllSize10\\datadream3p1_matrix.csv")
    datahatA <- as.matrix(datanum) ** as.matrix(datahat)
    datahatA[which(datahatA == 1)] <- 0
    ## GNIPLR need sample * gene. And is ready to use.
    
    
    ##.datahatA is simulated scRNA-seq data;
    ## datanum. is bulk data from DREAM3
    
    
    ## Save bukl data for Python input
    setwd("/home/lly/R/Boolean/Data/DREAM3result_GNIPLR/RData/bulkData")
    dataauto_bulk <- paste("bulkData_",dataname,filenum,"_Node",genenum,".txt", sep="")
    ## Python's input file is txt files.
    write.table(datanum, file = dataauto_bulk, sep = '\t')
    
    
    ## Save single-cell data for Python input    
    setwd("/home/lly/R/Boolean/Data/DREAM3result_GNIPLR/RData/scData")
    dataauto_sc <- paste("scData_",dataname,filenum,"_Node",genenum,".txt", sep="")
    ## Python's input file is txt files.
    write.table(datahatA, file = dataauto_sc, sep = '\t')
    
    
    ## print method
    print("GNIPLR method")
  }
  
  

}


## input data parameter
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))

## Run GNIPLRsc function, then Save bulk Data and scData
GNIPLRsc(fileEcoli, "Ecoli")
GNIPLRsc(fileYeast, "Yeast")





################################################################################
################################################################################
# 2023.2.16 Plot AUROC/ AUPR and so on of GNIPLR--------------------------------


rm(list = ls())

## def
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))


# ## set parameter
filelist <- fileEcoli
dataname <- "Ecoli"
# 
# filelist <- fileYeast
# dataname <- "Yeast"
# 
datafile = "scData"
# datafile = "bulkData"


## Read Expression data
pathname = "/home/lly/"
Boolean = "Data"

## Mac
# pathname = "/Users/lilingyu/E/PhD/"
# Boolean = "Rshell/Data"


## load function
source(paste(pathname,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))







DNIPLRperformance <- function(filelist, dataname, datafile){

  for (i in 1:length(filelist)) {
    # i <- 7
    
    file <- filelist[[i]]
    genenum <- file[1]
    filenum <- file[2]
    k <- file[3]
    geneselect <- file[4]
    
    
    ## path and data
    setwd(paste(pathname,'R/Boolean/',Boolean,'/DREAM3result_GNIPLR/RData/', sep = ""))
    # setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/RData/")
    dataauto_sc <- paste(datafile, "/", datafile, "_",dataname,filenum,"_Node",genenum,
                         "_lasso.txt", sep="")
    result <- read.table(dataauto_sc, sep = "\t")
    result[is.na(result[,3]),3] <- 0
    # list <- result[,-4]

    
    ## List to Adj
    tomaru2 <- result
    type_regulation <- tomaru2[,3]
    netINFO <- tomaru2[,-c(3,4)]
    
    ## covert to character
    netINFO$V1 <- as.character(netINFO$V1)
    netINFO$V2 <- as.character(netINFO$V2)
    
    adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
    
    ## give the gene name, here is number
    genename <- stringr::str_c("G",seq(1,genenum,1))
    
    SIGN = 1
    
    for (i in 1:dim(netINFO)[1]) {
      # i <- 1
      idxGENEsource <- match(netINFO[i,1],genename)
      idxGENEtarget <- match(netINFO[i,2],genename)
      if(SIGN==1){
        adj_ref[idxGENEsource,idxGENEtarget] <- type_regulation[i]
      }else{
        adj_ref[idxGENEsource,idxGENEtarget] <- 1
      }
    }
    
    
    # Perf on all TF
    weightMatrix0 <- adj_ref
    
    ## 2023.2.12   --- New
    ## edge link truthedges
    # geneselect <- LogBTF$links[which(cellnum == k)]
    # geneselect <- truthedges
    library(GENIE3)
    linkListNum <- min(getLinkList(adj_ref, reportMax = geneselect)[,3])
    weightMatrix0 <- adj_ref
    weightMatrix0[weightMatrix0 < linkListNum] = 0
    weightMatrix0[weightMatrix0 >= linkListNum] = 1
    
    
    
    ## 0-????
    SIGN = 0
    # SIGN = 1
    
    
    ## load ground-truth Adj matrix
    ## Gold standard GRN  
    
    datapath0 <- paste(pathname,"R/Boolean/",Boolean,"/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
    # datapath0 <- paste("D:/E/??ʿ/R_????/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
    adj_gold0 <- paste(datapath0," InSilicoSize",genenum,"-",dataname,filenum,"-adj .csv",sep="")
    adj_gold <- as.matrix(read.csv(file = adj_gold0))
    
    
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
    
    
    print("GNIPLR method")
    
  }
  
  ## output
  output <- cbind(meanAUCAccAll, ResultAll0)
  return(output)
  
}


## Run function and collect results

### Ecoli set parameter  -------------------------------------------------------
filelist <- fileEcoli
dataname <- "Ecoli"

AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()

Ecoli_sc <- DNIPLRperformance(fileEcoli, "Ecoli", "scData")



### Yeast set parameter  -------------------------------------------------------
filelist <- fileYeast
dataname <- "Yeast"

AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()

Yeast_sc <- DNIPLRperformance(fileYeast, "Yeast", "scData")


# Combinate Ecoli and Yeast Data ------------------------------------------

result <- rbind(Ecoli_sc, Yeast_sc)



# Ecoli_bulk <- DNIPLRperformance(fileEcoli, "Ecoli", "bulkData")
# Yeast_bulk <- DNIPLRperformance(fileYeast, "Yeast", "bulkData")




# GNIPLRresult <- cbind(rbind(Ecoli_sc, Yeast_sc), rbind(Ecoli_bulk, Yeast_bulk)[,-1])
# colnames(GNIPLRresult) <- c("Dataset", "scAUC", "bulkAUC")

setwd(paste(pathname,'R/Boolean/',Boolean,'/DREAM3result_GNIPLR/', sep = ""))
write.csv(result, file = "GNIPLRresult_DREAM3_sc.csv", row.names = F)







################################################################################
################################################################################
# Plot ROC curve of GNIPLR----------------------------------------------------------


rm(list = ls())

## def
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))


# ## set parameter
# filelist <- fileEcoli
# dataname <- "Ecoli"
# 
# filelist <- fileYeast
# dataname <- "Yeast"
# 
# datafile = "scData"
# datafile = "bulkData"


## Read Expression data
pathname = "/home/lly/"
Boolean = "Data"

## Mac
# pathname = "/Users/lilingyu/E/PhD/"
# Boolean = "Rshell/Data"


DNIPLRperformance <- function(filelist, dataname, datafile){
  
  aucvalue <- c()
  lab <- c()
  for (i in 1:length(filelist)) {
    # i <- 1
    
    file <- filelist[[i]]
    genenum <- file[1]
    filenum <- file[2]
    k <- file[3]
    geneselect <- file[4]
    
    
    ## path and data
    setwd(paste(pathname,'R/Boolean/',Boolean,'/DREAM3result_GNIPLR/RData/', sep = ""))
    # setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/RData/")
    dataauto_sc <- paste(datafile, "/", datafile, "_",dataname,filenum,"_Node",genenum,
                         "_lasso.txt", sep="")
    result <- read.table(dataauto_sc, sep = "\t")
    
    library(pROC)
    p_test <- result[,3]
    y.test <- result[,4]
    A_test <- data.frame(p_test, y.test)
    names(A_test)<- c("p", "outcome")
    # pdf(file = "ROC_independent46.pdf",width = 5,height = 5)
    p <- plot.roc(A_test$outcome, A_test$p, print.auc=T)
    aucvalue <- c(aucvalue, p$auc)
    lab <- c(lab, paste(dataname,filenum,"_Node",genenum,sep=""))
    
  }
  
  resultsc <- cbind(lab, aucvalue)
  return(resultsc)
}


## Run function and collect results
Ecoli_sc <- DNIPLRperformance(fileEcoli, "Ecoli", "scData")
Yeast_sc <- DNIPLRperformance(fileYeast, "Yeast", "scData")

Ecoli_bulk <- DNIPLRperformance(fileEcoli, "Ecoli", "bulkData")
Yeast_bulk <- DNIPLRperformance(fileYeast, "Yeast", "bulkData")


GNIPLRresult <- cbind(rbind(Ecoli_sc, Yeast_sc), rbind(Ecoli_bulk, Yeast_bulk)[,-1])
colnames(GNIPLRresult) <- c("Dataset", "scAUC", "bulkAUC")

setwd(paste(pathname,'R/Boolean/',Boolean,'/DREAM3result_GNIPLR/', sep = ""))
write.csv(GNIPLRresult, file = "GNIPLRresult_DREAM3.csv", row.names = F)








## Bellow is for single data
################################################################################
################################################################################
# For the single situation case -------------------------------------------

rm(list = ls())

## def
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))


## set parameter
filelist <- fileEcoli
dataname <- "Ecoli"

filelist <- fileYeast
dataname <- "Yeast"


aucvalue <- c()
lab <- c()
for (i in 1:length(filelist)) {
  # i <- 1
  
  file <- filelist[[i]]
  genenum <- file[1]
  filenum <- file[2]
  k <- file[3]
  geneselect <- file[4]
  
  
  ## path and data
  setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/RData/")
  dataauto_sc <- paste("scData/scData_",dataname,filenum,"_Node",genenum,
                       "_lasso.txt", sep="")
  result <- read.table(dataauto_sc, sep = "\t")
  
  library(pROC)
  p_test <- result[,3]
  y.test <- result[,4]
  A_test <- data.frame(p_test, y.test)
  names(A_test)<- c("p", "outcome")
  # pdf(file = "ROC_independent46.pdf",width = 5,height = 5)
  p <- plot.roc(A_test$outcome, A_test$p, print.auc=T)
  aucvalue <- c(aucvalue, p$auc)
  lab <- c(lab, paste(dataname,filenum,"_Node",genenum,sep=""))
  
}

resultsc <- cbind(lab, aucvalue)



