## 3023.2.11 according LMPPGlmPenalty.R in RShell by Lingyu LI


## 2023.2.15  Compute the AUROC  and other performance index of DeepSEM method

####################################################################################
####################################################################################
## 2023.2.14 List to Adj marix ---------------------------------

## input data
rm(list = ls())


## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'

##
set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
##
SIGN = 1

dataclass <- 'count_matrix'
# dataclass <- 'exper_clean'
genenum <- 100
celltyp <- 2
cellbum <- 10


## Notice that the gene node is brgin from 0, but not 1

setwd(paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', 
            sep = ""))
inter <- read.csv(file = "Interaction_cID_4SIGN.csv", header=T)
## Convert 2 to -1, Keep 1 is 1
inter[which(inter$type == 2), 3] <- -1

## List to Adj
tomaru2 <- inter
type_regulation <- tomaru2[,3]
netINFO <- tomaru2[,-3]
adj_ref <- matrix(0, nrow = genenum, ncol = genenum)

## give the gene name, here is number
# genename <- seq(1,genenum,1)
genename <- seq(0,(genenum-1),1)


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

## See how many edges
sum(abs(adj_ref))    # 256  

## Save adj matrix results
# write.csv(adj_ref, file = "Interaction_cID_4SIGNmatrix.csv", row.names = F)


## Matrix to List
library(tidyfst)
colnames(adj_ref) <- rownames(adj_ref) <- seq(1,100,1)
## with sign
# List <- mat_df(adj_ref)

## without sign
List <- mat_df(abs(adj_ref))

## delet duplication
List <- as.matrix(List[-which(List[,1] == List[,2]),])
dim(List)
List[1,2]

List <- apply(List, 2, as.numeric)

# write.csv(List, file = "Interaction_cID_4SIGNlist.csv", row.names = F)
## delete ""
# write.table(List, file = "Interaction_cID_4SIGNlist.txt", row.names = F, quote=F)


## for GNILPR method
setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/")
write.table(List, file = "Interaction_cID_4SIGNlist.txt", row.names = F, col.names = F, quote=F, sep = "\t")



####################################################################################
####################################################################################
# Input data and add gene label for GNIPLR -------------------------------------------

rm(list = ls())

## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50)


for (k in cellnum) {
  # k <- 10
  
  ## load data
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  Data = t(read.csv(file = dataauto, header=F))
  
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
  
  colnames(datanum) <- seq(1,p,1)
  
  # dataauto1 <-  "/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/"
  dataauto1 <-  paste(pathway, "R/Boolean/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  # write.table(datanum, file = paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, ".txt", 
  #                                 sep=""),sep="\t", row.names = F, quote = F)
  
  
  ## must row.names = F
  write.table(datanum, file = paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, ".txt", 
                                    sep=""),sep="\t", row.names = F)

  print("Data GNIPLR")

  
}


####################################################################################
####################################################################################
# 2023.2.18 Input Real Data5/16 and add gene label for GNIPLR -------------------------------------------

rm(list = ls())

## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'


## Parameter setting
# # dataclass <- 'count_matrix'
# dataclass <- 'exper_clean_'
# genenum <- 100
# celltyp <- 1
# # cellnum <- 10
# cellnum <- c(10,20,30,40,50)


## Parameter setting
# file <- c(5,2074)
file <- c(16,36003)
dataclass <- 'count_matrix'
# dataclass <- 'exper_clean_'

  
## load data
dataauto <- paste(pathway,'MATLAB/Boolean/GRISLI/SCODE-master/data', file[1], '/data.txt', sep = "")
Data = t(read.table(file = dataauto, header=F))

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
  library(BiTrinA)dui
  ## features*times - A n x m matrix comprising m raw measurements of n features
  datahatA <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  
  datahatA[datahatA>0] <-1 
  
  ## Ever, need to genereate simulated scData, Here we dont ues this
  # datahatA <- as.matrix(datanum) ** as.matrix(datahat)
  # datahatA[which(datahatA == 1)] <- 0
  
  datanum <- datahatA
  
}

colnames(datanum) <- seq(1,p,1)

dataauto1 <-  paste(pathway,'MATLAB/Boolean/GRISLI/SCODE-master/data', file[1], sep = "")

## must row.names = F
write.table(datanum, file = paste(dataauto1, "/Label_data.txt", 
                                  sep=""),sep="\t", row.names = F)
## for Data 16, dont use 0/1
# write.table(datanum, file = paste(dataauto1, "/Label0_data.txt", 
#                                   sep=""),sep="\t", row.names = F)

print("Data GNIPLR")
  


####################################################################################
####################################################################################
## 2023.2.18 Real Data5/16 and its Adj marix ---------------------------------

## input data
rm(list = ls())


## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'

##
set.seed(123)
##
SIGN = 1

## Parameter setting
# file <- c(5,2074)
file <- c(16,36003)
dataclass <- 'count_matrix'
# dataclass <- 'exper_clean_'


## input data
## Gold standard GRN 
datapath <- paste(pathway,'MATLAB/Boolean/GRISLI/SCODE-master/data', file[1], sep = "")
fA <- paste(datapath,"/A.txt",sep="") # file of gold standard network
adj_ref <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target
genenum <- dim(adj_ref)[1]

## Matrix to List
library(tidyfst)
colnames(adj_ref) <- rownames(adj_ref) <- seq(1,genenum,1)
## with sign
# List <- mat_df(adj_ref)

## without sign
List <- mat_df(abs(adj_ref))

## delet duplication
List <- as.matrix(List[-which(List[,1] == List[,2]),])
dim(List)
List[1,2]

List <- apply(List, 2, as.numeric)


## for GNILPR method
write.table(List, file = paste(datapath, "/Interaction_A_List.txt", sep = ""), row.names = F, col.names = F, quote=F, sep = "\t")


####################################################################################
####################################################################################
# Plot AUROC --------------------------------------------------------------



rm(list = ls())

## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50)


aucvalue <- c()
for (i in 1:length(cellnum)) {
  # i <- 1
  k <- cellnum[i]
  
  ## path and data
  dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  dataauto_sc <- paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, "_lasso.txt", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t")
  result[is.na(result[,3]),3] <- 0
  ## 0.52, 0.54, 0.58
  # result <- result[-is.na(result[,3]),]   
  
  library(pROC)
  p_test <- result[,3]
  y.test <- result[,4]
  A_test <- data.frame(p_test, y.test)
  names(A_test)<- c("p", "outcome")
  # pdf(file = "ROC_independent46.pdf",width = 5,height = 5)
  p <- plot.roc(A_test$outcome, A_test$p, print.auc=T)
  aucvalue <- c(aucvalue, p$auc)

}


## Run function and collect results

setwd(dataauto1)
write.csv(aucvalue, file = paste(dataclass, "GNIPLRresult_AUROC.csv", sep = ""), row.names = F)




####  Compute the AUROC  and other performance index of GNIPLR method


####################################################################################
####################################################################################
# 2022.2.16 GNIPLR method All node -- Use Performance Function  ------------------------------------------------

rm(list = ls())


## set pathway
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50)
truthedges <- 258

## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2

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
  
  
  ## path and data
  # dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  dataauto1 <- paste(pathway,"R/Boolean/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  dataauto_sc <- paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, "_lasso.txt", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t")
  result[is.na(result[,3]),3] <- 0
  # list <- result[,-4]
  
  
  
  
  ## List to Adj
  tomaru2 <- result
  type_regulation <- tomaru2[,3]
  netINFO <- tomaru2[,-c(3,4)]
  adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
  
  ## give the gene name, here is number
  genename <- seq(1,genenum,1)
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
  geneselect <- LogBTF$links[which(cellnum == k)]
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
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
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
  
  
  print("GNIPLR method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



# save
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_GNIPLR/",sep=""))
if(method == 2){
  write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 3){
  write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 1){
  write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}


####################################################################################
####################################################################################
# 2022.2.18 Real Data5/16 -- GNIPLR method -- Use Performance Function  ------------------------------------------------


rm(list = ls())


## set pathway
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
# genenum <- 100
# file <- c(5,2074)
genenum <- 948
file <- c(16,36003)
truthedges <- file[2]

## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2


## path and data
dataauto1 <- paste(pathway,'MATLAB/Boolean/GRISLI/SCODE-master/data', file[1], sep = "")
dataauto_sc <- paste(dataauto1, "/Label_data", "_lasso.txt", 
                     sep="")
result <- read.table(dataauto_sc, sep = "\t")
result[is.na(result[,3]),3] <- 0
# list <- result[,-4]


## List to Adj
tomaru2 <- result
type_regulation <- tomaru2[,3]
netINFO <- tomaru2[,-c(3,4)]
adj_ref <- matrix(0, nrow = genenum, ncol = genenum)

## give the gene name, here is number
genename <- seq(1,genenum,1)
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
geneselect <- truthedges
library(GENIE3)
linkListNum <- min(getLinkList(adj_ref, reportMax = geneselect)[,3])
weightMatrix0 <- adj_ref
weightMatrix0[weightMatrix0 < linkListNum] = 0
weightMatrix0[weightMatrix0 >= linkListNum] = 1



## 0-????
SIGN = 0
# SIGN = 1


AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()


## load ground-truth Adj matrix
## Gold standard GRN 
datapath <- paste(pathway,'MATLAB/Boolean/GRISLI/SCODE-master/data', file[1], sep = "")
fA <- paste(datapath,"/A.txt",sep="") # file of gold standard network
adj_gold <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target


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


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



# save
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_GNIPLR/",sep=""))
if(method == 2){
  write.csv(output, file = paste("ENet_output_RealData_",file[1],".csv",sep=""))
}else if(method == 3){
  write.csv(meanAUCAccAll, file = paste("Lasso_output_RealData_",file[1],".csv",sep=""))
}else if(method == 1){
  write.csv(meanAUCAccAll, file = paste("Ridge_output_RealData_",file[1],".csv",sep=""))
}


####################################################################################
####################################################################################
####  Compute the AUROC  and other performance index of DeepSEM method
## 2023.2.16 DeepSEM method --  Because only give the regulatory of All 10 Tf genes,  ----------------
## For all nodes --  So we test to comare the truth 258 edges and see the AUC values

rm(list = ls())


## set pathways
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50)
truthedges <- 258


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2

## input LogBTF information
## Dont select the regulatory link, according to LogBTF method
## We select the number of edges of 158, the truth regulatory edges

# setwd(paste(pathway,"R/Boolean/",dataLinux,"/DREAM3result_LogBTF/",sep=""))
# LogBTF <- read.csv(file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""), row.names = 1)




AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()




for (k in cellnum) {
  # k <- 10
  
  
  ## path and data
  # dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_DeepSEM/SERGIOData/", sep = "")
  dataauto1 <- paste(pathway,"Python/DeepSEM-master/out/", sep = "")
  dataauto_sc <- paste(dataauto1, "GRN_inference_result_", dataclass,genenum,"_",celltyp, "_", k, ".tsv", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t", header = T)
  
  ## delete "G"
  result$TF <- stringr::str_split_fixed(result$TF, "G", 2)[,2]
  result$Target <- stringr::str_split_fixed(result$Target, "G", 2)[,2]
  ## convert to number and form 0 to 1
  # as.numeric(result$Target)+1
  list <- result
  list[,1] <- as.numeric(result$TF)+1
  list[,2] <- as.numeric(result$Target)+1
  

  ## List to Adj
  tomaru2 <- list
  type_regulation <- tomaru2[,3]
  netINFO <- tomaru2[,-3]
  # unique(netINFO[,1])   # 只有TF  94 45 75 18  2 85 57 63 15 68
  # unique(netINFO[,2])   # 100
  adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
  
  ## give the gene name, here is number
  genename <- seq(1,genenum,1)
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
  sum(adj_ref != 0)    # 990
  
  ## Perf on all TF
  weightMatrix0 <- adj_ref
  
  ## 2023.2.12   --- New
  ## edge link
  # geneselect <- LogBTF$links[which(cellnum == k)]
  geneselect <- truthedges
  library(GENIE3)
  linkListNum <- min(getLinkList(adj_ref, reportMax = geneselect)[,3])
  weightMatrix0 <- adj_ref
  weightMatrix0[weightMatrix0 < linkListNum] = 0
  weightMatrix0[weightMatrix0 >= linkListNum] = 1
  
  
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  ## load ground-truth Adj matrix
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
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
  # meanAUCAcc <- AUC
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("link", "AUROC", "AUPR")
  # colnames(meanAUCAccAll) <- c("AUROC", "AUPR")
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Acc0 <- performance$Acc
  Recall0 <- performance$Recall
  Pre0 <- performance$Pre
  FPR0 <- performance$FPR
  Fmeasure0 <- performance$Fmeasure
  Result0 <- cbind(Acc0,Recall0,Pre0,FPR0,Fmeasure0)
  ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
  colnames(ResultAll0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")
  
  
  print("DeepSEM method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



## save
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_DeepSEM/",sep=""))
if(method == 2){
  write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 3){
  write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 1){
  write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}



####################################################################################
####################################################################################
## 2023.2.16 DeepSEM method --  Because only give the regulatory of All 10 Tf genes,  ----------------
## For all nodes --  So we test to comare the truth 258 edges and see the AUC values
## Different from the Last one, here we compare with LogBTF methds

rm(list = ls())


## set pathway
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 100
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50)


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2

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
  
  
  ## path and data
  # dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_DeepSEM/SERGIOData/", sep = "")
  dataauto1 <- paste(pathway,"Python/DeepSEM-master/out/", sep = "")
  dataauto_sc <- paste(dataauto1, "GRN_inference_result_", dataclass,genenum,"_",celltyp, "_", k, ".tsv", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t", header = T)
  
  ## delete "G"
  result$TF <- stringr::str_split_fixed(result$TF, "G", 2)[,2]
  result$Target <- stringr::str_split_fixed(result$Target, "G", 2)[,2]
  ## convert to number and form 0 to 1
  # as.numeric(result$Target)+1
  list <- result
  list[,1] <- as.numeric(result$TF)+1
  list[,2] <- as.numeric(result$Target)+1
  
  
  ## List to Adj
  tomaru2 <- list
  type_regulation <- tomaru2[,3]
  netINFO <- tomaru2[,-3]
  # unique(netINFO[,1])   # 只有TF  94 45 75 18  2 85 57 63 15 68
  # unique(netINFO[,2])   # 100
  adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
  
  ## give the gene name, here is number
  genename <- seq(1,genenum,1)
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
  sum(adj_ref != 0)    # 990
  
  ## Perf on all TF
  weightMatrix0 <- adj_ref
  
  ## 2023.2.12   --- New
  ## edge link
  # geneselect <- LogBTF$links[which(cellnum == k)]
  # library(GENIE3)
  # linkListNum <- min(getLinkList(adj_ref, reportMax = geneselect)[,3])
  # weightMatrix0 <- adj_ref
  # weightMatrix0[weightMatrix0 < linkListNum] = 0
  # weightMatrix0[weightMatrix0 >= linkListNum] = 1
  
  
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  ## load ground-truth Adj matrix
  ## input ground truth 258 edges
  
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
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
  # meanAUCAcc <- cbind(geneselect, AUC)
  meanAUCAcc <- AUC
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  # colnames(meanAUCAccAll) <- c("link", "AUROC", "AUPR")
  colnames(meanAUCAccAll) <- c("AUROC", "AUPR")
  
  
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
output



## save
# setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_DeepSEM/",sep=""))
# if(method == 2){
#   write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }else if(method == 3){
#   write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }else if(method == 1){
#   write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
# }




####################################################################################
####################################################################################
## 2023.2.16  List to Adj marix
# Using 20 node to generate SERGIO data -----------------------------------



## input data
rm(list = ls())


## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'

##
set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
##
SIGN = 1

dataclass <- 'count_matrix'
# dataclass <- 'exper_clean'
genenum <- 20
celltyp <- 1
cellbum <- 10



## Notice that the gene node is brgin from 0, but not 1

setwd(paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', 
            sep = ""))
inter <- read.csv(file = "Interaction_cID_4SIGN_20node.csv", header=T)

## Convert 2 to -1, Keep 1 is 1
inter[which(inter$type == 2), 3] <- -1

## List to Adj
tomaru2 <- inter
type_regulation <- tomaru2[,3]
netINFO <- tomaru2[,-3]
adj_ref <- matrix(0, nrow = genenum, ncol = genenum)

## give the gene name, here is number
# genename <- seq(1,100,1)
genename <- seq(0,(genenum-1),1)


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

## See how many edges
sum(abs(adj_ref))    # 36 

## Save adj matrix results
# write.csv(adj_ref, file = "Interaction_cID_4SIGNmatrix_20node.csv", row.names = F)


## Matrix to List
library(tidyfst)
colnames(adj_ref) <- rownames(adj_ref) <- seq(1,genenum,1)
## with sign
# List <- mat_df(adj_ref)

## without sign
List <- mat_df(abs(adj_ref))

## delet duplication
List <- as.matrix(List[-which(List[,1] == List[,2]),])
dim(List)
List[1,2]

List <- apply(List, 2, as.numeric)

# write.csv(List, file = "Interaction_cID_4SIGNlist.csv", row.names = F)
## delete ""
# write.table(List, file = "Interaction_cID_4SIGNlist.txt", row.names = F, quote=F)


## for GNILPR method
# setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/")
# write.table(List, file = "Interaction_cID_4SIGNlist_20node.txt", row.names = F, col.names = F, quote=F, sep = "\t")




####################################################################################
####################################################################################
## 2023.2.16 For 20 node Data --  List to Adj marix   --------------------------

## input data
rm(list = ls())


## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'

##
set.seed(123)
## 
method <- 2
## 0-???ų??Խ???
noDIAG = 1
##
SIGN = 1

dataclass <- 'count_matrix'
# dataclass <- 'exper_clean'
genenum <- 20
celltyp <- 1


## Notice that the gene node is brgin from 0, but not 1

setwd(paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', 
            sep = ""))
inter <- read.csv(file = "Interaction_cID_4SIGN_20node.csv", header=T)
## Convert 2 to -1, Keep 1 is 1
inter[which(inter$type == 2), 3] <- -1

## List to Adj
tomaru2 <- inter
type_regulation <- tomaru2[,3]
netINFO <- tomaru2[,-3]
adj_ref <- matrix(0, nrow = genenum, ncol = genenum)

## give the gene name, here is number
# genename <- seq(1,genenum,1)
genename <- seq(0,(genenum-1),1)


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

## See how many edges
sum(abs(adj_ref))    # 36 

## Save adj matrix results
# write.csv(adj_ref, file = "Interaction_cID_4SIGNmatrix.csv", row.names = F)


## Matrix to List
library(tidyfst)
colnames(adj_ref) <- rownames(adj_ref) <- seq(1,(genenum),1)
## with sign
# List <- mat_df(adj_ref)

## without sign
List <- mat_df(abs(adj_ref))

## delet duplication
List <- as.matrix(List[-which(List[,1] == List[,2]),])
dim(List)
List[1,2]

List <- apply(List, 2, as.numeric)

# write.csv(List, file = "Interaction_cID_4SIGNlist.csv", row.names = F)
## delete ""
# write.table(List, file = "Interaction_cID_4SIGNlist.txt", row.names = F, quote=F)


## for GNILPR method
# setwd("/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/")
setwd(paste(pathway, "R/Boolean/Data/DREAM3result_GNIPLR/SERGIOData/", sep = ""))
write.table(List, file = "Interaction_cID_4SIGNlist_20node.txt", row.names = F, col.names = F, quote=F, sep = "\t")



####################################################################################
####################################################################################
# 20 node Data -- Input data and add gene label for GNIPLR -------------------------------------------

rm(list = ls())

## set pathway
pathway <- '/home/lly/'
# pathway <- '/Users/lilingyu/E/PhD/'


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 20
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50,60)


for (k in cellnum) {
  # k <- 10
  
  ## load data
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  Data = t(read.csv(file = dataauto, header=F))
  
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
  
  colnames(datanum) <- seq(1,p,1)
  
  # dataauto1 <-  "/Users/lilingyu/E/PhD/R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/"
  dataauto1 <-  paste(pathway, "R/Boolean/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  # write.table(datanum, file = paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, ".txt", 
  #                                 sep=""),sep="\t", row.names = F, quote = F)
  
  
  ## must row.names = F
  write.table(datanum, file = paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, ".txt", 
                                    sep=""),sep="\t", row.names = F)
  
  print("Data GNIPLR")
  
  
}




# Then Run Python code  ---------------------------------------------------



####################################################################################
####################################################################################
# 2022.2.16 GNIPLR method 20 node -- Use Performance Function  -------------------

rm(list = ls())


## set pathway
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 20
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50,60)
# truthedges <- 258

## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2

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
  
  
  ## path and data
  # dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  dataauto1 <- paste(pathway,"R/Boolean/Data/DREAM3result_GNIPLR/SERGIOData/", sep = "")
  dataauto_sc <- paste(dataauto1, "Label_", dataclass,genenum,"_",celltyp, "_", k, "_lasso.txt", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t")
  result[is.na(result[,3]),3] <- 0
  # list <- result[,-4]
  
  
  
  
  ## List to Adj
  tomaru2 <- result
  type_regulation <- tomaru2[,3]
  netINFO <- tomaru2[,-c(3,4)]
  adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
  
  ## give the gene name, here is number
  genename <- seq(1,genenum,1)
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
  geneselect <- LogBTF$links[which(cellnum == k)]
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
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  setwd(datapath)
  adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix_20node.csv"))
  
  
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
output



# save
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_GNIPLR/",sep=""))
if(method == 2){
  write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 3){
  write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 1){
  write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}



####################################################################################
####################################################################################
####  Compute the AUROC  and other performance index of DeepSEM method
# ## 2023.2.16 DeepSEM method --  For 20 nodes --  ------------------------
## Because only give the regulatory of All 10 Tf genes,
## So we test to comare the truth 258 edges and see the AUC values

rm(list = ls())


## set pathways
pathway <- '/home/lly/'
dataLinux <- "Data"

# pathway <- '/Users/lilingyu/E/PhD/'
# dataLinux <- "RShell/Data"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 20
celltyp <- 1
# cellnum <- 10
cellnum <- c(10,20,30,40,50,60)
truthedges <- 36


## load function
source(paste(pathway,'R/Boolean/R/SINCERITIES functions/PermanceFunction.R', sep = ""))
method = 2

## input LogBTF information
## Dont select the regulatory link, according to LogBTF method
## We select the number of edges of 158, the truth regulatory edges

# setwd(paste(pathway,"R/Boolean/",dataLinux,"/DREAM3result_LogBTF/",sep=""))
# LogBTF <- read.csv(file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""), row.names = 1)




AUCALL <- c()
AUCALL0 <- c()
ResultAll <- c()
ResultAll0 <- c()
meanAUCAccAll <- c()




for (k in cellnum) {
  # k <- 10
  
  
  ## path and data
  # dataauto1 <- paste(pathway,"R/Boolean/Rshell/Data/DREAM3result_DeepSEM/SERGIOData/", sep = "")
  dataauto1 <- paste(pathway,"Python/DeepSEM-master/out/", sep = "")
  dataauto_sc <- paste(dataauto1, "GRN_inference_result_", dataclass,genenum,"_",celltyp, "_", k, ".tsv", 
                       sep="")
  result <- read.table(dataauto_sc, sep = "\t", header = T)
  
  ## delete "G"
  result$TF <- stringr::str_split_fixed(result$TF, "G", 2)[,2]
  result$Target <- stringr::str_split_fixed(result$Target, "G", 2)[,2]
  ## convert to number and form 0 to 1
  # as.numeric(result$Target)+1
  list <- result
  list[,1] <- as.numeric(result$TF)+1
  list[,2] <- as.numeric(result$Target)+1
  
  
  ## List to Adj
  tomaru2 <- list
  type_regulation <- tomaru2[,3]
  netINFO <- tomaru2[,-3]
  # unique(netINFO[,1])   # 只有TF  94 45 75 18  2 85 57 63 15 68
  # unique(netINFO[,2])   # 100
  adj_ref <- matrix(0, nrow = genenum, ncol = genenum)
  
  ## give the gene name, here is number
  genename <- seq(1,genenum,1)
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
  sum(adj_ref != 0)    # 990
  
  ## Perf on all TF
  weightMatrix0 <- adj_ref
  
  ## 2023.2.12   --- New
  ## edge link
  # geneselect <- LogBTF$links[which(cellnum == k)]
  geneselect <- truthedges
  library(GENIE3)
  linkListNum <- min(getLinkList(adj_ref, reportMax = geneselect)[,3])
  weightMatrix0 <- adj_ref
  weightMatrix0[weightMatrix0 < linkListNum] = 0
  weightMatrix0[weightMatrix0 >= linkListNum] = 1
  
  
  
  ## 0-????
  SIGN = 0
  # SIGN = 1
  
  
  ## load ground-truth Adj matrix
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  setwd(datapath)
  adj_gold <- as.matrix(read.csv(file = "Interaction_cID_4SIGNmatrix_20node.csv"))
  
  
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
  # meanAUCAcc <- AUC
  meanAUCAccAll <- as.matrix(rbind(meanAUCAccAll, meanAUCAcc))
  colnames(meanAUCAccAll) <- c("link", "AUROC", "AUPR")
  # colnames(meanAUCAccAll) <- c("AUROC", "AUPR")
  
  
  performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
  Acc0 <- performance$Acc
  Recall0 <- performance$Recall
  Pre0 <- performance$Pre
  FPR0 <- performance$FPR
  Fmeasure0 <- performance$Fmeasure
  Result0 <- cbind(Acc0,Recall0,Pre0,FPR0,Fmeasure0)
  ResultAll0 <- as.matrix(rbind(ResultAll0, Result0))
  colnames(ResultAll0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")
  
  
  print("DeepSEM method")
  
}


## output
output <- cbind(meanAUCAccAll, ResultAll0)
output



## save
setwd(paste(pathway,"R/Boolean/", dataLinux,"/DREAM3result_DeepSEM/",sep=""))
if(method == 2){
  write.csv(output, file = paste("ENet_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 3){
  write.csv(meanAUCAccAll, file = paste("Lasso_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}else if(method == 1){
  write.csv(meanAUCAccAll, file = paste("Ridge_output_SERGIO_",dataclass,genenum,"_",celltyp,".csv",sep=""))
}
