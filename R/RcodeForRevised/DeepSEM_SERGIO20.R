
## 2023.2.9 use Python in R
## Apply the DeepSEM method fo Prof. Zeng Jianyang

## 2023.2.14 Using SERGIO data by DeepSEM method
## Data progrecing

## 2023.2.16 Copy from SeepSEM_SERGIO.py, here we use 20 node SERGIO data by DeepSEM method



## Data progrecing

## clear workspace
rm(list = ls())
## set seed
set.seed(123)
## set parameter
pathway = "/home/lly/"
# pathway = "/Users/lilingyu/E/PhD/"


## Parameter setting
# dataclass <- 'count_matrix'
dataclass <- 'exper_clean_'
genenum <- 20
celltyp <- 1
# cellbum <- 10
cellnum <- c(10,20,30,40,50,60)



################################################################################
################################################################################


for (k in cellnum) {
  # k <- 10
  
  ## load data
  datapath <- paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = "")
  dataauto <- paste(datapath,dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  Data = t(read.csv(file = dataauto, header=F))
  # colnames(Data) <- seq(1,genenum,1)
  
  colnames(Data) <- stringr::str_c("G",seq(0,(genenum-1),1))
  rownames(Data) <- stringr::str_c("S",seq(1,dim(Data)[1],1))

  ## GNIPLR need sample * gene. And is ready to use.
  class(Data)
  
  ## Save  data for Python input
  setwd(paste(pathway,"Python/DeepSEM-master/SERGIO_data",sep = ""))
  dataauto <- paste("Data_",dataclass,genenum,"_",celltyp, "_", k, ".csv",sep="")
  ## Python's input file is txt files.
  write.csv(Data, file = dataauto, sep = ',', quote = F)
  # write.csv(Data, file = dataauto)
  
  ## print method
  print("DeepESM method")
}
  

## Obtain Ground-truth network 
setwd(paste(pathway,'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/', sep = ""))
## Link file
## all links
# fw <- read.table("Interaction_cID_4SIGNlist.txt",sep = "\t")
fw <- read.csv("Interaction_cID_4SIGN_20node.csv",sep = ",")
fw1 <- fw[,-3]
colnames(fw1) <- c("Gene1", "Gene2")
fw1$Gene1 <- stringr::str_c("G",fw1[,1])
fw1$Gene2 <- stringr::str_c("G",fw1[,2])

class(fw1)
## Save the csv files
setwd(paste(pathway,"Python/DeepSEM-master/SERGIO_data",sep = ""))
write.csv(fw1, file = "Interaction_cID_4SIGNlistDeep_20node.csv", row.names = F, quote = F)
  

################################################################################
################################################################################

# Performance  ------------------------------------------------------------




