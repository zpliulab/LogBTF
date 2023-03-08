
## 2023.2.9 use Python in R
## Apply the DeepSEM method fo Prof. Zeng Jianyang

## 2023.2.14 Using SERGIO data by DeepSEM method
## Data progrecing

## clear workspace
rm(list = ls())
## set seed
set.seed(123)
## set parameter
pathname = "/home/lly/"
# pathname = "/Users/lilingyu/E/PhD/"



################################################################################
################################################################################

## define function
DeepESMsc <- function(filelist, dataname){
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
    datapath <- paste(pathname,"R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
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
    class(datahatA)
    
    ##.datahatA is simulated scRNA-seq data;
    ## datanum. is bulk data from DREAM3
    
    
    ## Save bukl data for Python input
    # pathshort = paste(pathname,"Paper/Paper7/Bioinformatics/AddMethod/ncs曾坚阳/",sep = "")
    # setwd(paste(pathshort,"DeepSEM-master/DREAM3_data/bulkData",sep = ""))
    setwd(paste(pathname,"Python/DeepSEM-master/DREAM3_data/bulkData",sep = ""))
    dataauto_bulk <- paste("bulkData_",dataname,filenum,"_Node",genenum,".csv", sep="")
    ## Python's input file is txt files.
    write.csv(datanum, file = dataauto_bulk, sep = ',', quote = F)


    ## Save single-cell data for Python input
    # setwd(paste(pathshort,"DeepSEM-master/DREAM3_data/scData",sep = ""))
    setwd(paste(pathname,"Python/DeepSEM-master/DREAM3_data/scData",sep = ""))
    dataauto_sc <- paste("scData_",dataname,filenum,"_Node",genenum,".csv", sep="")
    ## Python's input file is txt files.
    write.csv(datahatA, file = dataauto_sc, sep = ',', quote = F)
    
    
    
    ## Obtain Ground-truth network
    filename = paste(pathname,'R/Boolean/Data/DREAM/DREAM3 in silico challenge/', sep = "")
    ## Link file
    setwd(paste(filename,'Size',genenum,'/DREAM3 gold standards/',sep = ""))
    fw <- read.table(paste('DREAM3GoldStandard_InSilicoSize',genenum,'_',dataname,filenum,'.txt', sep = ""),sep = "\t")
    fw1 <- fw[,-3]
    colnames(fw1) <- c("Gene1", "Gene2")
    class(fw1)
    ## Save the csv files
    # pathshort = paste(pathname,"Paper/Paper7/Bioinformatics/AddMethod/ncs曾坚阳/",sep = "")
    # setwd(paste(pathshort,"DeepSEM-master/DREAM3_data/label/Size",genenum,sep = ""))
    setwd(paste(pathname,"Python/DeepSEM-master/DREAM3_data/label/Size",genenum,sep = ""))
    dataauto_label <- paste("label_",dataname,filenum,"_Node",genenum,".csv", sep="")
    write.csv(fw1, file = dataauto_label, row.names = F, quote = F)
    
    
    ## print method
    print("DeepESM method")
  }
  
  
  
}


## input data parameter
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))

## Run GNIPLRsc function, then Save bulk Data and scData
DeepESMsc(fileEcoli, "Ecoli")
DeepESMsc(fileYeast, "Yeast")


