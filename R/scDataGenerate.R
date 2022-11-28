## 2022.11.13 plot the distribution bar figure before binary or not.


# Chose data ------------------------------------------------------------------

rm(list = ls())


## Begin. Ecoli  Yeast mute data to plots
## output scRNA-seq data and mute data to plots
library(BiTrinA)
DatascRNA <- function(dataname, l=21, genenum, filenum, k){
  
  ## load data
  datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size",genenum,"/Data without noise/",sep="")
  dataauto <- paste(datapath,"InSilicoSize",genenum,dataname,filenum,"-nonoise-trajectories.tsv",sep="")
  Data = as.matrix(read.table(file = dataauto, header=T))
  run <- dim(Data)[1]/l
  
  ## as numeric
  datanum <- apply(Data[c(((k-1)*l+1):(k*l)),-1],2, as.numeric)  
  n <- dim(datanum)[1]
  p <- dim(datanum)[2]
  datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
  
  ## 1 replace 0
  datahatA <- as.matrix(datanum) ** as.matrix(datahat)
  datahatA[which(datahatA == 1)] <- 0
  
  ## time 
  pseudotime <- seq(1,21,1)
  pseudotime <- pseudotime/max(pseudotime)
  datahatAtime <- cbind(pseudotime, datahatA)
  
  ## ground trouth -- Gold standard GRN 
  datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
  adj_gold0 <- paste(datapath0," InSilicoSize",genenum,dataname,filenum,"-adj .csv",sep="")
  adj_gold <- as.matrix(read.csv(file = adj_gold0))


  SIGN <- 0
  # no self-loof, direct
  if(SIGN == 1){
    adj_gold[row(adj_gold) == col(adj_gold)] <- 0
  }
  # undirect
  if(SIGN == 0){
    adj_gold[which(adj_gold != 0)] <- 1
  }
  
  ## add first col
  DataOld <- data.frame(cbind(as.matrix(colnames(datanum)), t(datanum)))
  DataNew <- data.frame(cbind(as.matrix(colnames(datanum)), t(datahatA)))
  colnames(DataNew) <- colnames(DataOld) <- c("label", rownames(datahat))
  
  ## Progreeding 
  data1 <- DataOld %>% melt(id='label')  
  data1$Type <- "Old"
  data2 <- DataNew %>% melt(id='label')  
  data2$Type <- "New"
  
  ## as mumeric for character
  data12 <- rbind(data1, data2)
  data21 <- as.data.frame(lapply(data12, as.numeric))
  data12$value <- data21$value
  
  return(list(datahatAtime, adj_gold, data12))
}

 
## Plot function
library(ggplot2)
library(latex2exp)
FrenPlot <- function(data, dataname, genenum, filenum){
  ggplot(data,aes(x=value)) +
    geom_histogram(aes(y=after_stat(count / sum(count)),
                       fill=Type),
                   bins = 150,
                   alpha = 0.9) + ## alpha - color deep
    scale_fill_manual(values = c("New"="#a3cd5b",
                                 "Old"="#8ea0cc"),
                      labels=c("New"="scRNA-seq",
                               "Old"="Bulk RNA-seq")) +
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(),
          legend.position = c(0.1,0.9),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank(),
          legend.justification = c(0,1)) +
    labs(x="Expression value", y="Frequency (%)") + 
    labs(title = paste(dataname, filenum, " with node ", genenum, sep="")) -> p
  
  return(p)
}


## parameters new
fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))
fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))
filelist <- c(fileEcoli, fileYeast)


# Plot --------------------------------------------------------------------
pEcoli <- list()
for (i in 1:length(fileEcoli)) {
  file <- fileEcoli[[i]]
  data101 <- DatascRNA(dataname = "-Ecoli", l=21, genenum=file[1], filenum=file[2], k=file[3])
  pEcoli[[i]] <- FrenPlot(data101[[3]], dataname = "Ecoli", genenum=file[1], filenum=file[2])
}

pYeast <- list()
for (i in 1:length(fileYeast)) {
  file <- fileYeast[[i]]
  data101 <- DatascRNA(dataname = "-Yeast", l=21, genenum=file[1], filenum=file[2], k=file[3])
  pYeast[[i]] <- FrenPlot(data101[[3]], dataname = "Yeast", genenum=file[1], filenum=file[2])
}


# Togather ----------------------------------------------------------------

# setwd('/home/lly/R/Boolean/Data/DREAM3result/Figuresc')
# pdf("scRNA-seq VS bulk Ecoli.pdf",width = 12, height = 6)
cowplot::plot_grid(plotlist = pEcoli, nrow = 2)
# dev.off()


# setwd('/home/lly/R/Boolean/Data/DREAM3result/Figuresc')
# pdf("scRNA-seq VS bulk Yeast.pdf",width = 12, height = 9)
cowplot::plot_grid(plotlist = pYeast, nrow = 3)
# dev.off()


######################################################################################################
######################################################################################################

## save scRNA-seq data and ground truth A 
for (i in 1:length(fileEcoli)) {
  file <- fileEcoli[[i]]
  data101 <- DatascRNA(dataname = "-Ecoli", l=21, genenum=file[1], filenum=file[2], k=file[3])
  datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAMsc/Ecoli",file[2],"-", file[1],"/",sep="")
  dataauto <- paste(datapath,"data.txt", sep="")
  write.table(data101[[1]], file = dataauto, row.names = F, sep = '\t', col.names = F)
  dataA <- paste(datapath,"A.txt", sep="")
  write.table(data101[[2]], file = dataA, row.names = F, sep = '\t', col.names = F)
}

for (i in 1:length(fileYeast)) {
  # i <- 1
  file <- fileYeast[[i]]
  data101 <- DatascRNA(dataname = "-Yeast", l=21, genenum=file[1], filenum=file[2], k=file[3])
  datapath <- paste("/home/lly/R/Boolean/Data/DREAM/DREAMsc/Yeast",file[2],"-", file[1],"/",sep="")
  dataauto <- paste(datapath,"data.txt", sep="")
  write.table(data101[[1]], file = dataauto, row.names = F, sep = '\t', col.names = F)
  dataA <- paste(datapath,"A.txt", sep="")
  write.table(data101[[2]], file = dataA, row.names = F, sep = '\t', col.names = F)
}
