

################################################################################################################
################################################################################################################

## 2022.11.28, SCODE code for simulated single-cell data, which is based on SCODE-master package

rm(list = ls())

## load function
source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
source('/home/lly/R/Boolean/R/MethodFuction.R')

## parameters
set.seed(123)

fileEcoli <- list(c(10,1,4,45),c(50,1,12,576),c(100,1,29,1708),
                  c(10,2,3,48),c(50,2,6,585), c(100,2,11,1996))

fileYeast <- list(c(10,1,1,53),c(50,1,10,415),c(100,1,19,1497),
                  c(10,2,4,59),c(50,2,20,557),c(100,2,4,1700),
                  c(10,3,2,40),c(50,3,16,417),c(100,3,42,1718))

## save results
SCODEsc(fileEcoli, "Ecoli")
SCODEsc(fileYeast, "Yeast")


################################################################################################################
################################################################################################################

## 2022.11.28, SCODE code for real scRNA-seq data, which is based on SCODE-master package



