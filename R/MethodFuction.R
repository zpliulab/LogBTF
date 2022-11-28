## 2022.11.28 MethodFunction.R
## It includes some methods in this work "LogBTF"
## logBTF-logBTFfunction; SCODE-SCODEsc function

##############################################################################################################################
##############################################################################################################################


# LogBTFmainfunction ------------------------------------------------------
## 2022.11.28 LogBTF function 

logBTFfunction <- function(method, nfold, numGENES, X_matrix, DISTANCE_matrix){
  # def method
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
  
  ## LOOCV settings
  foldid <- 1:nfold
  keep <- TRUE
  pred_lambda_min <- matrix(0, nrow = numGENES+1, ncol = numGENES)
  lambda_res <- vector()
  alpha_res <- vector()
  AUCall0 <- c()
  ptrainall0 <- c()
  
  library(glmnet)
  library(pROC)
  for (gi in 1:numGENES) {
    
    # gi <- 3
    AUCall <- c()
    ptrainall <- c()
    cverrorall <- c()
    
    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
    theta <- matrix(data=0,nrow = dim(X_matrix)[2]+1,ncol = length(alphas))
    
    
    # for (test in 1:length(alphas)) {
    test <- 1
    Y_vector <- DISTANCE_matrix[2:(num_time_points),gi]
    
    # if Y exist one 1/0, use noise 0/1 data.
    if(sum(Y_vector) == 0 | sum(Y_vector) == n){
      glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
      coef <- glm.fit$coefficients
      coef[is.na(coef)] <- 0
      pred_lambda_min[,gi] <- coef
      
      ptrain <- Y_vector    
      ptrainall <- cbind(ptrainall, ptrain)
      ptrainall0 <- cbind(ptrainall0, ptrainall)
    }else if(sum(Y_vector) == 1 | sum(Y_vector) == (n-1)){
      
      glm.fit <- glm(Y_vector~., xglm, family = "binomial", control = list(maxit = 100))
      coef <- glm.fit$coefficients
      coef[is.na(coef)] <- 0
      pred_lambda_min[,gi] <- coef
      
      ptrain <- sgn(as.matrix(xglm01) %*% coef[-1] + coef[1])
      ptrainall <- cbind(ptrainall, ptrain)
      ptrainall0 <- cbind(ptrainall0, ptrainall)
      
      aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
      auc <- aucplot$auc
      AUCall <- cbind(AUCall, auc)
      AUCall0 <- cbind(AUCall0, AUCall)
      
    }else{
      ## glmnet
      if(noDIAG==1){
        CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,
                                nfolds = nfold, foldid = foldid, keep = keep, grouped = FALSE)
      }else{
        CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],
                                nfolds = nfold, foldid = foldid, keep = keep, grouped = FALSE)
      }
      
      # plot(CV_results)
      lambda[test] <- CV_results$lambda.min
      cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
      cverrorall <- cbind(cverrorall, cvERROR)
      coef.CV_results <- coef(CV_results, s='lambda.min')
      beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
      theta[coef.CV_results@i+1,test] = coef.CV_results@x
      
      
      theta[1,test] <- lambda*theta[1,test]
      ptrain <- sgn(as.matrix(X_matrix) %*% theta[-1,test] + theta[1,test])
      ptrainall <- cbind(ptrainall, ptrain)
      
      aucplot <- plot.roc(Y_vector, as.numeric(ptrain), print.auc=T)
      auc <- aucplot$auc
      AUCall <- cbind(AUCall, auc)
      
      
      minIdx <- max(which(cvERROR==min(cvERROR)))
      lambda_res[gi] <- lambda[minIdx]
      alpha_res[gi] <- alphas[minIdx]
      pred_lambda_min[,gi] <- theta[,minIdx]
      
      ptrainall0 <- cbind(ptrainall0, as.matrix(ptrainall[,minIdx]))
      AUCall0 <- cbind(AUCall0, AUCall[,minIdx])
    }
    
    print(gi)
    
  }
  return(list(AUCall0, pred_lambda_min, ptrainall0))
}


##############################################################################################################################
##############################################################################################################################


# SCODEmainfunction ------------------------------------------------------
## 2022.11.28 SCODE function 


SCODEsc <- function(filelist, dataname){
  for (i in 1:length(filelist)) {
    
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
    datahat <- t(binarizeMatrix(t(datanum),method="kMeans")[,1:n])
    datahatA <- as.matrix(datanum) ** as.matrix(datahat)
    datahatA[which(datahatA == 1)] <- 0
    datanum <- t(datahatA)
    
    ## use MASS library to calculate pseudo inverse matrix.
    library(MASS)
    
    ## gene
    tfnum <- genenum
    pnum <- 4
    ## cell
    cnum <- 21
    
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
    pseudotime <- seq(1,21,1)
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
    
    #infer A
    B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
    for(i in 1:pnum){
      B[i,i] <- new_B[i]
    }
    invW <- ginv(W)
    A <- W %*% B %*% invW
    
    weightMatrix0 <- A
    
    ## Performance 
    SIGN = 0
    # SIGN = 1
    
    ## Gold standard GRN 
    datapath0 <- paste("/home/lly/R/Boolean/Data/DREAM/DREAM3 in silico challenge/Size", genenum, "/NetAdjMatrix/", sep="")
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
    
    
    source('/home/lly/R/Boolean/R/SINCERITIES functions/PermanceFunction.R')
    library(pracma)
    library(Matrix)
    AUCresult <- auc_from_ranks_TC_sign(adj_matrix, adj_gold, 1000)
    AUC <- cbind(AUCresult$AUROC,AUCresult$AUPR)
    meanAUCAcc <- cbind(geneselect, AUC)
    colnames(meanAUCAcc) <- c("link", "AUROC", "AUPR")
    
    
    performance <- performance_from_ranks_TC_sign(adj_matrix, adj_gold, 0)
    Result0 <- cbind(performance$Acc,performance$Recall,performance$Pre,performance$FPR,performance$Fmeasure)
    colnames(Result0) <- c("Acc0","Recall0","Pre0","FPR0","Fmeasure0")
    
    
    ## output
    output <- cbind(meanAUCAcc, Result0)
    setwd(paste("/home/lly/R/Boolean/Data/DREAM3result_SCODE01/",dataname,filenum,"Node",genenum, sep=""))
    
    
    if(method == 2){
      write.csv(output, file = "ENet_output.csv")
    }else if(method == 3){
      write.csv(meanAUCAccAll, file = "Lasso_output.csv")
    }else if(method == 1){
      write.csv(meanAUCAccAll, file = "Ridge_output.csv")
    }
    
    print("SCODE method")
    output
    
    
  }
  
}



##############################################################################################################################
##############################################################################################################################

# 2022.7.24  ------------------------------------------------------
## Modified GENIE3 method
## The same as that in PermanceFunction.R code, we cpopy it here.

myGENIE3 <- function(expr.matrix, K="sqrt", nb.trees=1000, input.idx=NULL, importance.measure="IncNodePurity", seed=NULL, trace=TRUE, ...) {
  
  ## input data
  # expr.matrix <- exprMatrix
  # K="sqrt"
  # nb.trees=1000
  # input.idx=NULL
  # importance.measure="IncNodePurity"
  # seed=NULL
  # trace=TRUE
  
  
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  # Check if nodesize parameter is in the input arguments
  args <- list(...)
  nodesize.in.args <- "nodesize" %in% names(args)
  # transpose expression matrix to (samples x genes)
  expr.matrix <- t(expr.matrix)
  # setup weight matrix
  num.samples <- dim(expr.matrix)[1]
  num.genes <- dim(expr.matrix)[2]
  gene.names <- colnames(expr.matrix)
  weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- gene.names
  # get number of input genes, names of input genes
  if (is.null(input.idx)) {
    input.gene.names <- gene.names
  } else {
    # input gene indices given as integers
    if (is.numeric(input.idx)) {
      input.gene.names <- gene.names[input.idx]
      # input gene indices given as names
    } else {
      input.gene.names <- input.idx
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) {
        for (missing.gene.name in missing.gene.names) {
          cat(paste("Gene ", missing.gene.name,
                    " was not in the expression matrix\n", sep=""))
        }
        stop("Aborting computation")
      }
    }
  }
  # compute importances for every target gene
  
  AUCall <- c()
  ptrainall <- c()
  
  for (target.gene.idx in seq(from=1, to=num.genes)) {
    # target.gene.idx <- 1
    if (trace) {
      cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
      flush.console()
    }
    target.gene.name <- gene.names[target.gene.idx]
    # remove target gene from input genes
    these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
    num.input.genes <- length(these.input.gene.names)
    x <- expr.matrix[,these.input.gene.names]
    y <- expr.matrix[,target.gene.name]
    # set mtry
    if (class(K) == "numeric") {
      mtry <- K
    } else if (K == "sqrt") {
      mtry <- round(sqrt(num.input.genes))
    } else if (K == "all") {
      mtry <- num.input.genes
    } else {
      stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
    }
    if (trace) {
      cat(paste("K = ", mtry,", ", nb.trees, " trees\n\n",
                sep=""))
      flush.console()
    }
    if (importance.measure == "%IncMSE") {
      if (nodesize.in.args) {
        rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE, ...)
      } else {
        # By default, grow fully developed trees
        rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE, nodesize=1, ...)
      }
      
    } else {
      # Normalize output
      y <- y / sd(y)
      if (nodesize.in.args) {
        rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=FALSE, ...)
      } else {
        # By default, grow fully developed trees
        rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=FALSE, nodesize=1, ...)
      }
    }
    
    ## predict
    pred = predict(rf,x)
    ptrain <- pred
    ptrainall <- cbind(ptrainall, ptrain)
    
    
    aucplot <- plot.roc(y, pred, print.auc=T)
    auc <- aucplot$auc
    AUCall <- cbind(AUCall, auc)
    
    
    im <- importance(rf)[,importance.measure]
    im.names <- names(im)
    weight.matrix[im.names, target.gene.name] <- im
  }
  return(list((weight.matrix / num.samples), ptrainall, AUCall))
}       


