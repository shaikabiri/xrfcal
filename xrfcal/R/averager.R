
quant <- function(X,p){
  return(as.numeric(stats::quantile(X,probs=p)))
}


averager <- function(X,Y){
  Xm <- stats::aggregate(X,list(Y[,1]),FUN=mean)
  Xm <- Xm[,-1]
  Ym <- stats::aggregate(Y,list(Y[,1]),FUN=mean)
  Ym <- Ym[,-1]  
  
  X1 <- matrix(ncol=ncol(Xm),nrow(Xm),data = 0)
  X3 <- matrix(ncol=ncol(Xm),nrow(Xm),data = 0)
  
  Y1 <- matrix(ncol=ncol(Ym),nrow(Ym),data = 0)
  Y3 <- matrix(ncol=ncol(Ym),nrow(Ym),data = 0)
  
  colnames(X1)<-colnames(Xm)
  colnames(X3)<-colnames(Xm)
  colnames(Y1)<-colnames(Ym)
  colnames(Y3)<-colnames(Ym)
  
  for (i in 1:ncol(X)){
    X1[,i] <- stats::aggregate(X[,i],list(Y[,1]),FUN = quant, p=0.25)[,2]
    X3[,i] <- stats::aggregate(X[,i],list(Y[,1]),FUN = quant, p=0.75)[,2]
  }
  
  for (i in 1:ncol(Y)){
    Y1[,i] <- stats::aggregate(Y[,i],list(Y[,1]),FUN = quant, p=0.25)[,2]
    Y3[,i] <- stats::aggregate(Y[,i],list(Y[,1]),FUN = quant, p=0.75)[,2]
  }
  
  X <- rbind(Xm,X1,X3)
  Y <- rbind(Ym,Y1,Y3)
  out <- list(X,Y)
  return(out)
}