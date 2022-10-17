
quant <- function(X,p){
  return(as.numeric(stats::quantile(X,probs=p)))
}

#' @title Reduce the dataset to a representative dataset
#'
#' @description This function takes counts matrix X and concentrations matrix Y and reduce it to a smaller representative dataset by taking the average, first and third quantile 
#' of counts for each unique concentration. This is useful when a unique concentration is associated with many count instances. 
#' @examples
#' #Replace zeros in X
#' dl <- apply(xrf$X, 2, agrmt::minnz)
#' nzX <- zCompositions::multRepl(xrf$X, dl = dl, label = 0)
#' 
#' #Reduce dataset
#' xrfRed <- averager(nzX,xrf$Y)
#' Xred <- xrfRed[[1]]
#' Yred <- xrfRed[[2]]
#' 
#' @param X A matrix of XRF counts. Zeros should be replaced beforehand. 
#' @param Y A matirx of reference concentrations. 
#' @return A list with two objects, reduced X, `obj[[1]]`, reduced Y, `obj[[2]]`
#' @export
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