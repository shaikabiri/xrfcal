#' Calibrate XRF counts to element concentrations
#' This function takes an XRF counts matrix and reference concentrations
#' count matrix and outputs alphas and betas. 
#' @param X Elemental counts matrix
#' @param Y Reference concentrations 
#' @param method 
#' @return A matrix of alphas and betas

#' @export

calib <- function(X,Y,method = "RF"){
  #align the two matrices 
  names <- vector(mode="character")

  
  for (j in colnames(X)){
    if (j %in% colnames(Y)){
      print(j)
      names<- append(names,j)
    }
  }
  
  X <- X[,names]
  Y <- Y[,names]
  
  #check if there are zeros if yes, replace zeros
  nzX <- X
  if (any(nzX==0))
  {
    dl <- apply(nzX, 2, agrmt::minnz)
    nzX <- zCompositions::multRepl(nzX,dl=dl,label = 0)
  }
  
  R2 <- matrix(ncol=(ncol(X)),nrow=(ncol(X)))
  
  #do an LRCE for each element as denominator and save the results
  for (i in 1:ncol(nzX))
  {
    R2i <- vector(length = 30)
    Xalr <- alr(nzX,i)
    Yalr <- alr(Y,i)
    for (j in 1:ncol(Xalr)){
      #model <- lmodel2::lmodel2(Yalr[,j]~Xalr[,j])
      model <- robslopes::RepeatedMedian(Xalr[,j],Yalr[,j])
      yhat <- Xalr[,j]*model$slope + model$intercept
      rsquare <- cor(Yalr[,j],yhat)^2
      R2i[j] <- rsquare}
    R2i <- append(R2i,0,after=i-1)
    R2[,i]<- R2i
  }
}