#' Predict new concentrations based on a trained model
#'
#' This function takes a trained calibration model and a matrix of elemental concentrations
#' and predicts elemental concentrations for the matrix.
#' 
#' @param mdl A model object created by xrfcal::calib()
#' @param newX An elemental counts matrix to make concentration predictions based on the trained model
#' @return A matrix of predicted elemental concentrations 
#' @export


#' @export

pred <- function(mdl, newX) {
  
  if (mdl$Method == "LRCE") {
    newX <- newX[, colnames(mdl$Y)]
  }
  
  #check if there are zeros if yes, replace zeros
  nzX <- newX
  
  if (any(nzX == 0))
  {
    nzX <- zCompositions::multRepl(nzX, dl = mdl$dl, label = 0)
  }
  
  
  
  if (mdl$Method == "LRCE") {
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))
    for (i in 1:ncol(mdl$Y)) {
      bdi <- which(colnames(mdl$Y)==mdl$BestDenom[i])
      
      Xt <- alr(mdl$X, bdi)
      Yt <- alr(mdl$Y, bdi)
      newXt <- alr(nzX, bdi)
      Yhatj <- matrix(nrow=nrow(newXt),ncol=ncol(mdl$Y)-1)
      
      
      for (j in 1:(ncol(mdl$Y) - 1)) {
        model <- stats::lm(Yt[,j] ~ as.matrix(Xt[,j]))
        Yhatj[, j] <-
          as.matrix(newXt[,j]) %*% model$coefficients[2] + model$coefficients[1]
      }
      Yhatj <- invalr(Yhatj,bdi)
      Yhat[,i] <- Yhatj[,i]
    } 
  }
  
  if (mdl$Method == "MLR") {
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))
    for (i in 1:ncol(mdl$Y)) {
      bdy <- which(colnames(mdl$Y)==mdl$BestDenom[i])
      bdx <- which(colnames(mdl$X)==mdl$BestDenom[i])
      Xt <- alr(mdl$X, bdx)
      Yt <- alr(mdl$Y, bdy)
      newXt <- alr(nzX, bdx)
      
      Yhatj <- matrix(nrow=nrow(newXt),ncol=ncol(mdl$Y)-1)
      
      
      for (j in 1:(ncol(mdl$Y) - 1)) {
        model <- stats::lm(Yt[,j] ~ as.matrix(Xt))
        Yhatj[, j] <-
          as.matrix(newXt) %*% model$coefficients[2:(ncol(Xt)+1)] + model$coefficients[1]
      }
      Yhatj <- invalr(Yhatj,bdy)
      Yhat[,i] <- Yhatj[,i]
    } 
  }
  
  if (mdl$Method == "Cubist"){
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))
    for (i in 1:ncol(mdl$Y)) {
      bdy <- which(colnames(mdl$Y)==mdl$BestDenom[i])
      bdx <- which(colnames(mdl$X)==mdl$BestDenom[i])
      Xt <- alr(mdl$X, bdx)
      Yt <- alr(mdl$Y, bdy)
      newXt <- alr(nzX, bdx)
      
      Yhatj <- matrix(nrow=nrow(newXt),ncol=ncol(mdl$Y)-1)
      
      
      for (j in 1:(ncol(mdl$Y) - 1)) {
        model <- Cubist::cubist(y = Yt[,j], x = Xt)
        Yhatj[, j] <- stat::predict(model,newXt)
      }
      Yhatj <- invalr(Yhatj,bdy)
      Yhat[,i] <- Yhatj[,i]
    } 
  }
  
  if (mdl$Method == "RF"){
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))
    for (i in 1:ncol(mdl$Y)) {
      bdy <- which(colnames(mdl$Y)==mdl$BestDenom[i])
      bdx <- which(colnames(mdl$X)==mdl$BestDenom[i])
      Xt <- alr(mdl$X, bdx)
      Yt <- alr(mdl$Y, bdy)
      newXt <- alr(nzX, bdx)
      
      Yhatj <- matrix(nrow=nrow(newXt),ncol=ncol(mdl$Y)-1)
      
      
      for (j in 1:(ncol(mdl$Y) - 1)) {
        model <-
          randomForest::randomForest(
            x = Xt,
            y = Yt[,j],
            mtry = mdl$mtry,
            ntree = mdl$ntree
          )
        Yhatj[, j] <- stats::predict(model, newXt)
      }
      Yhatj <- invalr(Yhatj,bdy)
      Yhat[,i] <- Yhatj[,i]
    } 
  }
  
  Yhat <- Yhat/rowSums(Yhat)
  Yhat <- as.data.frame(Yhat)
  colnames(Yhat) <- colnames(mdl$Y)
  
  return(Yhat)
}