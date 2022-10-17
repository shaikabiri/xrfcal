#' @title Predict new concentrations based on a trained calibration model
#'
#' @description This function takes a trained calibration model and a matrix of elemental counts
#' and predicts elemental concentrations for the matrix.
#' @details 
#' The procedure for predicting on new data in `pred` function is as followed:
#' \enumerate{
#'  \item take in the output of calibration function `xrfcal`.
#'  \item if there are zeroes in `newX`, replace zeroes with `multRepl` from `zCompositions`.
#'  \item for each element to be predicted, transform newX to alr with best denominator for that element.
#'  \item predict all elements for each best deonimnator and transform back to simplex but only get the element of interest as output.
#'  \item aggregate all element predictions and close them to 1.
#'  \item export predicted concentrations.
#' }
#' 
#' @examples 
#' #Remove count columns with more than 20 percent zeroes
#'
#' zeroes <- apply(xrf$X,2,function(x){return(length(which(x==0))/length(x))})
#' ind <- which(zeroes>0.2)
#' if (length(ind)>0)
#'   X <- xrf$X[,-ind]
#'
#' #Replace zeros in X
#' dl <- apply(X, 2, agrmt::minnz)
#' nzX <- zCompositions::multRepl(X, dl = dl, label = 0)
#' 
#' #Reduce dataset
#' xrfRed <- averager(nzX,xrf$Y)
#' Xred <- xrfRed[[1]]
#' Yred <- xrfRed[[2]]
#' 
#' #Split dataset into train and test sets
#' set.seed(123)
#' ind <- sample.int(ceiling(dim(Xred)[1]/3))
#' testX <- Xred[ind,]
#' testY <- Yred[ind,]
#' trainX <- Xred[-ind,]
#' trainY <- Yred[-ind,]
#' 
#' #Build a calibration model based on training data
#' mdl <- xrfcal(trainX,trainY, method = "RF", reduce=FALSE)
#' 
#' #Predict for testing data
#' Yhat <- pred(mdl,testX) 
#' 
#' #Calculate R2 and RMSE for test set
#' res <- caret::postResample(Yhat$K,testY$K)
#' 
#' #Plot predicted potassium concentrations vs actual concentrations for test set
#' plot(y=Yhat$K,x=testY$K,ylim = c(0,max(cbind(Yhat$K,testY$K))),
#'     xlim = c(0,max(cbind(Yhat$K,testY$K))),ylab="Calibrated Concentrations for K (mass fraction)",
#'     xlab ="Actual Concentrations for K (mass fraction) ")
#' lines(-1:1,-1:1)
#' text(0.2*max(cbind(Yhat$K,testY$K)),0.8*max(cbind(Yhat$K,testY$K)),
#'     paste("R2 =",signif(res[2],2),"\n", "RMSE =",signif(res[1],2)))
#' @references
#' Weltje, G. J., & Tjallingii, R. (2008). Calibration of XRF core scanners for quantitative geochemical logging of sediment cores: Theory and application. Earth and Planetary Science Letters, 274(3), 423-438. https://doi.org/10.1016/j.epsl.2008.07.054
#'  
#' Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5-32. https://doi.org/10.1023/A:1010933404324
#' 
#' Quinlan, J.R. (1993). Combining Instance-Based and Model-Based Learning. ICML.
#' 
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821 
#'
#' @param mdl A model object created by xrfcal::calib().
#' @param newX An elemental counts matrix in simplex space to make concentration predictions on.
#' @return A matrix of predicted elemental concentrations in simplex space.
#' @export


#' @export

pred <- function(mdl, newX) {
  
  
  
  newX <- newX[,colnames(mdl$X)]
  
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
  
  if (mdl$Method == "ENET") {
    #Yhat for all elements with different denoms
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))

    
    for (i in 1:ncol(mdl$Y)) {
      bdy <- which(colnames(mdl$Y)==mdl$BestDenom[i])
      bdx <- which(colnames(mdl$X)==mdl$BestDenom[i])
      Xt <- alr(mdl$X, bdx)
      Yt <- alr(mdl$Y, bdy)
      newXt <- alr(nzX, bdx)
      
      Yhatj <- matrix(nrow=nrow(newXt),ncol=ncol(mdl$Y)-1)
      
      
      for (j in 1:(ncol(mdl$Y) - 1)) {
        m <- glmnet::cv.glmnet(as.matrix(Xt),Yt[,j])
        model <- glmnet::glmnet(as.matrix(Xt), Yt[,j], nlambda = 25, alpha = 0, family = 'gaussian', lambda = m$lambda.min)
        Yhatj[,j] <- glmnet::predict.glmnet(model,as.matrix(newXt))
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
        Yhatj[, j] <- stats::predict(model,newXt)
      }
      Yhatj <- invalr(Yhatj,bdy)
      Yhat[,i] <- Yhatj[,i]
    } 
  }
  
  if (mdl$Method == "RF"){
    Yhat <- matrix(nrow=nrow(nzX),ncol=ncol(mdl$Y))
    
    
    for (i in 1:ncol(mdl$Y)) {
      print(i)
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