#' Calibrate XRF counts to element concentrations
#' This function takes an XRF counts matrix and reference concentrations 
#' and returns a model object with R2s, RMSEs and  a predictor function for new data.
#' @examples
#' #Using the example data in the package build a calibration model
#' mdl <- calib(xrf[[1]],xrf[[2]], method = "MLR", denom = "Auto")
#' 
#' 
#' @param X Elemental counts matrix
#' @param Y Reference concentrations matrix
#' @param method Regression model. Either "LRCE", "MLR", "Cubist" or "RF". RF and Cubist generally take longer time to train
#' @param dl (Optional) Detection limits to use for zero replacement algorithm. A vector with D elements. If not provided will be estimated. 
#' @param oalr A boolean parameter. If true the model results and predictor function will be in alr space.
#' @param mtry The number of variables for random forest to sample randomly at each split 
#' @param ntree Number of estimators for random forest
#' @return A model object that includes R2s, RMSEs and a predictor function, pred.
#' @export

calib <-
  function(X,
           Y,
           method = "RF",
           dl = 0,
           oalr = FALSE,
           mtry = ncol(X),
           ntree = 20) {
    
    if (!(TRUE %in% (method == c("RF","MLR","Cubist","LRCE"))))
    {
      stop(paste(method," is not a known regression method."))
    }
    
    if (nrow(X)!=nrow(Y))
    {
      stop("X and Y should have the same number of rows.")
    }
    
    if ((length(dl)!=ncol(X))&dl!=0)
    {
      stop("Detection limits should have the same length as number of elements in X.")
    }
    
    if ((length(dl)!=ncol(X))&dl!=0)
    {
      stop("Detection limits should have the same length as number of elements in X.")
    }
    
    if ((typeof(oalr)!="logical"))
    {
      stop("Oalr must be boolean")
    }
    
    
    #align the two matrices
    X = as.matrix(X)
    Y = as.matrix(Y)
    if (method == "LRCE") {
      names <- vector(mode = "character")
      for (j in colnames(X)) {
        if (j %in% colnames(Y)) {
          names <- append(names, j)
        }
      }
      
      X <- X[, names]
      Y <- Y[, names]
    }
    
    #check if there are zeros if yes, replace zeros
    nzX <- X
    
    if (any(nzX == 0))
    {
      if (dl == 0) {
        dl <- apply(nzX, 2, agrmt::minnz)
      }
      nzX <- zCompositions::multRepl(nzX, dl = dl, label = 0)
    }
    
    if (nrow(unique(Y))!=nrow(Y)){
      out <- averager(nzX, Y)
      nzX <- out[[1]]
      Y <- out[[2]]
    }
    set.seed(123)
    folds <- caret::createFolds(Y[, 1], k = 3)
    R2 <- matrix(ncol = (ncol(Y)), nrow = (ncol(Y)))
    RMSE <- matrix(ncol = (ncol(Y)), nrow = (ncol(Y)))
    #do an LRCE for each element as denominator and save the results
    for (i in 1:ncol(Y))
    {
      Xalr <- alr(nzX, match(colnames(Y)[i], colnames(nzX)))
      Yalr <- alr(Y, i)
      
      if (oalr) {
        R2k <- matrix(nrow = ncol(Y) - 1, ncol = length(folds))
        RMSEk <- matrix(nrow = ncol(Y) - 1, ncol = length(folds))
      } else {
        R2k <- matrix(nrow = ncol(Y), ncol = length(folds))
        RMSEk <- matrix(nrow = ncol(Y), ncol = length(folds))
      }
      
      for (k in 1:length(folds)) {
        Xtest <- Xalr[folds[[k]],]
        Xtrain <- Xalr[-folds[[k]],]
        Ytest <- Yalr[folds[[k]],]
        Ytrain <- Yalr[-folds[[k]],]
        
        preds <- matrix(ncol = ncol(Ytest), nrow = nrow(Ytest))
        
        for (j in 1:ncol(Yalr)) {
          if (method == "LRCE")
          {
            model <- stats::lm(Ytrain[, j] ~ as.matrix(Xtrain[, j]))
            pred <-
              as.matrix(Xtest[, j]) %*% model$coefficients[2] + model$coefficients[1]
          }
          
          
          if (method == "Cubist")
          {
            model <- Cubist::cubist(y = Ytrain[,j], x = Xtrain)
            pred <- stats::predict(model, Xtest)
          }
          
          if (method == "MLR")
          {
            model <- stats::lm(Ytrain[, j] ~ as.matrix(Xtrain))
            pred <-
              as.matrix(Xtest) %*% model$coefficients[2:(ncol(nzX))] + model$coefficients[1]
          }
          
          if (method == "RF")
          {
            model <-
              randomForest::randomForest(
                x = Xtrain,
                y = Ytrain[,j],
                mtry = mtry,
                ntree = ntree
              )
            pred <- stats::predict(model, Xtest)
          }
          
          preds[, j] <- pred
        }
        
        if (oalr) {
          for (l in 1:ncol(preds))
          {
            RMSEk[l, k] <- caret::postResample(preds[, l], Ytest[, l])[1]
            R2k[l, k] <- caret::postResample(preds[, l], Ytest[, l])[2]
          }
          
          
        } else {
          preds <- invalr(preds, j)
          Ytest <- invalr(Ytest, j)
          
          for (l in 1:ncol(preds))
          {
            RMSEk[l, k] <- caret::postResample(preds[, l], Ytest[, l])[1]
            R2k[l, k] <- caret::postResample(preds[, l], Ytest[, l])[2]
          }
        }
        
        R2i <- apply(R2k, 1, mean)
        RMSEi <- apply(RMSEk, 1, mean)
        
      }
      
      if (oalr) {
        R2i <- append(R2i, 0, after = i - 1)
        RMSEi <- append(RMSEi, 0, after = i - 1)
      }
      
      R2[, i] <- R2i
      RMSE[,  i] <- RMSEi
      
    }
  

mdl <- list()
mdl[["X"]] <- nzX
mdl[["Y"]] <- Y
mdl[["R2"]] <- R2
mdl[["RMSE"]] <- RMSE
mdl[["Method"]] <- method
mdl[["BestR2"]] <- apply(R2,1,max)
mdl[["BestRMSE"]] <- apply(RMSE,1,max)
mdl[["BestDenom"]] <- colnames(Y)[apply(R2,1,which.max)]
mdl[["dl"]] <- dl

if (method=="RF"){
  mdl[["mtry"]] <- mtry
  mdl[["ntree"]] <- ntree
}


return(mdl)
}