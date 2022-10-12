#' Calibrate XRF counts to element concentrations
#' This function takes an XRF counts matrix and reference concentrations 
#' and returns a model object with R2s, RMSEs and  a predictor function for new data.
#' @example
#' #Using the example data in the package build a calibration model
#' mdl <- calib(xrf[[1]],xrf[[2]], method = "MLR", denom = "Auto")
#' 
#' 
#' @param X Elemental counts matrix
#' @param Y Reference concentrations matrix
#' @param method Regression model. Either "LRCE", "MLR", "Cubist" or "RF". RF and Cubist generally take longer time to train
#' @param denom Denominator for the logratio transformation. Can take values "Auto", "All" or abbreviation of an element that is in references column names such as "Ca".
#' "Auto" chooses the best denominator for each element, "All" outputs the R2 and RMSEs for all the denominators and the predictor function will take extra argument 
#' for a denominator in this case. If the abbreviation for a specific denominator is used only the results for that specific denominator will be available in the model object output.
#' @param dl Detection limits to use for zero replacement algorithm. A vector with D elements. If not provided will be estimated. 
#' @param oalr A boolean parameter. If true the model results and predictor function will be in alr space.
#' @return A model object that includes R2s, RMSEs and a predictor function, pred.
#' @export

calib <-
  function(X,
           Y,
           method = "RF",
           denom = "Auto",
           dl = 0,
           oalr = FALSE) {
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
    
    
    out <- averager(nzX, Y)
    nzX <- out[[1]]
    Y <- out[[2]]
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
            model <- lm(Ytrain[, j] ~ as.matrix(Xtrain[, j]))
            pred <-
              as.matrix(Xtest[, j]) %*% model$coefficients[2] + model$coefficients[1]
          }
          
          
          if (method == "Cubist")
          {
            model <- Cubist::cubist(y = Ytrain[,j], x = Xtrain)
            pred <- predict(model, Xtest)
          }
          
          if (method == "MLR")
          {
            model <- lm(Ytrain[, j] ~ as.matrix(Xtrain))
            pred <-
              as.matrix(Xtest) %*% model$coefficients[2:ncol(nzX)] + model$coefficients[1]
          }
          
          if (method == "RF")
          {
            model <-
              randomForest::randomForest(
                x = Xtrain,
                y = Ytrain[,j],
                mtry = 30,
                ntree = 300
              )
            pred <- predict(model, Xtest)
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
mdl[["X"]] <- X
mdl[["Y"]] <- Y
mdl[["R2"]] <- R2
mdl[["RMSE"]] <- RMSE
mdl[["Method"]] <- method


return(mdl)
}