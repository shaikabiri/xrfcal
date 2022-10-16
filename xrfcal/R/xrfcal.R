#' @title Calibrate XRF counts to element concentrations
#' @description This function takes an XRF counts matrix and reference concentrations matrix
#' and returns a model object with R2s and RMSEs for each element with various denominators for alr and best R2, RMSE, and denominator for each element. Four methods are available to use, with logratio calibration equation (LRCE)
#' using a single input for regression for each element, and multiple linear regression (MLR), cubist (Cubist) and random forest (RF), taking multiple inputs. Random forest and cubist are more computationally expensive
#' to train but generally more powerful. 
#' @examples
#' #Using the example data in the package build a calibration model
#' mdl <- xrfcal(xrf$X,xrf$Y, method = "LRCE")
#' 
#' #Plot R2s for different denominators for potassium
#' barplot(mdl$R2[which(rownames(mdl$R2)=="K"),],names.arg = colnames(mdl$R2),
#'     xlab="Denominator",ylab=paste("Crossvalidated R² for potassium"))
#' 
#' #Print out best R2 and best denominator for each element
#' print(rbind(mdl$BestR2,mdl$BestDenom))
#' @details 
#' In order to linearlize and reduce diversion from normal distribution and to
#' make methods applicable to linear geometry, also applicable to compositional geometry as well,
#' it is argued based on work of Aitchison (1982) that elemental counts and concentrations should be transformed to an additive logratio space.
#' In the current package in addition to simple linear regression between each two elements in additive logratio space (LRCE), 
#' multiple regressions with cubist, random forest and multiple linear regression are also implemented. 
#' In addition to this for each element all possible denominators for alr transformation are considered and the best
#' denominator is chosen for prediction purposes.
#' 
#' The procedure for developing the model is as followed:
#' \enumerate{
#'  \item if "LRCE" is chosen as regression method, only use the mutual elements in X and Y. Otherwise, use all.
#'  \item if there are zeroes in X, use `multRepl` from `zCompositions` to replace zeroes. If detection limits are not provided
#'  estimate them as the minimum counts detected by XRF scanner for each element.
#'  \item if `reduce` is `TRUE`, reduce dataset to a representative dataset by taking mean, first and third percentile
#'  of counts for each unique Y and aggregate them to a new dataset. 
#'  \item split data in three folds.
#'  \item for each fold, transform training and test set to alr space for each element as denominator.
#'  \item for each element as denominator train the model, make predictions on test sets, transform back to simplex, and record average RMSE and R2 over three folds.
#'  \item export results for investigation and parameters for the prediction function `pred`.
#' } 
#' 
#' @references
#' Weltje, G. J., & Tjallingii, R. (2008). Calibration of XRF core scanners for quantitative geochemical logging of sediment cores: Theory and application. Earth and Planetary Science Letters, 274(3), 423-438. https://doi.org/10.1016/j.epsl.2008.07.054
#'  
#' Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5-32. https://doi.org/10.1023/A:1010933404324
#' 
#' Quinlan, J.R. (1993). Combining Instance-Based and Model-Based Learning. ICML.
#' 
#' Aitchison, J. (1982). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-177. http://www.jstor.org/stable/2345821 
#' 
#' @param X Elemental counts matrix in simplex space. Doesn't have to be closed.
#' @param Y Reference concentrations matrix in simplex space. Doesn't have to be closed.
#' @param method Regression method. Either "LRCE", "MLR", "Cubist" or "RF".
#' @param dl Detection limits to use for zero replacement algorithm. A vector with D elements. If not provided will be estimated. 
#' @param oalr A boolean parameter. If true the model results will be in alr space.
#' @param reduce A boolean parameter. If true X and Y will be reduced to a representative dataset. 
#' @param mtry The number of variables for random forest to sample randomly at each split.
#' @param ntree Number of estimators for random forest.
#' @return A model object with attributes R2, RMSE, BestR2, BestRMSE and BestDenom. Rest of attributes are used for prediction function.
#' @export

xrfcal <-
  function(X,
           Y,
           method = "RF",
           dl = 0,
           oalr = FALSE,
           reduce = TRUE,
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
      stop("oalr must be boolean")
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
    
    if (reduce){
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
              as.matrix(Xtest) %*% model$coefficients[2:(ncol(Xtest)+1)] + model$coefficients[1]
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
  
rownames(R2) <- colnames(Y)
colnames(R2) <- colnames(Y)
rownames(RMSE) <- colnames(Y)
colnames(RMSE) <- colnames(Y)

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