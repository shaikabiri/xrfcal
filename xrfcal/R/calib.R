#' Calibrate XRF counts to element concentrations
#' This function takes an XRF counts matrix and reference concentrations
#' count matrix and outputs alphas and betas. 
#' @param X Elemental counts matrix
#' @param Y Reference concentrations 
#' @param method 
#' @return A matrix of alphas and betas

#' @export

calib <- function(X,Y,method = "RF"){
  if (method == "LRCE"){
    #align the two matrices 
    names <- vector(mode="character")
    
    
    for (j in colnames(X)){
      if (j %in% colnames(Y)){
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
    
    out <- averager(nzX,Y)
    nzX <- out[[1]]
    Y <- out[[2]]
    set.seed(123)
    folds <- caret::createFolds(Y[,1],k=3)
    R2 <- matrix(ncol=(ncol(X)),nrow=(ncol(X)))
    #do an LRCE for each element as denominator and save the results
    for (i in 1:ncol(nzX))
    {
      R2i <- vector(length = 30)
      Xalr <- alr(nzX,i)
      Yalr <- alr(Y,i)
      
      for (j in 1:ncol(Xalr)){
        R2k = vector(length = length(folds))
        for (k in 1:length(folds)){
          Xtest <- Xalr[folds[[k]],]
          Xtrain <- Xalr[-folds[[k]],]
          Ytest <- Yalr[folds[[k]],j]
          Ytrain <- Yalr[-folds[[k]],j]

          
          if (method=="cubist")
          {
            model <- Cubist::cubist(y=Ytrain,x=Xtrain)
            pred <- predict(model,Xtest)
          }
          
          if (method=="linear")
          {
            model <- lm(Ytrain~as.matrix(Xtrain))
            pred <- as.matrix(Xtest)%*%model$coefficients[2:31] + model$coefficients[1]
          }
          
          if (method=="RF")
          {
            model <- randomForest::randomForest(x=Xtrain,y=Ytrain,mtry=30,ntree=300)
          }
          
          R2k[k] <- caret::postResample(pred,Ytest)[2]
        }
        R2i[j] <- mean(R2k)}
      R2i <- append(R2i,0,after=i-1)
      R2[,i]<- R2i
    }
  } 
  
  if (method == "RF"){
    
    out <- averager(X,Y)
    X <- out[[1]]
    Y <- out[[2]]
    
    
    set.seed(123)
    folds <- caret::createFolds(Y[,1],k=3)
    #do an LRCE for each element as denominator and save the results
    for (i in 1:ncol(X))
    {
      set.seed(123)
      
      repeat_cv <- caret::trainControl(method='repeatedcv',number=5,repeats=3)
      train_index <- caret::createDataPartition(y=as.matrix(X[,1]),p=0.7,list=FALSE)
      x_train <- X[train_index,]
      x_test <- X[-train_index,]
      y_train <- Y[train_index,i]
      y_test <- Y[-train_index,i]
      forest <- caret::train(x=X,y=X[,5],method='rf', trt = repeat_cv,metric="Rsquared")
      
    }
  }
  
  
  return(R2)
}