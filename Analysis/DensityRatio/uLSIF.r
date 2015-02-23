## uLSIF.r
library(MASS)
library(kernlab)

uLSIF <- function(xtrain, xtest, xreference=xtrain){
  ## Unconstrained least-squares importance fitting with leave-one-out cross validation
  ##
  ## Usage:
  ##      list(wtrain, wreference) = uLSIF(xtrain, xtest, xreference)
  ##
  ## Input:
  ##    xtrain : ntrain by d training input matrix (iid from density ptrain)
  ##    xtest  : ntest by d test input matrix (iid from density ptest)
  ## xreference: (OPTIONAL) d by nrefrence reference input matrix
  ##
  ## Output: list
  ##    wtrain    : b by length(lambda) estimates of normalized density ratio w=ptest/ptrain at xtrain
  ##    wreference: b by length(lambda) estimates of normalized density ratio w=ptest/ptrain at xreference
  ##
  ## (c) Takafumi Kanamori, Department of Computer Science and Mathematical Informatics, Nagoya University, Japan.
  ##     kanamori@is.nagoya-u.ac.jp,   http://www.math.cm.is.nagoya-u.ac.jp/~kanamori/software/LSIF/
  xdim <- nrow(xtrain); b <- min(100,nrow(xtest)); center <- matrix(xtest[sample(1:nte,b),],b,ncol(xtest))
  sigmai <- as.array(quantile((dist(center)))); lambdai=10^(seq(-3,1,l=9))
  LOOCV <- function(dat.train,dat.test,center,sigma=median(dist(center)),lambda=1){
    ntr <- nrow(dat.train); nte <- nrow(dat.test);
    Xtr <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(dat.train), as.matrix(center))
    Xte <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(dat.test), as.matrix(center))
    H <- t(Xtr) %*% Xtr/ntr; h <- colMeans(Xte); cvidx <- sample(1:nte,ntr,replace=T)
    tmpV <- ginv(H+(1-1/ntr)*lambda*diag(nrow(H))) %*% cbind(h,t(Xtr),t(Xte[cvidx,]))
    beta <- tmpV[,1]; V <- tmpV[,2:(ntr+1)]; W <- tmpV[,-(1:(ntr+1))]; dom <- ntr-colSums(V * t(Xtr))
    B <- beta + t(t(V) * as.vector((Xtr %*% beta)/dom)); B <- B*(ntr-1)/ntr
    B <- nte/(nte-1) * B - (ntr-1)/(ntr*(nte-1)) * (W+t(t(V)*(colSums(t(Xtr)*W)/dom)))
    A <- pmax(B,0); wtr <- colSums(A * t(Xtr)); wte <- colSums(A * t(Xte[cvidx,])); ans <- mean(wtr^2)/2-mean(wte)
    ans
  }
  LOOCV.lambda <- function(dat.train,dat.test,center,sigma=median(dist(center)),lambdai=10^(seq(-3,1,l=9))){
    res <- c()
    for (lambda in lambdai)
      res <- c(res, LOOCV(dat.train,dat.test,center,sigma=sigma,lambda=lambda))
    res
  }
  cv <- Inf
  for (sigma in sigmai){
    tl <- LOOCV.lambda(xtrain,xtest,center,sigma=sigma,lambdai=lambdai)
    if (cv > min(tl)){
      cv <- min(tl)
      best.para <- c(lambdai[which.min(tl)], sigma)
    }
  }
  alpha <- alpha.estimation(xtrain,xtest,center,sigma=best.para[2],lambda=best.para[1])
  trK <- kernelMatrix(rbfdot(sigma=1/(2*best.para[2]^2)), as.matrix(center), as.matrix(xtrain))
  refK <- kernelMatrix(rbfdot(sigma=1/(2*best.para[2]^2)), as.matrix(center), as.matrix(xreference))
  wtrain <- array(t(trK) %*% alpha)
  wreference <- array(t(refK) %*% alpha)
  list(wtrain=wtrain/sum(wtrain), wreference=wreference/sum(wreference))
}


## uLSIF estimator
alpha.estimation <- function(dat.train,dat.test,center,sigma=median(dist(center)),lambda=1){
  ntr <- nrow(dat.train)
  Xtr <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(dat.train), as.matrix(center))
  Xte <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(dat.test), as.matrix(center))
  H <- t(Xtr) %*% Xtr/ntr
  h <- colMeans(Xte)
  pmax(ginv(H+lambda*diag(nrow(H))) %*% h, 0)
}
