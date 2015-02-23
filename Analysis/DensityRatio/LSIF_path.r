## LSIF_path.r
library(MASS)
library(kernlab)

LSIF_path <- function(xtrain, xtest, xreference=xtrain){
  ## regularization path of least-squares importance fitting
  ##
  ## Usage:
  ##  list(lambda, wtrain, wreference) = LSIF_path(xtrain, xtest, xreference)
  ##
  ## Input:
  ##  xtrain    : ntrain by d training input matrix (iid from density ptrain)
  ##  xtest     : ntest by d test input matrix
  ##  xreference: (OPTIONAL) nreference by d reference input matrix
  ##
  ## Output: list
  ##  lambda : change points of regularization parameter
  ##  wtrain    : b by length(lambda) estimates of normalized density ratio w=ptest/ptrain at xtrain at each change point
  ##  wreference: b by length(lambda) estimates of normalized density ratio w=ptest/ptrain at xreferece at each change point
  ##
  ## (c) Takafumi Kanamori, Department of Computer Science and Mathematical Informatics, Nagoya University, Japan.
  ##     kanamori@is.nagoya-u.ac.jp,   http://www.math.cm.is.nagoya-u.ac.jp/~kanamori/software/LSIF/
  digit=10; xdim <- nrow(xtrain); b <- min(100,nrow(xtest))
  center <- matrix(xtest[sample(1:nte,b),],b,ncol(xtest))
  sigma <- median(dist(center))
  trK <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(center), as.matrix(xtrain))
  teK <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(center), as.matrix(xtest))
  refK <- kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), as.matrix(center), as.matrix(xreference))
  A <- trK %*% t(trK)/ncol(trK); addD <- min(c(eigen(A)$val,0)); A <- A + 3*abs(addD)*diag(1,b)
  d <- rowMeans(teK)
  H.mat <- function(A,J){
    if (length(J)==0){B <- A}else
    {
      J <- sort(J)
      E <- matrix(0,length(J),nrow(A))
      E[cbind(1:length(J),J)] <- 1
      B <- rbind(cbind(A,-t(E)), cbind(-E, matrix(0,length(J),length(J))))
    }
    B
  }
  lambda.list <- c(); g.list <- c(); h.list <- c(); J <- 1:b; lambda <- max(d)
  J <- J[J!=which.max(d)]
  while(lambda>0){
    Hmat <- H.mat(A,J) + 10^(-digit)*diag(1,nrow(H.mat(A,J)))
    gh <- solve(Hmat,cbind(c(d,array(0,length(J))),c(array(-1,b),array(0,length(J)))))
    gh <- round(gh,digit); g.list <- cbind(g.list, gh[1:b,1]); h.list <- cbind(h.list, gh[1:b,2])
    if(all(gh[,2]<=0)){
      lambda.list <- c(lambda.list,0); break
    }
    idx <- which(gh[,2]>0)
    tmplambda <- -gh[idx,1]/gh[idx,2]; tmplambda[tmplambda>=lambda] <- -Inf 
    tmpk <- idx[which.max(tmplambda)]; lambda <- max(c(tmplambda,0))
    if((length(lambda.list)>=1) && (lambda>rev(lambda.list)[1])){print(sigma);stop('! lambda numerical error !')}
    lambda.list <- c(lambda.list, lambda)
    if (tmpk<=b){
      J <- sort(c(J,tmpk))
    }else{
      J <- sort(J[J!=J[tmpk-b]])
    }
  }
  lambda.list <- rev(lambda.list)
  g.list <- as.matrix(g.list[,ncol(g.list):1]); h.list <- as.matrix(h.list[,ncol(h.list):1])
  alpha.list <- g.list+t(t(h.list)*lambda.list[1:(length(lambda.list))])
  wtrain <- t(trK) %*% alpha.list; wreference <- t(refK) %*% alpha.list
  ##  list(lambda=lambda.list, sol=round(alpha.list,digit), w=t(t(w)/colSums(w)))
  list(lambda=lambda.list,wtrain=t(t(wtrain)/colSums(wtrain)),wreference=t(t(wreference)/colSums(wreference)))
}


