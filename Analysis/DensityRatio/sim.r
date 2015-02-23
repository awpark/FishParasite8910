source('LSIF_path.r')
source('uLSIF.r')

#set.seed(2)
##
nte <- 200
ntr <- 40
xdim <- 2
xtrain <- matrix(rnorm(ntr*xdim,mean=1,sd=1/2),ntr,xdim)
xtest  <- matrix(rnorm(nte*xdim,mean=2,sd=1/4),nte,xdim)

p <- LSIF_path(xtrain, xtest)
up <- uLSIF(xtrain, xtest)

idx <- rev(order(rowMeans(p$wtrain)))[1:15]
weights <- p$wtrain[idx,]
lambdas <- p$l


plot(lambdas, weights[1,],ylim=c(0,max(weights)),xlim=c(0,max(lambdas)*1.1),type='l',lwd=2,xlab='lambda',ylab='weights',
     main='regularization paths of weights')
for (i in 2:length(idx)){
  lines(lambdas, weights[i,],col=i,lwd=2)
}
for (i in 1:length(idx))
  points(max(lambdas)*1.1,up$wtrain[idx[i]],col=i,cex=1.5,lwd=3)

