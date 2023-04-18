#' @title ModelSetup: Setup the index matrix of variables of candidate models
#' @description
#' ModelSetup is a function to generate an index matrix of variables of candidate models.
#' The sth row contains the indexes of variables used in the sth candidate model.
#'
#' @param y a vector consist of natural numbers.
#' @param X a numeric matrix.
#' @param modeltype a string chose from 'nested' or 'all'. 'nested' means
#' generating candidate models in nested from according to cor(y,x); 'all' means considering
#' all possible combinations.
#' @param intercept a bool. It means the first column of X is 1, or not.
#'
#' @return an index matrix of variables of candidate models. The sth row contains
#' the indexes of variables used in the sth candidate model.
#' @export
#'
#' @import magrittr
#' @import Rsolnp
#' @import mvtnorm
#' @importFrom stats cor
#' @importFrom utils combn
#'
#'
#' @examples
#' library(mvtnorm)
#' library(magrittr)
#' beta0=rep(0.1,5)
#' X=rmvnorm(100,mean=rep(0,5))
#' lambda=exp(X%*%beta0)
#' y=rpois(100,lambda)
#' ModelSetup(y,X,'nested',intercept=FALSE)
#'
#'
#' beta0=rep(0.1,5)
#' X=cbind(1,rmvnorm(100,mean=rep(0,4)))
#' lambda=exp(X%*%beta0)
#' y=rpois(100,lambda)
#' ModelSetup(y,X,'nested',intercept=TRUE)
#'
#'
#'

ModelSetup<-function(y,X,modeltype='nested',intercept=FALSE){
  if(modeltype=='all'){
    if(intercept){
      m=ncol(X)-1
      Index=matrix(ncol=m,nrow=0)
      for (i in 1:m) {
        index=combn(m,i)+1
        index=rbind(index,matrix(0,nrow=m-i,ncol=choose(m,i)))
        Index=rbind(Index,t(index))
      }
      Index=cbind(1,Index)

    }else{
      m=ncol(X)
      Index=matrix(ncol=m,nrow=0)
      for (i in 1:m) {
        index=combn(m,i)
        index=rbind(index,matrix(0,nrow=m-i,ncol=choose(m,i)))
        Index=rbind(Index,t(index))
      }
    }

  }else if(modeltype=='nested'){
    if(intercept){
      X2=cbind(y,X[,2:ncol(X)])
      orders=X2%>%cor%>%'['(1,2:ncol(X2))%>%abs%>%order(decreasing = T)%>%add(1)
      p=ncol(X)
      m=p-1
      Index=matrix(0,ncol=m,nrow=m)
      for (i in 1:m) {
        Index[i,1:i]=orders[1:i]
      }
      Index=cbind(matrix(rep(1,each=m),nrow=m),Index)
    }else{
      X2=cbind(y,X)
      orders=X2%>%cor%>%'['(1,2:ncol(X2))%>%abs%>%order(decreasing = T)
      m=ncol(X)
      Index=matrix(0,ncol=m,nrow=m)
      for (i in 1:m) {
        Index[i,1:i]=orders[1:i]
      }
    }
  }
  return(Index)
}
