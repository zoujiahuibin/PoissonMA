library(magrittr)
library(mvtnorm)
library(Rsolnp)



#setup candidate model set  
Modelsetup<-function(y,X,modeltype=1,intercept){
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
      orders=X2%>%cor%>%.[1,2:ncol(X2)]%>%abs%>%order(decreasing = T)%>%add(1)
      p=ncol(X)
      m=p-1
      Index=matrix(0,ncol=m,nrow=m) 
      for (i in 1:m) {
        Index[i,1:i]=orders[1:i]
      }
      Index%<>%cbind(matrix(rep(1,each=m),nrow=m),.)
    }else{
      X2=cbind(y,X)
      orders=X2%>%cor%>%.[1,2:ncol(X2)]%>%abs%>%order(decreasing = T)
      m=ncol(X)
      Index=matrix(0,ncol=m,nrow=m) 
      for (i in 1:m) {
        Index[i,1:i]=orders[1:i]
      }
    }
    }
  return(Index)
}




#Poisson model averaging
PoissonMA<-function(y,X,Index=NULL, modeltype,intercept=intercept){
  #setup candidate models
  if(is.null(Index)){
    Index=Modelsetup(y,X,modeltype=modeltype,intercept)
  }
  S=nrow(Index)
  p2=ncol(Index)
  n=length(y)

  
  # solve every candidate model
  modelsolver<-function(i){
    datai=data.frame(cbind(y,X[,Index[i,]]))
    set.seed(100)
    fit <- try(
      modeli<-glm(y~.-1,family = poisson(link = "log"),data=datai)#,control=list(maxit=maxit))
    )
    if("try-error" %in% class(fit)){
      print(paste0("error in solving candidate models s=",i))
      hatbetai=rep(0,p2)
      Bic=1e7
      Aic=1e7
      
    }else if(!modeli$converged){
      print(paste0("disconvergent in solving candidate model s=",i))
      hatbetai=rep(0,p2)
      hatbetai[c(Index[i,])]=modeli$coefficients
      hatbetai[is.na(hatbetai)]=0
      Bic=BIC(modeli)
      Aic=AIC(modeli)
    }else{
      hatbetai=rep(0,p2)
      hatbetai[c(Index[i,])]=modeli$coefficients
      hatbetai[is.na(hatbetai)]=0
      Bic=BIC(modeli)
      Aic=AIC(modeli)
    }
    return(list(hatbetai,Aic,Bic))
    
  }
  
  Aic=matrix(ncol=1,nrow = S)
  Bic=Aic
  Hatbeta=matrix(nrow = S,ncol=p2)
  for (i in 1:S) {
    parlresult=modelsolver(i)
    Hatbeta[i,]=parlresult[[1]]
    Aic[i]=parlresult[[2]]
    Bic[i]=parlresult[[3]]
  }
  rm(i,parlresult)
  print("finish step1\n")
  
  
  # calculate \sum_{i=1}^n y_i x_i^T \hatbeta_{(s)}^{(y_i-1)}, remark1
  Yxbetersolver<-function(s){
    yxbeta=0
    for (i in 1:n) {
      if(y[i]!=0){
        yi=y
        yi[i]=yi[i]-1
        datai=data.frame(cbind(y=yi,X[,Index[s,]]))
        fit<-try(
          modeli<-glm(y~.-1,family = poisson(link = "log"),data=datai)
        )
        if("try-error" %in% class(fit)){
          hatbetai=rep(0,length(Index[s,]))
        }else if(!modeli$converged){
          hatbetai=rep(0,length(modeli$coefficients))
          hatbetai[is.na(hatbetai)]=0
          print("没有收敛！")
        }else{
          hatbetai=modeli$coefficients
          hatbetai[is.na(hatbetai)]=0
        }
        yxbeta[i]=y[i]*X[i,Index[s,]]%*%hatbetai
      }else{
        yxbeta[i]=0
      }
    }
    Yxbeta=sum(yxbeta)
    return(Yxbeta)
  }
  
  
  Yxbeta=matrix(nrow = S,ncol=1)
  for(i in 1:S){
    Yxbeta[i]=Yxbetersolver(i)
  }
  rm(i)
  print("finish step2\n")

  
  
  # calcualte model-averaging weight
  solvew_fun_MA<-function(start0){
    ob_fun<-function(w){
      return((sum(exp(X[,1:(p2)]%*%t(Hatbeta)%*%w))-sum(Yxbeta*w)))
    }
    A=matrix(rep(1,S),1,S,byrow = TRUE)
    ret=solnp(start0,ob_fun,eqfun=sum,eqB=c(1),LB=rep(0,S),UB=rep(1,S))
    hatw=ret$par
    return(hatw)
  }
 
  
  # calculate KL loss 
  KLloss_semi<-function(w){
    #-\frac{1}{n_1}\sum_{i=1}^{n_1} \log\frac{\hat{\mu}^y_{test,i}e^{-\hat{\mu}_i}}{y_i !}
    part1=X%*%t(Hatbeta)%*%w
    klloss_semi=-mean(y*part1-exp(part1)-log(factorial(y)))
    return(klloss_semi)
  }
  
  # calculate MSE
  MSE<-function(w){
    part1=X%*%t(Hatbeta)%*%w
    return(mean((y-exp(part1))^2))
  }
  
  
  
  # calculate the scores of all methods 
  start0=rep(1,S)/S
  hatw=solvew_fun_MA(start0)
  KLlossMA=KLloss_semi(hatw)
  MSEMA=MSE(hatw)
  

  AICw=rep(0,S)
  AICselection=which.min(Aic)
  AICw[AICselection]=1
  KLlossAIC=KLloss_semi(AICw)
  MSEAIC=MSE(AICw)

  BICw=rep(0,S)
  BICselection=which.min(Bic)
  BICw[BICselection]=1
  KLlossBIC=KLloss_semi(BICw)
  MSEBIC=MSE(BICw)
  
  
  Aic2=t(matrix(rep(-Aic,S),ncol=S,byrow=TRUE)+c(Aic))
  SAICw=1/colSums(apply(Aic2/2, 2,exp))
  KLlossSAIC=KLloss_semi(SAICw)
  MSESAIC=MSE(SAICw)
  
  Bic2=t(matrix(rep(-Bic,S),ncol=S,byrow=TRUE)+c(Bic))
  SBICw=1/colSums(apply(Bic2/2, 2,exp)) 
  KLlossSBIC=KLloss_semi(SBICw)
  MSESBIC=MSE(SBICw)
  
  Fullw=rep(0,S)
  Fullw[S]=1
  KLlossFull=KLloss_semi(Fullw)
  MSEFull=MSE(Fullw)
  
  
  result=list(Index=Index,
       weight.MA=hatw, weight.SAIC=SAICw, weight.SBIC=SBICw, OptimalModel.AIC=AICselection, OptimalModel.BIC=BICselection,
       MSE.MA=MSEMA, MSE.AIC=MSEAIC, MSE.BIC=MSEBIC, MSE.SAIC=MSESAIC, MSE.SBIC=MSESBIC, MSE.FULL=MSEFull,
       KLloss.MA=KLlossMA,KLloss.AIC=KLlossAIC,KLloss.BIC=KLlossSBIC,KLloss.SAIC=KLlossSAIC,KLloss.SBIC=KLlossSBIC,KLloss.Full=KLlossFull)
  return(result)
}
 



beta02=c(0.1,0,0,-0.6,0.1,0.07,-0.07,0.3)
rho=0.8
n=700#用于评估损失的样本量

beta0=beta02[-length(beta02)]#最后一个变量不在候选模型中
m=length(beta02)-2#这是可变变量的个数
p=length(beta02)-1#参数的维数

XCov=matrix(rep(rho,p*p),p,p)+diag(1-rho,p)#X的协方差阵

set.seed(100)
X=rmvnorm(n,mean=rep(0,p),sigma=XCov)
X=cbind(rep(1,n),X)
lambda=exp(X%*%beta02) 
y=rpois(n, lambda)










Index=Modelsetup(y,X,modeltype='nested',intercept = T)

results=PoissonMA(y,X,Index=Index, modeltype='nested',intercept=T)
results

