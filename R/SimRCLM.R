SimRCLM <- function(I=100,J=10,beta,sigma){
  ###### generating subject property########
  Z<-NULL
  age<-runif(I,0.2,0.8)
  sex<-rbinom(I,1,0.5)
  education<-rmultinom(I,1,prob=c(0.1,0.45,0.4,0.05))
  education<-t(education)
  education<-education[,2:4]
  Z<-cbind(age, sex, education)
  
  lambda<-cbind(2*rnorm(I), 0.5*rnorm(I), rnorm(I), rnorm(I), rnorm(I))
  
  ###### generating product  property########
  W<-NULL
  age<-runif(J,0,1) ## the age of the product
  producer<-rmultinom(J,1,prob=c(0.2, 0.3, 0.4,0.1))
  producer<-t(producer)
  producer<-producer[,1:3]
  
  categ<-rmultinom(J,1,prob=c(0.25,0.35,0.15,0.25))
  categ<-t(categ)
  categ<-categ[,1:3]
  W<-cbind(age, producer, categ)
  
  gamma<-cbind(rnorm(J),0.5*rnorm(J),0.5*rnorm(J),0.5*rnorm(J),rnorm(J),rnorm(J),rnorm(J))
  ###### generating data  ##############3
  
  x0<-rep(1,I*J)
  x1<-rnorm(I*J,0,1)
  x2<-rnorm(I*J,0,1)
  x3<-rnorm(I*J,0,1)
  X<-cbind(x0,x1,x2,x3)
  
  eps<-sigma*rnorm(I*J,0,1)
  
  
  ########## coefficients setting #######
  
  SS1<-NULL
  for(i in 1:ncol(X)){
    SS1<-cbind(SS1,X[,i]*beta[i])
  }
  SS1<-apply(SS1,1,sum)
  
  SS2<-NULL
  for(i in 1:I){
    s2<-sum(Z[i,]*lambda[i,])
    SS2<-c(SS2,rep(s2,J))
  }
  
  SS3<-NULL
  for(j in 1:J){
    s3<-sum(W[j,]*gamma[j,])
    SS3<-c(SS3,s3)
  } 
  SS3<-rep(SS3,I)
  
  y<-SS1+SS2+SS3+eps
  
  ZZ<-NULL
  for(i in 1:I){
    z1<-rep(Z[i,],J)
    ZZ<-rbind(ZZ, t(matrix(z1, ncol=J)))
  }
  
  WW<-NULL
  for(j in 1:ncol(W))
    WW<-cbind(WW, rep(W[,j],I))
  
  D<-round(cbind(y,X,ZZ,WW),6)
  
  return(D)
}
