RCLM <- function(Data, iteration = 10000, warm = 100)
{
  DataNum=length(Data[,1]);
  DataP=length(Data[1,]);
  total=DataNum*DataP
  DataX=NULL;
  for(i in 1:DataNum)
  {
    for(j in 1:DataP)
      DataX[(i-1)*DataP+j]=Data[i,j];
  }
  
  .C("podm0",
     as.numeric(DataX),
     as.integer(total),
     as.integer(DataNum),
     as.integer(iteration),
     as.integer(warm)
  )
  path<-matrix(scan("aa00.path"), ncol=5, byrow=TRUE)
  coef <- matrix(scan("aa00.out"), ncol=5, byrow=TRUE)
  colnames(path) <- c("beta_0","beta_1","beta_2","beta_3","sigma^2")
  colnames(coef) <- c("beta_0","beta_1","beta_2","beta_3","sigma^2")
  result <- list()
  result$path = path
  result$coef = coef
  unlink("aa00.path")
  unlink("aa00.out")
  return(result)
}
