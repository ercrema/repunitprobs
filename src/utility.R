sim.settlement = function(K1=300,K2=300,mu=3,sigma=1,attempts=20000,rr=5)
{
  t1=round(rlnorm(K1,mu,sigma))
  t2=0
  x=0
  while((sum(t2)!=sum(t1))&(x<=attempts))
  {
    t2=round(runif(K2,min=sum(t1)/K1-rr,max=sum(t1)/K1+rr))
    x=x+1
  }
  PeriodsOut=c(rep("t1",length(t1)),rep("t2",length(t2)))
  return(list(t1=t1,t2=t2,t1t2=c(t1,t2),Periods=PeriodsOut))
}


biasedsampling<-function(x,r=0.2,b=0)
{
  
  index=sample(length(x$t1t2),size=c(length(x$t1t2)*r),replace=FALSE,prob=x$t1t2^b)
  
  t2=x$t1t2[index][which(x$Periods[index]=="t2")]
  t1=x$t1t2[index][which(x$Periods[index]=="t1")]
  
  pr=100*(sum(t2)-sum(t1))/sum(t1)
  return(pr)
}