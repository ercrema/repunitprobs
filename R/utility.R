### Utility Functions ####


#### Computes 95% Confidence Interval of the rate of change between two periods
marginSites<-function(t1,t2)
{
  p=t2/(t1+t2) #observed proportion
  n=t1+t2 #sample size
  se=sqrt(p*(1-p)/n) #standard error
  ci=1.96*se #95% confidence interval
  rate = p/(1-p) #observed rate of change
  lo = (p-ci)/(1-(p-ci))
  hi = (p+ci)/(1-(p+ci))
  return(list(obsRate=rate,ci=c(lo,hi)))
}

#### Computes 95% Confidence Interval of the rate of change between two periods via Boootstrap
marginHouses<-function(t1,t2,nsim=1000)
{
  rate = sum(t2)/sum(t1) #observed rate of change
  t1t2=c(t1,t2) #combine two vectors
  Periods<-c(rep("t1",length(t1)),rep("t2",length(t2))) #define period of each index
  bootStrapRate=numeric()
  
  for (s in 1:nsim)
  {
    index=sample(1:length(Periods),replace=TRUE)
    bootStrapRate[s]=sum(t1t2[index][which(Periods[index]=="t2")])/sum(t1t2[index][which(Periods[index]=="t1")])
  }
  lo=quantile(bootStrapRate,0.025)
  hi=quantile(bootStrapRate,0.975)
  return(list(obsRate=rate,ci=c(lo,hi)))
}


#### Simulates sampling bias (t = taphonomic bias; s = size bias)

biasedsampling<-function(t1t2,Periods,recoveryRate=0.2,bias,biasType=c("t","s","cf"))
{
  
  if (biasType=="t")
  {
    prob=rep(1,length(t1t2))
    prob[which(Periods=="t1")]=bias
    index=sample(length(t1t2),size=c(length(t1t2)*recoveryRate),replace=FALSE,prob=prob)
  }
  if (biasType=="s")
  {
    index=sample(length(t1t2),size=c(length(t1t2)*recoveryRate),replace=FALSE,prob=t1t2^bias)
  }
  
  t2=t1t2[index][which(Periods[index]=="t2")]
  t1=t1t2[index][which(Periods[index]=="t1")]
  PeriodsOut=c(rep("t1",length(t1)),rep("t2",length(t2)))
  return(list(t1=t1,t2=t2,t1t2=c(t1,t2),Periods=PeriodsOut))
}

#### Simulates group fission fusion 

fissionfusion<-function(nSettlements,t1Mean,t1Sd,fissionRate)
{
  t1 = round(rnorm(nSettlements,mean=t1Mean,sd=t1Sd))
  t2 = numeric()
  
  for (x in 1:length(t1))
  {
    tmp=t1[x]
    diff=tmp%%fissionRate
    offSpringSize=(tmp-diff)/fissionRate
    offSprings=rep(offSpringSize,fissionRate)
    offSprings[1]=offSprings[1]+diff
    t2=c(t2,offSprings)
  }
  PeriodsOut=c(rep("t1",length(t1)),rep("t2",length(t2)))
  return(list(t1=t1,t2=t2,t1t2=c(t1,t2),Periods=PeriodsOut))
}

#### Simulates Primate and Convex Patterns:

primateConvex<-function(N=300,mean1=3,sd1=1,attempts=20000,rr=5)
{
  t1=0
  while(all(t1==0))
  {
    t1=round(rlnorm(N,mean1,sd1))
  }
  t2=0
  x=0
  while((sum(t2)!=sum(t1))&(x<=attempts))
  {
    t2=round(runif(N,min=sum(t1)/N-rr,max=sum(t1)/N+rr))
    x=x+1
  }
  PeriodsOut=c(rep("t1",length(t1)),rep("t2",length(t2)))
  return(list(t1=t1,t2=t2,t1t2=c(t1,t2),Periods=PeriodsOut))
}



