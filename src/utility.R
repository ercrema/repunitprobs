#'@title Simulates two settlement patterns with different site-size distributions.
#'@description The simulation generates a log-normal distribution (time-period 1) and uniform distribution (time-period 2) of settlement sizes with user defined total number of residential units for each.
#'@param K1 Total number of residential units in period 1
#'@param K2 Total number of residential units in period 1
#'@param mu Mean of the distribution of settlement sizes in log-scale. 
#'@param sigma Standard deviation of the distribution of settlement sizes in log-scale. 
#'@param attempts Number of attempts to ensure outpu has used defined values of K1 and K2
#'@param rr Adjustement parameters for the uniform distribution. 

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

#'@title Simulates site-size biased sampling
#'@param x Output of the \code{sim.settlement} function.
#'@param r Sampling fraction. Default is 0.2
#'@param b Size-size sampling bias. Default is 0.
#'@details The function randomly samples sites with a site-bias factor b. When this is set to 0 there is no bias (i.e. all sites have the same chance of being selected), when b>0 larger settlements have higher chance of being selected. When b=1 the probability of a settlement being selected is directly proportional to its relative size. 

biasedsampling<-function(x,r=0.2,b=0)
{
  if (b<0|b>1) {b=0;warning("b must be between 0 and 1. running with b = 0")}
  
  index=sample(length(x$t1t2),size=c(length(x$t1t2)*r),replace=FALSE,prob=x$t1t2^b)
  
  t2=x$t1t2[index][which(x$Periods[index]=="t2")]
  t1=x$t1t2[index][which(x$Periods[index]=="t1")]
  
  pr=100*(sum(t2)-sum(t1))/sum(t1)
  return(pr)
}

#'@title Rescales a numerical vector to a user-defined range
#'@param x A vector of numeric values
#'@param M Target maximum value
#'@param m Target minimum value

rescale=function(x,M,m)
{
  return(((M-m)/(max(x)-min(x)))*(x-max(x))+M)
}


#'@title Simulating the effect of archaeological periodisation and modifiable temporal unit problem
#'@param n Number of events to be sampled. Default is 200.
#'@param years A vector of calendar years
#'@param p A vector of probabilities matching years.
#'@param nsim Number of Monte-Carlo simulations. Default is 1000.
#'@param breaks A vector breakpoint between archaeological periods.
#'@param bb Size of chronological blocks for aggregating simulated events.
#'@details. The function generates first n samples from the user provided model (parameters p and years), then assigns to each sample a membership to an archaeological phase defined by the parameter breaks. A simulation routine is then initialised. This consist of assigning randomly to each event a new time-stamp within each phase using a uniform probability distribution, and then aggregating dates by equally sized chronological blocks of size resolution. This process is repeated nsim times. 


mcunif=function(n=200,years,p,nsim=1000,breaks,resolution=100)
{
  tmp=sample(years,size=n,prob=p/sum(p),replace=T)
  tmp2 = cut(tmp,breaks=breaks,labels=letters[1:(length(breaks)-1)])
  tmp3 = rev(table(tmp2))
  bbSEQ = seq(max(breaks),min(breaks),by=-resolution)
  res = matrix(NA,nrow=length(bbSEQ)-1,ncol=nsim)
  
  for (s in 1:nsim) 
  {
    tmp4 = numeric()
    for (b in 1:(length(breaks)-1))
    {
      tmp4 = c(tmp4,runif(tmp3[b],min=breaks[b+1],max(breaks[b])))
    }
    tmp5 = cut(tmp4,breaks=bbSEQ,labels=1:(length(bbSEQ)-1))        
    res[,s]=(rev(table(tmp5)))
  }
  
  return(res)
}





