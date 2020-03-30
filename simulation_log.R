source("../src/utility.R")

# Simulation 1 (Periodisation) ####

## Set theorethical population
age = seq(300,800,1)
asymptote = 0.9
sl = .04
midpoint = 550
p = 0.1+asymptote / (1 + exp((midpoint - age) * sl))
d=data.frame(yr=rev(age),p=p)

## Sample sites
n=1000
set.seed(123)
ss = round(sample(d$yr,size=1000,prob=d$p/sum(d$p),replace=TRUE))

## Define time-block resolution and block mid points
resolution=50 #define time-block resolution
midPoints = seq(775,325,-50)

# Extract min-max from theorethical model
m=min((d$p/sum(d$p)*n*resolution))
M=max((d$p/sum(d$p)*n*resolution))

#Setup different periodisations:
LL = list(c(800,550,300),
          c(800,700,300),
          c(800,600,500,300),
          c(800,400,300),
          c(800,700,400,300),
          c(800,700,600,500,400,300))


simRes1 = vector("list",length=length(LL))

for (i in 1:length(simRes1))
{
  breaks=LL[[i]]
  simRes1[[i]]=mcsim(x=ss,nsim=1000,breaks=breaks,resolution=50)
}

save(simRes1,breaks,midPoints,d,n,LL,file="../R_Images/simRes1.RData")




# Simulation 2 (Duration) ####

## Parameters Setup 
edgec <- 500 # Begin 500 years earlier and end 500 years later to avoid edge effects
years <- c(-1750:-750) #Simulation Interval
simyears <- c((years[1]-edgec) : (years[length(years)]+edgec))
nmines <- 100 #Number of mines in a given yers
minepreallocation <- 10000

mDurationMax <- 200 #Maximum Duration of Sites
mDurationMin <- 10 #Minimum Duration of Sites
df <- data.frame(x=c(years[1],years[length(years)]), y=c(mDurationMax,mDurationMin))
mod <- lm(y~x, data=df) #

## Simulate
usemodel <- TRUE
nsim <- 100
simmat <- matrix(nrow=length(simyears), ncol=nsim)
set.seed(100)
pb <- txtProgressBar(min=1, max=nsim, style=3)


for (b in 1:nsim){
  setTxtProgressBar(pb, b)
  minedf <- data.frame(StartBCE=c(rep(simyears[1],nmines),rep(NA,(minepreallocation-nmines))), EndBCE=NA, Duration=NA, Mu=NA, Active=FALSE)
  minedf$Duration[!is.na(minedf$StartBCE)] <- rnbinom(nmines, mu=200, size=1)
  minedf$EndBCE <- minedf$StartBCE + minedf$Duration
  for (a in simyears){
    check1 <- minedf$StartBCE <= a & !is.na(minedf$StartBCE)
    check2 <- minedf$EndBCE >= a & !is.na(minedf$EndBCE)
    minedf[,"Active"] <- FALSE
    minedf[check1 & check2,"Active"] <- TRUE
    checksum1 <- 100-sum(minedf$Active)
    if (checksum1 > 0){
      myrows <- which(is.na(minedf$StartBCE))[1:checksum1]
      minedf$StartBCE[myrows] <- a
      if (usemodel){
        mu <- as.numeric(coefficients(mod)[1]) + (as.numeric(coefficients(mod)[2]) * a)
        if (mu > mDurationMax){
          mu <- mDurationMax
        } else if (mu < mDurationMin){
          mu <- mDurationMin
        }
      } else {
        mu <- 200
      }
      mu <- round(mu,0)
      minedf$Mu[myrows] <- mu
      minedf$Duration[myrows] <- rnbinom(checksum1, mu=mu, size=1)
      minedf$EndBCE[myrows] <- minedf$StartBCE[myrows] + minedf$Duration[myrows]
    }
  }
  minedf <- minedf[!is.na(minedf$StartBCE),]
  nrow(minedf)
  minedf$MidYear <- minedf$StartBCE+round(minedf$Duration/2,0)
  tmp <- density(minedf$MidYear, n=length(simyears), from=simyears[1],to=simyears[length(simyears)])$y
  simmat[,b] <- tmp
}
close(pb)
mediankd <- apply(simmat,1,median)
lastrun <- minedf

## Use the last simulation run and select some example mines in early, middle and late phases
preferreddates <- c(-1700,-1300,-900)
nsitesperphase <- 5
for (a in 1:length(preferreddates)){
  myrows <- row.names(lastrun[lastrun$StartBCE >= preferreddates[a],])
  if (a==1){
    examplemines <- lastrun[myrows[1:nsitesperphase],]
    examplemines$Set <- a
  } else {
    tmp <- lastrun[myrows[1:nsitesperphase],]
    tmp$Set <- a
    examplemines <- rbind(examplemines,tmp)
  }  
}

## Back-calibrate 5 random years within each duration and back calibrate to create 5 hypothetical radiocarbon dates (assigning a plausible measurement error of 3 years). Then run an OxCal span model per site on the results.
exampleerror <- 35
ndates <- 5
spansdf <- vector(mode="list", length=nrow(examplemines))
names(spansdf) <- row.names(examplemines)
pb <- txtProgressBar(min=1, max=nrow(examplemines), style=3)

for (a in 1:nrow(examplemines)){
  setTxtProgressBar(pb, a)
  datedyears <- round(runif(n=ndates,  min=BCADtoBP(examplemines[a,"EndBCE"]), max=BCADtoBP(examplemines[a,"StartBCE"])),0)
  errors <- rep(exampleerror,length(datedyears))
  Samples14C <- uncalibrate(datedyears, CRAerrors=30)[,c("rCRA","rError")]
  fn <- tempfile()
  oxcalSpanIntScript(ids=as.character(1:ndates), ages=Samples14C$rCRA, errors=Samples14C$rError, fn=fn, span=TRUE)
  myoxcal2 <- readChar(fn, file.info(fn)$size)
  myoxcalresf2 <- executeOxcalScript(myoxcal2)
  myoxcalres2 <- parseFullOxcalOutput(readOxcalOutput(myoxcalresf2))
  span2 <- extractSpan(myoxcalres2)
  spansdf[[a]] <- span2
}
close(pb)

## Pool and normalise the span probabilties for Early, Middle and Late phases.
maxduration <- max(unlist(lapply(spansdf,max)))
allspans <- data.frame(Years=seq(2.5,maxduration,5)) 
## Add new columns
for (d in 1: length(spansdf)){
  tmp <- merge(allspans[,"Years",drop=FALSE],spansdf[[d]], all.x=TRUE, by.x="Years", by.y="duration")
  tmp[is.na(tmp$prob),] <- 0
  allspans[,names(spansdf)[d]] <- tmp$prob
}

## Pool all
pooled <- allspans[,-1]
pooled <- cbind(allspans[,1], (rowSums(pooled)/ncol(pooled)))
## Pool Early
pooledEarly <- allspans[,-1]
pooledEarly <- pooledEarly[,1:5]
pooledEarly <- cbind(allspans[,1], (rowSums(pooledEarly)/ncol(pooledEarly)))
## Pool Middle
pooledMiddle <- allspans[,-1]
pooledMiddle <- pooledMiddle[,6:10]
pooledMiddle <- cbind(allspans[,1], (rowSums(pooledMiddle)/ncol(pooledMiddle)))
## Pool Late
pooledLate <- allspans[,-1]
pooledLate <- pooledLate[,11:15]
pooledLate <- cbind(allspans[,1], (rowSums(pooledLate)/ncol(pooledLate)))


save(simyears,mediankd,simmat,pooled,pooledEarly,pooledMiddle,pooledLate,file="../R_Images/simRes2.RData")





# Simulation 3 (Nucleation/Dispersion) ####
# Set Random Seed
set.seed(133)
# Simulate Settlements
simdata3 = sim.settlement(K1=1000,K2=1000)
# Simulate Sampling
nsim = 100
b = c(0,0.3,0.7)
r = c(0.1,0.3,0.7)
simRes3=expand.grid(r=r,b=b,nsim=1:nsim)
simRes3$pr=NA

for (i in 1:nrow(simRes3))
{
  simRes3$pr[i]=biasedsampling(simdata3,r=simRes3$r[i],b=simRes3$b[i])
}

save(b,r,simdata3,simRes3,file="../R_Images/simRes3.RData")
