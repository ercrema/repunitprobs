#Figure Log Files
source("../src/utility.R")
library(trapezoid)
library(reldist)
library(ggplot2)
library(ggpubr)
library(oxcAAR)
library(reldist)
library(rcarbon)

## Figure 1 ####
pdf(file = "./figure1.pdf",width = 12,height = 4)
layout(matrix(c(1,2,3,4,5,5),nrow=2,ncol=3))


### Panel a (Periodisation Example):
plot(runif(1))
plot(runif(1))



### Panel b (Duration Example):
#### Load data:
maukE <- read.csv("../data/mauk.csv", header=TRUE, stringsAsFactors=FALSE, encoding="UTF-8",na.strings=c("NA",""),strip.white=TRUE)
maukE <- maukE[maukE$SiteName == "Mauk E", ]
ids <- paste(maukE$LabID, " (", maukE$CRA, "±", maukE$Error,")", sep="")

#### Calibrate 
maukEcal <- calibrate(maukE$CRA, maukE$Error, verbose=FALSE, ids=ids)

#### Bayesian Modelling
maukEpref <- c("VERA-4448","VERA-4447","VERA-4452","VERA-4453", "VERA-4450", "VERA-4900", "VERA-4449")
maukEpref <- maukE[maukE$LabID %in% maukEpref,]
quickSetupOxcal()
fn1 <- "../oxcalscripts/maukE.oxcal"
oxcalSpanIntScript(ids=maukEpref$LabID, ages=maukEpref$CRA, errors=maukEpref$Error, fn=fn1, span=TRUE, interval=TRUE)
myoxcal1 <- readChar(fn1, file.info(fn1)$size)
myoxcalresf <- executeOxcalScript(myoxcal1)
myoxcalres <- parseFullOxcalOutput(readOxcalOutput(myoxcalresf))
spanE <- extractSpan(myoxcalres)
intervalE <- extractInterval(myoxcalres)

#### Plot
par(mar=c(4, 2, 0.2, 0.2)) #c(bottom, left, top, right)
multiplot(maukEcal, decreasing=FALSE, label=FALSE, cex.id=0.5, gapFactor=0.1, calendar="BCAD")
legend(x=-1250,y=0.03, legend="Dendro", col="red", lwd=3, bty="n", cex=1)
xmin <- par('usr')[1]
ymin <- par('usr')[3]
xoffset <- par('usr')[2]-xmin
yoffset <- par('usr')[4]-ymin
text(0.5*(xoffset)+xmin,0.9*(yoffset)+ymin, ids[1], cex=0.5)
text(0.55*(xoffset)+xmin,0.82*(yoffset)+ymin, ids[2], cex=0.5)
text(0.57*(xoffset)+xmin,0.74*(yoffset)+ymin, ids[3], cex=0.5)
text(0.57*(xoffset)+xmin,0.66*(yoffset)+ymin, ids[4], cex=0.5)
text(0.33*(xoffset)+xmin,0.56*(yoffset)+ymin, ids[5], cex=0.5)
text(0.33*(xoffset)+xmin,0.47*(yoffset)+ymin, ids[6], cex=0.5)
text(0.33*(xoffset)+xmin,0.39*(yoffset)+ymin, ids[7], cex=0.5)
text(0.33*(xoffset)+xmin,0.31*(yoffset)+ymin, ids[8], cex=0.5)
text(0.35*(xoffset)+xmin,0.22*(yoffset)+ymin, ids[9], cex=0.5)
text(0.35*(xoffset)+xmin,0.13*(yoffset)+ymin, ids[10], cex=0.5)
text(0.35*(xoffset)+xmin,0.06*(yoffset)+ymin, ids[11], cex=0.5)
polygon(c(-715,-715,-707,-707,-715), c(0,0.56*(yoffset)+ymin,0.56*(yoffset)+ymin,0,0), col=rgb(255,0,0,200,maxColorValue=255), border=NA)
mtext("C", side=1, font=2, cex=1.2, line=2.2, adj=-0.05)

par(mar=c(4, 2, 0.2, 0.2)) #c(bottom, left, top, right)
plot(spanE$duration,spanE$prob,type='l', col="white", xlim=c(0,750), xlab="Years")
polygon(c(spanE$duration,rev(spanE$duration),spanE$duration[1]), c(spanE$prob,rep(0,length(spanE$prob)),spanE$prob[1]), col="grey75", border="grey75")
lines(intervalE$duration,intervalE$prob, col="grey25", lty="dashed")
legend("topright", legend=c("span","interval"), col=c("grey75","grey25"), lwd=c(3,1), lty=c("solid","dashed"), bty="n", cex=1)
mtext("D", side=1, font=2, cex=1.2, line=2.2, adj=-0.05)

### Panel c (Nucleation/Dispersion Example):
#### Load data:
chiba=read.csv("../data/chibaJomon_fromCrema2013.csv")
chiba = subset(chiba,!is.na(chiba$SimplifiedPhase))
chiba$SimplifiedPhase=factor(chiba$SimplifiedPhase,levels=c("Otamadai-Katsuzaka","Kasori E (early)","Kasori E (late)","Shomyoji","Horinouchi","Kasori B","Soya","Angyo 1 & 2"))
x = table(chiba$SiteID,chiba$SimplifiedPhase)
par(mar=c(10,5,2,2))
plot(0,0,xlab="",ylab="",xlim=c(0,9),axes=FALSE,type='n',ylim=c(0,max(apply(x,2,sum))+20))
# gg = numeric(length=8)
for (i in 1:8)
{
  tmp=x[,i]
  tmp=tmp[which(tmp>0)]
#  gg[i]=gini(tmp)
  tmp=c(0,cumsum(tmp))
  for (k in 1:c(length(tmp)-1))
  {
    col="lightgrey"
    if (names(tmp[k+1])=='10529'){col='orange'} #Ariyoshi-Kita Site
    if (names(tmp[k+1])=='9825'){col='royalblue'} #Miyauchiidosaku Site
    rect(xleft=i-0.4,xright=i+0.4,ybottom=tmp[k],ytop=tmp[k+1],col=col,border='grey30')
  }
  text(x=i,y=tmp[k+1]+10,labels=paste0("n=",length(tmp)-1))
}
axis(2)
#axis(3,at=1:8,labels=round(gg,2))
mtext(side=2,line=3,"Number of Pithouses",cex=1)
#mtext(side=3,line=3,"Gini Coefficient (settlement size)",cex=1)
mtext(side=1,line=6.5,"Ceramic Phases",cex=1)
axis(1,at=1:8,labels=levels(chiba$SimplifiedPhase),las=2)
legend(6,300,legend=c("Ariyoshi-Kita","Aioi","Other Sites"),fill=c("orange","royalblue","grey"),cex=1)
mtext("E", side=1, font=2, cex=1.2, line=4.5, adj=-0.05)


dev.off()








## Figure 2 ####
pdf(file = "./figure2.pdf",width = 10,height = 3)

par(mfrow=c(1,3))
cx=1.2

# panel a
plot(0,0,xlim=c(800,300),ylim=c(0,10),type='n',xlab="BC",ylab="",axes=F,main="Phase Assignment Uncertainty",cex.main=cx)
axis(1,at=seq(800,300,-100))
rect(xleft = 800,xright=600,ybottom=7.5,ytop=9.5,col='darkgrey',border=NA) 
text(x=700,y=10,labels="Phase A",cex=cx)
rect(xleft = 660,xright=420,ybottom=4,ytop=6,border=NA,col='darkgrey') 
text(x=540,y=6.5,labels="Phase B",cex=cx)
rect(xleft = 420,xright=310,ybottom=0.5,ytop=2.5,border=NA,col='darkgrey') 
text(x=365,y=3,labels="Phase C",cex=cx)
legend("bottomleft",legend=c(expression(paste(P[x],"(Phase A)=0.05)")),expression(paste(P[x],"(Phase B)=0.8)")),expression(paste(P[x],"(Phase C)=0.15)"))),bty='n',cex=cx)

# panel b
rr = seq(660,420,-1)
plot(0,0,xlim=c(800,300),ylim=c(0,8),type='n',xlab="BC",ylab="",axes=F,main="Within-phase uncertainty",cex.main=cx)
scl = 5/max(dtrapezoid(-rr,-660,-550,-440,-420))
axis(1,at=seq(800,300,-100))

m = 1+dtrapezoid(-630,-660,-550,-440,-420)*scl
polygon(x=c(660,630,630,660),y=c(1,m,1,1),col='darkgrey',border=NA)
sum(dtrapezoid(-seq(660,630,-1),-660,-550,-440,-420))
text(x=700,y=5,labels=expression(paste(P[x],"(660-630)=0.02")),cex=cx)

m1 = 1+dtrapezoid(-600,-660,-550,-440,-420)*scl
m2 = 1+dtrapezoid(-540,-660,-550,-440,-420)*scl
polygon(x=c(600,550,540,540,600),y=c(m1,6,m2,1,1),col='darkgrey',border=NA)
sum(dtrapezoid(-seq(600,540,-1),-660,-550,-440,-420))
text(x=680,y=7,labels=expression(paste(P[x],"(600-540)=0.28")),cex=cx)

polygon(x=c(500,447,447,500),y=c(6,6,1,1),col='darkgrey',border=NA)
sum(dtrapezoid(-seq(500,447,-1),-660,-550,-440,-420))
text(x=400,y=7,labels=expression(paste(P[x],"(500-447)=0.31")),cex=cx)
lines(x=c(665,639),y=c(4.5,1.6))
lines(x=c(595,571),y=c(6.5,4.5))
lines(x=c(425,470),y=c(6.5,5.4))
polygon(x=c(660,550,440,420,660),y=c(1,6,6,1,1),col='NA')


# panel c
plot(0,0,xlim=c(800,300),ylim=c(0,8),type='n',xlab="BC",ylab="",axes=F,main="Phase boundary uncertainty",cex.main=cx)
axis(1,at=seq(800,300,-100))
polygon(x=c(680,550,430,390,680),y=c(1,6,6,1,1),border=rgb(0,0,0,0.5))
polygon(x=c(640,540,420,400,680),y=c(1,6,6,1,1),border=rgb(0,0,0,0.5))
polygon(x=c(670,555,470,396,680),y=c(1,6,6,1,1),border=rgb(0,0,0,0.5))
polygon(x=c(682,556,425,380,680),y=c(1,6,6,1,1),border=rgb(0,0,0,0.5))
polygon(x=c(650,556,423,330,680),y=c(1,6,6,1,1),border=rgb(0,0,0,0.5))
points(c(680,640,670,682,650),y=c(1,1,1,1,1),pch=20)
points(c(550,540,555,556,556),y=c(6,6,6,6,6),pch=20)
points(c(430,420,470,425,434),y=c(6,6,6,6,6),pch=20)
points(c(390,400,396,380,330),y=c(1,1,1,1,1),pch=20)
arrows(x0=680,x1=640,y=0.5,y1=0.5,length=0.02,angle=90,code=3)
text(660,0.1,labels="a",cex=cx)
arrows(x0=556,x1=540,y=6.5,y1=6.5,length=0.02,angle=90,code=3)
text(548,7,labels="b",cex=cx)
arrows(x0=470,x1=420,y=6.5,y1=6.5,length=0.02,angle=90,code=3)
text(445,7,labels="c",cex=cx)
arrows(x0=400,x1=330,y=0.5,y1=0.5,length=0.02,angle=90,code=3)
text(365,0.1,labels="d",cex=cx)

dev.off()



## Figure 3 ####
load("../R_Images/simRes1.RData")
pdf(file = "./figure3.pdf",width = 5.5,height = 8)
par(mfrow=c(3,2),mar=c(5,5,1,1))
for (i in 1:length(LL))
{
  breaks=LL[[i]]
  tmp=simRes1[[i]]
  avg = apply(tmp,1,mean,na.rm=TRUE)
  plot(0,0,type='n',xlab='BC',ylab='Number of Events',xlim=c(800,300),ylim=c(0,250),axes=FALSE)
  
  axis(1)
  axis(2)
  apply(tmp[,sample(1:1000,size=100)],2,lines,x=midPoints,col=rgb(0,0,0,0.05)) #
  lines(midPoints,avg,type='b',pch=20)
  lines(d$yr,(d$p/sum(d$p)*n*50),lwd=2,lty=2,col='darkred')
  legend("bottomright",legend=letters[i],bty='n',cex=2)
  
  for (b in 1:(length(breaks)))
  {
    col='lightgrey'
    if (as.logical(b%%2)){col='darkgrey'}
    rect(xleft=breaks[b],xright=breaks[b+1],ybottom=210,ytop=250,border=NA,col=col)
    text(x=breaks[b+1]+(breaks[b]-breaks[b+1])/2,y=230,labels=as.roman(b))
  }
  
}
dev.off()




## Figure 4 ####
load("../R_Images/simRes2.RData")
pdf(file = "./figure4.pdf",width = 9,height = 6)
layout(matrix(c(1,2,1,3,1,4),nrow=2,ncol=3))
par(mar=c(4, 4, 1, 1)) #c(bottom, left, top, right)
plot(simyears, mediankd, type="l",  xlim=c(-1750,-751), xlab="Years BCE", ylab="",axes=FALSE)
axis(1,at=seq(-1600,-800,+200),labels=seq(1600,800,-200))
axis(2)
box()

for (d in 1:ncol(simmat)){
  lines(simyears,simmat[,d], col=rgb(191,191,191,alpha=75,max=255))
}
abline(h=1/length(simyears), col="red", lty="dashed")
lines(simyears,mediankd, col="black")
legend("topleft", legend=c("actual mine density", "simulated densities","median simulated"), col=c("red", "grey75","black"), lwd=c(1,1,1), lty=c("dashed","solid","solid"), bty="n", cex=1)
mtext("A", side=1, font=2, cex=1.2, line=2, adj=-0.055)

par(mar=c(4, 3, 1, 1)) #c(bottom, left, top, right)
plot(pooledEarly[,1],pooledEarly[,2],type='l', col="white", xlim=c(0,750), xlab="Years", ylab="", ylim=c(0,0.018))
polygon(c(pooledEarly[,1],rev(pooledEarly[,1]),pooledEarly[1,1]), c(pooledEarly[,2],rep(0,length(pooledEarly[,2])),pooledEarly[1,2]), col="plum1", border=NA)
lines(pooled[,1],pooled[,2], col="grey25", lwd=2)
legend("topright", lwd=c(3,1), col=c("plum1","grey25"), legend=c("Early mines (c1700 BCE)","mean (EarlyMiddle/Late)"), bty="n", cex=0.7)
mtext("B", side=1, font=2, cex=1.2, line=2, adj=-0.15)
plot(pooledMiddle[,1],pooledMiddle[,2],type='l', col="white", xlim=c(0,750), xlab="Years", ylab="", ylim=c(0,0.018))
polygon(c(pooledMiddle[,1],rev(pooledMiddle[,1]),pooledMiddle[1,1]), c(pooledMiddle[,2],rep(0,length(pooledMiddle[,2])),pooledMiddle[1,2]), col="cadetblue", border=NA)
lines(pooled[,1],pooled[,2], col="grey25", lwd=2)
legend("topright", lwd=c(3,1), col=c("cadetblue","grey25"), legend=c("Middle mines (c1300 BCE)","mean (EarlyMiddle/Late)"), bty="n", cex=0.7)
mtext("C", side=1, font=2, cex=1.2, line=2, adj=-0.15)
plot(pooledLate[,1],pooledLate[,2],type='l', col="white", xlim=c(0,750), xlab="Years", ylab="", ylim=c(0,0.018))
polygon(c(pooledLate[,1],rev(pooledLate[,1]),pooledLate[1,1]), c(pooledLate[,2],rep(0,length(pooledLate[,2])),pooledLate[1,2]), col="palegreen", border=NA)
lines(pooled[,1],pooled[,2], col="grey25", lwd=2)
legend("topright", lwd=c(3,1), col=c("palegreen","grey25"), legend=c("Late mines (c900 BCE)","mean (EarlyMiddle/Late)"), bty="n", cex=0.7)
mtext("D", side=1, font=2, cex=1.2, line=2, adj=-0.15)
dev.off()

## Figure 5 ####
load("../R_Images/simRes3.RData")
pdf(file = "./figure5.pdf",width = 5,height = 4.5)
par(mar=c(6,4,3.3,1),bg="white")
layout(matrix(c(1,1,2,1),2,2),width=c(1,1),height=c(1,1))

plot(0,0,type='n',xlab=c("b (Sampling Bias)"),ylab="Observed Percentage Change",xlim=c(0.5,3.5),ylim=range(simRes3$pr),axes=FALSE)

colSeq=c("darkblue","darkorange","darkgrey")
alpha=0.2
colSeq2=c(rgb(0,0,0.54,alpha),rgb(1,0.55,0,alpha),rgb(0.66,0.66,0.66,alpha))
mids=c(0.75,1,1.25)
for (i in 1:length(b))
{
  for (j in 1:length(r))
  {
    bb=b[i]
    rr=r[j]
    y = subset(simRes3,b==bb&r==rr)$pr
    points(x=i-1+runif(100,mids[j]-0.05,mids[j]+0.05),y=y,pch=20,col=colSeq2[j])
    rect(ybottom=quantile(y,0.25),ytop=quantile(y,0.75),xleft=i-1+mids[j]-0.07,xright=i-1+mids[j]+0.07,border=colSeq[j])
    lines(x=c(i-1+mids[j]-0.07,i-1+mids[j]+0.07),y=c(median(y),median(y)),lwd=2,col=colSeq[j])
  }
}

axis(1,at=c(-2,1,2,3,4),labels=c(NA,0,0.3,0.7,NA))
axis(2)
abline(h=0,lty=2,lwd=1)
text(x=3,y=3,labels="True Percentage Change",cex=0.8)
legend("bottomleft",bty='n',legend=c("r=0.1","r=0.3","r=0.7"),col=c("darkblue","darkorange","darkgrey"),pch=20,title = "Sampling Fraction")


# Add Figure with Skewed and Dispersed settlement pattern
set.seed(122)
nucleated=simdata3$t1[sample(length(simdata3$t1),size=20)]
dispersed=simdata3$t2[sample(length(simdata3$t2),size=20)]
par(mar=c(3,2,3,1))
plot(0,0,type='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-0.1,1.1),axes=FALSE)

points(runif(20)-0.8,runif(20),cex=nucleated/30,pch=21,col=rgb(0,0,0,0.8),bg=rgb(1,0.2,0,0.8))
text(x=-0.30,y=-0.05,labels = "Nucleatation (1sr period)",cex=0.5)
points(runif(20)+0.8,runif(20),cex=dispersed/30,pch=21,col=rgb(0,0,0,0.8),bg=rgb(1,0.2,0,0.8))
text(x=1.30,y=-0.05,labels = "Dispersed Pattern (2nd period)",cex=0.5)
box(col='lightgrey')
dev.off()



















set.seed(133)
simdata = sim.settlement(K1=1000,K2=1000)
nsim = 100
b = c(0,0.3,0.7)
r = c(0.1,0.3,0.7)
params=expand.grid(r=r,b=b,nsim=1:nsim)
params$pr=NA

for (i in 1:nrow(params))
{
  params$pr[i]=biasedsampling(simdata,r=params$r[i],b=params$b[i])
}


pdf(file = "./figure5.pdf",width = 5,height = 4.5)
par(mar=c(6,4,3.3,1),bg="white")
layout(matrix(c(1,1,2,1),2,2),width=c(1,1),height=c(1,1))
plot(0,0,type='n',xlab=c("b (Sampling Bias)"),ylab="Observed Percentage Change",xlim=c(0.5,3.5),ylim=range(params$pr),axes=FALSE)

colSeq=c("darkblue","darkorange","darkgrey")
alpha=0.2
colSeq2=c(rgb(0,0,0.54,alpha),rgb(1,0.55,0,alpha),rgb(0.66,0.66,0.66,alpha))
mids=c(0.75,1,1.25)
for (i in 1:length(b))
{
  for (j in 1:length(r))
  {
    bb=b[i]
    rr=r[j]
    y = subset(params,b==bb&r==rr)$pr
    points(x=i-1+runif(100,mids[j]-0.05,mids[j]+0.05),y=y,pch=20,col=colSeq2[j])
    rect(ybottom=quantile(y,0.25),ytop=quantile(y,0.75),xleft=i-1+mids[j]-0.07,xright=i-1+mids[j]+0.07,border=colSeq[j])
    lines(x=c(i-1+mids[j]-0.07,i-1+mids[j]+0.07),y=c(median(y),median(y)),lwd=2,col=colSeq[j])
  }
}

axis(1,at=c(-2,1,2,3,4),labels=c(NA,0,0.3,0.7,NA))
axis(2)
abline(h=0,lty=2,lwd=1)
text(x=3,y=3,labels="True Percentage Change",cex=0.8)
legend("bottomleft",bty='n',legend=c("r=0.1","r=0.3","r=0.7"),col=c("darkblue","darkorange","darkgrey"),pch=20,title = "Sampling Fraction")


# Add Figure with Skewed and Dispersed settlement pattern
set.seed(122)
nucleated=simdata$t1[sample(length(simdata$t1),size=20)]
dispersed=simdata$t2[sample(length(simdata$t2),size=20)]
par(mar=c(3,2,3,1))
plot(0,0,type='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-0.1,1.1),axes=FALSE)

points(runif(20)-0.8,runif(20),cex=nucleated/30,pch=21,col=rgb(0,0,0,0.8),bg=rgb(1,0.2,0,0.8))
text(x=-0.30,y=-0.05,labels = "Nucleatation (1sr period)",cex=0.5)
points(runif(20)+0.8,runif(20),cex=dispersed/30,pch=21,col=rgb(0,0,0,0.8),bg=rgb(1,0.2,0,0.8))
text(x=1.30,y=-0.05,labels = "Dispersed Pattern (2nd period)",cex=0.5)
box(col='lightgrey')
dev.off()



