#Figure Log Files
source("../src/utility.R")
library(trapezoid)
library(reldist)
library(ggplot2)
library(ggpubr)


## Figure 1 ####

## Figure 2 ####
pdf(file = "./figure2.pdf",width = 7,height = 6)

chiba=read.csv("../data/chibaJomon_fromCrema2013.csv")
chiba = subset(chiba,!is.na(chiba$SimplifiedPhase))
chiba$SimplifiedPhase=factor(chiba$SimplifiedPhase,levels=c("Otamadai-Katsuzaka","Kasori E (early)","Kasori E (late)","Shomyoji","Horinouchi","Kasori B","Soya","Angyo 1 & 2"))
x = table(chiba$SiteID,chiba$SimplifiedPhase)

par(mar=c(9,5,5,3))
plot(0,0,xlab="",ylab="",xlim=c(0,9),axes=FALSE,type='n',ylim=c(0,max(apply(x,2,sum))+20))

gg = numeric(length=8)

for (i in 1:8)
{
  tmp=x[,i]
  tmp=tmp[which(tmp>0)]
  gg[i]=gini(tmp)
  tmp=c(0,cumsum(tmp))
  for (k in 1:c(length(tmp)-1))
  {
    col="lightgrey"
    if (names(tmp[k+1])=='10529'){col='orange'} #Ariyoshi-Kita Site
    if (names(tmp[k+1])=='9825'){col='royalblue'} #Miyauchiidosaku Site
    rect(xleft=i-0.4,xright=i+0.4,ybottom=tmp[k],ytop=tmp[k+1],col=col)
  }
  text(x=i,y=tmp[k+1]+10,labels=paste0("n=",length(tmp)-1))
}

axis(2)
axis(3,at=1:8,labels=round(gg,2))
mtext(side=2,line=3,"Number of Pithouses",cex=1)
mtext(side=3,line=3,"Gini Coefficient (settlement size)",cex=1)
mtext(side=1,line=6.5,"Ceramic Phases",cex=1)
axis(1,at=1:8,labels=levels(chiba$SimplifiedPhase),las=2)
legend(6.5,220,legend=c("Ariyoshi-Kita","Aioi","Other Sites"),fill=c("orange","royalblue","grey"))
dev.off()




## Figure 3 ####
pdf(file = "./figure3.pdf",width = 10,height = 3)

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



## Figure 4 ####

pdf(file = "./figure4.pdf",width = 5,height = 8)

age = seq(300,800,1)
asymptote = 0.9
sl = .04
midpoint = 550
p = 0.1+asymptote / (1 + exp((midpoint - age) * sl))
d=data.frame(yr=rev(age),p=p)

par(mfrow=c(3,2))
plot(d$yr,d$p,xlim=c(800,300),type='l',col='darkgrey',ylim=c(0,1),xlab="BC",axes=FALSE,ylab="Population Size")
axis(1,at=seq(800,300,-100))
legend("bottomright",legend=c("a"),bty='n',cex=2)

LL = list(c(800,550,300),c(800,700,300),c(800,400,300),c(800,700,400,300),c(800,700,600,500,400,300))
resolution=50
bbseq=seq(800-resolution/2,300+resolution/2,-resolution)
for (x in 1:length(LL))
{
  breaks=LL[[x]]
  tmp2=mcunif(n=1000,p=d$p,years=d$yr,resolution=resolution,breaks=breaks)
  
  plot(0,0, xlim=c(800,300),type='n',ylim=c(0,max(tmp2)*1.2),axes=FALSE,ylab="",xlab="BC")
  axis(side=1,at=seq(800,300,-100))
  
  for (b in 1:(length(breaks)))
  {
    col='lightgrey'
    if (as.logical(b%%2)){col='darkgrey'}
    rect(xleft=breaks[b],xright=breaks[b+1],ybottom=max(tmp2)*1,ytop=max(tmp2)*1.1,border=NA,col=col)
    text(x=breaks[b+1]+(breaks[b]-breaks[b+1])/2,y=max(tmp2)*1.05,labels=as.roman(b))
  }
  
  lines(bbseq,apply(tmp2,1,mean),type='b',pch=20)
  apply(tmp2[,1:100],2,lines,x=bbseq,col=rgb(0,0,0,0.05))
  axis(2)
  
  legend("bottomright",legend=letters[x+1],bty='n',cex=2)
  
}

dev.off()


## Figure 5 ####
pdf(file = "./figure5.pdf",width = 5,height = 5)

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

# Plot results
params$r=as.factor(params$r)
p <- ggboxplot(params, x = "b", y = "pr",
               color = "r", palette = "jco",
               add = "jitter",ggtheme=theme_grey(),alpha=0.9) + ylab("Percentage Change") + geom_hline(yintercept = 0, linetype="dashed") +theme(legend.position = 'top')

p
dev.off()
