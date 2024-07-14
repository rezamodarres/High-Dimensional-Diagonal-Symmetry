
library(datamicroarray)
library(MASS)
library("pamr")
options(digits=4)


########################################
load(file = "C:/Users/rmodarres/Desktop/Datasets/chin.Rdata")
str(chin)
ALLcomsam=chin$x
show(dim(ALLcomsam))
yy=chin$y
show(yy)
nx=n1=sum(yy=='positive')
ny=n2=sum(yy=='negative')

N=CAPN=nx+ny 
cat('pos',nx,'neg=',ny, fill=TRUE)
obslab=ifelse(yy=='positive',1,2)
xmat=ALLcomsam[obslab==1,]
show(dim(xmat))
 
p=px=dim(xmat)[2]
N=mm=dim(xmat)[1]
n1=nx=trunc(N/2)
n2=ny=N-n1
n1p=n1+1
xg1=xmat[1:n1,]
yg1=-xmat[n1p:N,]
mx=n1*(n1-1)/2
my=n2*(n2-1)/2
show(n1)
show(n2)
show(p)

########################
xdmat <- as.matrix(dist(xg1,method = "euclidean"))
ipdx=rep(0,mx)
count=1
for(i in 1:(nx-1)){for(j in (i+1):nx){ipdx[count]=xdmat[i,j]; count <- count+1 }}
 
show(ipdx)
 
ydmat= as.matrix(dist(yg1,method = "euclidean"))
ipdy=rep(0,my)
count=1
for(i in 1:(ny-1)){for(j in (i+1):ny){ipdy[count]=ydmat[i,j]; count <- count+1 }}
show(ipdy)
 
hist(ipdx,main='IPD of Chin Data (Images)',xlab='IPD',col = 'red')
hist(ipdy,main='IPD of Chin Data (Mirror Images)',xlab='IPD',col = 'green')
cat('mean(ipdx)=',mean(ipdx),fill=TRUE)
cat('sd(ipdx)=',sqrt(var(ipdx)),fill=TRUE)
cat('mean(ipdy)=',mean(ipdy),fill=TRUE)
cat('sd(ipdy)=',sqrt(var(ipdy)),fill=TRUE)
ipd=c(ipdx,ipdy)
ipdxbar=sum(ipdx)/mx
ipdybar=sum(ipdy)/my
 
ipdbars=c(ipdxbar,ipdybar)
show(ipdbars)
minipd=min(ipd)
maxipd=max(ipd)
deltaint=100
unit=(maxipd-minipd)/deltaint
delta=seq(minipd,maxipd,unit)
cdfx=cdfy =rep(0,deltaint)
for (i in 1:deltaint){
   cdfx[i]=sum(ipdx<=delta[i])/mx
   cdfy[i]=sum(ipdy<=delta[i])/my
   }
plot(cdfx,type="l",lty=2,  col="red", lwd=1,ylab="ECDF", xlab = "IPD")
title(main = "ECDFs of Chin Data")
lines(cdfy,lty=3,lwd=2,col="green")
legend("topleft",legend=c('Image','Mirror image'),  col=c("red",'green'), lty=c(1,2),lwd=2, ncol=1)
       
 