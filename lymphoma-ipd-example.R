
library(datamicroarray)
library(MASS)
library(spls)
library("pamr")
options(digits=4)


options(digits=5)

########################################
data(lymphoma)
lym=as.data.frame(lymphoma)
str(lym)
xxlym=lym
yylym=lym$y
show(yylym)
n1=sum(yylym==0)
n2=sum(yylym==1)
n3=sum(yylym==2)
N=dim(lym)[1]
p=px=dim(lym)[2]
cat('p=',p,'N=',N,'n1=',n1,'n2=',n2,'n3=',n3,fill=TRUE)
# toform group 1=2 use
samplex1=xxlym[yylym==0,]
samplex2=xxlym[yylym==1,]
samplex3=xxlym[yylym==2,]

xg=samplex1
N=dim(xg)[1]
n1=nx=N/2
n2=ny=N/2
n1p=n1+1
 mx=n1*(n1-1)/2
 my=n2*(n2-1)/2
p=px=dim(xg)[2]
xg1=xg[1:n1,]
yg1=-xg[n1p:N,]

  
########################
xmat <- as.matrix(dist(xg1,method = "euclidean"))
ipdx=rep(0,mx)
count=1
for(i in 1:(nx-1)){for(j in (i+1):nx){ipdx[count]=xmat[i,j]; count <- count+1 }}
 
show(ipdx)
 
ymat= as.matrix(dist(yg1,method = "euclidean"))
ipdy=rep(0,my)
count=1
for(i in 1:(ny-1)){for(j in (i+1):ny){ipdy[count]=ymat[i,j]; count <- count+1 }}
show(ipdy)
 
hist(ipdx,main='IPD of Lymphoma Data (Images)',xlab='IPD',col = 'red')
hist(ipdy,main='IPD of Lymphoma Data (Mirror Images)',xlab='IPD',col = 'green')
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
title(main = "ECDFs of Lymphoma Data")
lines(cdfy,lty=3,lwd=2,col="green")
legend("topleft",legend=c('Image','Mirror image'),  col=c("red",'green'), lty=c(1,2),lwd=2, ncol=1)
       


####################################################
