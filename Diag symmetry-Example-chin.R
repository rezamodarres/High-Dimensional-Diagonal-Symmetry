library(mvtnorm)
library(tidyverse)    # data manipulation and visualization
library(RColorBrewer) # customized coloring of plots
library(MASS)
library(matrixStats)
library(spls)
library(MatrixGenerics)
library(ISLR)         # contains example data set "Khan"
library(LaplacesDemon)
library(spls)

rm(list = ls()) #clear workspace
options(digits=3)

Rho0 <- function(ipdmat) 
{
myN=nrow(ipdmat)
R0=matrix(0,nrow=myN,ncol=myN)
con=1/((myN-2)*sqrt(px))
for (ix in 1:(myN-1))
{ 
   for (iy in (ix+1):myN)
   {   
      part1=ipdmat[ix,]
      p1=part1[-ix]
      iy1=iy-1
      Q1=p1[-iy1]
      part2=ipdmat[iy,]
      p2=part2[-ix]
      Q2=p2[-iy1]
     tt=mean(abs(Q1-Q2))
      R0[ix,iy]=con*mean(abs(Q1-Q2))
      R0[iy,ix]= R0[ix,iy]    
   }
}
return(R0)  
}
 
 

delta <- function(both) { 
  myN=nrow(both)
  myp=ncol(both)
  Newzbar=matrix(0,nrow=myN,ncol=myN)
  con=1/(myN-2)
  Mux=rowMeans2(both)
  sdx=sqrt(rowVars(both))
  for (ix in 1:(myN-1))
  { for (iy in (ix+1):myN)
  { part=0; RG=1:myN;  RG=RG[c(-ix,-iy)];   
  for (it in RG)  
  { dxz=sqrt(((Mux[ix]-Mux[it])^2+(sdx[ix]-sdx[it])^2)) 
  dyz=sqrt(((Mux[iy]-Mux[it])^2+(sdx[iy] -sdx[it])^2))
  part=part+abs(dxz-dyz)
  }
  Newzbar[ix,iy]<- con*part
  Newzbar[iy,ix]<-Newzbar[ix,iy]
  }  }
  return(Newzbar)  
}

 
##################################  Main *********************************
#set.seed(3166779)
 

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
n1=trunc(N/2)
n2=N-n1
n1p=n1+1
samplex1=xmat[1:n1,]
samplex2=-xmat[n1p:N,]
show(n1)
show(n2)
sim=100
numtest=5

BF_E=BG_E=BF_Rho=BG_Rho=BF_Delta=BG_Delta=zzu=mori=ASSDEL=rep(0,sim)
TSV=matrix(0,nrow=sim,ncol=numtest)
m1=n1*(n1-1)/2
m2=n2*(n2-1)/2
plist=c(500,1000)
mysample=xmat
for (xp in plist)
{
xmat=mysample[,1:xp]

N=dim(xmat)[1]
n1=trunc(N/2)
n2=N-n1
mm=N
n1p=n1+1

p=px=dim(xmat)[2]
samplex1=xmat[1:n1,]
samplex2=-xmat[n1p:N,]
#######################
con1=mm*(mm-1)/2
con2=con1^2
cat('xp=',xp,fill=TRUE)
############################### H0 #####################
 
# asymptotic delta-hat stat
samplex1=as.matrix(samplex1)
samplex2=as.matrix(samplex2)
Ssigx=cov(samplex1)
Ssigy=cov(samplex2)
sigF=sum(diag(Ssigx))/px
sigG=sum(diag(Ssigy))/px
#  cat('sigF=',sigF,'sigF=G=',sigG,fill=TRUE)
Amux= rowMeans2(samplex1)
AAmux=mean(Amux)
Amuy= rowMeans2(samplex2)
AAmuy=mean(Amuy)
sxp=mean(sqrt(rowVars(samplex1)))
syp=mean(sqrt(rowVars(samplex2)))
OBSASSDEL=(AAmux-AAmuy)^2 
mydist <- as.matrix(dist(xmat,method = "euclidean"))
avgx=sum(mydist[1:n1,1:n1])/m1
avgy=sum(mydist[(n1+1):N,(n1+1):N])/m2
avgxy=mean(mydist[1:n1,(n1+1):N]) 
OBSBF_E=2*avgxy-avgx-avgy

L0=delta(xmat)
avgx=sum(L0[1:n1,1:n1])/m1
avgy=sum(L0[(n1+1):N,(n1+1):N])/m2
avgxy=mean(L0[1:n1,(n1+1):N]) 
OBSBF_Delta=2*avgxy-avgx-avgy

# Y. Song
num1=num2=0

for (t1 in 1:(N-1)){
  { for (t2 in (t1+1):N)  
    xpt1=xmat[t1,]+xmat[t2,]
    xpt2=xmat[t1,]-xmat[t2,]
    sum1=sqrt(sum(xpt1*xpt1))
    diff1=sqrt(sum(xpt2*xpt2))
    h=sum1 - diff1
    num1=num1+h 
    num2=num2+h^2
  }}
Un=num1/con1
sig02=num2/con2
OBSzzu=Un/sqrt(sig02)

# G. J. Szckely T. F. Moori
sum1=0
for (t1 in 1:N){sum1=sum1+sqrt(sum(xmat[t1,]^2))}
OBSmori=num1/sum1

######################################### Ha   ############################

cat('OBSASSDEL=',OBSASSDEL,'OBS BF_E=',OBSBF_E,'OBSBF_Delta=',OBSBF_Delta,'OBSzzu=',OBSzzu,'OBSmori=',OBSmori,fill=TRUE)

#      hold=cbind(sBFe,sBFd,szzu,smori,Sdel)
for (it in 1:sim)
   {
  indx1=sample(1:N,n1)
  indx2=setdiff(1:N,indx1)
  samplex1=xmat[indx1,]
  samplex2=-xmat[indx2,]
  xmat=rbind(samplex1,samplex2)
  samplex1=as.matrix(samplex1)
  samplex2=as.matrix(samplex2)
  # asymptotic delta-hat stat
  Ssigx=cov(samplex1)
  Ssigy=cov(samplex2)
  sigF=sum(diag(Ssigx))/px
  sigG=sum(diag(Ssigy))/px
  #  cat('sigF=',sigF,'sigF=G=',sigG,fill=TRUE)
  Amux= rowMeans2(samplex1)
  AAmux=mean(Amux)
  Amuy= rowMeans2(samplex2)
  AAmuy=mean(Amuy)
  sxp=mean(sqrt(rowVars(samplex1)))
  syp=mean(sqrt(rowVars(samplex2)))
  ASSDEL[it]=(AAmux-AAmuy)^2 

  mydist <- as.matrix(dist(xmat,method = "euclidean"))
  avgx=sum(mydist[1:n1,1:n1])/m1
  avgy=sum(mydist[(n1+1):N,(n1+1):N])/m2
  avgxy=mean(mydist[1:n1,(n1+1):N]) 
  BF_E[it]=2*avgxy-avgx-avgy
  
  L0=delta(xmat)
  avgx=sum(L0[1:n1,1:n1])/m1
  avgy=sum(L0[(n1+1):N,(n1+1):N])/m2
  avgxy=mean(L0[1:n1,(n1+1):N]) 
  BF_Delta[it]=2*avgxy-avgx-avgy
  
  # Y. Song
  num1=num2=0
  for (t1 in 1:(N-1)){
    { for (t2 in (t1+1):N)  
      xpt1=xmat[t1,]+xmat[t2,]
    xpt2=xmat[t1,]-xmat[t2,]
    sum1=sqrt(sum(xpt1*xpt1))
    diff1=sqrt(sum(xpt2*xpt2))
    h=sum1 - diff1
    num1=num1+h 
    num2=num2+h^2
    }}
  Un=num1/con1
  sig02=num2/con2
  zzu[it]=Un/sqrt(sig02)
  
  # G. J. Szckely T. F. Moori
  sum1=0
  for (t1 in 1:N){sum1=sum1+sqrt(sum(xmat[t1,]^2))}
  mori[it]=num1/sum1
 
} 

     pASSDEL=sum(ASSDEL>=OBSASSDEL)/sim
     pBF_E=sum(BF_E>=OBSBF_E)/sim
     pBF_Delta=sum(BF_Delta>=OBSBF_Delta)/sim
     pzzu=sum(zzu>=OBSzzu)/sim
     Pmori =sum(mori>=OBSmori)/sim
 
 zz=cbind(px,pASSDEL,  pBF_E,pBF_Delta,pzzu,  Pmori)
  show(zz)
}  
 # plot(muxx, pBFe,      type="b",lty=1, lwd=1, pch=3,  bg='red',  col="red",  ylim = c(0, 1),xlab=expression(paste(mu)), ylab="Empirical Power")
#lines(muxx,     pBFd, type="b",lty=2, lwd=1, pch=21, bg='green',col="green"   )
#lines(muxx,     pzzu, type="b",lty=3, lwd=1, pch=4,  bg='blue', col="blue"    )
#lines(muxx,     pmori,type="b",lty=4, lwd=1, pch=2,  bg='black',col="black"    )
#lines(muxx,     PDEL, type="b",lty=5, lwd=1, pch=12, bg='brown',col="brown"    )
#legend("topleft",legend=c('R_Euc','R_delta','Z','T','A_D'), col=c("red","green","blue","black","brown"),
#lty=c(1,2,3,4,5),lwd=1.5, ncol=2,pch=c(3,21,4,2,12),cex=0.8)

  