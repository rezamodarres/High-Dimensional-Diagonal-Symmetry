library(mvtnorm)
library(tidyverse)    # data manipulation and visualization
library(RColorBrewer) # customized coloring of plots
library(MASS)
library(matrixStats)
library(spls)
library(MatrixGenerics)
library(ISLR)        
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
 
 

delta <- function(both) 
 
{ 
myN=nrow(both)
Newzbar=matrix(0,nrow=myN,ncol=myN)
con=1/(myN-2)
kk=1   
for (ix in 1:(myN-1))
{ xx=both[ix,];        mux=mean(xx);   vx=sqrt(var(xx))
  for (iy in (ix+1):myN)
    {   yy=both[iy,];      muy=mean(yy);   vy=sqrt(var(yy))
        part=0
        for (it in 1:myN)  if ((it!=ix) & (it!=iy)) {   
            zz=both[it,]
            muz=mean(zz)
            vz=sqrt(var(zz))
            dxz=sqrt(((mux-muz)^2+(vx-vz)^2)) 
            dyz=sqrt(((muy-muz)^2+(vy-vz)^2))
            part=part+abs(dxz-dyz)
         }
      Newzbar[ix,iy]<- con*part
      Newzbar[iy,ix]<-Newzbar[ix,iy]
      kk=kk+1
    }
}
return(Newzbar)  
}


##################################  Main *********************************

################################chin############################

load(file = "C:/Users/rmodarres/Desktop/Datasets/chin.Rdata")
str(chin)
ALLcomsam=chin$x
show(dim(ALLcomsam))
yy=chin$y
show(yy)
nx=n1=sum(yy=='positive')
ny=n2=sum(yy=='negative')
cat('pos',nx,'neg=',ny, fill=TRUE)
obslab=ifelse(yy=='positive',1,2)


xmat=ALLcomsam[obslab==1,]
show(dim(xmat))
p=px=dim(xmat)[2]
N=mm=dim(xmat)[1]
n1=trunc(N/2)
n2=N-n1
n1p=n1+1
sim=1000
numtest=5
con1=mm*(mm-1)/2
con2=con1^2
cat('xp=',xp,fill=TRUE)
BF_E=BG_E=BF_Rho=BG_Rho=BF_Delta=BG_Delta=zzu=mori=ASSDEL=rep(0,sim)
TSV=matrix(0,nrow=sim,ncol=numtest)

plist=c(1000,2000,4000,10000)
mysample=xmat
for (xp in plist)
{
xmat=mysample[,1:xp]
samplex1=xmat[1:n1,]
samplex2=-xmat[n1p:N,]
 

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
OBSASSDEL=(AAmux-AAmuy)^2+ (sxp-syp)^2


#  if(EST==1) {xbar=colMeans2(comsam); comsam=comsam-xbar}
#   if (ANG==1) for (cc in 1:N) comsam[,cc]=comsam[,cc]/sqrt(sum(comsam[,cc]^2))
mydist <- as.matrix(dist(xmat,method = "euclidean"))
avgx=mean(mydist[1:n1,1:n1])
avgy=mean(mydist[(n1+1):N,(n1+1):N])
avgxy=mean(mydist[1:n1,(n1+1):N]) 
OBSBF_E=2*avgxy-avgx-avgy

L0=delta(xmat)
avgx=mean(L0[1:n1,1:n1])
avgy=mean(L0[(n1+1):N,(n1+1):N])
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
  ASSDEL[it]=(AAmux-AAmuy)^2+ (sxp-syp)^2
  

  #  if(EST==1) {xbar=colMeans2(comsam); comsam=comsam-xbar}
  #   if (ANG==1) for (cc in 1:N) comsam[,cc]=comsam[,cc]/sqrt(sum(comsam[,cc]^2))
  mydist <- as.matrix(dist(xmat,method = "euclidean"))
  avgx=mean(mydist[1:n1,1:n1])
  avgy=mean(mydist[(n1+1):N,(n1+1):N])
  avgxy=mean(mydist[1:n1,(n1+1):N]) 
  BF_E[it]=2*avgxy-avgx-avgy
  
  L0=delta(xmat)
  avgx=mean(L0[1:n1,1:n1])
  avgy=mean(L0[(n1+1):N,(n1+1):N])
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
hist(ASSDEL)
hist(BF_E)
hist(BF_Delta)
hist(zzu)
hist(mori)
     pASSDEL=sum(ASSDEL>=OBSASSDEL)/sim
     pBF_E=sum(BF_E>=OBSBF_E)/sim
     pBF_Delta=sum(BF_Delta>=OBSBF_Delta)/sim
     pzzu=sum(zzu>=OBSzzu)/sim
     Pmori =sum(mori>=OBSmori)/sim
 
 zz=cbind(px,pASSDEL,  pBF_E,pBF_Delta,pzzu,  Pmori)
  show(zz)
}  
