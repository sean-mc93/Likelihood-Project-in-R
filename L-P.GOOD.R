###LIKELIHOOD PROJECT###
#
#

#Input Data
read.table("proj-lik.dat",header=TRUE)
#name data
dat <- read.table("proj-lik.dat",header=TRUE)

#defining data
xlow <- dat$xlow
xupp <- dat$xupp
#add indicator to data
d <- as.numeric(dat$xupp!=999)
dat$delta <- d
sum(d)

#No zeros==>0.0001
z <- c()
for(i in 1:94){
  if(xlow[i]==0){z[i]=0.001}else{z[i]=xlow[i]}
}
dat$xlowz <- z

#seperate data trt = 0/1
xlow.trt0 <- dat$xlow[1:46]
xupp.trt0 <- dat$xupp[1:46]
d.trt0    <- dat$delta[1:46]
xlowz.trt0<- dat$xlowz[1:46]
xlow.trt1 <- dat$xlow[47:94]
xupp.trt1 <- dat$xupp[47:94]
d.trt1    <- dat$delta[47:94]
xlowz.trt1<- dat$xlowz[47:96]

#datasets
xtrt0<-data.frame(xlow=xlow.trt0,xupp=xupp.trt0,delta=d.trt0)
xtrt1<-data.frame(xlow=xlow.trt1,xupp=xupp.trt1,delta=d.trt1)
xall<-data.frame(xlow=xlow,xupp=xupp,delta=d)

#Define n
N <- length(d)
#Count uncensored data
C <- table(dat$delta)
M <- sum(d)

#Hess function
hess<-function(f,y,x){
  
  ep<-0.0001
  eps<- ep*y
  n <- length(y)
  m <- matrix(0,ncol=n, nrow=n)
  for(i in 1:n){
    for( j in 1:n){
      y1<-y
      y1[i]<-y1[i]+eps[i]
      y1[j]<-y1[j]+eps[j]
      y2<-y
      y2[i]<-y2[i]+eps[i]
      y2[j]<-y2[j]-eps[j]
      y3<-y
      y3[i]<-y3[i]-eps[i]
      y3[j]<-y3[j]+eps[j]
      y4<-y
      y4[i]<-y4[i]-eps[i]
      y4[j]<-y4[j]-eps[j]
      m[i,j]<-(f(y1,x)-f(y2,x)-f(y3,x)+f(y4,x))/(4*eps[i]*eps[j])
    }
  }
  
  solve(m)
  
}

#EXPONENTIAL#

#Lower, Upper, and Start for optim exp
bl <- 0.0001
bu <- Inf
start <- 1

#(-)EXPONENTIAL LIKLIHOOD FUNCTION
llik.exp <- function(p,l,u,d){
 
  n <- length(l)
  m <- sum(d)
  
  return(
    
  p*sum(l*(1-d)) -
    
  sum(
    (log(exp((-l)*p)-exp((-u)*p)))*d
  )  
  
  )    
}


llik.exp2 <- function(p,xtrt0){
  
  n <- length(xtrt0$xlow)
  m <- sum(xtrt0$delta)
  
  return(
    
    p*sum(xtrt0$xlow*(1-xtrt0$delta)) -
      
      sum(
        (log(exp((-xtrt0$xlow)*p)-exp((-xtrt0$xupp)*p)))*xtrt0$delta
      )  
    
  )    
}

llik.exp3 <- function(p,x){
  
  n <- length(x$xlow)
  m <- sum(x$delta)
  
  return(
    
    p*sum(x$xlow*(1-x$delta)) -
      
      sum(
        (log(exp((-x$xlow)*p)-exp((-x$xupp)*p)))*x$delta
      )  
    
  )    
}

#OPTIM FUNCTIONS#
#
#OPTIM FUNCTION FOR EXP TRT=0
y.trt0 <- optim(start,llik.exp,method="L-BFGS-B",l=xlowz.trt0,u=xupp.trt0,d=d.trt0,lower=bl,upper=bu)
y.trt0

#OPTIM FUNCTION FOR EXP TRT=1
y.trt1 <- optim(start,llik.exp,method="L-BFGS-B",l=xlow.trt1,u=xupp.trt1,d=d.trt1,lower=bl,upper=bu)
y.trt1

#OPTIM FUNCTION FOR ALL DATA EXPONENTIAL
y <- optim(start,llik.exp2,method="L-BFGS-B",x=xall,lower=bl,upper=bu)
y

#PLOTS
#trt0
plot.ts(xtrt0,main="Time Series Plot trt=0")
qqnorm(xtrt0-y.trt0$par,main="QQ Plot ws1 on theta hat");lines(c(-10,10),c(-10,10)/sqrt(2),col=2,lwd=2)
n <- length(ws1)
plot(ws1[-1]-mean(ws1),ws1[-n]-mean(ws1),xlab="",ylab="",main="(c)")

dev.off()

help(qqnorm)

#HESS FUNCTIONS#
#
#hess trt0
vtrt0 <- hess(f=llik.exp2,y=y.trt0$par,x=xtrt0)
setrt0 <- sqrt(vtrt0)

#hess trt1
vtrt1 <-hess(f=llik.exp3,y=y.trt1$par,x=xtrt1)
setrt1 <- sqrt(vtrt1)

#hess all
vall <- hess(f=llik.exp3,y=y$par,x=xall)
seall <- sqrt(vall)

#Confidence Intervals
#trt=0 
CItrt0 <- c(y.trt0$par-(1.96)*setrt0,y.trt0$par+(1.96)*setrt0)
CItrt0

#trt=1 
CItrt1 <- c(y.trt1$par-(1.96)*setrt1,y.trt1$par+(1.96)*setrt1)
CItrt1

#ALL 
CI <- c(y$par-(1.96)*seall,y$par+(1.96)*seall)
CI


#WEIBULL#

#Lower, Upper and Start for optim weibull.
blw <- c(0.0001,0.0001)
buw <- c(Inf, Inf)
startw <- c(0.0001,0.0001)

#(-) Weibull liklihood
llik.weib <- function(p,l,u,d){
  
  n <- length(l)
  m <- sum(d)
  
  return(
    
    (p[2]^p[1])*sum((l^p[1])*(1-d)) -
      
      sum(
        (log(exp(-(p[2]*l)^p[1])-exp(-(p[2]*u)^p[1])))*d
      )
    
  )    
}

llik.weib(p=c(1,1),l=xlow.trt0,u=xupp.trt0,d=d.trt0)

llik.weib2 <- function(p,x){
  
  n <- length(x$xlow)
  m <- sum(x$xdelta)
  
  return(
    
    sum((1-x$delta)*((p[2]*x$xlow)^p[1])) 
    
    - sum(x$delta*(log(exp(-(p[2]*x$xlow)^p[1])-exp(-(p[2]*x$xupp)^p[1]))))
      
        
  )    
}

llik.weib2(p=c(1,1),x=xtrt0)


#OPTIM FUNCTIONS
#
#OPTIM FUNCTION FOR WEIBULL TRT=0
y.trt0.weib <- optim(startw,llik.weib,method="L-BFGS-B",l=xlow.trt0,u=xupp.trt0,d=d.trt0,lower=blw,upper=buw,hessian=T)
y.trt0.weib


#OPTIM FUNCTION FOR WEIBULL TRT=1
y.trt1.weib <- optim(startw,llik.weib,method="L-BFGS-B",l=xlow.trt1,u=xupp.trt1,d=d.trt1,lower=blw,upper=buw)
y.trt1.weib

#OPTIM FUNCTION FOR ALL DATA WEIBULL
y.weib <- optim(startw,llik.weib,method="L-BFGS-B",l=xlow,u=xupp,d=d,lower=blw,upper=buw)
y.weib


#HESS FUNCTIONS#
#
#Hess function trt0
vtrt0w <- hess(f=llik.weib2,y=y.trt0.weib$par,x=xtrt0)
setrt0w <- sqrt(diag(vtrt0w))

#Hess function trt1
vtrt1w <-hess(f=llik.weib2,y=y.trt1.weib$par,x=xtrt1)
setrt1w <- sqrt(diag(vtrt1w))
vtrt1w

#Hess function trt0
vallw <-hess(f=llik.weib2,y=y.weib$par,x=xall)
seallw <- sqrt(diag(vallw))


#Confidence Intervals
#trt=0 A
CItrt0weibA <- c(y.trt0.weib$par[1]-(1.96)*setrt0w[1],y.trt0.weib$par[1]+(1.96)*setrt0w[1])
CItrt0weibA
#trt=0 B
CItrt0weibB <- c(y.trt0.weib$par[2]-(1.96)*setrt0w[2],y.trt0.weib$par[2]+(1.96)*setrt0w[2])
CItrt0weibB
#trt=1 A
CItrt1weibA <- c(y.trt1.weib$par[1]-(1.96)*setrt1w[1],y.trt1.weib$par[1]+(1.96)*setrt1w[1])
CItrt1weibA
#trt=1 B
CItrt1weibB <- c(y.trt1.weib$par[2]-(1.96)*setrt1w[2],y.trt1.weib$par[2]+(1.96)*setrt1w[2])
CItrt1weibB
#ALL A
CIA <- c(y.weib$par[1]-(1.96)*seallw[1],y.weib$par[1]+(1.96)*seallw[1])
CIA
#ALL B
CIB <- c(y.weib$par[2]-(1.96)*seallw[2],y.weib$par[2]+(1.96)*seallw[2])
CIB
               
#---------------------------------------------------------------

#Liklihood ratio statistic
#Exponential

WE = -2*(
  
  llik.exp3(p=y.trt0$par,x=xtrt0) + llik.exp3(p=y.trt1$par,x=xtrt1) - llik.exp3(p=y$par,x=xall)
  
)

WW = -2*(
  
  llik.weib2(p=y.trt0.weib$par,x=xtrt0) + llik.weib2(p=y.trt1.weib$par,x=xtrt1) - llik.weib2(p=y.weib$par,x=xall)
  
)

WE
WW

#------------------------------------------------------------------

#Parametric Bootstrap
nn=100
MM=1000 ##number of replicate data sets to calculate

#trt 0 exp
#(1) Simulate Data

xsimtrt0=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimtrt0[i,]=rexp(nn,y.trt0$par) 

#(2) Calculate MLEs

mle.simtrt0=rep(0,MM) ##vector to store MLEs
for (i in 1:MM) mle.simtrt0[i]=max((mean(xsimtrt0[i,]))^(-1),0)

#(3) Confidence interval
cltrt0=quantile(mle.simtrt0,0.025)
cutrt0=quantile(mle.simtrt0,0.975)

cltrt0
cutrt0

#trt1 exp

#(1) Simulate Data

xsimtrt1=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimtrt1[i,]=rexp(nn,y.trt1$par) 

#(2) Calculate MLEs

mle.simtrt1=rep(0,MM) ##vector to store MLEs
for (i in 1:MM) mle.simtrt1[i]=max((mean(xsimtrt1[i,]))^(-1),0)

#(3) Confidence interval
cltrt1=quantile(mle.simtrt1,0.025)
cutrt1=quantile(mle.simtrt1,0.975)

cutrt1
cltrt1

#(1) Simulate Data

xsim=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsim[i,]=rexp(nn,y$par) 

#(2) Calculate MLEs

mle.sim=rep(0,MM) ##vector to store MLEs
for (i in 1:MM) mle.sim[i]=max((mean(xsim[i,]))^(-1),0)

#(3) Confidence interval
cl=quantile(mle.sim,0.025)
cu=quantile(mle.sim,0.975)

cl
cu

##NOTE it is interesting to look at the distribution of mle.sim:
hist(mle.sim)
qqnorm(mle.sim) ##NORMAL QQ PLOT


#Parametric Bootstrap Weibull
nn=100
MM=1000 ##number of replicate data sets to calculate

beta <- function(x,A){
  
  n=length(x)
  
  return(
    (n/sum(x^(A)))^(1/A)
  )
}


#(1) Simulate Data

xsimweibtrt0=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimweibtrt0[i,]=rweibull(nn,y.trt0.weib$par[1],1/(y.trt0.weib$par[2])) 

#(2) Calculate MLEs
mle.simweibtrt0=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simweibtrt0[i]=beta(x=xsimweibtrt0[i,],A=y.trt0.weib$par[1])}

#(3) Confidence interval
cltrt0w=quantile(mle.simweibtrt0,0.025)
cutrt0w=quantile(mle.simweibtrt0,0.975)

cltrt0w
cutrt0w




#(1) Simulate Data
xsimweibtrt1=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimweibtrt1[i,]=rweibull(nn,y.trt1.weib$par[1],1/(y.trt1.weib$par[2])) 

#(2) Calculate MLEs

mle.simweibtrt1=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simweibtrt1[i]=beta(x=xsimweibtrt1[i,],A=y.trt1.weib$par[1])}

#(3) Confidence interval
cltrt1w=quantile(mle.simweibtrt1,0.025)
cutrt1w=quantile(mle.simweibtrt1,0.975)

cltrt1w
cutrt1w



#(1) Simulate Data

xsimweib=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimweib[i,]=rweibull(nn,y.weib$par[1],1/(y.weib$par[2])) 

#(2) Calculate MLEs

mle.simweib=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simweib[i]=beta(x=xsimweib[i,],A=y.weib$par[1])}

#(3) Confidence interval
clw=quantile(mle.simweib,0.025)
cuw=quantile(mle.simweib,0.975)

clw
cuw


##NOTE it is interesting to look at the distribution of mle.sim:
hist(mle.simweib)
qqnorm(mle.simweib) ##NORMAL QQ PLOT


#alpha p boot test

#equation

f <- function(x,a,b){
  n=length(x)
  ( (n*a)+sum(x*a)-sum(((b*x)^a)*(log(b*x)) ))
}

fd <- function(x,a,b){
  n=length(x)
  (  n- sum(((b*x)^a)*((log(b*x)^2)) )  )
}

alpha <- function(x,a,b){
  n=length(x)
  return(
    a-( ( (n/a)+n*log(b)+sum(log(x)) -sum(((b*x)^a)*(log(b*x))) )/(  -(n/(a^2))- sum(((b*x)^a)*((log(b*x)^2)) )  ) )  
  )
  
}

#trt0


xsimalph0=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimalph0[i,]=rweibull(nn,y.trt0.weib$par[1],1/(y.trt0.weib$par[2])) 

#Newton method
a10 <- alpha(x=xsimalph0[i,], a =2 ,b=y.trt0.weib$par[2]) 
a20 <- alpha(x=xsimalph0[i,], a =a1 ,b=y.trt0.weib$par[2]) 
a30 <- alpha(x=xsimalph0[i,], a =a2 ,b=y.trt0.weib$par[2]) 
a40 <- alpha(x=xsimalph0[i,], a =a3 ,b=y.trt0.weib$par[2]) 
a40

#(2) Calculate MLEs

mle.simalph0=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simalph0[i]=alpha(x=xsimalph0[i,], a =a40 ,b=y.trt0.weib$par[2])}

#(3) Confidence interval
clA0=quantile(mle.simalph0,0.025)
cuA0=quantile(mle.simalph0,0.975)
clA0
cuA0


#trt1


xsimalph1=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimalph1[i,]=rweibull(nn,y.trt1.weib$par[1],1/(y.trt1.weib$par[2])) 


#Newton raphson
a11 <- alpha(x=xsimalph1[i,], a =2 ,b=y.trt1.weib$par[2]) 
a21 <- alpha(x=xsimalph1[i,], a =a1 ,b=y.trt1.weib$par[2]) 
a31 <- alpha(x=xsimalph1[i,], a =a2 ,b=y.trt1.weib$par[2]) 
a41 <- alpha(x=xsimalph1[i,], a =a3 ,b=y.trt1.weib$par[2]) 
a41

#(2) Calculate MLEs

mle.simalph1=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simalph1[i]=alpha(x=xsimalph1[i,], a =a41 ,b=y.trt1.weib$par[2])}

#(3) Confidence interval
clA1=quantile(mle.simalph1,0.025)
cuA1=quantile(mle.simalph1,0.975)
clA1
cuA1


#all
a1 <- alpha(x=xsimalph[i,], a =2 ,b=y.weib$par[2]) 
a2 <- alpha(x=xsimalph[i,], a =a1 ,b=y.weib$par[2]) 
a3 <- alpha(x=xsimalph[i,], a =a2 ,b=y.weib$par[2]) 
a4 <- alpha(x=xsimalph[i,], a =a3 ,b=y.weib$par[2]) 
a4

xsimalph=matrix(0,MM,nn) ##matrix to store data sets
for(i in 1:MM) xsimalph[i,]=rweibull(nn,y.weib$par[1],1/(y.weib$par[2])) 

#(2) Calculate MLEs

mle.simalph=rep(0,MM) ##vector to store MLEs
for(i in 1:MM){mle.simalph[i]=alpha(x=xsimalph[i,], a =a4 ,b=y.weib$par[2])}

#(3) Confidence interval
clA=quantile(mle.simalph,0.025)
cuA=quantile(mle.simalph,0.975)

clA

cuA









