
#### Function for method M3 when there are one Z, one confounding variable, and one 
#### instrumental variable.  
fm3=function(matcdata=mydata, capm=capm){
n=nrow(matcdata)/(capm+1)

myv=rep(rnorm(n), each=(capm+1))
myst= rep(1:n, each=(capm+1))
#outtr=coxph(Surv(myv, y)~tr+z+strata(myst), matcdata)
#store1=rbind(store1, outtr$coef)
##
outw=coxph(Surv(myv, y)~w+z+strata(myst), matcdata)
##
## Our proposed method 
###### estimation of the first likelihood
#
out300=glm(w~tstar+sv+z, family=binomial, data=matcdata)

 nlglik1=function(par){
  # prpi=pr(T=1|T^*, Z) # a logistic model is assumed 
  prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z))
## alpha.0=pr(W=1|T=0, Z)     
   alpha.0=1/(1+exp(-par[5])+exp(par[6]-par[5]) )
##   alpha.1=pr(W=0|T=1, Z)
   alpha.1=1/(1+exp(-par[6])+exp(par[5]-par[6]) )

#In the event that we get an extremely high or extremely low
#probability of treatment, set it to a specific value
prt1[prt1>0.99999]=0.99999
prt1[prt1<0.00001]=0.00001

# prw.1=pr(W=1|T^*, Z)
prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
prw1[prw1>0.99999]=0.99999
prw1[prw1<0.00001]=0.00001
#
#########################

neglk=- sum ((1-matcdata$y)*(log(prw1)*matcdata$w+(1-matcdata$w)*log(1-prw1)))
return(neglk)}
#We use the optim function to estimate the optimal values for the parameters from our first likelihood
out2=optim(c(out300$coef, 0, 0), nlglik1, method="L-BFGS-B", control=list(maxit=1500))#, hessian=T)
#out2=ucminf(c(out300$coef, 0, 0), nlglik1,  control=list(maxit=1500))#, hessian=T)

#
par=out2$par
alpha.0=1/(1+exp(-par[5])+exp(par[6]-par[5]) )# alpha.0=pr(W=1|T=0, Z)
alpha.1=1/(1+exp(-par[6])+exp(par[5]-par[6]) )# alpha.1=pr(W=0|T=1, Z)
prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z)) #=pr(T=1|T^*, SV)
prw1=alpha.0+(1-alpha.0-alpha.1)*prt1
###
##  Estimation of the beta parameters
###
##We use the nlglik2 function to create the log likelihood for the beta parameters, then use 
##optim function to determine the optimal parameters for this function

 nlglik2=function(beta){
exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)

exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0
################ conditional probability calculation 
term1=exp.g*exp(beta[2]*matcdata$z) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3
#
ee=-sum(matcdata$y*log(cond.prob))
return(ee)
}

out3=optim(as.numeric(outw$coef), nlglik2, method="L-BFGS-B", control=list(maxit=1500))
#out3=ucminf(as.numeric(outw$coef), nlglik2, hessian=3)
#  
beta=out3$par

############################# Standard error calculation 
##Next we calculate the standard errors of all of the parameters
##First we need to calculate the "middle term" of the asymptotic variance
##It is composed of the U equations
##We need the output from the first analysis

allpara=c(par, beta)

par=allpara[1:6]
beta= allpara[7:8]
#
alpha.0=1/(1+exp(-par[5])+exp(par[6]-par[5]) )# alpha.0=pr(W=1|T=0, Z)
alpha.1=1/(1+exp(-par[6])+exp(par[5]-par[6]) )# alpha.1=pr(W=0|T=1, Z)
prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z)) #=pr(T=1|T^*, Z)
prw1=alpha.0+(1-alpha.0-alpha.1)*prt1

#
exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)

#der.g.w1=deriv of log(exp.g.w1)= (exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0(1-prt1))
deriv.beta1.g.w1=(exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))
deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))
deriv.beta1.g=deriv.beta1.g.w1*matcdata$w+(1-matcdata$w)*deriv.beta1.g.w0
exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0

################ conditional probability calculation 
term1=exp.g*exp(beta[2]*matcdata$z) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3
#
##Here the estimates for the U equations are estimated
uvec=rep(0, 8);
uvec=matrix(NA,nrow=n,ncol=8,byrow=T)

term11=(
 (1-matcdata$y)*(1-alpha.0-alpha.1)*prt1*
(1-prt1)*(matcdata$w-prw1)/(prw1*(1-prw1))
)

term1.mat=matrix(term11, ncol=(capm+1), byrow=T)

uvec[,1]=apply(term1.mat,1,sum)

uvec21=term11*matcdata$tstar
uvec21.mat=matrix(uvec21,ncol=(capm+1),byrow=T)
uvec[,2]=apply(uvec21.mat,1,sum)
uvec31=term11*matcdata$sv
uvec31.mat=matrix(uvec31,ncol=(capm+1),byrow=T)
uvec[,3]=apply(uvec31.mat,1,sum)
uvec41=term11*matcdata$z
uvec41.mat=matrix(uvec41,ncol=(capm+1),byrow=T)
uvec[,4]=apply(uvec41.mat,1,sum)

term2= (1-matcdata$y)*(matcdata$w-prw1)/(prw1*(1-prw1))

uvec51=term2*( alpha.0*(1-alpha.0)-alpha.0*(1-alpha.0-alpha.1)*prt1 )
uvec51.mat=matrix(uvec51,ncol=(capm+1),byrow=T)
uvec[,5]=apply(uvec51.mat,1,sum)
uvec61= term2*( -alpha.0*alpha.1-alpha.1*(1-alpha.0-alpha.1)*prt1) 
uvec61.mat=matrix(uvec61,ncol=(capm+1),byrow=T)
uvec[,6]=apply(uvec61.mat,1,sum)


newterm1=matcdata$y-cond.prob
#am1=matrix(newterm1, ncol=(capm+1), byrow=T)

mult=(matcdata$y-cond.prob)*deriv.beta1.g
mult.mat=matrix(mult,ncol=(capm+1),byrow=T)
uvec[,7]=apply(mult.mat,1,sum)

uvec81=newterm1*matcdata$z
uvec81.mat=matrix(uvec81,ncol=(capm+1),byrow=T)
uvec[,8]=apply(uvec81.mat,1,sum)

##After calculating the uvec, we can now calculate the middle term 
##mid.term

mid.term=t(uvec)%*%uvec/n



f2=function(allpara){

par=allpara[1:6]
beta= allpara[7:8]
#
alpha.0=1/(1+exp(-par[5])+exp(par[6]-par[5]) )# alpha.0=pr(W=1|T=0, Z)
alpha.1=1/(1+exp(-par[6])+exp(par[5]-par[6]) )# alpha.1=pr(W=0|T=1, Z)
prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z)) #=pr(T=1|T^*, Z)
prw1=alpha.0+(1-alpha.0-alpha.1)*prt1

#
exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))/ prw1
exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))/(1-prw1)

#der.g.w1=deriv of log(exp.g.w1)= (exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0(1-prt1))
deriv.beta1.g.w1=(exp(beta[1])*(1-alpha.1)*prt1)/(exp(beta[1])*(1-alpha.1)*prt1 +alpha.0*(1-prt1))
deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-alpha.0)*(1-prt1))
deriv.beta1.g=deriv.beta1.g.w1*matcdata$w+(1-matcdata$w)*deriv.beta1.g.w0
exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0

################ conditional probability calculation 
term1=exp.g*exp(beta[2]*matcdata$z) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3
#
##Here the estimates for the U equations are estimated
uvec=rep(0, 8);
uvec=matrix(NA,nrow=n,ncol=8,byrow=T)

term11=(
 (1-matcdata$y)*(1-alpha.0-alpha.1)*prt1*
(1-prt1)*(matcdata$w-prw1)/(prw1*(1-prw1))
)

term1.mat=matrix(term11, ncol=(capm+1), byrow=T)

uvec[,1]=apply(term1.mat,1,sum)

uvec21=term11*matcdata$tstar
uvec21.mat=matrix(uvec21,ncol=(capm+1),byrow=T)
uvec[,2]=apply(uvec21.mat,1,sum)
uvec31=term11*matcdata$sv
uvec31.mat=matrix(uvec31,ncol=(capm+1),byrow=T)
uvec[,3]=apply(uvec31.mat,1,sum)
uvec41=term11*matcdata$z
uvec41.mat=matrix(uvec41,ncol=(capm+1),byrow=T)
uvec[,4]=apply(uvec41.mat,1,sum)

term2= (1-matcdata$y)*(matcdata$w-prw1)/(prw1*(1-prw1))

uvec51=term2*( alpha.0*(1-alpha.0)-alpha.0*(1-alpha.0-alpha.1)*prt1 )
uvec51.mat=matrix(uvec51,ncol=(capm+1),byrow=T)
uvec[,5]=apply(uvec51.mat,1,sum)
uvec61= term2*( -alpha.0*alpha.1-alpha.1*(1-alpha.0-alpha.1)*prt1) 
uvec61.mat=matrix(uvec61,ncol=(capm+1),byrow=T)
uvec[,6]=apply(uvec61.mat,1,sum)


newterm1=matcdata$y-cond.prob
#am1=matrix(newterm1, ncol=(capm+1), byrow=T)

mult=(matcdata$y-cond.prob)*deriv.beta1.g
mult.mat=matrix(mult,ncol=(capm+1),byrow=T)
uvec[,7]=apply(mult.mat,1,sum)

uvec81=newterm1*matcdata$z
uvec81.mat=matrix(uvec81,ncol=(capm+1),byrow=T)
uvec[,8]=apply(uvec81.mat,1,sum)

return(apply(uvec, 2, sum))
}
#
#
#
deriv.inv=ginv(jacobian(f2, allpara)/n)
var.cov=deriv.inv%*%mid.term%*%t(deriv.inv)/n
beta.sd=sqrt(diag(var.cov)[7:8])
#
########################################################################### 
m2=data.frame(as.numeric(summary(outw)$coeff[, 1]), as.numeric(summary(outw)$coeff[, 3]))
m3=data.frame(as.numeric(beta), as.numeric(beta.sd))
names(m2)<-c("coef", "se")
names(m3)<-c("coef", "se")

result=list(m2, m3)
names(result)<-c("M2", "M3")
return(result)}

#
#
#
#
#
#
#### Function for method M4 when there are one Z, one confounding variable, and one 
#### instrumental variable.  

fm4=function(matcdata=mydata, capm=capm){
n=nrow(matcdata)/(capm+1)

myv=rep(rnorm(n), each=(capm+1))
myst= rep(1:n, each=(capm+1))
#outtr=coxph(Surv(myv, y)~tr+z+strata(myst), matcdata)
#store1=rbind(store1, outtr$coef)
##
outw=coxph(Surv(myv, y)~w+z+strata(myst), matcdata)
##
## Our proposed method 
###### estimation of the first likelihood
#
out300=glm(w~tstar+sv+z, family=binomial, data=matcdata)

 nlglik=function(par){
  # prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1.0=1/(1+exp(-par[1]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z))
  prt1.1=1/(1+exp(-par[1]-par[7]-par[2]*matcdata$tstar-par[3]*matcdata$sv-par[4]*matcdata$z))
## alpha.0=pr(W=1|X=0, Z)     
   alpha.0=1/(1+exp(-par[5])+exp(par[6]-par[5]) )
##   alpha.1=pr(W=0|X=1, Z)
   alpha.1=1/(1+exp(-par[6])+exp(par[5]-par[6]) )


#In the event that we get an extremely high or extremely low
#probability of treatment, set it to a specific value
prt1.0[prt1.0>0.99999]=0.99999
prt1.0[prt1.0<0.00001]=0.00001
prt1.1[prt1.1>0.99999]=0.99999
prt1.1[prt1.1<0.00001]=0.00001

##
# prw1.0=pr(W=1|T^*, Z,Y=0)
# prw1.1=pr(W=1|T^*, Z,Y=1)
prw1.0=alpha.0+(1-alpha.0-alpha.1)*prt1.0
prw1.1=alpha.0+(1-alpha.0-alpha.1)*prt1.1

prw1.0[prw1.0>0.99999]=0.99999
prw1.0[prw1.0<0.00001]=0.00001
prw1.1[prw1.1>0.99999]=0.99999
prw1.1[prw1.1<0.00001]=0.00001

########################
#This section is used for the beta formulas
#g2=log(1-prt1.0+exp(par[7])*prt1.0)

exp.g2=(1-prt1.0+exp(par[7])*prt1.0)


################ conditional probability calculation 
term1=exp.g2*exp(par[8]*matcdata$z) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3

cond.prob[cond.prob<0.00001]=0.00001
#########################

neglk=- sum (matcdata$y*log(cond.prob)+
               (matcdata$y)*(log(prw1.1)*matcdata$w+(1-matcdata$w)*log(1-prw1.1))+
  (1-matcdata$y)*(log(prw1.0)*matcdata$w+(1-matcdata$w)*log(1-prw1.0)))

return(neglk)}

out2=optim(c(out300$coef, 0, 0, outw$coef), nlglik, method="L-BFGS-B", control=list(maxit=1500), hessian=T)

beta=out2$par[7:8]

hmat=out2$hessian
var.cov=solve(hmat)
beta.sd=sqrt(diag(var.cov)[7:8])

##
########################################################################### 
m2=data.frame(as.numeric(summary(outw)$coeff[, 1]), as.numeric(summary(outw)$coeff[, 3]))
m3=data.frame(as.numeric(beta), as.numeric(beta.sd))
names(m2)<-c("coef", "se")
names(m3)<-c("coef", "se")

result=list(m2, m3)
names(result)<-c("M2", "M4")
return(result)}






 


#Install packages necessary for analysis
require(MASS)
require(survival)
require(numDeriv)


##### Simulation of a population 
capn=80000
capm=2;
n=1000
######
##Simulate variables used in the analysis
sv=runif(capn, -1, 1)#-0.66 # stratification variable 
z=rnorm(capn, 0, 0.5)
# simulation of the instrumental variable tstar
#
#
myep=runif(capn, 0, 1)
#tstar=as.numeric(myep<0.3)
tstar=rnorm(capn, 0, 0.5)
# simulation of the actual treatment/exposure variable 
#True treatment variable tr
#Set probability for exposure for weak association
probx=1/(1+exp(-(-1+2*tstar+1*(sv)))) 
tr=rbinom(capn, 1, probx)
##w is the mismeasured form of tr
#Misclassification alpha1=alpha0=0.2
capbw=rbinom(capn, 1, 0.8) #capbw=pr(W=1|X=0) =1-alpha1
capbws=rbinom(capn, 1, 0.8)#capbws=pr(W=0|X=1)=1-alpha0
w=capbw*tr+(1-capbws)*(1-tr)
##simulation of y 
proby=1/(1+exp(-(-2-2*(sv)+1*tr+0.5*z)))
y=rbinom(capn,1,proby)
##Combine data for analysis
mydata=cbind(sv, tr, tstar, w, y, z)
############# Forming matched data ####
##First randomly select the id of cases from the data
index=rep(0, ((capm+1)*n))
idc=which(y==1)
samplec=sample(idc, n, replace=F)
##Next, for each case, find a matching control by comparing difference 
##of stratification variable for case and control, using tolerances of 
##0.01, 0.05, and 0.1

for( i1 in 1: n){
k1=(capm+1)*(i1-1)+1
ids=which(abs(sv-sv[samplec[i1]])<0.01 & y==0)
if(length(ids)>capm) {temp1=sample(ids, capm, replace=F)} else {
 ids=which(abs(sv-sv[samplec[i1]])<0.05 & y==0) 
 if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
ids=which(abs(sv-sv[samplec[i1]])<0.1 & y==0) 
 temp1=sample(ids, capm, replace=F)
}
} 
index[k1]=samplec[i1]
index[k1+(1:capm)]=temp1
} 
########################## the columns are sv, tr, tstar, w, y, z
##The subset of data that was selected in previous set will be matched data
matcdata=mydata[index, ] 
matcdata=data.frame(matcdata) 

### End of data simulation
##########################

out3000=fm3(matcdata=matcdata, capm=capm)


out4000=fm4(matcdata=matcdata, capm=capm)


#Using this function to speed up simulation
#library(doParallel) 
#no_cores <- detectCores() - 1 
#registerDoParallel(cores=no_cores) 
#cl <- makeCluster(no_cores, type="FORK") 
#store <- foreach(it = 1:5000, .combine="rbind") %dopar% myfn1(nt = it)
