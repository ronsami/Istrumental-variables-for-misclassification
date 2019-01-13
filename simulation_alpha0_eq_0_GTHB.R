#### The following function computes the estimates of method 3
#### along with its standard error. This function also returns estimates
#### for the naive method M2.

fm3_alpha0_eq_0=function(matcdata=mydata, capm=capm){

n=nrow(matcdata)/(capm+1)

myv=rep(rnorm(n), each=(capm+1))
myst= rep(1:n, each=(capm+1))


outw=coxph(Surv(myv, y)~w+myz+strata(myst),data=matcdata)

out300=glm(w~tstar1+tstar2+tstar3+tstar4+sv1_1+sv2_2+sv2_3+sv4_2+sv4_3+sv4_4+sv4_5+myz, family=binomial, data=matcdata)

nlglik1=function(par){
  # prpi=pr(T=1|T^*, Z) # a logistic model is assumed 
  prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
                -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2
                -par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
                -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4
                -par[12]*matcdata$sv4_5-par[13]*matcdata$myz))
  
  # prw.1=pr(W=1|T^*, Z)
  prt1[prt1>0.99999]=0.99999
  prt1[prt1<0.00001]=0.00001
  
  alpha.1=1/(1+exp(-par[14]))
  prw1=(1-alpha.1)*prt1
 
  prw1[prw1>0.99999]=0.99999
  prw1[prw1<0.00001]=0.00001
  ########################
  neglk=- sum ((1-matcdata$y)*(log(prw1)*matcdata$w+(1-matcdata$w)*log(1-prw1))) 
  return(neglk)}
#out2=optim(c(out300$coef, 0), nlglik1, method="L-BFGS-B", control=list(maxit=1500),hessian=T)

out2=ucminf(c(out300$coef, 0), nlglik1)#, method="L-BFGS-B", control=list(maxit=1500),hessian=T)


par=out2$par
alpha.1=1/(1+exp(-par[14]))# alpha.1=pr(W=0|T=1, Z)
prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
              -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2-par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
              -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4-par[12]*matcdata$sv4_5-par[13]*matcdata$myz))
prw1=(1-alpha.1)*prt1
###
###
###
nlglik2=function(beta){
  exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 )/ prw1
  exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-prt1))/(1-prw1)
  exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0
################ conditional probability calculation 
term1=exp.g*exp(beta[2]*matcdata$myz) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3
   ee=-sum(matcdata$y*log(cond.prob))
  return(ee)
  }
out3=optim(outw$coef, nlglik2, method="L-BFGS-B", control=list(maxit=1500))
beta=out3$par
########### Asymptotic standard error calculation 
allpara=c(par, beta)
#
par=allpara[1:14]
beta= allpara[15:16]
#
# alpha.0=pr(W=1|T=0, Z)
alpha.1=1/(1+exp(-par[14]))# alpha.1=pr(W=0|T=1, Z)
#alpha.1old=0.49/(1+exp(-par[13]))
prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
              -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2-par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
              -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4-par[12]*matcdata$sv4_5-par[13]*matcdata$myz))
prw1=(1-alpha.1)*prt1
#
exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 )/ prw1
exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-prt1))/(1-prw1)
exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0

#
deriv.beta1.g.w1=1#der.g.w1=derviv of log(exp.g.w1)
deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-prt1))
#der.g.w0=derviv of log(exp.g.w0)
deriv.beta1.g=deriv.beta1.g.w1*matcdata$w+(1-matcdata$w)*deriv.beta1.g.w0

################ conditional probability calculation 
term1=exp.g*exp(beta[2]*matcdata$myz) 
term2=matrix(term1, ncol=(capm+1), byrow=T)
term3=rep(apply(term2, 1, sum), each=(capm+1))
cond.prob=term1/term3
#
#
uvec=matrix(NA, nrow=n, ncol=16, byrow=T)
term11=(
    (1-matcdata$y)*(1-alpha.1)*prt1*
      (1-prt1)*(matcdata$w-prw1)/(prw1*(1-prw1))
  )
term1.mat=matrix(term11, ncol=(capm+1), byrow=T)  
uvec[, 1]=apply(term1.mat, 1, sum)
#
uvec21=term11*matcdata$tstar1
uvec21.mat=matrix(uvec21, ncol=(capm+1), byrow=T)  
uvec[, 2]=apply(uvec21.mat, 1, sum)
#
uvec31=term11*matcdata$tstar2
uvec31.mat=matrix(uvec31, ncol=(capm+1), byrow=T)  
uvec[, 3]=apply(uvec31.mat, 1, sum)
#
uvec41=term11*matcdata$tstar3
uvec41.mat=matrix(uvec41, ncol=(capm+1), byrow=T)  
uvec[, 4]=apply(uvec41.mat, 1, sum)
#
uvec51=term11*matcdata$tstar4
uvec51.mat=matrix(uvec51, ncol=(capm+1), byrow=T)  
uvec[, 5]=apply(uvec51.mat, 1, sum)
#
uvec61=term11*matcdata$sv1_1
uvec61.mat=matrix(uvec61, ncol=(capm+1), byrow=T)  
uvec[, 6]=apply(uvec61.mat, 1, sum)
#
uvec71=term11*matcdata$sv2_2
uvec71.mat=matrix(uvec71, ncol=(capm+1), byrow=T)  
uvec[, 7]=apply(uvec71.mat, 1, sum)
#
uvec81=term11*matcdata$sv2_3
uvec81.mat=matrix(uvec81, ncol=(capm+1), byrow=T)  
uvec[, 8]=apply(uvec81.mat, 1, sum)
#
uvec91=term11*matcdata$sv4_2
uvec91.mat=matrix(uvec91, ncol=(capm+1), byrow=T)  
uvec[, 9]=apply(uvec91.mat, 1, sum)
#
uvec101=term11*matcdata$sv4_3
uvec101.mat=matrix(uvec101, ncol=(capm+1), byrow=T)  
uvec[, 10]=apply(uvec101.mat, 1, sum)
#
uvec111=term11*matcdata$sv4_4
uvec111.mat=matrix(uvec111, ncol=(capm+1), byrow=T)  
uvec[, 11]=apply(uvec111.mat, 1, sum)
#
uvec121=term11*matcdata$sv4_5
uvec121.mat=matrix(uvec121, ncol=(capm+1), byrow=T)  
uvec[, 12]=apply(uvec121.mat, 1, sum)
#
uvec131=term11*matcdata$myz
uvec131.mat=matrix(uvec131, ncol=(capm+1), byrow=T)  
uvec[, 13]=apply(uvec131.mat, 1, sum)
#
term2= (1-matcdata$y)*(matcdata$w-prw1)/(prw1*(1-prw1))
#
uvec141=term2*( -alpha.1*(1-alpha.1)*prt1) 
uvec141.mat=matrix(uvec141, ncol=(capm+1), byrow=T)  
uvec[, 14]=apply(uvec141.mat, 1, sum)
#
uvec151=(matcdata$y-cond.prob)*deriv.beta1.g
uvec151.mat=matrix(uvec151, ncol=(capm+1), byrow=T)
uvec[, 15]=apply(uvec151.mat, 1, sum)
#
uvec161=(matcdata$y-cond.prob)*matcdata$myz
uvec161.mat=matrix(uvec161, ncol=(capm+1), byrow=T)
uvec[, 16]=apply(uvec161.mat, 1, sum)

mid.term=t(uvec)%*%uvec
mid.term=mid.term/n


##########################################################
###
f2=function(allpara){
  par=allpara[1:14]
  beta= allpara[15:16]
  #
  
  # alpha.0=pr(W=1|T=0, Z)
  alpha.1=1/(1+exp(-par[14]))# alpha.1=pr(W=0|T=1, Z)
  prt1=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
                -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2-par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
                -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4-par[12]*matcdata$sv4_5-par[13]*matcdata$myz))
  prw1=(1-alpha.1)*prt1
  #
  exp.g.w1= (exp(beta[1])*(1-alpha.1)*prt1 )/ prw1
  exp.g.w0= (exp(beta[1])*alpha.1*prt1+(1-prt1))/(1-prw1)
  exp.g=exp.g.w1*matcdata$w+(1-matcdata$w)*exp.g.w0
  
  ####
  deriv.beta1.g.w1=1
  deriv.beta1.g.w0=(exp(beta[1])*alpha.1*prt1)/(exp(beta[1])*alpha.1*prt1+(1-prt1))
  deriv.beta1.g=deriv.beta1.g.w1*matcdata$w+(1-matcdata$w)*deriv.beta1.g.w0
  
  ################ conditional probability calculation 
  term1=exp.g*exp(beta[2]*matcdata$myz) 
  term2=matrix(term1, ncol=(capm+1), byrow=T)
  term3=rep(apply(term2, 1, sum), each=(capm+1))
  cond.prob=term1/term3
  #
  #
  uvec=matrix(NA, nrow=n, ncol=16, byrow=T)
  term11=(
    (1-matcdata$y)*(1-alpha.1)*prt1*
      (1-prt1)*(matcdata$w-prw1)/(prw1*(1-prw1))
  )
  term1.mat=matrix(term11, ncol=(capm+1), byrow=T)  
  uvec[, 1]=apply(term1.mat, 1, sum)
  #
  uvec21=term11*matcdata$tstar1
  uvec21.mat=matrix(uvec21, ncol=(capm+1), byrow=T)  
  uvec[, 2]=apply(uvec21.mat, 1, sum)
  #
  uvec31=term11*matcdata$tstar2
  uvec31.mat=matrix(uvec31, ncol=(capm+1), byrow=T)  
  uvec[, 3]=apply(uvec31.mat, 1, sum)
  #
  uvec41=term11*matcdata$tstar3
  uvec41.mat=matrix(uvec41, ncol=(capm+1), byrow=T)  
  uvec[, 4]=apply(uvec41.mat, 1, sum)
  #
  uvec51=term11*matcdata$tstar4
  uvec51.mat=matrix(uvec51, ncol=(capm+1), byrow=T)  
  uvec[, 5]=apply(uvec51.mat, 1, sum)
  #
  uvec61=term11*matcdata$sv1_1
  uvec61.mat=matrix(uvec61, ncol=(capm+1), byrow=T)  
  uvec[, 6]=apply(uvec61.mat, 1, sum)
  #
  uvec71=term11*matcdata$sv2_2
  uvec71.mat=matrix(uvec71, ncol=(capm+1), byrow=T)  
  uvec[, 7]=apply(uvec71.mat, 1, sum)
  #
  uvec81=term11*matcdata$sv2_3
  uvec81.mat=matrix(uvec81, ncol=(capm+1), byrow=T)  
  uvec[, 8]=apply(uvec81.mat, 1, sum)
  #
  uvec91=term11*matcdata$sv4_2
  uvec91.mat=matrix(uvec91, ncol=(capm+1), byrow=T)  
  uvec[, 9]=apply(uvec91.mat, 1, sum)
  #
  uvec101=term11*matcdata$sv4_3
  uvec101.mat=matrix(uvec101, ncol=(capm+1), byrow=T)  
  uvec[, 10]=apply(uvec101.mat, 1, sum)
  #
  uvec111=term11*matcdata$sv4_4
  uvec111.mat=matrix(uvec111, ncol=(capm+1), byrow=T)  
  uvec[, 11]=apply(uvec111.mat, 1, sum)
  #
  uvec121=term11*matcdata$sv4_5
  uvec121.mat=matrix(uvec121, ncol=(capm+1), byrow=T)  
  uvec[, 12]=apply(uvec121.mat, 1, sum)
  #
  uvec131=term11*matcdata$myz
  uvec131.mat=matrix(uvec131, ncol=(capm+1), byrow=T)  
  uvec[, 13]=apply(uvec131.mat, 1, sum)
  #
  term2= (1-matcdata$y)*(matcdata$w-prw1)/(prw1*(1-prw1))
  #
  uvec141=term2*( -alpha.1*(1-alpha.1)*prt1) 
  uvec141.mat=matrix(uvec141, ncol=(capm+1), byrow=T)  
  uvec[, 14]=apply(uvec141.mat, 1, sum)
  #
  uvec151=(matcdata$y-cond.prob)*deriv.beta1.g
  uvec151.mat=matrix(uvec151, ncol=(capm+1), byrow=T)
  uvec[, 15]=apply(uvec151.mat, 1, sum)
  #
  uvec161=(matcdata$y-cond.prob)*matcdata$myz
  uvec161.mat=matrix(uvec161, ncol=(capm+1), byrow=T)
  uvec[, 16]=apply(uvec161.mat, 1, sum)
  
return(apply(uvec, 2, sum))
}
out100=jacobian(f2, allpara)
dmati=ginv(out100/n)
var.cov=dmati%*%mid.term%*%t(dmati)/n
beta.sd=sqrt(diag(var.cov)[15:16])

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
#
#
# 
#### The following function computes the estimates of method 4
#### along with its standard error. This function also returns estimates
#### for the naive method M2.

fm4_alpha0_eq_0=function(matcdata=matcdata, capm=capm){
  
  n=nrow(matcdata)/(capm+1)
  
  myv=rep(rnorm(n), each=(capm+1))
  myst= rep(1:n, each=(capm+1))
  
  
  outw=coxph(Surv(myv, y)~w+myz+strata(myst),data=matcdata)
  
  out300=glm(w~tstar1+tstar2+tstar3+tstar4+sv1_1+sv2_2+sv2_3+sv4_2+sv4_3+sv4_4+sv4_5+myz, family=binomial, data=matcdata)
  
 nlglik1=function(par){
  # prpi=pr(T=1|T^*, Z) # a logistic model is assumed 
  prt1.0=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
                  -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2-par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
                  -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4-par[12]*matcdata$sv4_5-par[13]*matcdata$myz))
  prt1.1=1/(1+exp(-par[1]-par[2]*matcdata$tstar1-par[3]*matcdata$tstar2-par[4]*matcdata$tstar3
                  -par[5]*matcdata$tstar4-par[6]*matcdata$sv1_1-par[7]*matcdata$sv2_2-par[8]*matcdata$sv2_3-par[9]*matcdata$sv4_2
                  -par[10]*matcdata$sv4_3-par[11]*matcdata$sv4_4-par[12]*matcdata$sv4_5-par[13]*matcdata$myz-par[15]))
  
  ## we shall never allow the distribution of W to depend on Z. Conditional on the TR, W is independent of Z
  alpha.1=1/(1+exp(-par[14]))
  
  
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1.0[prt1.0>0.99999]=0.99999
  prt1.0[prt1.0<0.00001]=0.00001
  prt1.1[prt1.1>0.99999]=0.99999
  prt1.1[prt1.1<0.00001]=0.00001
  
  ##
  # prw1.0=pr(W=1|T^*, Z,Y=0)
  # prw1.1=pr(W=1|T^*, Z,Y=1)
  prw1.0=(1-alpha.1)*prt1.0
  prw1.1=(1-alpha.1)*prt1.1
  
  prw1.0[prw1.0>0.99999]=0.99999
  prw1.0[prw1.0<0.00001]=0.00001
  prw1.1[prw1.1>0.99999]=0.99999
  prw1.1[prw1.1<0.00001]=0.00001
  ########################
  exp.g2=(1-prt1.0+exp(par[15])*prt1.0)
  ################ conditional probability calculation 
  term1=exp.g2*exp(par[16]*matcdata$myz) 
  term2=matrix(term1, ncol=(capm+1), byrow=T)
  term3=rep(apply(term2, 1, sum), each=(capm+1))
  cond.prob=term1/term3
  #
 
  neglk=- sum(matcdata$y*log(cond.prob)+
                (matcdata$y)*(log(prw1.1)*matcdata$w+(1-matcdata$w)*log(1-prw1.1))+
                (1-matcdata$y)*(log(prw1.0)*matcdata$w+(1-matcdata$w)*log(1-prw1.0))
              ) 
      return(neglk)}


out2=optim(c(out300$coef, 0, outw$coef), nlglik1, method="L-BFGS-B",
           control=list(maxit=1500),hessian=T)
#out2=ucminf(runif(16, -0.75, 0.75), nlglik1,hessian=3)
beta=out2$par[15:16]

hmat=out2$hessian
var.cov=solve(hmat)
beta.sd=sqrt(diag(var.cov)[15:16])

m2=data.frame(as.numeric(summary(outw)$coeff[, 1]), as.numeric(summary(outw)$coeff[, 3]))
m4=data.frame(as.numeric(beta), as.numeric(beta.sd))
names(m2)<-c("coef", "se")
names(m4)<-c("coef", "se")

result=list(m2, m4)
names(result)<-c("M2", "M4")

return(result)
}

######## Data simulation 

data.par=c( -3.3939,  
            0.9875,
            -0.3194,  0.1382,
            -0.2133,  1.2873,  0.3640,  1.6599,  0.6690,  0.9097,  1.0858, 
            0.2768,-0.1940,  0.6904, -0.1070,  0.3325)

coef.y=c(-3.20565, 0.19619, -0.02983, 0.36578, 0.09228, 0.31728, 0.56554, 0.54982  )
###
###
##### Simulation of a population 
capn=42933
capm=2;
n=1000
#data generation
#Parameters estimated from first likelihood in real data analysis
#Need gamma parameters from actual dataset to conduct simulation.

tstar1=rnorm(capn);
tstar2=rnorm(capn);
tstar3=rnorm(capn);
tstar4=rnorm(capn);
sv1_1=rbinom(capn, 1, 0.5)
sv2=as.numeric(t(rmultinom(capn, 1, c(0.3, 0.4, 0.3)))%*%c(1, 2, 3))
sv2_2=sv2_3=rep(0, capn)
sv2_2[sv2==2]<-1
sv2_3[sv2==3]<-1
sv4=as.numeric(t(rmultinom(capn, 1, c(0.2,.2, 0.2, 0.2, 0.2)))%*%c(1, 2, 3, 4, 5))
sv4_2=sv4_3=sv4_4=sv4_5=rep(0, capn)
sv4_2[sv4==2]<-1
sv4_3[sv4==3]<-1
sv4_4[sv4==4]<-1
sv4_5[sv4==5]<-1
myz=rnorm(capn)
mydata=data.frame(tstar1, tstar2, tstar3, tstar4, sv1_1, sv2_2, sv2_3,
                  sv4_2, sv4_3, sv4_4, sv4_5, myz)

#Generate treatment data
probtr=1/(1+exp(-(data.par[1]+data.par[2]*mydata$tstar1+data.par[3]*mydata$tstar2
                  +data.par[4]*mydata$tstar3+data.par[5]*mydata$tstar4
                  +data.par[6]*mydata$sv1_1+data.par[7]*mydata$sv2_2
                  +data.par[8]*mydata$sv2_3+data.par[9]*mydata$sv4_2
                  +data.par[10]*mydata$sv4_3+data.par[11]*mydata$sv4_4
                  +data.par[12]*mydata$sv4_5))) 
tr=rbinom(capn,1,probtr)
mydata=cbind(mydata,tr)

#Generate W based on misclassification probability alpha.1 estimated from real data
alpha.1=1/(1+exp(-data.par[13]))
capbw=rbinom(capn,1,1-alpha.1) 

#capbws=rbinom(capn, 1, 0.8)#capbws=pr(W=0|X=1)=1-alpha0
w=capbw*tr

#w=capw*tr+(1-capw)*(1-tr)
mydata=cbind(mydata,w)


#Generate y
proby=1/(1+exp(-(coef.y[1]+data.par[14]*mydata$tr+
                   data.par[16]*mydata$myz+
                   coef.y[2]*mydata$sv1_1+coef.y[3]*mydata$sv2_2+
                   coef.y[4]*mydata$sv2_3
                 +coef.y[5]*mydata$sv4_2+coef.y[6]*mydata$sv4_3
                 +coef.y[7]*mydata$sv4_4+coef.y[8]*mydata$sv4_5)))
y=rbinom(capn,1,proby)

mydata=cbind(mydata,y)
index=rep(0, ((capm+1)*n))
idc=which(mydata$y==1)
samplec=sample(idc, n, replace=F)
for( i1 in 1: n){
  #seed 25, i1=374 no matching criteria
  #print(i1)
  k1=(capm+1)*(i1-1)+1
  ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<0.01 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<0.01 
            &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<0.01 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<0.01 &
              abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<0.01&mydata$y==0&
              abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
              abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0)
  if(length(ids)>capm) {temp1=sample(ids, capm, replace=F)} else {
    ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<0.05 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<0.05 
              &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<0.05 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<0.05 &
                abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<0.05 & mydata$y==0&
                abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0) 
    if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
      ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<0.1 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<0.1 
                &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<0.1 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<0.1 &
                  abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<0.1 & mydata$y==0&
                  abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                  abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0) 
      if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
        ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<0.5 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<0.5 
                  &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<0.5 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<0.5 &
                    abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<0.5  & mydata$y==0&
                    abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                    abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0)
        if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
          ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<0.75 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<0.75 
                    &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<0.75 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<0.75 &
                      abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<0.75 & mydata$y==0&
                      abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                      abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0)
          if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
            ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<1 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<1 
                      &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<1 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<1 &
                        abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<1  & mydata$y==0&
                        abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                        abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0)         
            if(length(ids)>capm) temp1=sample(ids, capm, replace=F) else {
              ids=which(abs(mydata$sv1_1-mydata$sv1_1[samplec[i1]])<1.25 &abs(mydata$sv2_2-mydata$sv2_2[samplec[i1]])<1.25 
                        &abs(mydata$sv2_3-mydata$sv2_3[samplec[i1]])<1.25 & abs(mydata$sv4_2-mydata$sv4_2[samplec[i1]])<1.25 &
                          abs(mydata$sv4_3-mydata$sv4_3[samplec[i1]])<1.25  & mydata$y==0&
                          abs(mydata$sv4_4-mydata$sv4_4[samplec[i1]])<0.01&mydata$y==0&
                          abs(mydata$sv4_5-mydata$sv4_5[samplec[i1]])<0.01&mydata$y==0)       
            }
          }
        }
      }
    }
  } 
  index[k1]=samplec[i1]
  index[k1+(1:capm)]=temp1
} 
###

matcdata=mydata[index, ] 
matcdata=data.frame(matcdata) 
# End of data simulation
#### Beginning of analysis 

####### Required packages 
require(MASS)
require(survival)
require(ucminf)
require(numDeriv)
##########################
# For computing the estimates under M3
out3000=fm3_alpha0_eq_0(matcdata=matcdata, capm=2)
# For computing the estimates under M4
out4000=fm4_alpha0_eq_0(matcdata=matcdata, capm=2)


