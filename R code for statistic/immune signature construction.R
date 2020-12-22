library(glmnet)
library(survival)

rm(list=ls(all=TRUE))
dat1=read.csv("Resec1.csv",header=TRUE)
dat2=read.csv("Resec2.csv",header=TRUE)

###############immunomarker
##Ratio
set.seed(10)
y=Surv(dat1$OS_time,dat1$OS_statue) 
x=as.matrix(dat1[2:7])
fit_cv=cv.glmnet(x,y,family="cox",nfolds=10)
plot(fit_cv,cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
#abline(v=log(fit_cv$lambda.1se),col="darkgray",lwd=2)
abline(v=log(fit_cv$lambda.min),col="darkgray",lwd=2)
axis=(3,cex.lab=1.25)

model.final=fit_cv$glmnet.fit
plot(model.final,xvar=c("lambda"),label=TRUE,cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=log(fit_cv$lambda.min),col="darkgray",lwd=2)

#beta value
c=round(coef(model.final,s=fit_cv$lambda.min),8)
cc=summary(c)
#selected variables of immune features
LassoM.coef=coef(model.final,s=fit_cv$lambda.min)
VV=names(LassoM.coef[as.vector(LassoM.coef[,1]!=0),])
VV

#####immune signature construction
attach(dat1)
dat1$PI_im=get(VV[1])*cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat1$PI_im=dat1$PI_im+ve
}

attach(dat2)
dat2$PI_im=get(VV[1])*cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat2$PI_im=dat2$PI_im+ve
}

######C-index for immune signature
library(Hmisc)
r=rcorrcens(Surv(dat1$OS_time,dat1$OS_statue==1,type="right")~(dat1$PI_im))
1-r
library(Hmisc)
r=rcorrcens(Surv(dat2$OS_time,dat2$OS_statue==1,type="right")~(dat2$PI_im))
1-r

########KM curve for immune signature
library(survivalROC)
cutoff=36
dat1$OS_statue=ifelse(dat1$OS_statue==1,1,0)
roc=survivalROC(Stime=dat1$OS_time,status=dat1$OS_statue,marker=dat1$PI_im,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]

dat1$group_im=ifelse(dat1$PI_im<=cut.op,2,1)
s<-survdiff(Surv(dat1$OS_time,dat1$OS_statue)~group_im,data=dat1,rho=1)

dat2$group_im=ifelse(dat2$PI_im<=cut.op,2,1)
s<-survdiff(Surv(dat2$OS_time,dat2$OS_statue)~group_im,data=dat2,rho=1)

km<-survfit(Surv(dat1$OS_time,dat1$OS_statue)~group_im,data=dat1,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec1 dataset",col=c('deepskyblue3','firebrick1'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(80,1.05,c("Immune high-level","Immune low-level"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

km<-survfit(Surv(dat2$OS_time,dat2$OS_statue)~group_im,data=dat2,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec2 dataset",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(80,1.05,c("Immune high-level","Immune low-level"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)


write.csv(dat1,'dat1.csv')
write.csv(dat2,'dat2.csv')

