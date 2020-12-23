library(glmnet)
library(survival)

rm(list=ls(all=TRUE))
dat1=read.csv("dat1.csv",header=TRUE)
dat2=read.csv("dat2.csv",header=TRUE)
dat3=read.csv("Resec3.csv",header=TRUE)
dat4=read.csv("Resec4.csv",header=TRUE)
dat5=read.csv("TCIA.csv",header=TRUE)
dat6=read.csv("IMU1.csv",header=TRUE)
dat7=read.csv("IMU2.csv",header=TRUE)
dat8=read.csv("IMU3.csv",header=TRUE)
dat9=read.csv("Chem1.csv",header=TRUE)
dat10=read.csv("Chem2.csv",header=TRUE)
dat11=read.csv("TKI1.csv",header=TRUE)
dat12=read.csv("TKI2.csv",header=TRUE)


################  iEDI construction
set.seed(101)
y=dat1$group_im 
y=ifelse(y==1,2,1)
x=as.matrix(dat[109:153])
fit_cv=cv.glmnet(x,y,family="binomial",nfolds=10)
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
#selected variables
LassoM.coef=coef(model.final,s=fit_cv$lambda.min)
VV=names(LassoM.coef[as.vector(LassoM.coef[,1]!=0),])
VV
 
#####Resec1 cohort
 attach(dat1)
dat1$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat1$PI_5=dat1$PI_5+ve
}

#####Resec2 cohort
attach(dat2)
dat2$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat2$PI_5=dat2$PI_5+ve
}

#####Resec3 cohort
attach(dat3)
dat3$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat3$PI_5=dat3$PI_5+ve
}

#####Resec4 cohort
attach(dat4)
dat4$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat4$PI_5=dat4$PI_5+ve
}

#####TCIA cohort
attach(dat5)
dat5$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat5$PI_5=dat5$PI_5+ve
}

#####IMU1 cohort
attach(dat6)
dat6$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat6$PI_5=dat6$PI_5+ve
}

#####IMU2 cohort
attach(dat7)
dat7$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat7$PI_5=dat7$PI_5+ve
}

#####IMU3 cohort
attach(dat8)
dat8$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat8$PI_5=dat8$PI_5+ve
}

#####Chem1 cohort
attach(dat9)
dat9$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat9$PI_5=dat9$PI_5+ve
}

#####Chem2 cohort
attach(dat10)
dat10$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat10$PI_5=dat10$PI_5+ve
}

#####TKI1 cohort
attach(dat11)
dat11$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat11$PI_5=dat11$PI_5+ve
}

#####TKI2 cohort
attach(dat12)
dat12$PI_5=cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat12$PI_5=dat12$PI_5+ve
}

#################################  KM curve  ###############################
#########  Resect cohorts

#####Resec1 cohort
library(survivalROC)
cutoff=36
dat1$OS_statue=ifelse(dat1$OS_statue==1,1,0)
roc=survivalROC(Stime=dat1$OS_time,status=dat1$OS_statue,marker=dat1$PI_5,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]
cut.op

dat1$group_5=ifelse(dat1$PI_5<=cut.op,2,1)
s<-survdiff(Surv(dat1$OS_time,dat1$OS_statue)~group_5,data=dat1,rho=1)

km<-survfit(Surv(dat1$OS_time,dat1$OS_statue)~group_5,data=dat1,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec1 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(80,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####Resec2 cohort
dat2$group_5=ifelse(dat2$PI_5<=-0.2905835,2,1)
s<-survdiff(Surv(dat2$OS_time,dat2$OS_statue)~group_5,data=dat2,rho=1)

km<-survfit(Surv(dat2$OS_time,dat2$OS_statue)~group_5,data=dat2,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec2 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(58,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####Resec3 cohort
dat3$group_5=ifelse(dat3$PI_5<=-0.2905835,2,1)
s<-survdiff(Surv(dat3$OS_time,dat3$OS_statue)~group_5,data=dat3,rho=1)

km<-survfit(Surv(dat3$OS_time,dat3$OS_statue)~group_5,data=dat3,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec3 cohort",col=c('deepskyblue3','firebrick1'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(82,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####Resec4 cohort
dat4$group_5=ifelse(dat4$PI_5<=-0.2905835,2,1)
s<-survdiff(Surv(dat4$OS_time,dat4$OS_statue)~group_5,data=dat4,rho=1)

km<-survfit(Surv(dat4$OS_time,dat4$OS_statue)~group_5,data=dat4,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Resec4 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(58,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####TCIA cohort
dat5$group_5=ifelse(dat5$PI_5<=-0.2905835,2,1)
s<-survdiff(Surv(dat5$OS_time,dat5$OS_statue)~group_5,data=dat5,rho=1)

km<-survfit(Surv(dat5$OS_time,dat5$OS_statue)~group_5,data=dat5,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="TCIA cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(58,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#########  Immunotherapy cohorts
#####IMU1 cohort
library(survivalROC)
cutoff=10
dat6$status_PFS=ifelse(dat6$status_PFS==1,1,0)
roc=survivalROC(Stime=dat6$PFS,status=dat6$status_PFS,marker=dat6$PI_5,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]
cut.op

dat6$group_5=ifelse(dat6$PI_5<=cut.op,2,1)
s<-survdiff(Surv(dat6$PFS,dat6$status_PFS)~group_5,data=dat6,rho=1)

km<-survfit(Surv(dat6$PFS,dat6$status_PFS)~group_5,data=dat6,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-Free Survival",xlab="Time Since Surgery (months)",main="IMU1 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(20,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####IMU2 cohort
dat7$group_5=ifelse(dat7$PI_5<=-0.5178462,2,1)
s<-survdiff(Surv(dat7$PFS,dat7$status_PFS)~group_5,data=dat7,rho=1)

km<-survfit(Surv(dat7$PFS,dat7$status_PFS)~group_5,data=dat7,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-Free Survival",xlab="Time Since Surgery (months)",main="IMU2 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(20,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####IMU3 cohort
dat8$group_5=ifelse(dat8$PI_5<=-0.5178462,2,1)
s<-survdiff(Surv(dat8$PFS,dat8$status_PFS)~group_5,data=dat8,rho=1)

km<-survfit(Surv(dat8$PFS,dat8$status_PFS)~group_5,data=dat8,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-Free Survival",xlab="Time Since Surgery (months)",main="IMU3 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(20,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#########  Chemotherapy cohorts
#####Chem1 cohort
library(survivalROC)
cutoff=10
dat9$status_PFS=ifelse(dat9$status_PFS==1,1,0)
roc=survivalROC(Stime=dat9$PFS,status=dat9$status_PFS,marker=dat9$PI_5,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]
cut.op

dat9$group_5=ifelse(dat9$PI_5<=-0.5237359,2,1)
s<-survdiff(Surv(dat9$PFS,dat9$status_PFS)~group_5,data=dat9,rho=1)

km<-survfit(Surv(dat9$PFS,dat9$status_PFS)~group_5,data=dat9,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-Free Survival",xlab="Time Since Surgery (months)",main="Chem1 cohort",col=c('deepskyblue3','firebrick1'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(36,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####Chem2 cohort
dat10$group_5=ifelse(dat10$PI_5<=-0.5237359,2,1)
s<-survdiff(Surv(dat10$PFS,dat10$status_PFS)~group_5,data=dat10,rho=1)

km<-survfit(Surv(dat10$PFS,dat10$status_PFS)~group_5,data=dat10,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-Free Survival",xlab="Time Since Surgery (months)",main="Chem2 cohort",col=c('deepskyblue3','firebrick1'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96),labels=c(0,12,24,36,48,60,72,84,96),cex.axis=1.25)
legend(31,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#########  EGFR-TKI therapy cohorts
#####TKI1 cohort
library(survivalROC)
cutoff=10
dat11$status_PFS=ifelse(dat11$status_PFS==1,1,0)
roc=survivalROC(Stime=dat11$PFS,status=dat11$status_PFS,marker=dat11$PI_5,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]
cut.op

dat11$group_5=ifelse(dat11$PI_5<=-0.057,2,1)
s=survdiff(Surv(dat11$PFS,dat11$status_PFS)~group_5,data=dat11,rho=0)

km<-survfit(Surv(dat11$PFS,dat11$status_PFS)~group_5,data=dat11,conf.int=TRUE)
par(mar=c(9,6,3,2))
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progression-free Survival",xlab="Time Since Surgery (months)",main="TKI1 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(93,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)

#####TKI2 cohort
dat12$group_5=ifelse(dat12$PI_5<=-0.057,2,1)
s<-survdiff(Surv(dat12$PFS,dat12$status_PFS)~group_5,data=dat12,rho=0)

km<-survfit(Surv(dat12$PFS,dat12$status_PFS)~group_5,data=dat12,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Progrssion-free Survival",xlab="Time Since Surgery (months)",main="TKI2 cohort",col=c('firebrick1','deepskyblue3'),cex.lab=1.25,cex.main=1.25,yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1.25)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108),labels=c(0,12,24,36,48,60,72,84,96,108),cex.axis=1.25)
legend(50,1.05,c("High iEDI-score","Low iEDI-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','firebrick1'),box.lty=0)



