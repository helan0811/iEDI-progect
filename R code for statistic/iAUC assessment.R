###########################iAUC 
###########Resect cohorts
library(MASS)
library(risksetROC)
library(survival)
##Resec1 cohort
surv.prob=unique(survfit(Surv(dat1$OS_time,dat1$OS_statue)~1)$surv)
fit0=coxph(Surv(dat1$OS_time,dat1$OS_statue)~dat1$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat1$OS_time[dat1$OS_statue==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat1$OS_time,dat1$OS_statue,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

##Resec2 cohort
surv.prob=unique(survfit(Surv(dat2$OS_time,dat2$OS_statue)~1)$surv)
fit0=coxph(Surv(dat2$OS_time,dat2$OS_statue)~dat2$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat2$OS_time[dat2$OS_statue==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat2$OS_time,dat2$OS_statue,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

##Resec3 cohort
surv.prob=unique(survfit(Surv(dat3$OS_time,dat3$OS_statue)~1)$surv)
fit0=coxph(Surv(dat3$OS_time,dat3$OS_statue)~dat3$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat3$OS_time[dat3$OS_statue==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat3$OS_time,dat3$OS_statue,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=365)
iAUC
AUC3=AUC
iAUC3=iAUC
utimes3=utimes

##Resec4 cohort
surv.prob=unique(survfit(Surv(dat4$OS_time,dat4$OS_statue)~1)$surv)
fit0=coxph(Surv(dat4$OS_time,dat4$OS_statue)~dat4$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat4$OS_time[dat4$OS_statue==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat4$OS_time,dat4$OS_statue,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC4=AUC
iAUC4=iAUC
utimes4=utimes

##TCIA cohort
surv.prob=unique(survfit(Surv(dat5$OS_time,dat5$OS_statue)~1)$surv)
fit0=coxph(Surv(dat5$OS_time,dat5$OS_statue)~dat5$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat5$OS_time[dat5$OS_statue==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat5$OS_time,dat5$OS_statue,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=365)
iAUC
AUC5=AUC
iAUC5=iAUC
utimes5=utimes

plot(utimes1,AUC1,type="l",xlim=c(0,60),ylim=c(0.4,1.0),bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Overall survival',col='dodgerblue3',main="Resect cohorts",pch=7)
lines(utimes2,AUC2,type="l",ylim=c(0.4,1.0),col='seagreen',pch=7)
lines(utimes3,AUC3,type="l",ylim=c(0.4,1.0),col='darkorange3',pch=7)
lines(utimes4,AUC4,type="l",ylim=c(0.4,1.0),col='orangered2',pch=7)
lines(utimes5,AUC5,type="l",ylim=c(0.4,1.0),col='darkorchid2',pch=7)
axis(1,seq(0,60,10),seq(0,60,10))
axis(2,seq(0.4,1.0,0.1),seq(0.4,1.0,0.1))
legend(35,1.0,c("Resec1  iAUC=0.755","Resec2  iAUC=0.684","Resec3  iAUC=0.655","Resec4  iAUC=0.707","TCIA      iAUC=0.687"),cex=c(1,1),lty=1,col=c('dodgerblue3','seagreen','darkorange3','orangered2','darkorchid2'),box.lty=0)

#####C-index for forestplot
r.mean=c(0.790,0.701,0.692,0.705,0.745)
r.lower=c(0.780,0.684,0.687,0.678,0.721)
r.upper=c(0.800,0.718,0.697,0.732,0.769)
#labeltext <- cbind(c("Resec1", "Resec2", "Resec3", "Resec4","TCIA"),c(rep(myspace,10)))
labeltext=cbind(c("", "", "", "",""),c(rep(myspace,5)))
forestplot.surv(labeltext=labeltext ,mean=r.mean,lower=r.lower,upper=r.upper,zero=0.5,boxsize=NULL,graphwidth=grid::unit(0.5,"inches"),x.ticks=seq(0.2,1.0,0.2),boxsize=c(0.1,0.1,0.1,0.1,0.1),lwd.ci=2,align=c("l"),col=meta.colors(box=c('dodgerblue3','seagreen','darkorange3','orangered2','darkorchid2'),line=c('dodgerblue3','seagreen','darkorange3','orangered2','darkorchid2')),xlab="C-index")

forestplot.surv(labeltext=labeltext ,mean=r.mean,lower=r.lower,upper=r.upper,zero=0.5,graph.pos=3,graphwidth=grid::unit(2,"inches"),boxsize=c(0.1,0.1,0.1,0.1,0.1),lwd.ci=2,align=c("l"),col=meta.colors(box=c('dodgerblue3','seagreen','darkorange3','orangered2','darkorchid2'),line=c('dodgerblue3','seagreen','darkorange3','orangered2','darkorchid2'),zero="darkred"),xlab="C-index")

##IMU1 cohort
surv.prob=unique(survfit(Surv(dat6$PFS,dat6$status_PFS)~1)$surv)
fit0=coxph(Surv(dat6$PFS,dat6$status_PFS)~dat6$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat6$PFS[dat6$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat6$PFS,dat6$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

##IMU2 cohort
surv.prob=unique(survfit(Surv(dat7$PFS,dat7$status_PFS)~1)$surv)
fit0=coxph(Surv(dat7$PFS,dat7$status_PFS)~dat7$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat7$PFS[dat7$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat7$PFS,dat7$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

##IMU3 cohort
surv.prob=unique(survfit(Surv(dat8$PFS,dat8$status_PFS)~1)$surv)
fit0=coxph(Surv(dat8$PFS,dat8$status_PFS)~dat8$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat8$PFS[dat8$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat8$PFS,dat8$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC3=AUC
iAUC3=iAUC
utimes3=utimes

plot(utimes1,AUC1,type="l",xlim=c(0,12),ylim=c(0.4,1.0),xlab='Times(Months)',bty="n",xaxt="n",yaxt="n",ylab='Time-dependent AUC for Progression-free survival',col='dodgerblue3',main="Immunotherapy cohorts",pch=7)
lines(utimes2,AUC2,type="l",ylim=c(0.5,1.0),col='seagreen',pch=7)
lines(utimes3,AUC3,type="l",ylim=c(0.5,1.0),col='orangered2',pch=7)
axis(1,seq(0,12,2),seq(0,12,2))
axis(2,seq(0.4,1.0,0.1),seq(0.4,1.0,0.1))
legend(7.5,1.0,c("IMU1  iAUC=0.636","IMU2  iAUC=0.698","IMU3  iAUC=0.593"),cex=c(1,1),lty=1,col=c('dodgerblue3','seagreen','orangered2'),box.lty=0)


##Chem1 cohort
surv.prob=unique(survfit(Surv(dat9$PFS,dat9$status_PFS)~1)$surv)
fit0=coxph(Surv(dat9$PFS,dat9$status_PFS)~dat9$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat9$PFS[dat9$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat9$PFS,dat9$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

##Chem2 cohort
surv.prob=unique(survfit(Surv(dat10$PFS,dat10$status_PFS)~1)$surv)
fit0=coxph(Surv(dat10$PFS,dat10$status_PFS)~dat10$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat10$PFS[dat10$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat10$PFS,dat10$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

plot(utimes1,AUC1,type="l",xlim=c(0,12),ylim=c(0.4,1.0),bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Progression-free survival',col='dodgerblue3',main="Chemotherapy cohorts",pch=7)
lines(utimes2,AUC2,type="l",ylim=c(0.4,1.0),col='brown2',pch=7)
axis(1,seq(0,12,2),seq(0,12,2))
axis(2,seq(0.4,1.0,0.1),seq(0.4,1.0,0.1))
legend(7.5,1.0,c("Chem1  iAUC=0.647","Chem2  iAUC=0.605"),cex=c(1,1),lty=1,col=c('dodgerblue3','brown2'),box.lty=0)

##TKI1 cohort
surv.prob=unique(survfit(Surv(dat11$PFS,dat11$status_PFS)~1)$surv)
fit0=coxph(Surv(dat11$PFS,dat11$status_PFS)~dat11$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat11$PFS[dat11$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat11$PFS,dat11$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

##TKI2 cohort
surv.prob=unique(survfit(Surv(dat12$PFS,dat12$status_PFS)~1)$surv)
fit0=coxph(Surv(dat12$PFS,dat12$status_PFS)~dat12$PI_5)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat12$PFS[dat12$status_PFS==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat12$PFS,dat12$status_PFS,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=10)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

plot(utimes1,AUC1,type="l",xlim=c(0,12),ylim=c(0.4,1.0),bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Progression-free survival',col='dodgerblue3',main="EGFR-TKI therapy cohorts",pch=7)
lines(utimes2,AUC2,type="l",ylim=c(0.4,1.0),col='brown2',pch=7)
axis(1,seq(0,12,2),seq(0,12,2))
axis(2,seq(0.4,1.0,0.1),seq(0.4,1.0,0.1))
legend(7.5,1.0,c("TKI1  iAUC=0.564","TKI2  iAUC=0.572"),cex=c(1,1),lty=1,col=c('dodgerblue3','brown2'),box.lty=0)
