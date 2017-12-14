####                Final Exam BIOS612    ###############
#########################################################
### Question 1
#########################################################
## rats dataset in the nlme library
rm(list = ls())
library(nlme)
?BodyWeight
# Analyse the first group of rats i.e. the rats on diet 1
data("BodyWeight")

# Dataset
ratsdata<-data.frame(y = BodyWeight$weight[BodyWeight$Diet==1],
                     tj =  BodyWeight$Time[BodyWeight$Diet==1],
                     id = BodyWeight$Rat[BodyWeight$Diet==1])
summary(ratsdata)
# 

################################################################
library(R2jags)
m <- 8
n <- 11
Y=structure(
  .Data=ratsdata$y,
  .Dim=c(8, 11))
Tj = unique(ratsdata$tj)
jags.data <- list(Y=Y, x=Tj)
#jags.params<-c("beta.zero", "beta.one", "tau", "b", "tau.zero")
jags.params<-c("beta.zero", "beta.one","sigma.b", "b", "sigma.e")
jags.inits <- function(){
  list(beta.zero=rnorm(1,0,0.01), beta.one=rnorm(1,0,0.01), sigma.b=0.01,
       b=rep(rnorm(1,0,0.001),8), sigma.e=0.01)
}

cat(
  "model{
  for(i in 1:8){
    for(j in 1:11){
      Y[i,j] ~ dnorm(mu[i,j], tau)
      mu[i,j] <- beta.zero + b[i] + beta.one*(x[j]-33.54545)
    }
    b[i] ~ dnorm(0, tau.zero)
  }
  beta.zero ~ dnorm(260, 100)
  beta.one ~ dnorm(0, 10)

  tau ~ dgamma(1,0.0260)     # for measurement error
	sigma.e <- 1.0/sqrt(tau)

  tau.zero ~ dgamma(0.1, 0.5)  # for the random intercept
  sigma.b <- 1.0/sqrt(tau.zero)

}" , fill=TRUE, file="rats.txt")

jagsfit<-jags.model("rats.txt", data = jags.data, n.chains = 3, n.adapt = 2000)
samps <- coda.samples(jagsfit, jags.params, n.iter = 15000)
burn.in <- 1000
summary(window(samps, start = burn.in))
#plot(samps)

#par(mar=c(1,1,1,1))
traceplot(samps)


print(summary(samps))



##########################################################
###############################################################################
### Question 3
############################################################################### 
# Randomized trial to determine the efficacy of a treatment on arthritis

# Outcome is self-assessment of arthritis 0|1  = poor/good; binary variable
# 294 patients
# TXAURA treatment 0/1 == placebo/auranofin
# time; 0,1,2 == baseline, 1 month, 3 month
# 

setwd("~/Desktop/Semester03/Bios 612 - Advanced Generalised Linear Models/final_assignment/")
arthritis<-read.table("arthritis.txt", header = F)
head(arthritis)

names(arthritis)<-c("ID","ARTHRIT","TXAURA","TIME")
summary(arthritis)

# Question: The Effect of aurarofin treatment on an arthritis 
# Write a report This is an example of a clustered binary data using GEE. Interpretation is 
# much easier in this setting.

# Table of number of participants vs times of measurements

table(table(arthritis$ID))  # 294 participants with 3 measurements at baseline, 1, 3 months
#table(arthritis$TIME)


table(arthritis$ARTHRIT)
table(arthritis$TXAURA)
xtabs(~arthritis$TXAURA+arthritis$ARTHRIT)
## the crude OR
(310*156)/(285*131)

# lets check the model estimates
summary(glm(cbind(ARTHRIT,1-ARTHRIT) ~ TXAURA, data = arthritis,
            family="binomial"))
exp(0.2587)  # this is the same as the crude estimate

# lets look at the empirical correlation between visits or measurements
visit <- NULL
id<-arthritis$ID
y<-arthritis$ARTHRIT

#
uid<-unique(arthritis$ID)
for(i in uid) {
  visit <- c(visit,1:length(id[id==i]))}
table(visit)
#


new.dat <- data.frame(cbind(y,id,visit))
new.dat <- new.dat[order(id),]
wide <- reshape(new.dat, v.names="y", idvar="id", timevar="visit",direction="wide")
cmat <- nmat <- matrix(0,3,3)
rmat <- wide[,-1]
#
for(i in 1:3) { for(j in 1:3) {
  njk <- sum(!is.na(rmat[,i]*rmat[,j]))
  sjk <- sum(rmat[,i]*rmat[,j],na.rm=T)/njk
  cmat[i,j] <- sjk
  nmat[i,j] <- njk
} }
vvec <- diag(cmat)  # the Variance of the covariance matrix
round(sqrt(vvec),2)

round(sqrt(cmat),2) # the covariance matrix

cormat <- cmat/(outer(sqrt(vvec),sqrt(vvec)))
print(round(cormat,2))  # the correlation matrix




################################################################
# Some set-ups based on wakefield and plots
##################
#rm(list = ls())
library(xtable)
alldat <- read.table("arthritis.txt",header=F,na.strings=".")
head(alldat)
dim(alldat)

# names(arthritis)<-c("ID","ARTHRIT","TXAURA","TIME")
names(alldat)=c("id","y","dose","time")
nrow(alldat)

# first plot averages over time
y <- dosemat <- idmat <- timemat <- matrix(0,nrow=294,ncol=3)
sum(alldat$dose==1)
sum(alldat$dose==0)
dose <- NULL
y0 <- matrix(0,nrow=147,ncol=3)
time0 <- matrix(0,nrow=147,ncol=3)
dose0 <- matrix(0,nrow=147,ncol=3)
id0 <- matrix(0,nrow=147,ncol=3)
y1 <- matrix(0,nrow=147,ncol=3)
time1 <- matrix(0,nrow=147,ncol=3)
dose1 <- matrix(0,nrow=147,ncol=3)
id1 <- matrix(0,nrow=147,ncol=3)
count0 <- count1 <- 1
p00 <- p01 <- p10 <- p11 <- rep(0,3)
counter <- 1
n0 <- n1 <- rep(0,3)
for (i in 1:294){
  dose[i] <- alldat[counter,3]
  cat("dose counter: ",dose[i],counter,"\n")
  for (j in 1:3){
    y[i,j] <- alldat[counter,2]
    
    if (dose[i]==0){
      y0[count0,j] <- y[i,j]
      if (is.na(y[i,j]) == F) {n0[j] <- n0[j]+1}
      if(j==3){count0 <-count0+1}
    }
    if (dose[i]==1){
      y1[count1,j] <- y[i,j]
      if (is.na(y[i,j]) == F) {n1[j] <- n1[j]+1}
      if(j==3){count1 <- count1+1}
      
    }
    if (is.na(y[i,j])==F ){
      if (dose[i]==0){
        if (y[i,j]==0) {p00[j]<-p00[j]+1}
        if (y[i,j]==1) {p01[j]<-p01[j]+1}
      }
      if (dose[i]==1){
        if (y[i,j]==0) {p10[j]<-p10[j]+1}
        if (y[i,j]==1) {p11[j]<-p11[j]+1}
        
      }
    }
    counter <- counter+1
  }
}
alldat0 <- as.table(id0,dose0,y0,time0)
cor0 <- cor(y0,use="pairwise")
print(xtable(cor0))
cor1 <- cor(y1,use="pairwise")
# 
# Creates Tab 9.1
#
print(xtable(cor1))
cov(y0,use="pairwise")
cov(y1,use="pairwise")
for (j in 1:3){
  p00[j] <- p00[j]/(p00[j]+p01[j])
  p01[j] <- 1-p00[j]
  p10[j] <- p10[j]/(p10[j]+p11[j])
  p11[j] <- 1-p10[j]
}
# variance of logits:
vlogit0 <- 1/(n0*p01*(1-p01))
vlogit1 <- 1/(n1*p11*(1-p11))
cat("SEs of logits for placebo: ",sqrt(vlogit0),"\n")
cat("SEs of logits for aurarofin: ",sqrt(vlogit1),"\n")
#
# Fig 9.1
#

plot(seq(1,3),p01,pch=19,ylim=c(min(p01,p11),max(p01,p11)),
     ylab="Probability of arthritis score",xlab="Measurement Occasion",xaxp=c(1,4,3))
points(seq(1,3),p11,pch=21)
legend("topleft",legend=c("auranofin","placebo"),bty=2,pch=c(21,19))
title("Probability of arthritis over time")
#placebo/auranofin

##############################################################
### plots to visualize the logits of the fits
# Fitted logits versus empirical
#

# GEE, standard errors via sandwich estimation
#names(arthritis)
library(geepack)

# For small number of measurements, it is okay to consider an independent covariance Wi
summary(geese(ARTHRIT ~ TXAURA,corstr="exchangeable", id=ID,family="binomial", 
              data = arthritis))

# Just checking here. This may not be useful
summary(geese(ARTHRIT ~ as.factor(TXAURA) + TIME + as.factor(TXAURA)*TIME,corstr="independence", id=ID,family="binomial", 
              data = arthritis))
summary(geese(ARTHRIT ~ as.factor(TXAURA)+ TIME + as.factor(TXAURA)*TIME,corstr="exchangeable", id=ID,family="binomial", 
              data = arthritis))

arthritis$TIME2<-arthritis$TIME^2
modge<-geese(ARTHRIT ~ TIME + TIME2 + TXAURA*TIME + TXAURA*TIME2,corstr="exchangeable", 
             id=ID,family="binomial",   data = arthritis)
summary(modge)


alldat$time2<-alldat$time*alldat$time
modge<-geese(y ~ time + time2 + dose:time + dose:time2,corstr="exchangeable",
             id=id, family="binomial",   data = alldat)

summary(modge)
# interpretations as odds ratios


modge$beta
gamma0 <- modge$beta[1]
gamma1 <- modge$beta[2]
gamma2 <- modge$beta[3]
gamma3 <- modge$beta[4]
gamma4 <- modge$beta[5]
t4 <- seq(1,3,1)
exp(gamma0+gamma1*t4+gamma2*t4^2)/(1+exp(gamma0+gamma1*t4+gamma2*t4^2)) # the model without int / no dose
exp(gamma0+(gamma1+gamma3)*t4+(gamma2+gamma4)*t4^2)/
  (1+exp(gamma0+(gamma1+gamma3)*t4+(gamma2+gamma4)*t4^2)) # the model with interaction 
tseq <- seq(.99,3.02,.01)
# Fig 9.8
#
plot(gamma0+gamma1*tseq+gamma2*tseq^2~tseq,type="l",xlab="Measurement Occasion",
     ylab="Logit of arthritis score", ylim=c(-2.5,0.1),xaxp=c(1,4,3) )
# , pch=".", lwd=0.01
lines(gamma0+(gamma1+gamma3)*tseq+(gamma2+gamma4)*tseq^2~tseq,type="l",
      col="grey")



#logit0 <- log(p01/(1-p01))
#logit1 <- log(p11/(1-p11))

logit0<-gamma0+gamma1*t4+gamma2*t4^2  # cheating
logit1<-gamma0+(gamma1+gamma3)*t4+(gamma2+gamma4)*t4^2

shift <- .03
points(logit0~seq(1+shift,3+shift,1),pch=19)
points(logit1~seq(1-shift,3-shift,1),pch=21)
for(j in 1:3){
  lines(x=c(j+shift,j+shift),y=c(logit0[j]-1.96*sqrt(vlogit0[j]),
                                 logit0[j]+1.96*sqrt(vlogit0[j])),lwd=2)
  lines(x=c(j-shift,j-shift),y=c(logit1[j]-1.96*sqrt(vlogit1[j]),
                                 logit1[j]+1.96*sqrt(vlogit1[j])),lwd=2,col="grey")
}
legend("bottomleft",legend=c("auranofin","placebo"),bty=2,pch=c(21,19))

c(exp(-2.5), exp(-0.5)) # The odds of arthritis in much lower in users than non-users.
#
#############################################################################
### Question 4
#############################################################################

# some summaries
# 30 patients of leprosy
# Assign 1 of 2 treatments | treatment A or treatdment B | or a placebo labelled C
# Baseline measures from all the six sites
# After 6 months, a second measure was recorded
# The outcome variable is the count of leprosy bacilli out of the six sites


rm(list = ls())
library(geepack)
## read in the data and set it up for analysis
setwd("~/Desktop/Semester03/Bios 612 - Advanced Generalised Linear Models/final_assignment/")
leprosy <- read.table("Leprosy.txt", header=T)
head(leprosy)

leprosy$Drug <- as.factor(leprosy$Drug)
leprosy$Subject <- seq(1,nrow(leprosy))
lleprosy <- reshape(leprosy, direction="long", idvar="Subject",
                    varying=list(c("PreBC", "PosBC")), v.names=c("Count"), timevar="time",
                    times=c(0,1))
lleprosy <- lleprosy[order(lleprosy$Subject),]
trt1<-I(lleprosy$Drug=="A")
trt2<-I(lleprosy$Drug=="B")
head(lleprosy)


## (a) use GEE to fit the log linear model
rhs <- Count ~ time + time:trt1 + time:trt2
gee1 <- geese(rhs, id=Subject, data=lleprosy, corstr="exchangeable", family="poisson")
summary (gee1)
# The estimates and SEs
round(cbind(gee1$beta, sqrt(diag(gee1$vbeta))),3)
# The estimates and 95% confidence intervals
exp(cbind(gee1$beta, gee1$beta + t(matrix(c(-1,1),nrow=2)%*%matrix(sqrt(diag(gee1$vbeta)),nrow=1)) ))


## (b) lme with random intercepts and condition poisson likelihood
library(lme4)

lme1<-glmer( Count~time+time:trt1+time:trt2 + (1|Subject), family=poisson, data=lleprosy)
summary(lme1)


# (c) Use bayesian inference for the analysis

lleprosy$T1<-trt1
lleprosy$T2<-trt2
head(lleprosy)


# Models 
#
cat("model
    {
    ## Specify likelihood
    for (i in 1:n){
    for (j in 1:k){
    Y[i,j] ~ dpois(mu[i,j])
    log(mu[i,j]) <- beta0+beta1*time[i,j]+beta2*time[i,j]*T1[i,j]+beta3*time[i,j]*T2[i,j]
    }
    }
    beta0 ~ dnorm(0,10.04) # ~dunif(-10000,10000)#
    beta1 ~ dnorm(0,1.04) # ~dunif(-10000,10000)#
    beta2 ~ dnorm(0,10.04) # ~dunif(-10000,10000)#
    beta3 ~ dnorm(0,10.04) # ~dunif(-10000,10000)#
    tau ~ dgamma(.5,.01)
    }", file="lep.txt")


# n is the number of unique subjects and k is the two reps

# Initial estimates     
initials<-list(beta0=1, beta1=0.05, beta2=0.02, beta3=0.04)

# Data
datlep<-list(n = length(unique(lleprosy$Subject)),
             k = 2,
             time = matrix(lleprosy$time, ncol = 2, byrow = T),
             Y = matrix(lleprosy$Count, ncol = 2, byrow = T), 
             T1 = matrix(lleprosy$T1, ncol = 2, byrow = T), 
             T2 = matrix(lleprosy$T2, ncol = 2, byrow = T))

params<-c("beta0","beta1","beta2","beta3")


lep.m1<-jags.model("lep.txt",data=datlep,n.chains = 2, n.adapt = 10000,inits=initials)

samps <- coda.samples(lep.m1, params, n.iter = 50000)
#
burn.in <- 1000
summary(window(samps, start = burn.in,frequency=45))

par(mar=c(1,1,1,1))
#
traceplot(samps)
#
plot(samps)

density(samps)

# (d) Summarize the results in words for non-statisticians

# solved





