rm(list = ls())

########################################################


## Question II
setwd("~/Desktop/Semester03/Bios 612 - Advanced Generalised Linear Models/assignment04/")
cowdata<-read.table("CowData.txt", header = F, sep = "")
head(cowdata)

names(cowdata)<-c("diet","cow","week","protein")
cowdata<-data.frame(diet=cowdata$diet, cow=cowdata$cow,
                    week=cowdata$week, protein= cowdata$protein)

misscowdata<-cowdata

##### Data summary
#
summary(cowdata)
head(cowdata)

require(tableone)
# the number of cows per diet type
c(length(unique(cowdata[cowdata$diet==1,"cow"])),
  length(unique(cowdata[cowdata$diet==2,"cow"])),
  length(unique(cowdata[cowdata$diet==3,"cow"])))

# the number of week per diet type
c(length(unique(cowdata[cowdata$diet==1,"week"])),
  length(unique(cowdata[cowdata$diet==2,"week"])),
  length(unique(cowdata[cowdata$diet==3,"week"])))

#
##### Subjects and obs/subject
#
cat( paste("Number of cows =", length(table(cowdata$cow)),"\n\n") )
cat( paste("Number of measurements =", length(cowdata$cow),"\n\n") )
cat( "Number of subjects with a given observations-per-subject:\n\n" )
print( table( table( cowdata$cow ) ) )  # the only important part of the plot 
cat("\n\n")
#
#####

#####
##### mean plots
#####
#

##### protein content with time for the different diet exposures
#
## lets do this for each type of diet  barley/mixture/lupin

cowdata.barley<-cowdata[cowdata$diet==1,]
cowdata.mixture<-cowdata[cowdata$diet==2,]
cowdata.lupin<-cowdata[cowdata$diet==3,]


uid <- unique( cowdata.barley$cow )
plot(cowdata.barley$week, cowdata.barley$protein,
     xlab="Weeks since calving",
     ylab="Protein measurement", pch="o", col="green", lwd=0.05)
for( j in 1:length(uid) ){
  lines( cowdata.barley$week[ cowdata.barley$cow==uid[j] ], 
         cowdata.barley$protein[ cowdata.barley$cow==uid[j] ], col=j, lwd=1 )
}
title("barley diet")

###########
uid <- unique( cowdata.mixture$cow )
plot(cowdata.mixture$week, cowdata.mixture$protein,
     xlab="Weeks since calving",
     ylab="Protein measurement", pch="o", col="green", lwd=0.05)
for( j in 1:length(uid) ){
  lines( cowdata.mixture$week[ cowdata.mixture$cow==uid[j] ], 
         cowdata.mixture$protein[ cowdata.mixture$cow==uid[j] ],col=j, lwd=1 )
}
title("mixture diet")

###########
uid <- unique( cowdata.lupin$cow )
plot(cowdata.lupin$week, cowdata.lupin$protein,
     xlab="Weeks since calving",
     ylab="Protein measurement", pch="o", col="green", lwd=0.05)
for( j in 1:length(uid) ){
  lines( cowdata.lupin$week[ cowdata.lupin$cow==uid[j] ], 
         cowdata.lupin$protein[ cowdata.lupin$cow==uid[j] ] ,col=j, lwd=1)
}
title("lupin diet")
#

## overall mean plots
plot( cowdata$week, cowdata$protein, pch=".", xlab="Weeks since calving",
      ylab="Protein measurement" )
#lines( smooth.spline( cowdata$week, cowdata$protein, df=7 ) )
lines( smooth.spline( cowdata.barley$week, cowdata.barley$protein, df=7 ), col=2 )
lines( smooth.spline( cowdata.mixture$week, cowdata.mixture$protein, df=7 ), col=3 )
lines( smooth.spline( cowdata.lupin$week, cowdata.lupin$protein, df=7 ), col=4 )
legend(x = 2, y = 3, legend = c("barley", "mixture", "lupin"), col = c(2,3,4), lty = -1, pch = 1, bg = "gray90")
title("mean protein content by diet")



#####
##### covariance summaries
##### 
#
require(stats); require(graphics)
require(splines)  
summary(cowdata$week)  
fit <- lm( protein ~ as.factor(diet)+ns(week, knots = 3), data=cowdata )
resids <- cowdata$protein - fitted( fit )
#
#################################################################
#### testing | used as reference matrix here to index the positions with values
# here we created two vectors of values Rij and NA
matprotein<-matrix(rep(rep(NA,19),79), nrow = 79,byrow = T)
matRij<-matrix(rep(rep(0,19),79), nrow = 79,byrow = T)

for(i in 1:79){
  for(j in 1:19){
    if(j %in% misscowdata[misscowdata$cow==i,"week"])
    {
      matprotein[i,j]<-misscowdata[misscowdata$cow==i&misscowdata$week==j,"protein"]
      matRij[i,j]<-1
    }
  }
}


##### look at the covariance structure

nobs <- length( cowdata$protein )
nsubjects <- length( table( cowdata$cow ) )
rmat <- matrix( NA, nsubjects, 19 )
ycat <- c(1:19)
nj <- unlist( lapply( split( cowdata$cow, cowdata$cow ), length ) )
for( j in 1:19 ){
  legal <- ( cowdata$week >= ycat[j]-0.5 )&( cowdata$week < ycat[j]+0.5 )
  jtime <- cowdata$week + 0.01*rnorm(nobs)
  t0 <- unlist( lapply(
    split( abs(jtime - ycat[j]) , cowdata$cow ), min ) )
  tj <- rep( t0, nj )
  keep <- ( abs( jtime - ycat[j] )==tj )&( legal )
  yj <- rep( NA, nobs )
  yj[keep] <- resids[keep]
  yj <- unlist( lapply( split( yj, cowdata$cow), min, na.rm=T ) )
  rmat[ , j ] <- yj
}
#
dimnames( rmat ) <- list( NULL, paste("week",c(1:19)))
pairs( rmat )

##### covariance matrix
#
cmat <- matrix( 0, 19, 19 )
nmat <- matrix( 0, 19, 19 )
#
for( j in 1:1 ){
  for( k in j:1 ){
    njk <- sum( !is.na( rmat[,j]*rmat[,k] ) )
    sjk <- sum( rmat[,j]*rmat[,k], na.rm=T )/njk
    cmat[j,k] <- sjk
    nmat[j,k] <- njk
  }
}
print( round( cmat, 2 ) )
vvec <- diag(cmat)
cormat <- cmat/( outer( sqrt(vvec), sqrt(vvec) ) )
print( round( cormat, 2 ) )
print( nmat )
#
##### Variogram estimation
#
source("variogram.q")
#
out <- lda.variogram( id=cowdata$cow, y=resids, x=cowdata$week )
dr <- out$delta.y
dt <- out$delta.x
#
var.est <- var( resids )
#
plot( dt, dr, pch=".", ylim=c(0, 1.2*var.est) )
lines( smooth.spline( dt, dr, df=5 ), lwd=3 )
abline( h=var.est, lty=2, lwd=2 )
title("Protein residual variogram")

#
########### missing cowdata  #############################################

summary(misscowdata)
#misscowdata$Rij<-ifelse(is.na(misscowdata$protein),0,1)

# here we created two vectors of values Rij and NA
matprotein<-matrix(rep(rep(NA,19),79), nrow = 79,byrow = T)
matRij<-matrix(rep(rep(0,19),79), nrow = 79,byrow = T)
#dim(matfoo)
for(i in 1:79){
  for(j in 1:19){
    if(j %in% misscowdata[misscowdata$cow==i,"week"])
    {
      matprotein[i,j]<-misscowdata[misscowdata$cow==i&misscowdata$week==j,"protein"]
      matRij[i,j]<-1
    }
  }
}


# lets reconsruct the dataset
proteinvec<-c()  # unstacked
for(i in 1:79){
  proteinvec<-c(proteinvec,matprotein[i,])
}
Rijvec<-c()
for(i in 1:79){
  Rijvec<-c(Rijvec,matRij[i,])
}

diet<-c()
for(i in 1:79){ # creating a vector for the unique diet variables
  diet<-c(diet,unique(misscowdata[misscowdata$cow==i,"diet"]))
}

## dataframe of the re-unstacked values
remissdata<-data.frame(cow=sort(rep(1:79,19)), week=rep(1:19,79), diet =rep(diet, each=19) ,
                       protein= proteinvec, Rij=Rijvec )

# A matrix for all the previous protein measures
prevmatrix<-matrix(data = rep(NA, (79*19)), nrow = 79, byrow = T)
prevmatrix[,2:19]<-matprotein[,1:18]
prevproteinvec<-c()  # 
for(i in 1:79){
  prevproteinvec<-c(prevproteinvec,prevmatrix[i,])
}
# lets add the vector to the dataframe
remissdata$PrevProtein<-prevproteinvec

summary(fit<-glm(Rij~PrevProtein+as.factor(diet), family=binomial,  data = remissdata))
exp(fit$coef)

exp(confint(fit))



## overall mean plots
###################################################
## over 19 weeks
# collect index of subs with 19 measurement
vec19<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==19){
    vec19<-c(vec19,i)
  }
}
uid <- unique( vec19 )
plot(remissdata[remissdata$cow%in%vec19,"week"], remissdata[remissdata$cow%in%vec19,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green")
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
lines( smooth.spline(remissdata[remissdata$cow%in%vec19,"week"], remissdata[remissdata$cow%in%vec19,"protein"], df=7 ), lwd=2 )
title("protein intake for 19 weeks")



# collect index of subs with 18 measurement
vec18<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==18){
    vec18<-c(vec18,i)
  }
}
uid <- unique( vec18 )
plot(remissdata[remissdata$cow%in%vec18,"week"], remissdata[remissdata$cow%in%vec18,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green" )
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata18<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec18,"week"], 
                     missdata[missdata$cow%in%vec18,"protein"], df=5 ), lwd=2 )
title("protein intake for 18 weeks")


# collect index of subs with 17 measurement
vec17<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==17){
    vec17<-c(vec17,i)
  }
}
uid <- unique( vec17 )
plot(remissdata[remissdata$cow%in%vec17,"week"], remissdata[remissdata$cow%in%vec17,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green" )
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata17<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec17,"week"], missdata[missdata$cow%in%vec17,"protein"], df=5 ), lwd=2 )
title("protein intake for 17 weeks")


# collect index of subs with 16 measurement
vec16<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==16){
    vec16<-c(vec16,i)
  }
}
uid <- unique( vec16 )
plot(remissdata[remissdata$cow%in%vec16,"week"], remissdata[remissdata$cow%in%vec16,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green")
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata16<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec16,"week"], 
                     missdata[missdata$cow%in%vec16,"protein"], df=5 ), lwd=1.5 )
title("protein intake for 16 weeks")

############################
# collect index of subs with 15 measurement
vec15<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==15){
    vec15<-c(vec15,i)
  }
}
uid <- unique( vec15 )
plot(remissdata[remissdata$cow%in%vec15,"week"], remissdata[remissdata$cow%in%vec15,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green" )
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green")
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata15<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec15,"week"], 
                     missdata[missdata$cow%in%vec15,"protein"], df=5 ), lwd=1.5)
title("protein intake for 15 weeks")


# collect index of subs with 14 measurement
vec14<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==14){
    vec14<-c(vec14,i)
  }
}
uid <- unique( vec14 )
plot(remissdata[remissdata$cow%in%vec14,"week"], remissdata[remissdata$cow%in%vec14,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green")
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata14<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec14,"week"],
                     missdata[missdata$cow%in%vec14,"protein"], df=4 ), lwd=1.5 )
title("protein intake for 14 weeks")

###############################
# collect index of subs with 13 measurement
vec13<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==13){
    vec13<-c(vec13,i)
  }
}
uid <- unique( vec13 )
plot(remissdata[remissdata$cow%in%vec13,"week"], remissdata[remissdata$cow%in%vec13,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green")
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=0.8, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata13<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec13,"week"],
                     missdata[missdata$cow%in%vec13,"protein"], df=4 ), lwd=1.5 )
title("protein intake for 13 weeks")


#################
# collect index of subs with 12 measurement
vec12<-c()
for(i in 1:79){
  if( sum(remissdata[remissdata$cow==i&sum(remissdata$Rij),"Rij"])==12){
    vec12<-c(vec12,i)
  }
}
uid <- unique( vec12 )
plot(remissdata[remissdata$cow%in%vec12,"week"], remissdata[remissdata$cow%in%vec12,"protein"],
     xlab="Weeks since calving",
     ylab="Response", pch="o", lwd=1, col="green")
for( j in 1:length(uid) ){
  lines( remissdata$week[ remissdata$cow==uid[j] ], 
         remissdata$protein[ remissdata$cow==uid[j] ], lwd=1, col="green" )
}
missval<-!is.na(remissdata$protein); sum(missval)
missdata12<-missdata<-remissdata[missval,]
lines( smooth.spline(missdata[missdata$cow%in%vec12,"week"],
                     missdata[missdata$cow%in%vec12,"protein"], df=4 ), lwd=2 )
title("protein intake for 12 weeks")

#################################################################


#### The Response profiles

plot(remissdata$week, remissdata$protein,
     xlab="Time (Weeks since calving)",
     ylab="Response", pch=".", lwd=0.01, col="green")
lines( smooth.spline(remissdata[remissdata$cow%in%vec19,"week"], 
                     remissdata[remissdata$cow%in%vec19,"protein"], df=6 ), lwd=1.2, col=2 )
lines( smooth.spline(missdata18[missdata18$cow%in%vec18,"week"], 
                     missdata18[missdata18$cow%in%vec18,"protein"], df=6 ), lwd=1.2, col=3  )
lines( smooth.spline(missdata17[missdata17$cow%in%vec17,"week"], 
                     missdata17[missdata17$cow%in%vec17,"protein"], df=6 ), lwd=1.2, col=4  )
lines( smooth.spline(missdata[missdata$cow%in%vec16,"week"], 
                     missdata[missdata$cow%in%vec16,"protein"], df=6 ), lwd=1.2, col=5 )
lines( smooth.spline(missdata[missdata$cow%in%vec15,"week"], 
                     missdata[missdata$cow%in%vec15,"protein"], df=6 ), lwd=1.2, col=6 )
lines( smooth.spline(missdata[missdata$cow%in%vec14,"week"], 
                     missdata[missdata$cow%in%vec14,"protein"], df=6 ), lwd=1.2, col=7 )
lines( smooth.spline(missdata[missdata$cow%in%vec13,"week"], 
                     missdata[missdata$cow%in%vec13,"protein"], df=6 ), lwd=1.2, col=8 )
lines( smooth.spline(missdata12[missdata12$cow%in%vec12,"week"], 
                     missdata12[missdata12$cow%in%vec12,"protein"], df=6 ), lwd=1.2, col=9 )
legend("topright", legend = c("19 weeks", "18 weeks", "17 weeks", "16 weeks", "15 weeks",
                              "14 weeks", "13 weeks", "12 weeks"),
       col = c(2,3,4,5,6,7,8,9), lty = -1, pch = 1, bg = "gray90", cex = 0.75)
title("mean response by weeks since calving")



## fitting the model and residual plots
require(nlme)
require(splines)
summary(remissdata)
summary(cowdata)
lme1 <- lme( protein~as.factor(diet)+ns(week, knots = 3) ,
             random=~1|cow, method = "REML",
             #random = reStruct( ~ 1 | cow, pdClass="pdSymm", REML=F),
             correlation = corAR1( form = ~ 1 | cow ), data = cowdata )
summary(lme1)
# 

kk<-lme1$coefficients$fixed
se<-c(0.04469435,0.04930606,0.04938486,0.07021451,0.04373151)
ci<-cbind(ll= (kk-1.96*se),up= kk+1.96*se)
ci

# observed versus fitted values by Subject
plot(lme1, protein ~ fitted(.) | diet, abline = c(0,1))

# standardized residuals versus fitted values by gender
plot(lme1, resid(., type = "p") ~ fitted(.) | diet, abline = 0)


## normal plots of random effects
qqnorm(lme1, ~ranef(.))
plot(lme1)

##############################################################################################
#  Question I Jags
##############################################################################################
#  
##########################################################################################
#  boi and b1i are independent

cat("model{ 
    for (i in 1:2000) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1*(time[i] - mean(time[]))
    
    # missing data likelihood
    miss[i] ~ dbern(p[i])
    logit(p[i]) <- a + b*(y[i])
    }
    
    # Priors
    a <- -0.5
    b <- -0.01
    

    beta0 ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)
    tau ~ dgamma(0.001, 0.001)
    }", file="MAR.txt")



# Bundle data
require(R2jags)
dat<-list(y = stacked.bayes$y, 
          time = stacked.bayes$time,
          miss = stacked.bayes$miss )

params<-c("mu", "beta0", "beta1")


modMAR<-jags(data = dat,
              parameters.to.save = params,
              model.file = "MAR.txt",
              n.iter = 5000,
              n.burnin = 1000,
              n.thin = 1)

print(modMAR)
traceplot(outj)
res<-modMAR$BUGSoutput
res$median


jagsfit<-jags.model("MAR.txt", data = dat, n.chains = 3, n.adapt = 1000)
samps <- coda.samples(jagsfit, params, n.iter = 15000)
burn.in <- 1000
summary(window(samps, start = burn.in))
plot(samps)
#################
# 
# # Bayesian analysis using JAGS also works

#Define model
cat("model {
    # Likelihood
    for (i in 1:n_obs) {
    mu[i] <- beta0[subjects[i]]+beta1[subjects[i]]*x[i]
    #mu[i] <- beta0+beta1*x[i]
    y[i]~dnorm(mu[i], tau_res)

    # missing data likelihood
    miss[i] ~ dbern(p[i])
    logit(p[i]) <- a + b*(y[i])
    }

    # Priors
    a <- -0.5
    b <- -0.01

    mu_int~dnorm(0, 0.001)            # Mean hyperparameter for random intercepts
    sigma_int~dunif(0, 100)           # SD hyperparameter for random intercepts
    tau_int <- 1/(sigma_int*sigma_int)
    mu_slope~dnorm(0, 0.001)          # Mean hyperparameter for random slopes
    sigma_slope~dunif(0, 100)         # SD hyperparameter for slopes
    tau_slope <- 1/(sigma_slope*sigma_slope)
    for (i in 1:n_subj) {
    beta0[i]~dnorm(mu_int, tau_int)      # Random intercepts
    beta1[i]~dnorm(mu_slope, tau_slope)  # Random slopes
    }
    sigma_res~dunif(0, 100) # Residual standard deviation
    tau_res <- 1/(sigma_res*sigma_res)   # Residual precision

    }", fill=TRUE, file="lme_model2.txt")

# Bundle data

dat<-list(y = stacked.bayes$y,
          time = stacked.bayes$time,
          miss = stacked.bayes$miss )
#attach(dat) # centered age already
jags_data <- list(y=as.numeric(dat$y), subjects=as.numeric(stacked.bayes$id),
                  x=as.numeric(dat$time), n_subj=length(unique(stacked.bayes$id)), n_obs=2000 )

n_subjects<-200
# Inits function
inits <- function() {
  list(beta0=rnorm(n_subjects, 0, 2), beta1=rnorm(n_subjects, 0, 2), mu_int=rnorm(1, 0, 1),
       sigma_int=rlnorm(1), mu_slope=rnorm(1, 0, 1), sigma_slope=rlnorm(1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("beta0", "beta1", "mu_int", "sigma_int", "mu_slope", "sigma_slope", "sigma_res")

# MCMC settings
ni <- 11000; nb <- 1000; nt <- 10; nc <- 3

# Start Gibbs sampling
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model2.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
traceplot(outj)
# -- you can see that the chains are mixing much better when we use the appropriate model for the data
# -- see Gelman's `folk theorem of statistical computing' (http://andrewgelman.com/2008/05/13/the_folk_theore/):
# when you have computational problems, often there's a problem with your model

print(outj, dig=3)

out <- outj$BUGSoutput
out$median
#
#
tt0<-out$median$beta0; median(tt0)
hist(tt0)
# 101.2785

tt1<-out$median$beta1; median(tt1)
hist(tt1)
# -4.91573

jagsfit<-jags.model("lme_model2.txt", data = jags_data, n.chains = 3, n.adapt = 1000)
samps <- coda.samples(jagsfit, params, n.iter = 15000)
burn.in <- 1000
summary(window(samps, start = burn.in))
plot(samps)
########################################################








