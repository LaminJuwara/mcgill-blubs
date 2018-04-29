library(foreign) 

download.file("https://meps.ahrq.gov/mepsweb/data_files/pufs/h163ssp.zip",
              temp <- tempfile())
unzipped_file = unzip(temp)
h163 = read.xport(unzipped_file)
unlink(temp)

my.h163<-data.frame(TOTEXP13=h163$TOTEXP13,DIABDX=h163$DIABDX,
                    AGE13X=h163$AGE13X,REGION13=h163$REGION13,
                    SEX=h163$SEX,RACEV1X=h163$RACEV1X,MARRY13X=h163$MARRY13X,
                    EDRECODE=h163$EDRECODE,POVCAT13=h163$POVCAT13,FAMINC13=h163$FAMINC13,
                    HIBPDX=h163$HIBPDX,STRKDX=h163$STRKDX,CHOLDX=h163$CHOLDX,
                    ASTHDX=h163$ASTHDX,ARTHDX=h163$ARTHDX,
                    BMINDX53=h163$BMINDX53,ADSMOK42=h163$ADSMOK42,
                    EMPST31=h163$EMPST31,INSCOV13=h163$INSCOV13)

dim(my.h163)
## Outcome TOTEXP13~N(mu,sigma) 
## The treatment DIABDX is a binary variable eg. propensity score model with logistic
## 17 other potential confounders

### Causal Adjustment Analysis
str(my.h163)
summary(my.h163$TOTEXP13)  # Outcome TOTEXP13 USD, no NAs

(table(my.h163$DIABDX))  # the Z variable 
# We will exclude the inapplicable individual age <18
subdata<-my.h163[my.h163$DIABDX!=-1,]
subdata$DIABDX[subdata$DIABDX<0]<-NA  # All other unknowns are coded NAs
subdata$DIABDX<-ifelse(subdata$DIABDX==1,1,0)
table(subdata$DIABDX)

subdata$SEX<-as.factor(subdata$SEX)   # Male = 1 and Female = 2

subdata$REGION13[subdata$REGION13<0]<-NA
subdata$REGION13<- as.factor(subdata$REGION13) 

subdata$AGE13X[subdata$AGE13X<0]<-NA  # Age since 2013. Less than 2013 is -1
subdata$AGE13X<-as.numeric(subdata$AGE13X)

table(subdata$RACEV1X)
#subdata$RACEV1X<-ifelse(subdata$RACEV1X==1,"White","Not White")
subdata$RACEV1X<-as.factor(subdata$RACEV1X)
table(subdata$RACEV1X)

subdata$MARRY13X[subdata$MARRY13X<0]<-NA
subdata$EDRECODE[subdata$EDRECODE<0]<-NA
subdata$FAMINC13<-subdata$FAMINC13/1000   ## family income in 1000s
subdata$HIBPDX[subdata$HIBPDX<0]<-NA
subdata$STRKDX[subdata$STRKDX<0]<-NA
subdata$CHOLDX[subdata$CHOLDX<0]<-NA
subdata$ARTHDX[subdata$ARTHDX<0]<-NA
subdata$ASTHDX[subdata$ASTHDX<0]<-NA
subdata$BMINDX53[subdata$BMINDX53<0]<-NA
subdata$ADSMOK42[subdata$ADSMOK42<0]<-NA
subdata$EMPST31[subdata$EMPST31<0]<-NA
summary(subdata)


require(tableone)  # lest summarize by the outcome diabdx
str(subdata)
table(subdata$DIABDX)
Vars<-c( #"TOTEXP13",
         "AGE13X","REGION13","SEX","RACEV1X", "MARRY13X","EDRECODE",
         "POVCAT13","FAMINC13","HIBPDX","STRKDX","CHOLDX","ASTHDX","ARTHDX",
         "BMINDX53","ADSMOK42","EMPST31","INSCOV13")

tableOne<-CreateTableOne(vars = Vars, strata = "DIABDX",data = subdata, test = TRUE)
print(tableOne,smd = T)

#dim(subdata)
#summary(subdata)

# reformat the datafram into numeric type
subdata$AGE13X<-as.numeric(subdata$AGE13X)
subdata$SEX<-as.numeric(subdata$SEX) # reset
subdata$RACEV1X<-as.numeric(subdata$RACEV1X)
subdata$REGION13<- as.numeric(subdata$REGION13) 
newdata<-subdata


### Dealing with missing values
# % of missings
mean(is.na(subdata))*100
colMeans(is.na(subdata))*100  # column % of missing values
summary(newdata) 


## obtain 100 imputations of the dataset
library(mice)
tempdatamice <- mice(newdata,m = 100,maxit=2,meth='pmm',seed=4490, printFlag=F)
save(tempdatamice,file = "tempmice.RData") # save a copy of imputed dataset

completeData<-complete(tempdatamice,1) # visualize a complete dataset
summary(completeData)

#############
# Propensity score based causal adjustment methods

# lets consider complete case analysis
subdata<-subdata[complete.cases(subdata),]
# Checking Confounder balance eg. age

Diabetes<-subdata$DIABDX
age0 <- subdata$AGE13X[Diabetes==0]
age1 <- subdata$AGE13X[Diabetes==1]
ecdf0 <- ecdf(age0)
ecdf1 <- ecdf(age1)
par(mar=c(5,5,3,3));par(mfrow=c(1,1))
plot(ecdf0, verticals=TRUE, do.points=FALSE,main='Empirical cdfs for Age',col='blue')
plot(ecdf1, verticals=TRUE, do.points=FALSE, add=TRUE, col='red')
legend(20,0.9,c('Diabetes = 0', 'Diabetes = 1'),col=c('blue','red'),lty=1)
# we observe clear age inbalance in diabetes patients as the empirical cdf do not
# overlap

# lets build the propensity score model:

varsPS<-names(subdata)[4:19]
PS.form = formula(paste("DIABDX~AGE13X+",paste(varsPS,collapse="+")))
ps.mod <- glm(PS.form, data=subdata,family="binomial")
ps.lr <- predict(ps.mod,type="response")
summary(ps.lr)

ps0<-ps.lr[Diabetes==0]
ps1<-ps.lr[Diabetes==1]
#quints <- c(0,quantile(ps.lr,seq(0.28,1,.16)))
quints <- c(0,quantile(ps.lr,seq(0.30,1,.16)))
rbind(table(cut(ps.lr[Diabetes==0],quints)),
      table(cut(ps.lr[Diabetes==1],quints))) # see the distribution of cases

#
par(mar=c(2,2,2,0));par(mfrow=c(1,1))
hist(ps0, col=rgb(0,0,1,0.5), breaks=seq(0,1,by=0.05), ylim=c(0,15000),
     main="Propensity Score overlap", xlab="PS")
hist(ps1, col=rgb(1,0,0,0.5), breaks=seq(0,1,by=0.05), add=T);box()
legend(0,0.9,c('Diabetes = 0', 'Diabetes = 1'),col=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)),lty=1)

boxplot(ps0,ps1,ylab="PS",xlab="Treatment Group",names=c('Diabetes = 0', 'Diabetes = 1'),
        col=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)),ylim=c(-0.002,1));abline(h=quints,col="red")


# lets look at quintile specific ecdf
Pcat<-as.numeric(cut(ps.lr,quints,include.lowest=T))
par(mar=c(5,4,3,3),mfrow=c(3,2),oma=c(0,0,2,0))
for(k in 1:5){
  age0 <- subdata$AGE13X[Diabetes==0 & Pcat==k]
  age1 <- subdata$AGE13X[Diabetes==1 & Pcat==k]
  ecdf0 <- ecdf(age0)
  ecdf1 <- ecdf(age1)
  plot(ecdf0, verticals=TRUE, do.points=FALSE,
       main=substitute(paste('Quintile ',k),list(k=k)),col='blue')
  plot(ecdf1, verticals=TRUE, do.points=FALSE, add=TRUE, col='red')
}
plot(age0,type='n',ylim=range(0,1),axes=FALSE)
title("ECDFs for Age by PS quintile",outer = TRUE)
legend(30,0.75,c('Diabetes = 0', 'Diabetes = 1'),col=c('blue','red'),lty=1)
# Balance is not achieved in all strata

#### lets obtain stratum specific SMD estimates
W<-(Diabetes==0)/(1-ps.lr)+(Diabetes==1)/ps.lr
smd.mat<-ExtractSmd(tableOne)
vars = c("AGE13X",varsPS)
for(k in 1:5){
  subQ<-subdata[Pcat == k,]
  tabQs <- CreateTableOne(vars = c("AGE13X",varsPS), strata = "DIABDX", data = subQ, test = FALSE)
  smd.mat<-cbind(smd.mat,ExtractSmd(tabQs))
}
colnames(smd.mat)<-c('Original','Q1','Q2','Q3','Q4','Q5')
round(smd.mat,4)
## We see a reasonable balance in each strata. Increasing the number
## of strata would not necessarilly improve the balance


####################### Estimating ATE ###################################
## Using PS matching
library(Matching)
subdata$ps.lr<-ps.lr
ps.lr.match <- Match(Tr=subdata$DIABDX,X=subdata$ps.lr,estimand="ATE",ties=FALSE)
matched.samp <- subdata[c(ps.lr.match$index.control,ps.lr.match$index.treated),]
table(table(c(ps.lr.match$index.control, ps.lr.match$index.treated)))

tabMatched <- CreateTableOne(vars = vars, strata = "DIABDX",data = matched.samp, test = FALSE)
MatchBalance(PS.form,data=subdata,match.out=ps.lr.match)

age0 <- matched.samp$AGE13X[matched.samp$DIABDX==0]
age1 <- matched.samp$AGE13X[matched.samp$DIABDX==1]
ecdf0 <- ecdf(age0)
ecdf1 <- ecdf(age1)
par(mfrow=c(1,1))
par(mar=c(5,4,3,2))
plot(ecdf0, verticals=TRUE, do.points=FALSE,main='Empirical cdfs for Age in matched sample',col='blue')
plot(ecdf1, verticals=TRUE, do.points=FALSE, add=TRUE, col='red')
legend(20,0.9,c('Diabetes = 0', 'Diabetes = 1'),col=c('blue','red'),lty=1)
## We see a clear evidence of overlap in the strata for age ecdf


## PS for IPW method
library(survey);library(Hmisc);library(htmlwidgets)
ps.lr.weight <- subdata$DIABDX/ps.lr + (1-subdata$DIABDX)/(1-ps.lr)
subdata.IPW.lr <- svydesign(ids=~0, data=subdata, weights=ps.lr.weight)
tabIPW <- svyCreateTableOne(vars = vars, strata = "DIABDX",data = subdata.IPW.lr, test = FALSE)
print(tabIPW, smd = TRUE)
## smd in age is problematic

temp0 <- Ecdf(subdata$AGE13X[Diabetes==0],weights=ps.lr.weight[Diabetes==0],
              pl=FALSE)
temp1 <- Ecdf(subdata$AGE13X[Diabetes==1],weights=ps.lr.weight[Diabetes==1],pl=FALSE)
par(mar=c(2,2,2,0))
plot(temp0$x,temp0$y,ylab="ECDF(Age)",xlab="Age",
     main='Empirical cdfs for Age in weighted sample',col='blue',type="s",lwd=1)
lines(temp1$x,temp1$y,col="red",lwd=1,type='s')
# the estimate gets worse when we use IPW
# Here the matched samples show good balance but not the IPW samples



# Using outcome regression estimators
subdata.alldiab <- subdata
subdata.alldiab$DIABDX <- 1
subdata.nodiab <- subdata
subdata.nodiab$DIABDX <- 0
vars.list<-names(subdata)[3:19]
# no interaction model
out.form = formula(paste("TOTEXP13~DIABDX+",paste(vars.list,collapse="+")))
mod1.lmX <- lm(out.form,data=subdata)
APO.lmX.1 <- mean(predict(mod1.lmX,subdata.alldiab))
APO.lmX.0 <- mean(predict(mod1.lmX,subdata.nodiab))
APO.lmX.1 - APO.lmX.0

# with interactions
vars.list.int<-c(vars.list,paste("DIABDX:",paste(vars.list,collapse=":")))
out.form.int = formula(paste("TOTEXP13~DIABDX+",paste(vars.list.int,collapse="+")))
mod1.lmX.int <- lm(out.form.int,data=subdata)
APO.lmX.1.int <- mean(predict(mod1.lmX.int,subdata.alldiab))
APO.lmX.0.int <- mean(predict(mod1.lmX.int,subdata.nodiab))
APO.lmX.1.int - APO.lmX.0.int
# we obtain better estimates when interactions are added into the Out model


# Using PS matching to compute ATE
ps.lr.match <- Match(Tr=subdata$DIABDX,
                     X=subdata$ps.lr,estimand="ATE",ties=FALSE)
matched.samp <- subdata[c(ps.lr.match$index.control,ps.lr.match$index.treated),]
dim(matched.samp)
mean(matched.samp$TOTEXP13[matched.samp$DIABDX == 1]) -
  mean(matched.samp$TOTEXP13[matched.samp$DIABDX == 0])


## PS regression
library(splines)
mod1.PSlm1 <- lm(TOTEXP13~DIABDX+ps.lr,data=subdata)
APO.PSlm1.1 <- mean(predict(mod1.PSlm1,subdata.alldiab))
APO.PSlm1.0 <- mean(predict(mod1.PSlm1,subdata.nodiab))
APO.PSlm1.1 - APO.PSlm1.0

mod1.PSlm2 <- lm(TOTEXP13~DIABDX+ps.lr+I(ps.lr^2),
                 data=subdata)
APO.PSlm2.1 <- mean(predict(mod1.PSlm2,subdata.alldiab))
APO.PSlm2.0 <- mean(predict(mod1.PSlm2,subdata.nodiab))
APO.PSlm2.1 - APO.PSlm2.0

mod1.PSlm3 <- lm(TOTEXP13~DIABDX+bs(ps.lr,df=4),
                 data=subdata)
APO.PSlm3.1 <- mean(predict(mod1.PSlm3,subdata.alldiab))
APO.PSlm3.0 <- mean(predict(mod1.PSlm3,subdata.nodiab))
APO.PSlm3.1 - APO.PSlm3.0
## Obtain some plots for these models
patient1<-predict(mod1.PSlm3,subdata.alldiab)
patient0<-predict(mod1.PSlm3,subdata.nodiab)
ll=sort(sample(1:dim(subdata)[1],50))
plot(ll,patient1[ll], xlim = c(1,27000), pch="", lwd=1,ylim=c(0,13000),
     ylab = "Total Expenditure",xlab = "Participants")
lines(smooth.spline(ll,patient1[ll],df=12), col=2)
lines(smooth.spline(ll,patient0[ll],df=12), col=3)
title("Total Expenditure: Diabetic vs Non-Diabetic",cex=0.3)
legend("topright",legend = c("Diabetic","Non-Diabetic"),
       col = 2:3,pch = 5)

# IPW estimator
ps.lr.weight <- Diabetes/ps.lr + (1-Diabetes)/(1-ps.lr)
mean(Diabetes*Y*ps.lr.weight) -mean((1-Diabetes)*Y*ps.lr.weight)
#coef(lm(Y ~ Diabetes, weights = ps.lr.weight))[2]


## See outcome_adap_lasso.R file implementation of outcome adaptive lasso

############################################################################
########### Using the multiple generated datasets #########################

library(splines); 
library(Matching)
library(lqa) # version 1.0-3
library(MASS) # version 3.3.1
##############
### Define function for outcome adaptive lasso, ATE estimates, and the wAMD,
expit = function(x){ 
  pr = ( exp(x) / (1+exp(x)) ) 
  return(pr)
}
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}
create_weights = function(fp,fA,fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

####
micedata<-tempdatamice #storedata dataset (list of dataframes)
nreps<-100
res.mat<-matrix(0,nrow=nreps,ncol=4)

# Outcome Regression estimator
for (irep in 1:nreps) { 
  subdata<-complete(tempdatamice,irep)
  Diabetes<-subdata$DIABDX
  # lets build the propensity score model
  varsPS<-names(subdata)[4:19]
  PS.form = formula(paste("DIABDX~AGE13X+",paste(varsPS,collapse="+")))
  ps.mod <- glm(PS.form, data=subdata,family="binomial")
  ps.lr <- predict(ps.mod,type="response")
  summary(ps.lr)
  ps0<-ps.lr[Diabetes==0]
  ps1<-ps.lr[Diabetes==1]
  quints <- c(0,quantile(ps.lr,seq(0.30,1,.16)))
  
  subdata.alldiab <- subdata
  subdata.alldiab$DIABDX <- 1
  subdata.nodiab <- subdata
  subdata.nodiab$DIABDX <- 0
  vars.list<-names(subdata)[3:19]
  # outcome model with interaction terms
  vars.list.int<-c(vars.list,paste("DIABDX:",paste(vars.list,collapse=":")))
  out.form.int = formula(paste("TOTEXP13~DIABDX+",paste(vars.list.int,collapse="+")))
  mod1.lmX.int <- lm(out.form.int,data=subdata)
  APO.lmX.1.int <- mean(predict(mod1.lmX.int,subdata.alldiab))
  APO.lmX.0.int <- mean(predict(mod1.lmX.int,subdata.nodiab))
  res.mat[irep,1]<-APO.lmX.1.int - APO.lmX.0.int
}

#### Ps matching is too slow, so we only consider a complete case analysis

## Inverse Probability Weighting method
for (irep in 1:nreps) {  
  subdata<-complete(tempdatamice,irep)
  Diabetes<-subdata$DIABDX # lets build the propensity score model
  varsPS<-names(subdata)[4:19]
  PS.form = formula(paste("DIABDX~AGE13X+",paste(varsPS,collapse="+")))
  ps.mod <- glm(PS.form, data=subdata,family="binomial")
  ps.lr <- predict(ps.mod,type="response")
  
  subdata.alldiab <- subdata
  subdata.alldiab$DIABDX <- 1
  subdata.nodiab <- subdata
  subdata.nodiab$DIABDX <- 0
  vars.list<-names(subdata)[3:19]
  subdata$ps.lr<-ps.lr
  Y<-subdata$TOTEXP13
  
  # IPW estimator for ATE
  ps.lr.weight <- Diabetes/ps.lr + (1-Diabetes)/(1-ps.lr)
  res.mat[irep,3]<-mean(Diabetes*Y*ps.lr.weight) -mean((1-Diabetes)*Y*ps.lr.weight)
  #res.mat[irep,3]<-coef(lm(Y ~ Diabetes, weights = ps.lr.weight))[2]
}

## Propensity Score regression to compute ATE
for (irep in 1:nreps) { 
  subdata<-complete(tempdatamice,irep)
  Diabetes<-subdata$DIABDX
  # lets build the propensity score model
  varsPS<-names(subdata)[4:19]
  PS.form = formula(paste("DIABDX~AGE13X+",paste(varsPS,collapse="+")))
  ps.mod <- glm(PS.form, data=subdata,family="binomial")
  ps.lr <- predict(ps.mod,type="response")
  
  subdata$ps.lr<-ps.lr
  
  ## PS model is included a linear term in the regression model
  mod1.PSlm3 <- lm(TOTEXP13~DIABDX+bs(ps.lr,df=4), data=subdata)
  APO.PSlm3.1 <- mean(predict(mod1.PSlm3,subdata.alldiab))
  APO.PSlm3.0 <- mean(predict(mod1.PSlm3,subdata.nodiab))
  res.mat[irep,4]<-APO.PSlm3.1 - APO.PSlm3.0
}


for(irep in 1:nreps){ ## Outcome adaptive lasso-- ATE
  Data<-complete(tempdatamice,irep)
  Data<-Data[sample(1:dim(Data)[1],6000),]
  names(Data)<-c("Y","A",names(subdata)[3:19])
  n<-dim(Data)[1]
  var.list<-names(subdata)[3:19];  bA = 0
  
  # set vector of possible lambda's to try
  lambda_vec =  c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25,0.5)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  betaXY = coef(lm.Y)[var.list] 
  ###### Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  ##
  #####  Run outcome adaptive lasso for each lambda value 
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
    
    ### run outcome-adaptive lasso model with appropriate penalty
    logit_oal = lqa( w.full.form, data=Data, penalty=oal_pen, family=binomial(logit) )
    # generate propensity score
    Data[,paste("f.pA",lil,sep="")] = predict(logit_oal)$mu.new
    # save propensity score coefficients
    coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
    # create inverse probability of treatment weights
    Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
  }
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  res.mat[irep,2]<-ATE[tt]
}

par(mar=c(5,5,3,3))
nvec<-c("OR","OA-lass0","Ps-IPW","Ps-Reg")
boxplot(res.mat, ylab = "ATE", names=nvec, main = "ATE by estimator", 
        col = "gray", ylim=c(2500,4500))
abline(h = median(res.mat[,1]),col=2)
median(res.mat[,2])

DirATE<-median(res.mat); nsamp<-dim(subdata)[1]
# Bias and Variance considerations
A.T.E<-round(apply(res.mat,2,median),2)
Bias<-(apply(res.mat,2,mean)-DirATE)* sqrt(nsamp) #Bias
Variance<-apply(res.mat,2,var)*nsamp #Variance
MSE<-apply((res.mat-DirATE)^2,2,mean)#*nsamp #MSE
labs<-c("ATE","Bias", "Variance", "MSE")
specificationname<-nvec
data.frame(specificationname,A.T.E,Bias,Variance,MSE)

Variance^(1/2)


