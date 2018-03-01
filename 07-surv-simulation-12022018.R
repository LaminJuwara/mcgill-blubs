###########################################################################
############     Simulating censored survival data   ######################
rm(list = ls())
library(survival)

n = 1000
beta1 = 1.5
beta2 = -2
x1 = rnorm(n,mean = 3,sd = .2)
x2 = rnorm(n,mean = 3.5, sd = .6)
lambdaT = .002  # baseline hazard
lambdaC = .012  # hazard of censoring

t = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) # true event time
C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
time = pmin(t,C)  #observed time is min of censored and true time
event<-ifelse(time==t,1,0)
survdata<-data.frame(time=time, event=event, x1 = x1, x2 = x2)
survdata[1:20,]
table(event)

summary(event)


##################  simulating r repetitions of the survival datasets ########
r<-100   # the number of repititions
simdata<-list()
for(i in 1:r){
  x1 = rnorm(n,mean = 3,sd = .2)
  x2 = rnorm(n,mean = 3.5, sd = .6)
  # baseline hazard and censoring hazard are controlled above
  t = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) # true event time
  C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
  time = pmin(t,C)  #observed time is min of censored and true time
  event<-ifelse(time==t,1,0)
  survdata<-data.frame(time=time, event=event, x1 = x1, x2 = x2)
  simdata[[i]]<-survdata 
  }

head(simdata[[1]]) # preview of the first dataset | 1000 datasets in the list
head(simdata[[r]]) # the last dataset

####################################################################################
## Estimates for the unpooled values from the cox proportion regression model

unpoolcox<-function(r,simdata){
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(i in 1:r){
    survdata<-simdata[[i]]
    ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
    beta0vec[i]<-ress$coefficients[1,1]
    SEbeta0vec[i]<-ress$coefficients[1,3]
    beta1vec[i]<-ress$coefficients[2,1]
    SEbeta1vec[i]<-ress$coefficients[2,3]
  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec), median(beta1vec) ),
                     SE = c(median(SEbeta0vec), median(SEbeta1vec) ) )
  return(result)
}

#unpoolcox(r,simdata)


##############################################################################################
## Now lets establish that this holds for the reduced form: conditional logistic regression

unpoolclogit<-function(r,simdata){
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event;   
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength)),
                         event = c(survdata$event,rep(NA,rowlength)),
                         x1 = c(survdata$x1,rep(NA,rowlength)),
                         x2 = c(survdata$x2,rep(NA,rowlength)),
                         risksetid = c(survdata$risksetid,rep(NA,rowlength)) )
    
    # create a list of the event ids starting with the first event that took place
    #for(i in 1:length(time[event==1])){
    for(i in 1:sum(event)){
      eventid[i]<-which(time==sort(time[event==1], decreasing = FALSE)[i])
    }
    summary(eventid)
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]])

      if(length(l1)>=1 ){
        resurvdata$risksetid[eventid[i]]<-i  # risksetid for the case | this is in the new dataset
        index<-sample(x = l1,size = 1, replace = T)  # randomly select one of the controls to match on; note this could be a future case
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-1 case-control setting that is used in the next line
        resurvdata[(1000+i),]<-survdata[index,] # copied control in the new dataset
        resurvdata$event[(1000+i)]<-0  # event is fixed as a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-1 matched sets
    survdata<-NA
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]

    
    ress<-summary(clogit(formula = event~x1+x2+strata(risksetid), data = survdata))
    beta0vec[z]<-ress$coefficients[1,1]
    SEbeta0vec[z]<-ress$coefficients[1,3]
    beta1vec[z]<-ress$coefficients[2,1]
    SEbeta1vec[z]<-ress$coefficients[2,3]

  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  return(result)
}

#unpoolclogit(r, simdata) # this estimate is practically the same as the estimate
                          # for the unpooled estimates from clogit model



###################################################################################
### pools of size 2 in a 1-1 case/control matched setting 
###################################################################################


s=2  # the pool size
# r is the number of repition from the simulated dataset
# simdata is a list of datasets of r repetitions

pool2<-function(s,r,simdata){ # Need to add a simulated data component
  # pools of size 2
  poolsize = s
  r = r
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event   
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength)),
                           event = c(survdata$event,rep(NA,rowlength)),
                           x1 = c(survdata$x1,rep(NA,rowlength)),
                           x2 = c(survdata$x2,rep(NA,rowlength)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength)) )
    
    # create a list of the event ids starting with the first event that is taking place
    #for(i in 1:length(time[event==1])){
    for(i in 1:sum(event)){
      eventid[i]<-which(time==sort(time[event==1], decreasing = FALSE)[i])
    }
    summary(eventid)
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]])
      
      if(length(l1)>=1 ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = 1, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-1 case-control setting
        resurvdata[(1000+i),]<-survdata[index,]
        resurvdata$event[(1000+i)]<-0        # the selected person could be a future case. make sure it is a control
      }
      # if the number of elements in the riskset is less than 1, then there must be no element in the riskset
      
    }
    #head(resurvdata)  # A dataset of 1-1 matched sets
    survdata<-NA
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)

    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count-1),2) # identifier for picking a new riskset to create a pool e.g. case 1&2 form pool 1, 3&4 form pool 2, etc
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|survdata$risksetid==(poolmarker[i]+1)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),] 
    
    # case counts and control counts
    control.count<-floor(length(which(survdata$event==0))/poolsize)*poolsize              # the controls
    
    # 
    survdata=NA
    poolsurvdata$sumx1<-NA   # summing the x1 variables
    for(i in 1:case.count){
      poolsurvdata$sumx1[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx1[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    poolsurvdata$sumx2<-NA   # summing the x2 variables
    for(i in 1:case.count){
      poolsurvdata$sumx2[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx2[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    # now lets try the conditional from the clogit function
    
    ll=dim(poolsurvdata)[1]/poolsize
    # first we need a subset of the dataset with sums
    survdata<-data.frame(event=rep(NA,ll), 
                         sumx1=rep(NA,ll ),
                         sumx2=rep(NA,ll ),
                         id=rep(NA,ll )  )
    index0<-seq(1,ll-1,2)
    index1<-seq(2,ll,2)
    for(j in 1:(ll/poolsize)){
      survdata$event[index0[j]]=0  # fix the event = 0 and then ...
      survdata$sumx1[index0[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index0[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$id[index0[j]]<-j
      survdata$event[index1[j]]<-1 # fix the event = 1 and then ...
      survdata$sumx1[index1[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index1[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$id[index1[j]]<-j
    }
    
    
    Rg=case.count/control.count
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id),
                         #offset = rep(log(Rg),length(survsumx1)),
                         data = survdata))
    #ssp2$coefficients
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]

  }
  
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  return(result)
}

#pool2(s,r,simdata)

####  compiling the results

p0cox<-unpoolcox(r, simdata)
p0clogit<-unpoolclogit(r, simdata)
p2=pool2(s,r,simdata)
data.frame(parameter=c("beta1","SE(beta 1)","beta2","SE(beta 2)"),
           Fixed = c(beta1,"",beta2,""),
           unpooled.cox=c(p0cox$estimate[1],p0cox$SE[1] ,p0cox$estimate[2],p0cox$SE[2] ),
           unpooled.cond.reg = c(p0clogit$estimate[1],p0clogit$SE[1] ,
                                 p0clogit$estimate[2],p0clogit$SE[2]),
           pool2 = c(p2$estimate[1],p2$SE[1] ,p2$estimate[2],p2$SE[2])
           )


######################################################################################################
#### lets try different pool sizes for 1:m matched sets # this is likely to produce less bias estimates
######################################################################################################

multi.cc.pool<-function(s,r,simdata){
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r

  
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    m = floor(sum(event==0)/sum(event==1)) # determine the number of eligible controls for each case
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    
    # create a list of the event ids starting with the first event that is taking place
    for(i in 1:sum(event)){
      eventid[i]<-which(time==sort(time[event==1], decreasing = FALSE)[i])
    }
    #summary(eventid)
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-m matched set
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)
    
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count-1),2) # identifier for picking a new riskset to create a pool e.g. case 1and2 form pool 1, 3and4 form pool 2, etc
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|survdata$risksetid==(poolmarker[i]+1)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    poolsurvdata$sumx1<-NA   # summing the x1 variables
    for(i in 1:case.count){
      poolsurvdata$sumx1[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx1[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    poolsurvdata$sumx2<-NA   # summing the x2 variables
    for(i in 1:case.count){
      poolsurvdata$sumx2[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx2[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    # now lets try the conditional from the clogit function
    
    ll=dim(poolsurvdata)[1]/poolsize
    # first we need a subset of the dataset with sums
    survdata<-data.frame(event=rep(NA,ll), 
                         sumx1=rep(NA,ll ),
                         sumx2=rep(NA,ll ),
                         id=rep(NA,ll )  )
    index0<-seq(1,ll-1,2)
    index1<-seq(2,ll,2)
    for(j in 1:(ll/poolsize)){
      survdata$event[index0[j]]=0  # fix the event = 0 and then ...
      survdata$sumx1[index0[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index0[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$id[index0[j]]<-j
      survdata$event[index1[j]]<-1 # fix the event = 1 and then ...
      survdata$sumx1[index1[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index1[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$id[index1[j]]<-j
    }
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    ssp2$coefficients
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    
  }

  result<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  return(result)
}

# multi.cc.pool(s,r, simdata)
mmp2<-multi.cc.pool(s,r, simdata)
mmp2



#########
# We need to rewrite this code to accomdate a 1:m case-control settings and 
# a generic pool size of >2


##############################################################################
##############################################################################

###############  Simulation using the SMART dataset ##########################
                
##############################################################################
require(survival)
require(hdnom)
smartdata<-(smarto)

head(smartdata)
#setwd(/home/laminj/Desktop/directed_reading_course/simulations/")
#write.csv(smartdata, file = "smartdataset")

# some summary statistics
summary(smartdata)


# The biomakers to consider in this analysis 
#  Markers of atherosclerosis

# HOMOC - Homocyste; 463 missing values
# GLUT - Glutamine; 19 missing values
# CREAT - Creatinine clearance, in mL/min; 17 missing values
# ALBUMIN - Albumin (no, micro, macro); 28 missing values
# IMT - Intima media thickness, in mm; 98 missing values
# STENOSIS - Carotid artery stenosis > 50%; 93 missing values


# Tableone for all the variable
require(tableone)
str(smartdata)
names(smartdata)
vars1<-c("TEVENT","SEX","AGE","DIABETES","CEREBRAL","CARDIAC","AAA","PERIPH",
        "STENOSIS","SYSTBP","DIASTBP","SYSTH","DIASTH","LENGTH","WEIGHT","BMI","CHOL","HDL",
        "LDL","TRIG","HOMOC","GLUT","CREAT","IMT","ALBUMIN","SMOKING","PACKYRS","ALCOHOL")
Table1<-CreateTableOne(vars = vars1, strata = "EVENT", data = smartdata, test = TRUE)
print(Table1, smd=TRUE)

# first, lets select a subset that is useful for our analysis
subsmart<- smartdata[,c("TEVENT","EVENT","SEX","AGE","DIABETES",
                        "STENOSIS","HOMOC","GLUT","CREAT","IMT")]
#subsmart$EVENT<-as.factor(subsmart$EVENT)
subsmart$SEX<-as.factor(subsmart$SEX)
str(subsmart)
summary(subsmart)
# lets check the amount of missing data; we will exclude HOMOC-- 463 missing

names(subsmart)
subsmart<-subsmart[,c("TEVENT","EVENT","HOMOC","CREAT","IMT")]
summary(subsmart)
#   
hist(subsmart$TEVENT)
hist(subsmart$TEVENT[subsmart$EVENT==1])

plot(density(subsmart$TEVENT) )
#scatter.smooth(subsmart$HOMOC)
#plot(subsmart$CREAT)

pairs(~TEVENT+HOMOC+CREAT+IMT, data = subsmart) # A pairwise plot of the continuous variables
# The dataset doesn't look too bad


############### Simalation calculations ###########################
### lets subset this dataset and try it with the available code
data1<-subsmart[,c("TEVENT","EVENT","CREAT","IMT")]
names(data1)<-c("time","event","x1","x2")
dataf<-list()
dataf[[1]]<-data1[!is.na(data1$x1&data1$x2),]
length(dataf[1])


## lets compile all the functions
simdata<-dataf
r<-length(dataf[1])
sum(dataf[[1]]$event)  # this is the number of events


####################################################################################
## Estimates for the unpooled values from the cox proportion regression model
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  # 

unpoolcox<-function(r,simdata){
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(i in 1:r){
    survdata<-simdata[[i]]
    ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
    beta0vec[i]<-ress$coefficients[1,1]
    SEbeta0vec[i]<-ress$coefficients[1,3]
    beta1vec[i]<-ress$coefficients[2,1]
    SEbeta1vec[i]<-ress$coefficients[2,3]
  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec), median(beta1vec) ),
                     SE = c(median(SEbeta0vec), median(SEbeta1vec) ) )
  return(result)
}
p0cox<-unpoolcox(r, simdata) 
p0cox
# Equivalently
survdata<-simdata[[1]]
ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
ress$coefficients
#unpoolcox(r,simdata)


##############################################################################################
## Now lets establish that this holds for the reduced form:
#  conditional logistic regression

unpoolclogit<-function(r,simdata){  # Testing this model
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(z in 1:r){
    survdata<-simdata[[z]]
    lendata<-dim(survdata)[1]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event;   
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength)),
                           event = c(survdata$event,rep(NA,rowlength)),
                           x1 = c(survdata$x1,rep(NA,rowlength)),
                           x2 = c(survdata$x2,rep(NA,rowlength)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    # a list of the event index
    tempfoo<-lfoo$ix
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    # resurvdata[eventid,"time"]  # list of times starting with the smallest time, verified
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]])
      if(length(l1)>=1 ){
        resurvdata$risksetid[eventid[i]]<-i  # risksetid for the case | this is in the new dataset
        index<-sample(x = l1,size = 1, replace = T)  # randomly select one of the controls to match on; note this could be a future case
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-1 case-control setting that is used in the next line
        resurvdata[(lendata+i),]<-survdata[index,] # copied control in the new dataset
        resurvdata$event[(lendata+i)]<-0  # event is fixed as a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-1 matched sets
    survdata<-NA
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    
    ress<-summary(clogit(formula = event~x1+x2+strata(risksetid), data = survdata))
    beta0vec[z]<-ress$coefficients[1,1]
    SEbeta0vec[z]<-ress$coefficients[1,3]
    beta1vec[z]<-ress$coefficients[2,1]
    SEbeta1vec[z]<-ress$coefficients[2,3]
    
  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  return(result)
}

unpoolclogit(r, simdata) # slightly different results

# lets compile the results here
p0clogit<-unpoolclogit(r,simdata)
data.frame(parameter=c("beta1","SE(beta 1)","beta2","SE(beta 2)"),
           unpooled.cox=c(p0cox$estimate[1],p0cox$SE[1] ,p0cox$estimate[2],p0cox$SE[2] ),
           unpooled.cond.reg = c(p0clogit$estimate[1],p0clogit$SE[1] ,
                                 p0clogit$estimate[2],p0clogit$SE[2]))

## Reworking the equivalent conditional logistic form   #################################

## Lets do it for pools of size 2 in a 1:1 case-control setting


s=2  # the pool size
# r is the number of repition from the simulated dataset
# simdata is a list of datasets of r repititions

pool2<-function(s,r,simdata){ # Need to add a simulated data component
  # pools of size 2
  poolsize = s
  r = r
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(z in 1:r){
    survdata<-simdata[[z]]
    lendata<-dim(survdata)[1]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event   
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength)),
                           event = c(survdata$event,rep(NA,rowlength)),
                           x1 = c(survdata$x1,rep(NA,rowlength)),
                           x2 = c(survdata$x2,rep(NA,rowlength)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    # a list of the event index
    tempfoo<-lfoo$ix
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    # resurvdata[eventid,"time"]  # list of times starting with the smallest time, verified
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]])
      
      if(length(l1)>=1 ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = 1, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-1 case-control setting
        resurvdata[(lendata+i),]<-survdata[index,]
        resurvdata$event[(lendata+i)]<-0        # the selected person could be a future case. make sure it is a control
      }
      # if the number of elements in the riskset is less than 1, then there must be no element in the riskset
      
    }
    #head(resurvdata)  # A dataset of 1-1 matched sets
    survdata<-NA
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)
    
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count-1),2) # identifier for picking a new riskset to create a pool e.g. case 1&2 form pool 1, 3&4 form pool 2, etc
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|survdata$risksetid==(poolmarker[i]+1)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),] 
    
    # case counts and control counts
    control.count<-floor(length(which(survdata$event==0))/poolsize)*poolsize              # the controls
    
    # 
    survdata=NA
    poolsurvdata$sumx1<-NA   # summing the x1 variables
    for(i in 1:case.count){
      poolsurvdata$sumx1[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx1[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x1[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    poolsurvdata$sumx2<-NA   # summing the x2 variables
    for(i in 1:case.count){
      poolsurvdata$sumx2[poolsurvdata$event==1&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==1&poolsurvdata$poolid==i])
      poolsurvdata$sumx2[poolsurvdata$event==0&poolsurvdata$poolid==i]<-sum(poolsurvdata$x2[poolsurvdata$event==0&poolsurvdata$poolid==i])
    }
    
    # now lets try the conditional from the clogit function
    
    ll=dim(poolsurvdata)[1]/poolsize
    # first we need a subset of the dataset with sums
    survdata<-data.frame(event=rep(NA,ll), 
                         sumx1=rep(NA,ll ),
                         sumx2=rep(NA,ll ),
                         id=rep(NA,ll )  )
    index0<-seq(1,ll-1,2)
    index1<-seq(2,ll,2)
    for(j in 1:(ll/poolsize)){
      survdata$event[index0[j]]=0  # fix the event = 0 and then ...
      survdata$sumx1[index0[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index0[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==0&poolsurvdata$poolid==j)][1]
      survdata$id[index0[j]]<-j
      survdata$event[index1[j]]<-1 # fix the event = 1 and then ...
      survdata$sumx1[index1[j]]<-poolsurvdata$sumx1[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$sumx2[index1[j]]<-poolsurvdata$sumx2[(poolsurvdata$event==1&poolsurvdata$poolid==j)][1]
      survdata$id[index1[j]]<-j
    }
    
    
    Rg=case.count/control.count # This bias term is not useful in a cond logistic setting
    
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id),
                         #offset = rep(log(Rg),length(survsumx1)),
                         data = survdata))
    #ssp2$coefficients
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    
  }
  
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  return(result)
}
 pool2(s,r,simdata)

# lets compile all the results
#p0clogit<-unpoolclogit(r,simdata)  # computed above
p2clogit<- pool2(s,r,simdata)

data.frame(parameter=c("beta1","SE(beta 1)","beta2","SE(beta 2)"),
           unpooled.cox=c(p0cox$estimate[1],p0cox$SE[1] ,p0cox$estimate[2],p0cox$SE[2] ),
           unpooled.cond.reg = c(p0clogit$estimate[1],p0clogit$SE[1] ,
                                 p0clogit$estimate[2],p0clogit$SE[2]),
           pool2.cond.reg = c(p2clogit$estimate[1],p2clogit$SE[1] ,
                              p2clogit$estimate[2],p2clogit$SE[2])
           )

######  The estimates are comparable for the different pool sizes





