
#plotting
#install.packages("ggplot2")
library(ggplot2)

#fitting simulated distributions
install.packages("fitdistrplus")
library(fitdistrplus)

install.packages("logspline")
library(logspline)

#hospital processes often involve coordinating care for patients. This simulation
#tests whether a sufficient cost reduction can be achieved with behavioral health
#patients under a new staffing and coordination plan

#simulation patients
members<-1000

#creating patient comorbidities as a shifted poisson. Every patient has atleast 
#one condition

T<-5.1                                # pre-truncation mean of Poisson
U<-runif(members)                     # the uniform sample
t = -log(1 - U*(1 - exp(-T)))         # the "first" event-times
T1<-(T - t)                           # the set of (T-t)


Comorbidities <- rpois(members,T1)+1  # the final shifted Poisson sample

hist(Comorbidities)

#Get actual cost information
costdf<- read.csv('mmbrexp.depr.csv', stringsAsFactors = F)
costdf$year<- as.numeric(format(as.Date(costdf$Begin_Date,format = "%Y-%m-%d"),"%Y"))

#### GET PRESENT YEAR COST DATA SUMMARY #

#split cost by year  
options(scipen=999)
costdf %>%
  group_by(FC_Member_ID, year) %>% 
  summarize(Paid= sum(Paid_Total)) %>% 
  mutate(label= paste("Actual Member Cost:",year)) %>%
  ggplot(aes(x = Paid)) +
  geom_histogram(color = "white") +
  facet_wrap(~label) +h
  xlab("Cost Buckets") +
  theme_bw()

#summarize cost data, grouping by member and year 
cost.PY<- costdf[costdf$year==as.numeric(format(as.Date(Sys.Date(),format = "%Y-%m-%d"),"%Y"))-1,]  
cost.PY<-cost.PY %>% group_by(Member_ID, year) %>% summarize(Paid= sum(Paid_Total)) 

#Summary statistics about the distribution  
descdist(cost.PY$Paid)

#check lognormal
cost.PY$Paid_log<-log(cost.PY$Paid)
hist(cost.PY$Paid_log)

#check normal
fit.norm<- fitdist(cost.PY$Paid,"norm")
plot(fit.norm)
  
#check beta
beta<-cost.PY$Paid/max(cost.PY$Paid)
fit.beta<-fitdist(beta,"beta")

summary(cost.PY$Paid)
  
#check 0 values for the year
qualitycheck<-merge(cost.PY[cost.PY$Paid==0,],costdf)
  
write.csv(qualitycheck,"qualitycheck.csv")

#Simulate cost of treatment in a year
  
#Oddly the data fits a chi sq distribution with less than 1 degree of freedom best
#using a non central chisq with non centrality parameter
#The non-central rchisq is computed as a Poisson mixture central of chi-squares (Johnson et al, 1995, p.436).
  
  
#generate costs as a chisq with .3 degrees of freedom 
costs<-rchisq(members,.5,ncp=0)*40000
  max(costs)
  hist(costs)
  plot(density(costs))
  lines(density(cost.PY$Paid),col="red") #actual costs 

#### SET ENROLLMENT LEVELS #
#Medicaid enrollment needs to be renewed by patients - it is automatic for Medicare
#This leads to substantial churn in the member population

#simulate member months in Plan    
#monthly attrition rate
attr<-0.07                                  #~ 42% retention annually as per Medicaid enrollment dashboad 
enrollment<- rgeom(members,attr)            #assuming a geometric distribution

#visualize enrollment
enrollment  
summary(enrollment)
hist(enrollment)  
plot(density(enrollment))
  
paste("Percent Retention after 1 year: ",
length(enrollment[enrollment>=12])/members)

#Number of interactions with nurse care manager per member per month: 
#Visits: Poisson ~ 10 visits over 12 months
#Calls: Poisson ~ 10.8 visits over 12 months
#SOURCE: TEAM CARE Cost Paper
  
#NMvisits.mthly<- rpois(members,10)/12
#summary(NMvisits.monthly)
  
NMcalls.monthly<- rpois(members,10.8)/12
  
#Number of calls as a Negative Binomial
#modelled as a poisson where the rate = mu or mean of neg binomial, and the dispersion/size/number of successes 
#is allowed to vary to fit the datapoint we have from the study. Essentially added on a gamma distribution term
#for error correction by mixing poissons    
NMvisits.monthly<- rnbinom(members,mu=10,size =11)/12
summary(NMvisits.monthly)

#patient visits by an NM are a function of length of enrollment
totalvisits<-enrollment*NMvisits.monthly
hist(NMvisits.monthly, breaks=50)

#99% over 1 visit in 12 month period according to the study
length(totalvisits[totalvisits>=1])/members

#Costs of interactions with NM
costNMvisit<-54                   #$ for upto 10 visits
costNMcall<- 30                   #$ per call

#set years of follow up
years<-2
months<-years*12

#reformat time to not exceed follow up time
enrollment.period<- apply(cbind(enrollment,months),1,FUN=min)
totalvisits<-round(NMvisits.monthly*enrollment.period,0)
   
#cost incurred on NM visits
totalvisits*costNMvisit
   
# non intervention costs over enrollment period 
costs.ni<-costs*(1+(enrollment.period-12)/12)
   
#additional costs associated with intervention
costs.i<-totalvisits*costNMvisit
   
#percentage cost saving/improvement estimated per member 
curve(dbeta(x, 2, 40))                  # beta distribution for simulating saving percentages somewhere between 2 and 15%
qbeta(.5, 2, 40)                        # median is 4% 

sav.perc<-rbeta(members,2,40)   
plot(density(sav.perc))
   
#savings
savings<-sav.perc*costs.ni -costs.i
mean(savings)

plot(ecdf(savings),breaks=50)
   
#function to wrap all of the steps within a loop
  
runs<-500
attr<-0.07 #~ 42% retention annually as per Medicaid enrollment dashboad 
  
#Costs of interactions with NM
costNMvisit<-54 #$ for upto 10 visits
costNMcall<- 30 #$ per call
  
#Cost of interactions with physician
#$100 per member - replace with fixed cost for psychologist
#$140 per hour 3-4 weeks per month
  
#set years of follow up
years<-2
months<-years*12
  
runloop<- function(runs, members) {
  costMX<-replicate(runs,rchisq(members,.5,ncp=0)*40000)
  enrollMX<- replicate(runs,rgeom(members,attr))
  NMvisitsMX<- replicate(runs,rnbinom(members,mu=10,size =11)/12)
  NMcallssMX<-replicate(runs,rpois(members,10.8)/12)
  newEnrollMX<- pmin(enrollMX,replicate(runs,months))
  totVisitsMX<-newEnrollMX*NMvisitsMX
  totCallsMX<-newEnrollMX*NMcallssMX
  costIntMX<-totCallsMX*costNMcall + totVisitsMX*costNMvisit 
  costNIntMX<-((newEnrollMX-12)/12)*costMX+costMX
  savingsMX<-replicate(runs,rbeta(members,2,40))
  netsaveMX<-costNIntMX*savingsMX -costIntMX
  netsaving<-data.frame(Mean=apply(netsaveMX,2,mean))
  netsaving<-netsaving[order(netsaving$Mean),]
  return(netsaving)
    
  }

plot(x<-ecdf(netsaving<-runloop(runs<-1000,members<-12000)),breaks=50 )
  
y<-x(sort(netsaving))

vertices<- data.frame(cbind(x=c(netsaving,rev(netsaving)),y=abs(c(y,rep(0,length(y))))))
  
  
#polygon(vertices$x,vertices$y,col="light gray")
  
#95% of all possible outcomes
plot(x,breaks=50, main=paste("Average Net Income Distribution: ", runs, " simulations"), sub=paste(members," members with 95% of all possible outcomes"),xlab="Net Income Per Member",ylab="Percent Outcomes") + polygon(c(vertices$x[vertices$y>.05 &vertices$y<=0.95],rev(vertices$x[vertices$y>.05 &vertices$y<=0.95])),c(y[y>.05 & y<=0.95],rep(0,length(y[y>.05 & y<=0.95]))),col="light blue")
   

