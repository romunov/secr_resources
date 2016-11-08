# rm(list=ls())
# From https://sites.google.com/site/wild8390/software/simulate/simulation-for-study-design
require(RMark)

#IMBEDDING SIMULATION MODEL IN LOOP TO INVESTIGATE SAMPLE SIZE AND OTHER INPUTS

##FIRST DEFINE SOME FUNCTIONS NEEDED
## FUNCTION TO CREATE CAPTURE HISTORY CHARACTER STRINGS 

pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}

#Simulating data over n_sim_reps (ideally 1000 times)
sim.data<-function(N,p,k)
{
  #simulate capture histories
  y<-array(dim=c(N,k))
  ind<-array(dim=N)
  
  for(i in 1:N)
  {
    y[i,1:k]<-rbinom(k,1,p)
    ind[i]<-sum(y[i,])>0
  }
  #capture history data frame 
  capt.hist<-data.frame(ch=pasty(y[,1:k]),ind=ind)
  capt.hist<-subset(capt.hist,ind==T,select=c(ch))
  #end of function
  return(capt.hist)
}


####FUNCTION TO SIMULATE CAPTURE HISTORIES UNDER SPECIFIED INPUTS

#simulate capture histories from assumed model and estimated parameter values
#simulate data under homogeneous p 
#set up to simulate for specified inputs

sample_sim<-function(n_sim_reps,N,p,k)
{
  #WE CAN DEFINE SOME THINGS OUTSIDE THE LOOP FOR EXAMPLE WE WILL USE THE SAME TIME CONSTANT MODEL EACH SIMU;LLATION
  #Define parameters
  p.dot=list(formula=~1,share=TRUE)
  #SET UP AN 
  #set up an empty data frame to store all the simulation results
  output.data<-data.frame(k=numeric(0),rep=numeric(0),N.hat=numeric(0))
  
  for(r in 1:n_sim_reps)
  {
    #
    #function to create capture histories from simulated data
    capt.hist<-sim.data(N=N,p=p,k=k)
    
    #additional parameters to supress output each simulation
    m0<-mark(capt.hist,model="Closed",model.parameters=list(p=p.dot),silent=T,output=F)
    #pull off just the estimate of N (you can select other parameters if you want)
    N.hat<-get.real(m0,"N")
    #put in data frame
    new.data<-data.frame(k=k,rep=r,N.hat=N.hat)
    #append to output data
    output.data<-rbind(output.data,new.data)
  }
  #create  summay statistics
  #empirical mean, sd, and cv
  N.avg<-mean(output.data$N.hat)
  N.sd<-sd(output.data$N.hat)
  N.cv<-N.sd/N.avg
  #return simulated data and summary stats in list objects
  N.summary<-list(N.avg=N.avg,N.sd=N.sd,N.cv=N.cv)
  
  output=list(N.summary=N.summary,N.output=output.data)
}





#call simulation with specific inputs and save as object
sim.results<-sample_sim(n_sim_reps=100,N=100,p=0.2,k=8)
#suppose you just want cv for N
sim.results$N.summary$N.cv

#you could run this for a series of inputs for k etc -- or put in another loop
#create a data frame to hold the results
summary.data<-data.frame(N=numeric(0),p=numeric(0),k=numeric(0),cv=numeric(0))
for (k in 3:8)
{
  
  sim.results<-sample_sim(n_sim_reps=10,N=100,p=0.2,k)
  new.data<-data.frame(N=100,p=0.2,k=k,cv=sim.results$N.summary$N.cv)
  #append to output data
  summary.data<-rbind(summary.data,new.data)
  
}
#print summary data
summary.data
with(summary.data,plot(k,cv))

##############

#comparison to function based on approximation/ expected value
#seen in week 6
cmr_sample<-function(N,p,k)
{
  ndot<-N*p*k
  var<-N/(1/(1-p)^k-k/(1-p)+k-1)
  cv<-sqrt(var)/N
  return(cv)
}
approx.data<-data.frame(N=numeric(0),p=numeric(0),k=numeric(0),cv=numeric(0))
for (k in 3:8)
{
  cmr_sample(N=100,p=0.2,k=2)
  new.data<-data.frame(N=100,p=0.2,k=k,cv=cmr_sample(N=100,p=0.2,k=k))
  approx.data<-rbind(approx.data,new.data)
  
}
approx.data
with(approx.data,plot(k,cv))

############SIMULATE DATA WITH HETEROGENEITY IN P 
#########################################

############# heterogeneity in p among individuals

require(gtools) #needed for built in logit/ inv logit functions

#Simulating data over n_sim_reps (ideally 1000 times)
#average p with sigma logit scale random effect
sim.het.data<-function(N,p.avg,k,sigma)
{
  #simulate capture histories
  y<-array(dim=c(N,k))
  ind<-array(dim=N)
  
  #n random effects
  p<-inv.logit(rnorm(N,logit(p.avg),sigma))
  
  for(i in 1:N)
  {
    y[i,1:k]<-rbinom(k,1,p[i])
    ind[i]<-sum(y[i,])>0
  }
  #capture history data frame 
  capt.hist<-data.frame(ch=pasty(y[,1:k]),ind=ind)
  capt.hist<-subset(capt.hist,ind==T,select=c(ch))
  #end of function
  return(capt.hist)
}

#NOW W'ELL ESTIMATE N AS IF P IS CONSTANT 
sample_sim.het<-function(n_sim_reps,N,p.avg,k,sigma)
{
  #WE CAN DEFINE SOME THINGS OUTSIDE THE LOOP FOR EXAMPLE WE WILL USE THE SAME TIME CONSTANT MODEL EACH SIMU;LLATION
  #Define parameters
  p.dot=list(formula=~1,share=TRUE)
  #SET UP AN 
  #set up an empty data frame to store all the simulation results
  output.data<-data.frame(k=numeric(0),rep=numeric(0),N.hat=numeric(0))
  
  for(r in 1:n_sim_reps)
  {
    #parametric bootstrap 
    #function to create capture histories from simulated data
    capt.hist<-sim.het.data(N=N,p.avg=p.avg,k=k,sigma=sigma)
    
    #additional parameters to supress output each simulation
    m0<-mark(capt.hist,model="Closed",model.parameters=list(p=p.dot),silent=T,output=F)
    #pull off just the estimate of N (you can select other parameters if you want)
    N.hat<-get.real(m0,"N")
    #put in data frame
    new.data<-data.frame(k=k,rep=r,N.hat=N.hat)
    #append to output data
    output.data<-rbind(output.data,new.data)
  }
  #create  summay statistics
  #empirical mean, sd, and cv
  N.avg<-mean(output.data$N.hat)
  N.sd<-sd(output.data$N.hat)
  N.cv<-N.sd/N.avg
  #return simulated data and summary stats in list objects
  N.summary<-list(N.avg=N.avg,N.sd=N.sd,N.cv=N.cv)
  
  output=list(N.summary=N.summary,N.output=output.data)
}


#results for specific inputs
sim.results<-sample_sim.het(n_sim_reps=10,N=100,p.avg=0.2,k=8,sigma=0.8)
#summary results -- pay attention to the average N estimate!!
sim.results$N.summary


