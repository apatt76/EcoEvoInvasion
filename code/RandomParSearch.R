
library(here)
source(here("code/SolverFunctions.R"))


# Goal:
# To do a random parameter search with different native species scenarios (with theta=m)
# Also add some random parameter search cases with theta>m and theta<m

setID<-34981
set.seed(setID)
ntrials=10000
modellistFeb16par5<-list()
modellistNoEvoFeb16par5<-list()



# ORIGINAL PARAMETERS from previous submission
# (not using again, just have as reference)
parameters0<-c(m=runif(1,0,5), n=2, b=runif(1,0,5), u=1, A=runif(1,0,5), H=runif(1,0,1), c1=runif(1,0,1),
               delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
               tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.25,
               rf=1, rv=1)
#initial conditions- invasive pop should be small, native pop should be bigger- test for carrying capacity alone
# invasive trait should be bigger, but only by a little
state0=c(Nf=0.1, Nv=1, xbar=0.4,ybar=0.1) 


# This code runs for a lot of reps, expect it to take hours to a day or two, depending on machine speed
# Requires manually selecting parameters=parameters# to run that scenario 
start=Sys.time()
for(i in 1:ntrials) {
  # Okay, so first I want to randomly sample, and note what values i'm using
  
  # Parameter definition
  # Set 1: same but lower A AND lower n for this and all subsequent
  parameters1<-c(m=runif(1,0,5), n=1, b=runif(1,0,5), u=1, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.25,
                 rf=1, rv=1)
  
  parameters1[which(names(parameters1)=="thetaC")]<-parameters1[1]
  parameters1[which(names(parameters1)=="thetaH")]<-parameters1[2]
  
  # Set 2: Effect of m: decrease m native (to 0.5)
  parameters2<-c(m=runif(1,0,5), n=0.5, b=runif(1,0,5), u=1, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.25,
                 rf=1, rv=1)
  
  parameters2[which(names(parameters2)=="thetaC")]<-parameters2[1]
  parameters2[which(names(parameters2)=="thetaH")]<-parameters2[2]
  
  # Set 3: Make flatter: increase tau native (to 1) and decrease b native (to 0.5)
  parameters3<-c(m=runif(1,0,5), n=1, b=runif(1,0,5), u=0.5, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=1, sigmax=runif(1,0,2), sigmay=0.25,
                 rf=1, rv=1)
  
  parameters3[which(names(parameters3)=="thetaC")]<-parameters3[1]
  parameters3[which(names(parameters3)=="thetaH")]<-parameters3[2]
  
  # Set 4: Consider variation: increase sigma (to 0.8)
  parameters4<-c(m=runif(1,0,5), n=1, b=runif(1,0,5), u=1, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.8,
                 rf=1, rv=1)
  
  parameters4[which(names(parameters4)=="thetaC")]<-parameters4[1]
  parameters4[which(names(parameters4)=="thetaH")]<-parameters4[2]
  
  # Set 5: conspecific optimum < heterospecific optimum
  parameters5<-c(m=runif(1,0,5), n=1, b=runif(1,0,5), u=1, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.25,
                 rf=1, rv=1)
  
  parameters5[which(names(parameters5)=="thetaC")]<-(parameters5[1]-runif(1,0,parameters5[1]))
  parameters5[which(names(parameters5)=="thetaH")]<-(parameters5[2]-runif(1,0,parameters5[1]))
  
  # Set 6: conspecific optimum > heterospecific optimum
  
  parameters6<-c(m=runif(1,0,5), n=1, b=runif(1,0,5), u=1, A=runif(1,0,1), H=runif(1,0,1), c1=runif(1,0,1),
                 delta=runif(1,0,0.5), thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                 tauInter1=runif(1,0,2), tauInter2=0.412715, sigmax=runif(1,0,2), sigmay=0.25,
                 rf=1, rv=1)
  
  parameters6[which(names(parameters6)=="thetaC")]<-(parameters6[1]+runif(1,0,parameters6[1]))
  parameters6[which(names(parameters6)=="thetaH")]<-(parameters6[2]+runif(1,0,parameters6[1]))
  
  #############################################
  # Parameter choice
  parameters=parameters5
  
  #############################################
  # Solve for initial conditions based on steady state for native species!
  res=solveAdaptive(parameters=parameters,simlengthmin = 500,func=GMspecial)
  
  resNoEvo=solveAdaptive(parameters=parameters,simlengthmin = 500,func=GMspecialNoEvo)
  
  output<-c(parameters,res[[1]],res[[2]][5,])
  outputNoEvo<-c(parameters,resNoEvo[[1]],resNoEvo[[2]][5,])
  
  modellistFeb16par5[[i]]<-output
  modellistNoEvoFeb16par5[[i]]<-outputNoEvo
  
  
  print(paste("trial",i,"done"))
  print(paste("time elapsed:",Sys.time()-start))
  
}

totRunTime=Sys.time()-start
totRunTime

saveRDS(modellistFeb16par5,here("output/modellistFeb16par5.R"))
saveRDS(modellistNoEvoFeb16par5,here("output/modellistNoEvoFeb16par5.R"))

