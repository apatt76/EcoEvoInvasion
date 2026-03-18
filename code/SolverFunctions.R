# Contains a list of functions for solving differential equations for the eco-evo invasion system
# Functions also include functions to run batch solving for different parameter sweeps,
# functions to manipulate data related to batch solving,
# and functions related to visualization of the results

library(deSolve)


GMspecial<-function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    
    c=A
    
    # calculate alphas
    alphavv=1+c1-delta-(1-2*delta)*(tauIntra1/sqrt(tauIntra1^2+sigmax^2))*exp(-(xbar-thetaC)^2/(2*tauIntra1^2+2*sigmax^2))
    
    star1Exp= -(xbar^2+ybar^2+m^2-2*b*m*sigmay^2-2*xbar*(ybar+m-b*sigmay^2)-b^2*sigmay^2*sigmax^2-b^2*sigmay^2*tauInter1^2+2*ybar*(m+b*(sigmax^2+tauInter1^2)))/(2*(sigmay^2+sigmax^2+tauInter1^2)) #E1
    star2Exp= -(xbar^2+ybar^2-2*ybar*n+n^2+2*ybar*sigmax^2*u-2*n*sigmax^2*u-tauInter2^2*sigmax^2*u^2-sigmay^2*sigmax^2*u^2+2*xbar*(-ybar+n+(tauInter2^2+sigmay^2)*u))/(2*(tauInter2^2+sigmay^2+sigmax^2)) #E2 
    star1Coeff= sqrt(tauInter1^2/(sigmay^2+sigmax^2+tauInter1^2)) #C1
    star2Coeff= sqrt(tauInter2^2/(sigmay^2+sigmax^2+tauInter2^2)) #C2
    
    alphavf=c+A*(star1Coeff*exp(star1Exp)-star2Coeff*exp(star2Exp))
    alphafv=c+A*(star2Coeff*exp(star2Exp)-star1Coeff*exp(star1Exp))
    alphaff=1+c1-delta-(1-2*delta)*(tauIntra2/sqrt(tauIntra2^2+sigmay^2))*exp(-(ybar-thetaH)^2/(2*tauIntra2^2+2*sigmay^2))
    
    # calculate partials
    derivCoef1xbar= -(2*xbar-2*(ybar+m-b*sigmay^2))/(2*(sigmay^2+sigmax^2+tauInter1^2)) #D11
    derivCoef1ybar= -(-2*xbar+2*ybar+2*(m+b*(sigmax^2+tauInter1^2)))/(2*(sigmay^2+sigmax^2+tauInter1^2)) #D12
    derivCoef2xbar= -(2*xbar+2*(-ybar+n+(tauInter2^2+sigmay^2)*u))/(2*(sigmay^2+sigmax^2+tauInter2^2)) #D21
    derivCoef2ybar= -(-2*xbar+2*ybar-2*n+2*sigmax^2*u)/(2*(sigmay^2+sigmax^2+tauInter2^2)) #D22
    
    partialx=rf*(-Nv*(derivCoef1xbar*star1Coeff*exp(star1Exp)-derivCoef2xbar*star2Coeff*exp(star2Exp))-Nf*(-1)*(1-2*delta)*(tauIntra1/(sqrt(tauIntra1^2+sigmax^2)))*(-1/(2*tauIntra1^2+2*sigmax^2))*2*(xbar-thetaC)*exp(-(xbar-thetaC)^2/(2*sigmax^2+2*tauIntra1^2)))
    partialy=rv*(-Nf*(derivCoef2ybar*star2Coeff*exp(star2Exp)-derivCoef1ybar*star1Coeff*exp(star1Exp))-Nv*(-1)*(1-2*delta)*(tauIntra2/(sqrt(tauIntra2^2+sigmay^2)))*(-1/(2*tauIntra2^2+2*sigmay^2))*2*(ybar-thetaH)*exp(-(ybar-thetaH)^2/(2*sigmay^2+2*tauIntra2^2)))
    
    
    # rates of change
    dNf=rf*Nf*(1-alphafv*Nv-alphaff*Nf)
    dNv=rv*Nv*(1-alphavf*Nf-alphavv*Nv)
    dxbar= H*sigmax^2*partialx
    dybar= H*sigmay^2*partialy
    
    # zero population provision- no trait changes
    if(Nf==0){
      dxbar=0
    }
    if(Nv==0){
      dybar=0
    }
    
    # return the rates of change as a list
    list(c(dNf,dNv,dxbar,dybar))
    
  }) # end the with environment statement
}

########################################################################
########################################################################
########################################################################
# Same function specialize for no evolution
########################################################################
########################################################################
########################################################################


GMspecialNoEvo<-function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    
    c=A
    
    # calculate alphas
    alphavv=1+c1-delta-(1-2*delta)*(tauIntra1/sqrt(tauIntra1^2+sigmax^2))*exp(-(xbar-thetaC)^2/(2*tauIntra1^2+2*sigmax^2))
    
    star1Exp= -(xbar^2+ybar^2+m^2-2*b*m*sigmay^2-2*xbar*(ybar+m-b*sigmay^2)-b^2*sigmay^2*sigmax^2-b^2*sigmay^2*tauInter1^2+2*ybar*(m+b*(sigmax^2+tauInter1^2)))/(2*(sigmay^2+sigmax^2+tauInter1^2))
    star2Exp= -(xbar^2+ybar^2-2*ybar*n+n^2+2*ybar*sigmax^2*u-2*n*sigmax^2*u-tauInter2^2*sigmax^2*u^2-sigmay^2*sigmax^2*u^2+2*xbar*(-ybar+n+(tauInter2^2+sigmay^2)*u))/(2*(tauInter2^2+sigmay^2+sigmax^2)) 
    star1Coeff= sqrt(tauInter1^2/(sigmay^2+sigmax^2+tauInter1^2))
    star2Coeff= sqrt(tauInter2^2/(sigmay^2+sigmax^2+tauInter2^2))
    
    alphavf=c+A*(star1Coeff*exp(star1Exp)-star2Coeff*exp(star2Exp))
    alphafv=c+A*(star2Coeff*exp(star2Exp)-star1Coeff*exp(star1Exp))
    alphaff=1+c1-delta-(1-2*delta)*(tauIntra2/sqrt(tauIntra2^2+sigmay^2))*exp(-(ybar-thetaH)^2/(2*tauIntra2^2+2*sigmay^2))
    
    # rates of change
    dNf=rf*Nf*(1-alphafv*Nv-alphaff*Nf)
    dNv=rv*Nv*(1-alphavf*Nf-alphavv*Nv)
    
    # return the rates of change as a list
    list(c(dNf,dNv,0,0))
    
  }) # end the with environment statement
}

########################################################################
########################################################################
########################################################################
# Same function but specialized for no population dynamics
########################################################################
########################################################################
########################################################################


GMspecialNoEco<-function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    
    c=A
    
    star1Exp= -(xbar^2+ybar^2+m^2-2*b*m*sigmay^2-2*xbar*(ybar+m-b*sigmay^2)-b^2*sigmay^2*sigmax^2-b^2*sigmay^2*tauInter1^2+2*ybar*(m+b*(sigmax^2+tauInter1^2)))/(2*(sigmay^2+sigmax^2+tauInter1^2))
    star2Exp= -(xbar^2+ybar^2-2*ybar*n+n^2+2*ybar*sigmax^2*u-2*n*sigmax^2*u-tauInter2^2*sigmax^2*u^2-sigmay^2*sigmax^2*u^2+2*xbar*(-ybar+n+(tauInter2^2+sigmay^2)*u))/(2*(tauInter2^2+sigmay^2+sigmax^2)) 
    star1Coeff= sqrt(tauInter1^2/(sigmay^2+sigmax^2+tauInter1^2))
    star2Coeff= sqrt(tauInter2^2/(sigmay^2+sigmax^2+tauInter2^2))
    
    # calculate partials
    derivCoef1xbar= -(2*xbar-2*(ybar+m-b*sigmay^2))/(2*(sigmay^2+sigmax^2+tauInter1^2))
    derivCoef1ybar= -(-2*xbar+2*ybar+2*(m+b*(sigmax^2+tauInter1^2)))/(2*(sigmay^2+sigmax^2+tauInter1^2))
    derivCoef2xbar= -(2*xbar+2*(-ybar+n+(tauInter2^2+sigmay^2)*u))/(2*(sigmay^2+sigmax^2+tauInter2^2))
    derivCoef2ybar= -(-2*xbar+2*ybar-2*n+2*sigmax^2*u)/(2*(sigmay^2+sigmax^2+tauInter2^2))
    
    partialx=rf*(-Nv*(derivCoef1xbar*star1Coeff*exp(star1Exp)-derivCoef2xbar*star2Coeff*exp(star2Exp))-Nf*(-1)*(1-2*delta)*(tauIntra1/(sqrt(tauIntra1^2+sigmax^2)))*(-1/(2*tauIntra1^2+2*sigmax^2))*2*(xbar-thetaC)*exp(-(xbar-thetaC)^2/(2*sigmax^2+2*tauIntra1^2)))
    partialy=rv*(-Nf*(derivCoef2ybar*star2Coeff*exp(star2Exp)-derivCoef1ybar*star1Coeff*exp(star1Exp))-Nv*(-1)*(1-2*delta)*(tauIntra2/(sqrt(tauIntra2^2+sigmay^2)))*(-1/(2*tauIntra2^2+2*sigmay^2))*2*(ybar-thetaH)*exp(-(ybar-thetaH)^2/(2*sigmay^2+2*tauIntra2^2)))
    
    
    # rates of change
    dxbar= H*sigmax^2*partialx
    dybar= H*sigmay^2*partialy
    
    
    # return the rates of change as a list
    list(c(0,0,dxbar,dybar))
    
  }) # end the with environment statement
}

########################################################################
########################################################################
########################################################################
# runGM using the lsoda solver
########################################################################
########################################################################
########################################################################

# try to keep non-negative
# note: it will set the state variable to zero for the NEXT time step, leaving the negative value
# recorded for the CURRENT time step, we can just erase these negative values and replace them with
# zeros at the end of the sim, because the negative values are NEVER used to calculate a next time step
eventfun <- function(t, y, parms){
  with(as.list(y), {
    # eliminate functionally extinct populations
    y1<-y[1:2]
    y1[y1<1e-6]<-0
    # eliminate negative trait values
    y2<-y[3:4]
    y2[y2<0]<-0
    return(c(y1,y2))
  })
}

# plotting functions for alpha values
ridgeDif1<-function(x,y,parameters) {
  with(as.list(parameters), {
    c=A
    return(z=c+A*(exp(-(x-(m+y))^2/(2*tauInter1^2))*exp(-b*y)-(exp(-(y-(n+x))^2/(2*tauInter2^2))*exp(-b*x))))
  })
}

ridgeDif2<-function(x,y,parameters) {
  with(as.list(parameters), {
    c=A
    return(z=c+A*(exp(-(y-(n+x))^2/(2*tauInter1^2))*exp(-b*x)-(exp(-(x-(m+y))^2/(2*tauInter2^2))*exp(-b*y))))
  })
} # this is the one for native/resistant species


plotDif1<-function(parameters,toxvec=1,resvec=1,col="blue",xmin=0,xmax=5,ymin=0,ymax=5,zmin=0,zmax=10,by=0.25,main="",xlab="",ylab="",zlab="",theta=0,phi=15,type="persp") {
    x<-seq(xmin,xmax,by)
    y<-seq(ymin,ymax,by)
    z<-outer(x,y,ridgeDif1,parameters=parameters)
    #persp(x,y,z,col=col,main=main,xlab=xlab,ylab=ylab,zlab=zlab)
    if(type=="persp"){
      persp(x,y,z,col=col,main=main,xlab=xlab,ylab=ylab,zlab=zlab,xlim=c(xmin,xmax),ylim=c(xmin,xmax),zlim=c(zmin,zmax),ticktype = "simple",theta=theta,phi=phi)
    }
    if(type=="heat"){
      image(seq(xmin,xmax,by=by),seq(ymin,ymax,by=by),zlim=c(zmin,zmax),z,xlab=xlab,ylab=ylab,main=main,col = hcl.colors(11, "Inferno", rev = T))
      legend((ymax-0.5),(ymax/1.5),c(zmin,"","","","","","","","","",zmax),fill=hcl.colors(11, "Inferno", rev = T),xpd = T,cex=0.6)
      points(toxvec[1],resvec[1],pch=1)
      points(toxvec,resvec,type="l")
    }
    if(type=="heatonly"){
      par(mar=c(5.1, 4.1, 2.1, 4.1),xpd=T)
      image(seq(xmin,xmax,by=by),seq(ymin,ymax,by=by),zlim=c(zmin,zmax),z,xlab=xlab,ylab=ylab,main=main,col = hcl.colors(11, "Inferno", rev = T))
      legend((ymax+0.5),(ymax/1.5),c(zmin,"","","","","","","","","",zmax),fill=hcl.colors(11, "Inferno", rev = T),xpd = T,cex=0.6)
    }
}


runGMSolver<-function(state,parameters,tmax=10000, deltat=0.1, repsiz=500, plotout=F, report=T, save=F, outputOnly=F, func=GMspecial,extinctlim=1e-06,method="lsoda",plotoutadd=F) {
  
  
  times <- seq(0, tmax, by = deltat)
  
  out<-ode(y=state,times=times,func=func,parms=parameters,events=list(func=eventfun,time=times),method=method)
  
  # replacing negative values with zero because ode still returns them even when replacing with zeros (see above note)
  # first populations
  out[which(out[,2]<1e-6),2]<-0
  out[which(out[,3]<1e-6),3]<-0
  # then traits
  out[which(out[,4]<0),4]<-0
  out[which(out[,5]<0),5]<-0
  
  if(save) {
    pdf(paste0("timeseries",i,namer,".pdf"))
    par(mfrow=c(2,2),mar=c(4,4,4,2))
    plot(out[,1],out[,4],type="l",main="Output",ylab="Average trait species 1",xlab="Time")
    plot(out[,1],out[,5],type="l", ylab="Average trait species 2",xlab="Time")
    plot(out[,1],out[,2],type="l",col="red",ylab="Population species 1",xlab="Time")
    plot(out[,1],out[,3],type="l",col="blue",ylab="Population species 2",xlab="Time")
    dev.off()
    pdf(paste0("alphaplot",i,namer,".pdf"))
    par(mfrow=c(1,1))
    dev.off()
  }
  
  if(plotout) {
    if(plotoutadd){
      
      out2<-ode(y=state,times=times,func=GMspecialNoEvo,parms=parameters,events=list(func=eventfun,time=times),method=method)
      # replacing negative values with zero because ode still returns them even when replacing with zeros (see above note)
      # first populations
      out2[which(out2[,2]<1e-6),2]<-0
      out2[which(out2[,3]<1e-6),3]<-0
      # then traits
      out2[which(out2[,4]<0),4]<-0
      out2[which(out2[,5]<0),5]<-0
      
      idx=seq(1,length(out[,1]),by=10)
      plot(out[,1],out[,4],type="l",main="",ylab="Average trait",xlab="Time",col="deeppink3",ylim=c(0,6.5),lwd=1.5)
      lines(out[,1],out[,5],col="darkgreen",lwd=1.5,lty=1)
      lines(out2[idx,1],out2[idx,4],col="deeppink3",lty=3,lwd=1.5)
      lines(out2[idx,1],out2[idx,5],col="darkgreen",lty=3,lwd=1.5)
      
      plot(out[,1],out[,2],type="l",col="deeppink3",ylab="Population",xlab="Time",lwd=1.5,ylim=c(0,1))
      lines(out[,1],out[,3],col="darkgreen",lwd=1.5,lty=1)
      lines(out2[idx,1],out2[idx,2],col="deeppink3",lty=3,lwd=1.5)
      lines(out2[idx,1],out2[idx,3],col="darkgreen",lty=3,lwd=1.5)
  
      
    } else {
      par(mfrow=c(2,2),mar=c(4,4,4,2))
      plot(out[,1],out[,4],type="l",main="Output",ylab="Average trait",xlab="Time",col="deeppink3")
      lines(out[,1],out[,5],type="l",col="darkgreen")
      plot(out[,1],out[,2],type="l",col="deeppink3",ylab="Population",xlab="Time")
      lines(out[,1],out[,3],type="l",col="darkgreen")
      par(mfrow=c(1,1))
     
    }
    

  }
  

  
  if(outputOnly){
    return(out)
  }
  
  if(report){
    
    trunc<-out[which(out[,1]>tmax-(repsiz)),]
    # trait1
    trait1mean2<-mean(trunc[,4]) # NEED these parentheses in the brackets! Nov 9
    trait1var2<-var(trunc[,4])
    trait1min2<-min(trunc[,4])
    trait1max2<-max(trunc[,4])
    trait1fin<-trunc[repsiz/deltat,4]
    # trait2
    trait2mean2<-mean(trunc[,5])
    trait2var2<-var(trunc[,5])
    trait2min2<-min(trunc[,5])
    trait2max2<-max(trunc[,5])
    trait2fin<-trunc[repsiz/deltat,5]
    # population 1
    pop1mean2<-mean(trunc[,2])
    pop1var2<-var(trunc[,2])
    pop1min2<-min(trunc[,2])
    pop1max2<-max(trunc[,2])
    pop1fin<-trunc[repsiz/deltat,2]
    # population 2
    pop2mean2<-mean(trunc[,3])
    pop2var2<-var(trunc[,3])
    pop2min2<-min(trunc[,3])
    pop2max2<-max(trunc[,3])
    pop2fin<-trunc[repsiz/deltat,3]
    
    result1<-as.data.frame(c(trait1mean2,trait1var2,trait1min2,trait1max2,trait1fin),rownames=c("mean","var","min","max","fin"))
    colnames(result1)<-"trait1"
    result1$trait2<-c(trait2mean2,trait2var2,trait2min2,trait2max2,trait2fin)
    result1$pop1<-c(pop1mean2,pop1var2,pop1min2,pop1max2,pop1fin)
    result1$pop2<-c(pop2mean2,pop2var2,pop2min2,pop2max2,pop2fin)
    
    return(result1)
    
  }
  
  
}



########################################################################
# sweepSolver updating to use solveAdaptive
########################################################################


sweepSolverModGM3<-function(parameters,param1="m",param2="n",param3=NA,start1=0, start2=0, end1=4, end2=4, by1=0.1, by2=0.1,param3set=NA, simlengthmin=10000, simlengthmax=50000,func=GMspecial,symmetric=F, invT0=5, invP0=0.05) {
  
  focal1=(start1+end1)/2
  focal2=(start2+end2)/2
  
  while(start1<0){
    start1=start1+by1
    end1=end1+by1
  }
  while(start2<0){
    start2=start2+by2
    end2=end2+by2
  }
  traitvec1<-seq(start1,end1,by=by1)
  traitvec2<-seq(start2,end2,by=by2) 
  
  # setup data container
  biglist<-list()
  
  for(j in 1:length(traitvec1)) {
    
    # setup data container
    smalllist<-list()
    
    for(k in 1:length(traitvec2)) { # iterating over sp2 (native)
      
      # set initial values
      # state<-c(Nf=0.1, Nv=0.9,xbar=0.5, ybar=0.1)
      
      parameters[which(names(parameters)==param1)]<-traitvec1[j]
      parameters[which(names(parameters)==param2)]<-traitvec2[k]
      if(!is.na(param3)){
        parameters[which(names(parameters)==param3)]<-param3set
      }
      if(symmetric) {
        parameters["n"]=parameters["m"]
        parameters["u"]=parameters["b"]
        parameters["tauInter2"]=parameters["tauInter1"]
        parameters["sigmay"]=parameters["sigmax"]
        parameters["rf"]=parameters["rv"]
      }
      

      #############################################
      # update to use adaptive solver
      
      data1=solveAdaptive(parameters,simlengthmin = 500,func=func, invT0=invT0, invP0=invP0)
      
      # record data
      
      smalllist[[k]]<-data1
      
      print(paste("finished round",j, " ",k))
      
    }
    
    biglist[[j]]<-smalllist #small list has the native row/column
    # biglist looks like [[invasive parameter]][[native parameter]]
    
  }
  
  biglist[[length(traitvec1)+1]]<-c(start1,end1,by1,start2,end2,by2)
  biglist[[length(traitvec1)+2]]<-c(param1,param2)
  biglist[[length(traitvec1)+3]]<-c(focal1,focal2)
  return(biglist)
  
}




########################################################################
# standardized solver- apply to PCA runs, grid search, parameter search
########################################################################

solveAdaptive=function(parameters, simlengthmin=1000 ,simlengthmax=50000, simlengthstep=2500,func=GMspecial, plotout=F, report=T, save=F,repsiz=500, extinctlim=1e-06, invT0=5, invP0=0.05){
  

  #############################################
  # checks 2: solve for native species alone in order to set initial conditions
  # Note: we need to parameters first!
  state<-c(Nf=0, Nv=1, xbar=0,ybar=1) 
  
  stateSolve<-tryCatch({
    namer="withEvo"
    runGMSolver(state,parameters,tmax=simlengthmin, deltat=0.1, repsiz=100, plotout=F, report=T, save=F,func=GMspecial,extinctlim=extinctlim)
  },error=function(error_message){
    message(error_message)
    return(as.data.frame(matrix(rep(-1,20),nrow=5)))
  })
  
  # adaptive solving for native species
  # add the adaptive time number modification
  stab=F
  stab2=F
  while((stab&stab2)==F) {
    
    # test the data to see if it is stabilized
    
    # check populations
    if(stateSolve[1,3]<0.0001&stateSolve[1,4]<0.0001){
      stab=T # both extinct
    } else if(stateSolve[1,3]<0.0001){ # else if pop 1 mean
      stab=T # species 1 extinct
    } else if(stateSolve[1,4]<0.0001){ # else if pop 2 mean
      stab=T # species 2 extinct
    } else if(stateSolve[2,3]<0.001&stateSolve[2,4]<0.0001& # else if pop 1 var and pop 2 var #small variance
              # pop 1 mean minus pop 1 final
              abs(stateSolve[1,3]-stateSolve[5,3])<0.01& #mean and final value are close
              # pop 2 mean minus pop 2 final
              abs(stateSolve[1,4]-stateSolve[5,4])<0.01){
      stab=T #stabilized
    } else if (nSampCur>=simlengthmax) {
      # end simulations- we are up against the maximum sim length!
      stab=T
    }
    
    # check traits
    if(stateSolve[1,1]<0.0001&stateSolve[1,2]<0.0001){
      stab2=T # both trait extinct
    } else if(stateSolve[1,1]<0.0001){ # else if pop 1 mean
      stab2=T # trait 1 extinct
    } else if(stateSolve[1,2]<0.0001){# else if pop 2 mean
      stab2=T # trait 2 extinct
    } else if(stateSolve[2,1]<0.001&stateSolve[2,2]<0.0001& # else if pop 1 var and pop 2 var #small variance
              # pop 1 mean minus pop 1 final
              abs(stateSolve[1,1]-stateSolve[5,1])<0.01& #mean and final value are close
              # pop 2 mean minus pop 2 final
              abs(stateSolve[1,2]-stateSolve[5,2])<0.01){
      stab2=T #stabilized
    } else if (nSampCur>=simlengthmax) {
      # end simulations- we are up against the maximum sim length!
      stab2=T
    }
    
    # see if we need to run sim based on population AND traits
    if((stab&stab2)==F) {
      # run sim for longer- another simlengthstep time steps
      # use longer deltat if farther along (less likely to have very bumpy behavior)
      state<-c(Nf=stateSolve$pop1[5], Nv=stateSolve$pop2[5],xbar=stateSolve$trait1[5], ybar=stateSolve$trait2[5])
      if(nSampCur<7000){
        dtcur=0.1
      } else if(nSampCur<15000){
        dtcur=0.5
      } else {
        dtcur=1
      }
      stateSolve<-runGMSolver(state=state,parameters=parameters,tmax=simlengthstep,deltat=dtcur,report=T,plotout=F,func=func)
      nSampCur<-nSampCur+simlengthstep
    }
    
  }
  
  # use the solution from the native species only for the initial condition
  state<-c(Nf=(stateSolve[5,4]*invP0), Nv=stateSolve[5,4], xbar=(stateSolve[5,2]*invT0),ybar=stateSolve[5,2]) 
  state_ini=state
  
  # run model 
  
  data1<-tryCatch({
    namer="newTryCatch"
    runGMSolver(state=state, parameters=parameters, tmax = simlengthmin, deltat=0.1, report=T, plotout=plotout, func=func)
  },error=function(error_message){
    message(error_message)
    return(as.data.frame(matrix(rep(-1,20),nrow=5)))
  })
  nSampCur<-simlengthmin
  
  # add the adaptive time number modification
  stab=F
  stab2=F
  while((stab&stab2)==F) {
    
    # new code!
    # test the data to see if it is stabilized
    
    # check populations
    if(data1[1,3]<0.0001&data1[1,4]<0.0001){
      stab=T # both extinct
    } else if(data1[1,3]<0.0001){ # else if pop 1 mean
      stab=T # species 1 extinct
    } else if(data1[1,4]<0.0001){ # else if pop 2 mean
      stab=T # species 2 extinct
    } else if(data1[2,3]<0.001&data1[2,4]<0.0001& # else if pop 1 var and pop 2 var #small variance
              # pop 1 mean minus pop 1 final
              abs(data1[1,3]-data1[5,3])<0.01& #mean and final value are close
              # pop 2 mean minus pop 2 final
              abs(data1[1,4]-data1[5,4])<0.01){
      stab=T #stabilized
    } else if (nSampCur>=simlengthmax) {
      # end simulations- we are up against the maximum sim length!
      stab=T
    }
    
    # check traits
    if(data1[1,1]<0.0001&data1[1,2]<0.0001){
      stab2=T # both extinct
    } else if(data1[1,1]<0.0001){ # else if pop 1 mean
      stab2=T # species 1 extinct
    } else if(data1[1,2]<0.0001){# else if pop 2 mean
      stab2=T # species 2 extinct
    } else if(data1[2,1]<0.001&data1[2,2]<0.0001& # else if pop 1 var and pop 2 var #small variance
              # pop 1 mean minus pop 1 final
              abs(data1[1,1]-data1[5,1])<0.01& #mean and final value are close
              # pop 2 mean minus pop 2 final
              abs(data1[1,2]-data1[5,2])<0.01){
      stab2=T #stabilized
    } else if (nSampCur>=simlengthmax) {
      # end simulations- we are up against the maximum sim length!
      stab2=T
    }
    
    # see if we need to run sim based on population AND traits
    if((stab&stab2)==F) {
      # run sim for longer- another 10,000 time steps
      state<-c(Nf=data1$pop1[5], Nv=data1$pop2[5],xbar=data1$trait1[5], ybar=data1$trait2[5])
      if(nSampCur<7000){
        dtcur=0.1
      } else if(nSampCur<15000){
        dtcur=0.5
      } else {
        dtcur=1
      }
      data1<-runGMSolver(state=state,parameters=parameters,tmax=simlengthstep,deltat=dtcur,report=T,plotout=plotout,func=func)
      nSampCur<-nSampCur+10000
    }
    
  }
  
  return(list(state_ini,data1))
  
}


#########################
# Modified Jun 22 2021
# trait categories are conditional on population values
# takes in data from sweepGM 
# outputs a plot, option to save or display
# option to output data too
plotGMB<-function(biglist,save=T,givedat=F,plotname="defaultplot",plotlines=F,outline1=NA,outline2=NA,nokey=F,lwd=4){
  
  dim1<-length(biglist)-2
  dim2<-length(biglist[[1]])
  
  conditionMat1<-matrix(NA,nrow=dim1,ncol=dim2)
  
  for(j in 1:dim1) {
    
    for(k in 1:dim2) {
      
      # if trait 1 mean and trait 2 mean are about 0, and both species are NOT EXTINCT
      if(biglist[[j]][[k]][1,1]<0.0001&biglist[[j]][[k]][1,2]<0.0001&biglist[[j]][[k]][1,3]>0.0001&biglist[[j]][[k]][1,4]>0.0001){
        conditionMat1[j,k]=2 # both zero trait value
      }
      # else if trait 1 mean is zero and species 2 is extinct
      else if(biglist[[j]][[k]][1,1]<0.0001&biglist[[j]][[k]][1,4]<0.0001){
        conditionMat1[j,k]=8 # species 1 zero
      }
      # else if trait 2 mean is zero and species 1 is extinct
      else if(biglist[[j]][[k]][1,2]<0.0001&biglist[[j]][[k]][1,3]<0.0001){
        conditionMat1[j,k]=1 # species 2 zero
      }
      # else if trait 1 mean is zero and trait 2 is non zero
      else if(biglist[[j]][[k]][1,1]<0.0001&biglist[[j]][[k]][1,2]>0.0001){
        conditionMat1[j,k]=3 # species 2 zero
      }
      # else if trait 2 mean is zero and trait 1 is non zero
      else if(biglist[[j]][[k]][1,2]<0.0001&biglist[[j]][[k]][1,1]>0.0001){
        conditionMat1[j,k]=4 # species 2 zero
      }
      # else if trait 1 mean is non zero and species 2 is extinct
      else if(biglist[[j]][[k]][1,1]>0.0001&biglist[[j]][[k]][1,4]<0.0001){
        conditionMat1[j,k]=7 # species 2 zero
      }
      # else if trait 2 mean is non zero and species 1 is extinct
      else if(biglist[[j]][[k]][1,2]>0.0001&biglist[[j]][[k]][1,3]<0.0001){
        conditionMat1[j,k]=6 # species 2 zero
      }
      
      # else if trait 1 var and trait 2 var
      else if(biglist[[j]][[k]][2,1]<0.001&biglist[[j]][[k]][2,2]<0.0001& #small variance
              # trait 1 mean minus trait 1 final
              abs(biglist[[j]][[k]][1,1]-biglist[[j]][[k]][5,1])<0.01& #mean and final value are close
              # trait 2 mean minus trait 2 final
              abs(biglist[[j]][[k]][1,2]-biglist[[j]][[k]][5,2])<0.01){
        conditionMat1[j,k]=5 #stabilized
      }
      # else if trait 1 var and trait 2 var
      else if(biglist[[j]][[k]][2,1]>=0.001&biglist[[j]][[k]][2,2]>=0.001& #larger variance
              # trait 1 mean minus trait 1 final
              abs(biglist[[j]][[k]][1,1]-biglist[[j]][[k]][5,1])>=0.01& # mean and final value are farther
              # trait 2 mean minus trait 2 final
              abs(biglist[[j]][[k]][1,2]-biglist[[j]][[k]][5,2])>=0.01){
        conditionMat1[j,k]=9 #oscillating or not stabilized
      }
      else {
        conditionMat1[j,k]=10 # not quite stabilized
      }
      
    }
    
  }
  
  
  
  # plot the data (depending on if data is to be saved or displayed)
  
  start1<-biglist[[length(biglist)-1]][1]
  end1<-biglist[[length(biglist)-1]][2]
  by1<-biglist[[length(biglist)-1]][3]
  start2<-biglist[[length(biglist)-1]][4]
  end2<-biglist[[length(biglist)-1]][5]
  by2<-biglist[[length(biglist)-1]][6]
  param1<-biglist[[length(biglist)]][1]
  param2<-biglist[[length(biglist)]][2]
  if(plotlines){
    x1<-findOutlines(outline1)
    x2<-findOutlines(outline2)
  }
  
  
  if(save) {
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,8))
    
    pdf(paste0(plotname,"1.pdf"))
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(1,10),conditionMat1,xlab=paste0(param1,"invasive"),ylab=paste0(param2," native"),main="Traits",col = c(hcl.colors(8, "Green-Brown", rev = F),"red","orangered"))
    if(!nokey){
      legend(5.3,3,c("Inv ext Nat 0","Both 0 not ext","inv 0 nat pos","inv pos nat 0","stabilized","inv ext nat pos","inv pos nat ext","inv 0 nat ext","oscillating","not quite stabilized"),fill=c(hcl.colors(8, "Green-Brown", rev = F),"red","orangered"),xpd = T,cex=0.6)
    }
      
    # plot outlines
    if(plotlines){
      
      xdim<-dim(x1[[1]])[1]
      ydim<-dim(x1[[1]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x1[[1]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1+by1/2),((j-1)*by2-by2/2),lwd=lwd)
          }
        }
      }
      
      xdim<-dim(x1[[2]])[1]
      ydim<-dim(x1[[2]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x1[[2]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1-by1/2),((j-1)*by2+by2/2),lwd=lwd)
          }
        }
      }
      ########
      
      xdim<-dim(x2[[1]])[1]
      ydim<-dim(x2[[1]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x2[[1]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1+by1/2),((j-1)*by2-by2/2),lwd=lwd,col="grey70")
          }
        }
      }

      xdim<-dim(x2[[2]])[1]
      ydim<-dim(x2[[2]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x2[[2]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1-by1/2),((j-1)*by2+by2/2),lwd=lwd,col="grey70")
          }
        }
      }
    }
    dev.off()
    
  } else {
    
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,8))
    
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(1,10),conditionMat1,xlab=paste0(param1,"invasive"),ylab=paste0(param2,"native"),main="Traits",col = c(hcl.colors(8, "Green-Brown", rev = F),"red","orangered"))
    if(!nokey){
      legend(5.3,3,c("Inv ext Nat 0","Both 0 not ext","inv 0 nat pos","inv pos nat 0","stabilized","inv ext nat pos","inv pos nat ext","inv 0 nat ext","oscillating","not quite stabilized"),fill=c(hcl.colors(8, "Green-Brown", rev = F),"red","orangered"),xpd = T,cex=0.6)
    }    # plot outlines
    if(plotlines){
      
      xdim<-dim(x1[[1]])[1]
      ydim<-dim(x1[[1]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x1[[1]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1+by1/2),((j-1)*by2-by2/2),lwd=lwd)
          }
        }
      }
      
      xdim<-dim(x1[[2]])[1]
      ydim<-dim(x1[[2]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x1[[2]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1-by1/2),((j-1)*by2+by2/2),lwd=lwd)
          }
        }
      }
      ########
      
      xdim<-dim(x2[[1]])[1]
      ydim<-dim(x2[[1]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x2[[1]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1+by1/2),((j-1)*by2-by2/2),lwd=lwd,col="grey70")
          }
        }
      }
      
      xdim<-dim(x2[[2]])[1]
      ydim<-dim(x2[[2]])[2]
      for(i in 1:xdim) {
        for(j in 1:ydim) {
          if(x2[[2]][i,j]==1) {
            segments(((i-1)*by1-by1/2),((j-1)*by2-by2/2),((i-1)*by1-by1/2),((j-1)*by2+by2/2),lwd=lwd,col="grey70")
          }
        }
      }
    }

  }
  
  # output data
  if(givedat) {
    return(list(conditionMat1))
  }
  
}

# takes in an object with 2 matrices of the numbers to plot
# returns a grids of which line segments should be drawn to outline sections of the same color
findOutlines<-function(x) {
  
  # start with traits
  xdim1<-dim(x[[1]])[1]
  ydim1<-dim(x[[1]])[2]
  
  # first do the dividers across
  outlineMatHoriz1<-matrix(1,nrow=xdim1,ncol=ydim1+1)
  
  for(i in 1:(xdim1)) {
    for(j in 1:(ydim1-1)) {
      if(x[[1]][i,j]==x[[1]][i,j+1]) {
        outlineMatHoriz1[i,j+1]<-0
      }
    }
  }
  
  # then to the dividers down
  outlineMatVert1<-matrix(1,nrow=xdim1+1,ncol=ydim1)
  
  for(i in 1:(xdim1-1)) {
    for(j in 1:(ydim1)) {
      if(x[[1]][i,j]==x[[1]][i+1,j]) {
        outlineMatVert1[i+1,j]<-0
      }
    }
  }
  
  
  return(list(outlineMatHoriz1,outlineMatVert1))
  
}

########################################
########################################
# Code to plot the equilibrium values 
########################################
########################################

plotGMEQ<-function(biglist,save=T,givedat=F,plotname="defaultplot"){
  # just plot final values
  
  # get dimensions
  dim1<-length(biglist)-2
  dim2<-length(biglist[[1]])
  
  finValMatT1<-matrix(NA,nrow=dim1,ncol=dim2)
  finValMatT2<-matrix(NA,nrow=dim1,ncol=dim2)
  finValMatP1<-matrix(NA,nrow=dim1,ncol=dim2)
  finValMatP2<-matrix(NA,nrow=dim1,ncol=dim2)
  meanValMatT1<-matrix(NA,nrow=dim1,ncol=dim2)
  meanValMatT2<-matrix(NA,nrow=dim1,ncol=dim2)
  meanValMatP1<-matrix(NA,nrow=dim1,ncol=dim2)
  meanValMatP2<-matrix(NA,nrow=dim1,ncol=dim2)
  
  for(j in 1:dim1) {
    
    for(k in 1:dim2) {
      
      finValMatT1[j,k]<-biglist[[j]][[k]][5,1]
      finValMatT2[j,k]<-biglist[[j]][[k]][5,2]
      finValMatP1[j,k]<-biglist[[j]][[k]][5,3]
      finValMatP2[j,k]<-biglist[[j]][[k]][5,4]
      meanValMatT1<-biglist[[j]][[k]][1,1]
      meanValMatT2<-biglist[[j]][[k]][1,2]
      meanValMatP1<-biglist[[j]][[k]][1,3]
      meanValMatP2<-biglist[[j]][[k]][1,4]
      
    }
    
  }
  
  # now for plotting the data
  
  start1<-biglist[[length(biglist)-1]][1]
  end1<-biglist[[length(biglist)-1]][2]
  by1<-biglist[[length(biglist)-1]][3]
  start2<-biglist[[length(biglist)-1]][4]
  end2<-biglist[[length(biglist)-1]][5]
  by2<-biglist[[length(biglist)-1]][6]
  param1<-biglist[[length(biglist)]][1]
  param2<-biglist[[length(biglist)]][2]
  
  if(save) {
    
    pdf(paste0(plotname,"1.pdf"))
    
    par(mfrow=c(2,2),mar=c(5.1,4.1,4.1,8))
    
    loc1<-end2+(end2-start2)/10
    loc2<-end1-(end1-start1)/10
    
    # Instead of just using min and max, we want to establish a range that doesn't change as unpredictably
    # I will make the botton of the range 0 for all cases, and 
    minT1<-0
    maxT1<-max(max(finValMatT1),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minT1,maxT1),finValMatT1,xlab=paste0(param1," invasive"),ylab=paste0(param1," native"),main="Equilibrium traits invasive")
    legend(loc1,loc2,c(round(minT1,2),"","","","","","","","","","",round(maxT1,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minT2<-0
    maxT2<-max(max(finValMatT2),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minT2,maxT2),finValMatT2,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium traits native")
    legend(loc1,loc2,c(round(minT2,2),"","","","","","","","","","",round(maxT2,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minP1<-0
    maxP1<-max(max(finValMatP1),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minP1,maxP1),finValMatP1,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium population invasive")
    legend(loc1,loc2,c(round(minP1,2),"","","","","","","","","","",round(maxP1,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minP2<-0
    maxP2<-max(max(finValMatP2),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minP2,maxP2),finValMatP2,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium population native")
    legend(loc1,loc2,c(round(minP2,2),"","","","","","","","","","",round(maxP2,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    
    dev.off()
    
  } else {
    
    par(mfrow=c(2,2),mar=c(5.1,4.1,4.1,8))
    
    loc1<-end2+(end2-start2)/10
    loc2<-end1-(end1-start1)/10
    
    minT1<-0
    maxT1<-max(max(finValMatT1),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minT1,maxT1),finValMatT1,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium traits invasive")
    legend(loc1,loc2,c(round(minT1,2),"","","","","","","","","","",round(maxT1,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minT2<-0
    maxT2<-max(max(finValMatT2),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minT2,maxT2),finValMatT2,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium traits native")
    legend(loc1,loc2,c(round(minT2,2),"","","","","","","","","","",round(maxT2,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minP1<-0
    maxP1<-max(max(finValMatP1),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minP1,maxP1),finValMatP1,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium population native")
    legend(loc1,loc2,c(round(minP1,2),"","","","","","","","","","",round(maxP1,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
    minP2<-0
    maxP2<-max(max(finValMatP2),1)
    image(seq(start1,end1,by=by1),seq(start2,end2,by=by2),zlim=c(minP2,maxP2),finValMatP2,xlab=paste0(param1,"invasive"),ylab=paste0(param1,"native"),main="Equilibrium population invasive")
    legend(loc1,loc2,c(round(minP2,2),"","","","","","","","","","",round(maxP2,2)),fill=hcl.colors(12, "YlOrRd", rev = TRUE),xpd = T,cex=0.5)
    
  }
  
  # output data
  if(givedat) {
    return(list(finValMatT1,finValMatT2,finValMatP1,finValMatP2,meanValMatT1,meanValMatT2,meanValMatP1,meanValMatP2))
  }
  
}

#### Extract parameters from one of the lists returned ####
extractPars<-function(x){
  parout<-c(m=x$m, n=x$n, b=x$b, u=x$u, A=x$A, H=x$H, c1=x$c1, 
            delta=x$delta, thetaC=x$thetaC, thetaH=x$thetaH, tauIntra1=x$tauIntra1, tauIntra2=x$tauIntra2,
            tauInter1=x$tauInter1, tauInter2=x$tauInter2, sigmax=x$sigmax, sigmay=x$sigmay,
            rf=x$rf, rv=x$rv)
  return(parout)
}


########################################################################


sweepSolverModGM3<-function(parameters,param1="m",param2="n",param3=NA,start1=0, start2=0, end1=4, end2=4, by1=0.1, by2=0.1,param3set=NA, simlengthmin=10000, simlengthmax=50000,func=GMspecial,symmetric=F, invT0=5, invP0=0.05) {
  
  focal1=(start1+end1)/2
  focal2=(start2+end2)/2
  
  while(start1<0){
    start1=start1+by1
    end1=end1+by1
  }
  while(start2<0){
    start2=start2+by2
    end2=end2+by2
  }
  traitvec1<-seq(start1,end1,by=by1)
  traitvec2<-seq(start2,end2,by=by2) 
  
  # setup data container
  biglist<-list()
  
  for(j in 1:length(traitvec1)) {
    
    # setup data container
    smalllist<-list()
    
    for(k in 1:length(traitvec2)) { # iterating over sp2 (native)
      
      parameters[which(names(parameters)==param1)]<-traitvec1[j]
      parameters[which(names(parameters)==param2)]<-traitvec2[k]
      if(!is.na(param3)){
        parameters[which(names(parameters)==param3)]<-param3set
      }
      if(symmetric) {
        parameters["n"]=parameters["m"]
        parameters["u"]=parameters["b"]
        parameters["tauInter2"]=parameters["tauInter1"]
        parameters["sigmay"]=parameters["sigmax"]
        parameters["rf"]=parameters["rv"]
      }
      

      #############################################
      # updated to use adaptive solver function
      ridgeDif1()
      data1=ridge(parameters,simlengthmin = 500,func=func, invT0=invT0, invP0=invP0)
      data1=solveAdaptive(parameters,simlengthmin = 500,func=func, invT0=invT0, invP0=invP0)
      
      # record data
      
      smalllist[[k]]<-data1
      
      print(paste("finished round",j, " ",k))
      
    }
    
    biglist[[j]]<-smalllist #small list has the native row/column
    # biglist looks like [[invasive parameter]][[native parameter]]
    
  }
  
  biglist[[length(traitvec1)+1]]<-c(start1,end1,by1,start2,end2,by2)
  biglist[[length(traitvec1)+2]]<-c(param1,param2)
  biglist[[length(traitvec1)+3]]<-c(focal1,focal2)
  return(biglist)
  
}


NatCG_traits=function(gridResult,gridResultNoEvo,parameter){
  focal1=gridResult[[length(gridResult)]][1]
  focal2=gridResult[[length(gridResult)]][2]
  
  start1=gridResult[[length(gridResult)-2]][1]
  end1=gridResult[[length(gridResult)-2]][2]
  by1=gridResult[[length(gridResult)-2]][3]
  start2=gridResult[[length(gridResult)-2]][4]
  end2=gridResult[[length(gridResult)-2]][5]
  by2=gridResult[[length(gridResult)-2]][6]
  gridSize=length(seq(start1,end1,by=by1))
  
  # get which
  values1=seq(start1,end1,by=by1)
  values2=seq(start2,end2,by=by2)
  # this gives us which square is the focal square
  focal_which1=which(values1==focal1) 
  focal_which2=which(values2==focal2)
  if(length(focal_which1)==0|length(focal_which2)==0){
    focal_which1=which(abs(values1-focal1)<1e-5)
    focal_which2=which(abs(values2-focal2)<1e-5)
  }
  
  
  incLev=1.1
  decLev=0.9
  matrix1=matrix(NA,nrow=gridSize,ncol=gridSize)
  for(i in 1:gridSize){
    for(j in 1:gridSize){
      finalNat=gridResult[[i]][[j]][[2]][5,4] # ok so we have [[i]] for invasive and [[j]] for native
      finalNatNoEvo=gridResultNoEvo[[i]][[j]][[2]][5,4]
      if(is.na(finalNat)|is.nan(finalNat)|is.infinite(finalNat)|finalNat<0|is.na(finalNatNoEvo)|is.nan(finalNatNoEvo)|is.infinite(finalNatNoEvo)|finalNatNoEvo<0){
        cat(finalNat)
        cat("\n")
        cat(finalNatNoEvo)
        cat("\n")
        matrix1[i,j]=4 # E for error
      } else if(finalNatNoEvo==0&finalNat==0){
        matrix1[i,j]=3 # both are zero -> same
      } else if(finalNatNoEvo==0&finalNat!=0){
        matrix1[i,j]=1 # only no evo is zero -> increase
      } else if(finalNat/finalNatNoEvo>=incLev){
        matrix1[i,j]=1
      } else if(finalNat/finalNatNoEvo<=decLev){
        matrix1[i,j]=2
      } else {
        matrix1[i,j]=3
        # matrix1 looks like [invasive parameter, native parameter]
      }
    }
  }
  
  # each gridResult object comes from the sweepSolverModGM3 function
  # note that when using image(), the columns of the matrix will be plotted as rows of the image
  
  par(xpd=F)
  if(parameter=="sigma") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(sigma[1]),ylab=expression(sigma[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="tau"){
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(tau[1]),ylab=expression(tau[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="b") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(b[1]),ylab=expression(b[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="m") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(m[1]),ylab=expression(m[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=paste(parameter, "invasive"),ylab=paste(parameter, native),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  }
  
  par(xpd=T)
  text(x=1,y=0,round(start1,digits=1)) # x-axis label- should be invasive = 1
  text(x=gridSize,y=0,round(end1,digits=1)) # x-axis label- should be invasive 1
  text(x=0.3,y=1,round(start2,digits=1),adj=1) # y-axis label- should be native 2
  text(x=0.3,y=gridSize,round(end2,digits=1),adj=1) # y-axis label- should be native 2
  #text(x=-0.5,y=gridSize*focal2/(start2+end2)/2,round(focal2))
  #par(xpd=F)
  # add gridlines
  par(xpd=F)
  for (i in 1:gridSize) {
    # Horizontal lines
    abline(h = i - 0.5, col = "black", lwd = 0.5)
    # Vertical lines
    abline(v = i - 0.5, col = "black", lwd = 0.5)
  }
  
  # highlight the "focal" cell
  # check: x should be the native parameter (focal_which2)
  segments(x0=focal_which1-0.5,x1=focal_which1+0.5,y0=focal_which2-0.5,y1=focal_which2-0.5,lwd=3)
  segments(x0=focal_which1-0.5,x1=focal_which1+0.5,y0=focal_which2+0.5,y1=focal_which2+0.5,lwd=3)
  segments(x0=focal_which1-0.5,x1=focal_which1-0.5,y0=focal_which2-0.5,y1=focal_which2+0.5,lwd=3)
  segments(x0=focal_which1+0.5,x1=focal_which1+0.5,y0=focal_which2-0.5,y1=focal_which2+0.5,lwd=3)
  
  
}

NatCG_alphas=function(gridResult,gridResultNoEvo,parameter){
  focal1=gridResult[[length(gridResult)]][1]
  focal2=gridResult[[length(gridResult)]][2]
  
  start1=gridResult[[length(gridResult)-2]][1]
  end1=gridResult[[length(gridResult)-2]][2]
  by1=gridResult[[length(gridResult)-2]][3]
  start2=gridResult[[length(gridResult)-2]][4]
  end2=gridResult[[length(gridResult)-2]][5]
  by2=gridResult[[length(gridResult)-2]][6]
  gridSize=length(seq(start1,end1,by=by1))
  
  # get which
  values1=seq(start1,end1,by=by1)
  values2=seq(start2,end2,by=by2)
  # this gives us the which square is the focal square
  focal_which1=which(values1==focal1) 
  focal_which2=which(values2==focal2)
  if(length(focal_which1)==0|length(focal_which2)==0){
    focal_which1=which(abs(values1-focal1)<1e-5)
    focal_which2=which(abs(values2-focal2)<1e-5)
  }
  
  
  incLev=1.1
  decLev=0.9
  matrix1=matrix(NA,nrow=gridSize,ncol=gridSize)
  for(i in 1:gridSize){
    for(j in 1:gridSize){
      traitNat=gridResult[[i]][[j]][[2]][5,2] 
      traitInv=gridResult[[i]][[j]][[2]][5,1] 
      
      inv_alpha=ridgeDif1(traitInv,traitNat,parameters)
      nat_alpha=ridgeDif1(traitInv,traitNat,parameters)
      # how to recover parameters
      
      
      
    }
  }
  

  # each gridResult object comes from the sweepSolverModGM3 function
  # when using image(), the columns of the matrix will be plotted as rows of the image

  par(xpd=F)
  if(parameter=="sigma") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(sigma[1]),ylab=expression(sigma[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="tau"){
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(tau[1]),ylab=expression(tau[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="b") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(b[1]),ylab=expression(b[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="m") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(m[1]),ylab=expression(m[2]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=paste(parameter, "invasive"),ylab=paste(parameter, native),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  }
  
  par(xpd=T)
  text(x=1,y=0,round(start1,digits=1)) # x-axis label- should be invasive = 1
  text(x=gridSize,y=0,round(end1,digits=1)) # x-axis label- should be invasive 1
  text(x=0.3,y=1,round(start2,digits=1),adj=1) # y-axis label- should be native 2
  text(x=0.3,y=gridSize,round(end2,digits=1),adj=1) # y-axis label- should be native 2
  #text(x=-0.5,y=gridSize*focal2/(start2+end2)/2,round(focal2))
  #par(xpd=F)
  # add gridlines
  par(xpd=F)
  for (i in 1:gridSize) {
    # Horizontal lines
    abline(h = i - 0.5, col = "black", lwd = 0.5)
    # Vertical lines
    abline(v = i - 0.5, col = "black", lwd = 0.5)
  }
  
  # highlight the "focal" cell
  # check: x should be the native parameter (focal_which2)
  segments(x0=focal_which1-0.5,x1=focal_which1+0.5,y0=focal_which2-0.5,y1=focal_which2-0.5,lwd=3)
  segments(x0=focal_which1-0.5,x1=focal_which1+0.5,y0=focal_which2+0.5,y1=focal_which2+0.5,lwd=3)
  segments(x0=focal_which1-0.5,x1=focal_which1-0.5,y0=focal_which2-0.5,y1=focal_which2+0.5,lwd=3)
  segments(x0=focal_which1+0.5,x1=focal_which1+0.5,y0=focal_which2-0.5,y1=focal_which2+0.5,lwd=3)
  

  
  
}

ICSearch=function(pcpoint,traitmult=NA,popmult=NA, ntrials=100, nsplit=1, splitsec=1){
  
  if(is.na(traitmult)){
    traitmult=c(0.2,0.5,1,2,4)
  }
  if(is.na(popmult)){
    popmult=c(0.0001,0.001,0.01,0.1)
  }
  
  
  
  y_all=1:gridLength
  if(nsplit>1){
    categories=cut(y_all,breaks=nsplit,labels=F)
    y_select=y_all[which(categories==splitsec)]
  } else {
    y_select=y_all
  }
  
  # set up lists
  modellistFig2<-list()
  modellistFig2NoEvo<-list()
  # simulation loop
  for(j in y_select){# j is y 9-13
    for(k in 1:gridLength){ # k is x 16
      smallList=list()
      smallListNoEvo=list()
      for(i in 1:ntrials){
        
        parStart_i=getFullPCcoords(c(xseq[k],yseq[j]),commMatAxes[[2]]$rotation,scale=commMatAxes[[2]]$scale,center=commMatAxes[[2]]$center,maxTime=120)
        if(length(parStart_i)==1){ # if parameter list cant find positive parameters in two minutes, stop
          smallList[[i]]<-NA
          smallListNoEvo[[i]]<-NA
        } else {
          parametersFig1A<-c(m=parStart_i[1], n=1, b=parStart_i[2], u=1, A=parStart_i[3], H=parStart_i[4], c1=parStart_i[5],
                             delta=parStart_i[6], thetaC=0, thetaH=0, tauIntra1=0.412715, tauIntra2=0.412715,
                             tauInter1=parStart_i[7], tauInter2=0.412715, sigmax=parStart_i[8], sigmay=0.25,
                             rf=1, rv=1)
          
          
          
          # set optima for single species to be same as m, n
          parametersFig1A[which(names(parametersFig1A)=="thetaC")]<-parametersFig1A[which(names(parametersFig1A)=="m")]
          parametersFig1A[which(names(parametersFig1A)=="thetaH")]<-parametersFig1A[which(names(parametersFig1A)=="n")]
          
          parameters=parametersFig1A
          
          
          # use solve adaptive here
          res=solveAdaptive(parameters,simlengthmin = 500,func=GMspecial)
          resNoEvo=solveAdaptive(parameters,simlengthmin = 500,func=GMspecialNoEvo)
          
          
          # get output ready
          output<-c(parameters,res[[1]],res[[2]][5,])
          outputNoEvo<-c(parameters,resNoEvo[[1]],resNoEvo[[2]][5,])
          
          smallList[[i]]<-output
          smallListNoEvo[[i]]<-outputNoEvo
          
          if(length(smallList)==10&length(which(is.na(smallList)))==10){
            fail=T # if we have 10 in a row that come back NA, just skip the rest
          } else {
            # nothing
          }
        }
        
        
        cat(paste("j is",j,"\n"))
        cat(paste("k is",k,"\n"))
        cat(paste("trial",i,"done","\n"))
        
      }
      
      modellistFig2[[length(modellistFig2)+1]]=smallList
      modellistFig2NoEvo[[length(modellistFig2NoEvo)+1]]=smallListNoEvo
    }
  }
  
  return(list(modellistFig2,modellistFig2NoEvo))
  
}
