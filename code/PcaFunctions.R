# This is a list of functions for doing a structured parameter search by drawing from PCA
# Also includes functions for visualizing PCA plots

########################################################################
########################################################################
########################################################################
########################################################################
#Abundance
#- 1) Who is more abundant
# definition: invasive more abundant: inv/nat>=1.5
# definition: similar abundance: 0.67<inv/nat<1.5
# definition: native more abundant: inv/nat<=0.67



#' Title Run Scenario 1:Who is more abundant
#'
#' @param modList list of full model results with evolution
#' @param modListNoEvo list of full model results without evolution
#' @param natLev level of inv/nat below which native is considered more abundant
#' @param invLev level of inv/nat above which invasive is considered more abundant
#'
#' @return list of outcomes, with each element corresponding to a certain outcome, see code for details
#' @export
#'
#' @examples
runScen1=function(modList,modListNoEvo,natLev=0.67,invLev=1.5){
  
  outlier1<-0
  outlier2<-0
  listerA<-numeric(0) # Native more abundant -> Native more abundant
  listerB<-numeric(0) # Native more abundant -> Invasive more abundant
  listerC<-numeric(0) # Invasive more abundant -> Invasive more abundant
  listerD<-numeric(0) # Invasive more abundant -> Native more abundant
  listerE<-numeric(0) # Native more abundant -> similar abundance
  listerF<-numeric(0) # Invasive more abundant -> similar abundance
  listerG<-numeric(0) # similar abundance -> similar abundance
  listerH<-numeric(0) # similar abundance -> Native more abundant
  listerI<-numeric(0) # similar abundance -> Invasive more abundant
  samp1<-numeric(0)
  changeEvo<-matrix(0, nrow=2,ncol=3)
  
  for(i in 1:length(modList)){
    if(is.na(modListNoEvo[[i]][25][[1]])|is.nan(modListNoEvo[[i]][25][[1]])|is.infinite(modListNoEvo[[i]][25][[1]])|is.na(modListNoEvo[[i]][26][[1]])|is.nan(modListNoEvo[[i]][26][[1]])|is.infinite(modListNoEvo[[i]][26][[1]])){
      outlier1=outlier1+1
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<=natLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]<=natLev){ # Native more abundant -> Native more abundant
      # native abundant                                     GOING TO native abundant
      changeEvo[1,1]<-changeEvo[1,1]+1
      listerA=c(listerA,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<=natLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]>=invLev){ # Native more abundant -> Invasive more abundant
      # native abundant                                     GOING TO invasive abundant
      changeEvo[2,1]<-changeEvo[2,1]+1
      listerB=c(listerB,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>=invLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]>=invLev){ # Invasive more abundant -> Invasive more abundant
      # invasive abundant                                     GOING TO invasive abundant
      changeEvo[1,2]<-changeEvo[1,2]+1
      listerC=c(listerC,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>=invLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]<=natLev){ # Invasive more abundant -> Native more abundant
      # invasive abundant                                     GOING TO native abundant
      changeEvo[2,2]<-changeEvo[2,2]+1
      listerD=c(listerD,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<=natLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]>0.8&modList[[i]][25][[1]]/modList[[i]][26][[1]]<1.2){ # Native more abundant -> similar abundance
      # native abundant                                     GOING TO similar abundant
      changeEvo[1,3]<-changeEvo[1,3]+1
      listerE=c(listerE,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>=invLev&modList[[i]][25][[1]]/modList[[i]][26][[1]]>0.8&modList[[i]][25][[1]]/modList[[i]][26][[1]]<1.2){ # Invasive more abundant -> similar abundance
      # invasive abundant                                     GOING TO similar abundant
      changeEvo[2,3]<-changeEvo[2,3]+1
      listerF=c(listerF,i)
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>0.8&modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<1.2&modList[[i]][25][[1]]/modList[[i]][26][[1]]>0.8&modList[[i]][25][[1]]/modList[[i]][26][[1]]<1.2){ # similar abundance -> similar abundance
      # similar abundant                                     GOING TO similar abundant
      listerG=c(listerG,i) #inv wins to inv wins
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>0.8&modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<1.2&modList[[i]][25][[1]]/modList[[i]][26][[1]]<=natLev){ # similar abundance -> Native more abundant
      # similar abundant                                     GOING TO native abundant
      listerH=c(listerH,i) #coexist to coexist
    } else if(modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]>0.8&modListNoEvo[[i]][25][[1]]/modListNoEvo[[i]][26][[1]]<1.2&modList[[i]][25][[1]]/modList[[i]][26][[1]]>=invLev){ # similar abundance -> Invasive more abundant
      # similar abundant                                     GOING TO invasive abundant
      listerI=c(listerI,i) #nat wins to nat wins
    } else {
      outlier2=outlier2+1
    }
  }
  
  return(list(listerA,listerB,listerC,listerD,listerE,listerF,listerG,listerH,listerI,outlier1,outlier2))
  
}



########################################################################
########################################################################
########################################################################
########################################################################
#- 2) Native species new equilibrium w/ invader w/ evolution vs w/o evolution
#- native species equilibrium higher: evo/noevo=1.1
#- native species equilibrium lower: evo/noevo=0.91
#- native specie equilibrium same: 0.91<evo/noevo<1.1
#- native species equilibrium totally extinct with and without

runScen2=function(modList,modListNoEvo,incLev=1.1,decLev=0.91){
  
  outlier1<-0
  outlier2<-0
  listerA<-numeric(0) # native species equilibrium higher
  listerB<-numeric(0) # native species equilibrium lower
  listerC<-numeric(0) # native species equilibrium same- non-extinct
  listerD<-numeric(0) # native species equilibrium same- extinct
  listerE<-numeric(0) #
  listerF<-numeric(0) # 
  listerG<-numeric(0) # 
  listerH<-numeric(0) # 
  listerI<-numeric(0) 
  
  samp1<-numeric(0)
  changeEvo<-matrix(0, nrow=2,ncol=3)
  for(i in 1:length(modList)){
    if(is.na(modListNoEvo[[i]][25][[1]])|is.nan(modListNoEvo[[i]][25][[1]])|is.infinite(modListNoEvo[[i]][25][[1]])|is.na(modListNoEvo[[i]][26][[1]])|is.nan(modListNoEvo[[i]][26][[1]])|is.infinite(modListNoEvo[[i]][26][[1]])){
      outlier1=outlier1+1
    } else if(modList[[i]][26][[1]]==0&modListNoEvo[[i]][26][[1]]==0) {
      listerD=c(listerD,i)
    } else if(modList[[i]][26][[1]]/modListNoEvo[[i]][26][[1]]>=incLev){ # Native gets more abundant
      changeEvo[1,1]<-changeEvo[1,1]+1
      listerA=c(listerA,i)
    } else if(modList[[i]][26][[1]]/modListNoEvo[[i]][26][[1]]<=decLev){ # Native gets less abundant
      changeEvo[2,1]<-changeEvo[2,1]+1
      listerB=c(listerB,i)
    } else if(modList[[i]][26][[1]]/modListNoEvo[[i]][26][[1]]<incLev&modList[[i]][26][[1]]/modListNoEvo[[i]][26][[1]]>decLev){ # Native abundance is same
      changeEvo[1,2]<-changeEvo[1,2]+1
      listerC=c(listerC,i)
    } else {
      outlier2=outlier2+1
    }
  }
  
  return(list(listerA,listerB,listerC,listerD,outlier1,outlier2))
  
}





########################################################################
########################################################################
########################################################################
########################################################################
#- 3) Invasion threshold
#- invasive species equilibrium higher: evo/noevo=1.1
#- invasive species equilibrium lower: evo/noevo=0.91
#- invasive species equilibrium same: 0.91<evo/noevo<1.1

runScen3=function(modList,modListNoEvo,incLev=1.1,decLev=0.91){
  
  outlier1<-0
  outlier2<-0
  listerA<-numeric(0) # invasive species equilibrium higher
  listerB<-numeric(0) # invasive species equilibrium lower
  listerC<-numeric(0) # invasive species equilibrium same- non-extinct
  listerD<-numeric(0) # invasive species equilibrium same- extinct
  listerE<-numeric(0) #
  listerF<-numeric(0) # 
  listerG<-numeric(0) # 
  listerH<-numeric(0) # 
  listerI<-numeric(0) #
  
  samp1<-numeric(0)
  changeEvo<-matrix(0, nrow=2,ncol=3)
  for(i in 1:length(modList)){
    if(is.na(modListNoEvo[[i]][25][[1]])|is.nan(modListNoEvo[[i]][25][[1]])|is.infinite(modListNoEvo[[i]][25][[1]])|is.na(modListNoEvo[[i]][26][[1]])|is.nan(modListNoEvo[[i]][26][[1]])|is.infinite(modListNoEvo[[i]][26][[1]])){
      outlier1=outlier1+1
    } else if(modList[[i]][25][[1]]==0&modListNoEvo[[i]][25][[1]]==0) {
      listerD=c(listerD,i)
    } else if(modList[[i]][25][[1]]/modListNoEvo[[i]][25][[1]]>=incLev){ # Invasive gets more abundant
      changeEvo[1,1]<-changeEvo[1,1]+1
      listerA=c(listerA,i)
    } else if(modList[[i]][25][[1]]/modListNoEvo[[i]][25][[1]]<=decLev){ # Invasive gets less abundant
      changeEvo[2,1]<-changeEvo[2,1]+1
      listerB=c(listerB,i)
    } else if(modList[[i]][25][[1]]/modListNoEvo[[i]][25][[1]]<incLev&modList[[i]][25][[1]]/modListNoEvo[[i]][25][[1]]>decLev){ # Invasive abundance is same
      changeEvo[1,2]<-changeEvo[1,2]+1
      listerC=c(listerC,i)
    } else {
      outlier2=outlier2+1
    }
  }
  return(list(listerA,listerB,listerC,listerD,outlier1,outlier2))
}




########################################################################
########################################################################
########################################################################
########################################################################
#- 4) New Sept 25: Total population increase or decrease
#- total pop higher: evo/noevo=1.1
#- total pop lower: evo/noevo=0.91
#- total pop same: 0.91<evo/noevo<1.1


runScen4=function(modList,modListNoEvo,incLev=1.1,decLev=0.91){
  
  outlier1<-0
  outlier2<-0
  listerA<-numeric(0) # total pop higher
  listerB<-numeric(0) # total pop lower
  listerC<-numeric(0) # total pop same
  listerD<-numeric(0) # 
  listerE<-numeric(0) #
  listerF<-numeric(0) # 
  listerG<-numeric(0) # 
  listerH<-numeric(0) # 
  listerI<-numeric(0) #
  
  samp1<-numeric(0)
  changeEvo<-matrix(0, nrow=2,ncol=2)
  for(i in 1:length(modList)){
    if(is.na(modListNoEvo[[i]][25][[1]])|is.nan(modListNoEvo[[i]][25][[1]])|is.infinite(modListNoEvo[[i]][25][[1]])|is.na(modListNoEvo[[i]][26][[1]])|is.nan(modListNoEvo[[i]][26][[1]])|is.infinite(modListNoEvo[[i]][26][[1]])){
      outlier1=outlier1+1
    } else if((modList[[i]][25][[1]]+modList[[i]][26][[1]])/(modListNoEvo[[i]][25][[1]]+modListNoEvo[[i]][26][[1]])>=incLev){ # Invasive gets more abundant
      changeEvo[1,1]<-changeEvo[1,1]+1
      listerA=c(listerA,i)
    } else if((modList[[i]][25][[1]]+modList[[i]][26][[1]])/(modListNoEvo[[i]][25][[1]]+modListNoEvo[[i]][26][[1]])<=decLev){ # Invasive gets less abundant
      changeEvo[2,1]<-changeEvo[2,1]+1
      listerB=c(listerB,i)
    } else if((modList[[i]][25][[1]]+modList[[i]][26][[1]])/(modListNoEvo[[i]][25][[1]]+modListNoEvo[[i]][26][[1]])<incLev&(modList[[i]][25][[1]]+modList[[i]][26][[1]])/(modListNoEvo[[i]][25][[1]]+modListNoEvo[[i]][26][[1]])>decLev){ # Invasive abundance is same
      changeEvo[1,2]<-changeEvo[1,2]+1
      listerC=c(listerC,i)
    } else {
      outlier2=outlier2+1
    }
  }
  
  return(list(listerA,listerB,listerC,outlier1,outlier2))
}


makePcaMat=function(modList,scenRes){
  set.seed(934)
  maxLetter=length(scenRes)-2
  
  # start new data frame
  if(length(scenRes[[1]])==0){
    newFrame=as.data.frame(rbind(extractPars(modList[[scenRes[[2]][1]]]),extractPars(modList[[scenRes[[2]][2]]])))
    newFrame$Sample[1]<-scenRes[[2]][1]
    newFrame$Sample[2]<-scenRes[[2]][2]
    newFrame$StbCat[1]<-"B"
    newFrame$StbCat[2]<-"B"
    startnum=3
    Bstart=2
  } else if(length(scenRes[[1]])==1){
    newFrame=as.data.frame(rbind(extractPars(modList[[scenRes[[1]][1]]]),extractPars(modList[[scenRes[[2]][1]]])))
    newFrame$Sample[1]<-scenRes[[1]][1]
    newFrame$Sample[2]<-scenRes[[2]][1]
    newFrame$StbCat[1]<-"A"
    newFrame$StbCat[2]<-"B"
    startnum=3
    Bstart=1
  } else {
    newFrame=as.data.frame(rbind(extractPars(modList[[scenRes[[1]][1]]]),extractPars(modList[[scenRes[[1]][2]]])))
    newFrame$Sample[1]<-scenRes[[1]][1]
    newFrame$Sample[2]<-scenRes[[1]][2]
    newFrame$StbCat[1]<-"A"
    newFrame$StbCat[2]<-"A"
    Bstart=0
    if(length(scenRes[[1]])>2){
      for(i in 3:length(scenRes[[1]])){
        
        newFrame<-rbind(newFrame,c(extractPars(modList[[scenRes[[1]][i]]]),0,0))
        newFrame$Sample[i]<-scenRes[[1]][i]
        newFrame$StbCat[i]<-"A"
      }
      startnum<-nrow(newFrame)
    }
  }
  
  if(length(scenRes[[2]])>Bstart){
    for(i in (1+Bstart):length(scenRes[[2]])){
      
      newFrame<-rbind(newFrame,c(extractPars(modList[[scenRes[[2]][i]]]),0,0))
      newFrame$Sample[i+startnum]<-scenRes[[2]][i]
      newFrame$StbCat[i+startnum]<-"B"
      
    }
    
    startnum<-nrow(newFrame)
  }
  
  lettersAll=c("A","B","C","D","E","F","G","H","I")
  
  for(k in 3:maxLetter){
    if(length(scenRes[[k]])>0){
      for(i in 1:length(scenRes[[k]])){
        
        newFrame<-rbind(newFrame,c(extractPars(modList[[scenRes[[k]][i]]]),0,0))
        newFrame$Sample[i+startnum]<-scenRes[[k]][i]
        newFrame$StbCat[i+startnum]<-lettersAll[k]
        
      }
      
      startnum<-nrow(newFrame)
    }
  }
  
  
  return(newFrame)
  
}

getCommMat=function(newFrame){
  newFrameCom<-newFrame[,1:18]
  CommMat<-as.matrix(newFrameCom)
  return(CommMat)
}


pcaPlotChange1=function(newFrame){
  newFrame2=newFrame[-which(newFrame$StbCat=="A"),]
  newFrame2=newFrame2[-which(newFrame2$StbCat=="C"),]
  if(length(which(newFrame2$StbCat=="G"))>0){
    newFrame2=newFrame2[-which(newFrame2$StbCat=="G"),]
  }


  newFrameCom2<-newFrame2[,1:18]

  CommMat2<-as.matrix(newFrameCom2)
  CommMat2.pca<-prcomp(CommMat2[,c(1,3,5:8,13,15)],center = T,scale. = T) # remove fixed parameters
  return(list(CommMat2,CommMat2.pca,newFrame2))
}


########################################################################
########################################################################
########################################################################
# add second plot that focuses only on those that CHANGE
# Remove C, D
pcaPlotChange2=function(newFrame){
  newFrame2=newFrame[-which(newFrame$StbCat=="C"),]
  newFrame2=newFrame2[-which(newFrame2$StbCat=="D"),]
  
  newFrameCom2<-newFrame2[,1:18]
  
  CommMat2<-as.matrix(newFrameCom2)
  CommMat2.pca<-prcomp(CommMat2[,c(1,3,5:8,13,15)],center = T,scale. = T) # remove fixed parameters
  return(list(CommMat2,CommMat2.pca,newFrame2))
}



########################################################################
########################################################################
########################################################################
# add second plot that focuses only on those that CHANGE
# Remove C
pcaPlotChange4=function(newFrame){
  newFrame2=newFrame[-which(newFrame$StbCat=="C"),]
  
  
  newFrameCom2<-newFrame2[,1:18]
  
  CommMat2<-as.matrix(newFrameCom2)
  CommMat2.pca<-prcomp(CommMat2[,c(1,3,5:8,13,15)],center = T,scale. = T) # remove fixed parameters
  return(list(CommMat2,CommMat2.pca,newFrame2))
}

# function that adds an ellipse to a plot
addEllipse=function(x,y,col="black",lty=1,lwd=1,verbose=F){
  el1=confidence_ellipse(as.data.frame(cbind(x,y)),x="x",y="y")
  lines(el1,type="l",col=col,lty=lty,lwd=lwd)
  if(verbose){
    return(el1)
  }
}

# function that adds all the ellipses
pcaEllipse=function(matrixDat,categories,col,lty=1,lwd=1,verbose=F){
  allEllipses=list()
  for(i in 1:length(unique(categories))){
    xi=matrixDat[which(categories==unique(categories)[i]),1]
    yi=matrixDat[which(categories==unique(categories)[i]),2]
    if(length(xi)>1){
      allEllipses[[i]]=addEllipse(xi,yi,col=col[i],lty=lty,lwd=lwd,verbose=T)
    }
  }
  if(verbose){
    return(allEllipses)
  }
}

pcaEllipse2=function(commMatAll,col,lty=1,lwd=1,verbose=F){
  matrixDat=commMatAll[[2]]$x[,1:2]
  categories=commMatAll[[3]]$StbCat
  allEllipses=list()
  for(i in 1:length(unique(categories))){
    xi=matrixDat[which(categories==unique(categories)[i]),1]
    yi=matrixDat[which(categories==unique(categories)[i]),2]
    if(length(xi)>1){
      allEllipses[[i]]=addEllipse(xi,yi,col=col[i],lty=lty,lwd=lwd,verbose=T)
    }
  }
  if(verbose){
    return(allEllipses)
  }
}

# add ellipse, on different axes
pcaEllipseDif=function(commMatPoints,commMatAxes,col,lty=1,lwd=1,verbose=F){
  
  meanScaler=matrix(rep(commMatAxes[[2]]$center,nrow(commMatPoints[[1]])),byrow=T,nrow=nrow(commMatPoints[[1]]))
  sdScaler=matrix(rep(commMatAxes[[2]]$scale,nrow(commMatPoints[[1]])),byrow=T,nrow=nrow(commMatPoints[[1]]))
  
  matrixDat=((commMatPoints[[1]][,c(1,3,5:8,13,15)]-meanScaler)/sdScaler)%*%commMatAxes[[2]]$rotation
  
  categories=commMatPoints[[3]]$StbCat
  allEllipses=list()
  for(i in 1:length(unique(categories))){
    xi=matrixDat[which(categories==unique(categories)[i]),1]
    yi=matrixDat[which(categories==unique(categories)[i]),2]
    if(length(xi)>1){
      allEllipses[[i]]=addEllipse(xi,yi,col=col[i],lty=lty,lwd=lwd,verbose=T)
    }
  }
  if(verbose){
    return(allEllipses)
  }
}


# plots PCA on axes but not with new axes
plotOnAxes2=function(commMatAll){
  coordinates1=commMatAll[[2]]$x[,1:2]
  categories=commMatAll[[3]]
  plot(coordinates1,pch=categories$StbCat,col=rainbow(5)[as.factor(categories$StbCat)],
       ylim=c(min(min(coordinates1[,2])-1.5,-2),max(max(coordinates1[,2])+1.5,2)),xlim=c(min(min(coordinates1[,1])-1.5,-2),max(max(coordinates1[,1])+1.5,2)),axes=F)
  abline(h=0)
  abline(v=0)
}

# plots PCA on axes with NEW axes
plotOnAxesDif=function(commMatPoints,commMatAxes,type="p",pbuff=1.5,ylim_min=-2,ylim_max=2,xlim_min=-2,xlim_max=2,shrinkaxes=NA,ylabopt="PC2",xlabopt="PC1",cex=1){
  allPoints=commMatPoints[[1]][,c(1,3,5:8,13,15)]
  allAxes=commMatAxes[[2]]$rotation
  meanScaler=matrix(rep(commMatAxes[[2]]$center,nrow(commMatPoints[[1]])),byrow=T,nrow=nrow(commMatPoints[[1]]))
  sdScaler=matrix(rep(commMatAxes[[2]]$scale,nrow(commMatPoints[[1]])),byrow=T,nrow=nrow(commMatPoints[[1]]))
  
  coordinates1=((allPoints-meanScaler)/sdScaler)%*%allAxes
  categories=commMatPoints[[3]]
  ylim1=c(min(min(coordinates1[,2])-pbuff,ylim_min),max(max(coordinates1[,2])+pbuff,ylim_max))
  xlim1=c(min(min(coordinates1[,1])-pbuff,xlim_min),max(max(coordinates1[,1])+pbuff,xlim_max))
  par(mgp = c(0.8, 0.7, 0))
  plot(coordinates1,pch=categories$StbCat,col=rainbow(5)[as.factor(categories$StbCat)],
       ylim=ylim1,xlim=xlim1,axes=F,type=type,ylab=ylabopt,xlab=xlabopt,cex.lab=cex)
  if(is.na(shrinkaxes)){
    abline(h=0)
    abline(v=0)
  } else {
    segments(0, ylim1[1]*shrinkaxes, 0, ylim1[2]*shrinkaxes, col = "black", lwd = 1, lty = 1)
    segments(xlim1[1]*shrinkaxes, 0, xlim1[2]*shrinkaxes, 0, col = "black", lwd = 1, lty = 1)
  }
  mgp = c(3, 1, 0)
}

addArrows=function(rotation,mult=2,lwd=1,tsp=1,labels=T,fancy=F,fancylabels=NA,cex=1){
  for(i in 1:nrow(rotation)){
    arrows(0,0,rotation[i,1]*mult,rotation[i,2]*mult,length=0.1,lwd=lwd)
    # add text close by
    if(labels){
      
      if(fancy){ # fancy labels
        text(rotation[i,1]*mult*tsp,rotation[i,2]*mult*tsp,fancylabels[i],cex=cex)
      } else { # normal labels
        text(rotation[i,1]*mult*tsp,rotation[i,2]*mult*tsp,rownames(rotation)[i],cex=cex)
      }
    }
  }
}


tryAllAxes=function(rotation,scenList=scenList,type="p",fancy=F,fancylabels=NA){
  par(mfrow=c(4,4))
  par(mar=c(0,0,0,0))
  for(i in 1:length(scenList)){
    plotOnAxesDif(scenList[[i]],rotation,type=type)
    pcaEllipseDif(scenList[[i]],rotation,col=rainbow(5))
    addArrows(rotation[[2]]$rotation,fancy=fancy,fancylabels=fancylabels)
  }
  par(mfrow=c(1,1))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}


screeplot=function(pcaCommunity){
  par(mfrow=c(1,1))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  eigenvalues1=pcaCommunity[[2]]$sdev^2
  plot(eigenvalues1/sum(eigenvalues1),xlab="PC number",ylab="percent variance explained")
  return("done")
}


# draw a point inside one ellipse
drawPoint=function(ellipseCoord){
  ellipseCoord[nrow(ellipseCoord),]=ellipseCoord[1,]
  polygon=st_polygon(list(as.matrix(ellipseCoord)))
  polygon_sf <- st_sfc(polygon)
  bbox=st_bbox(polygon_sf)
  goodPoint=F
  while(!goodPoint){
    x=runif(1,bbox["xmin"],bbox["xmax"])
    y=runif(1,bbox["ymin"],bbox["ymax"])
    point=st_point(c(x,y))
    
    if(st_within(st_sfc(point),polygon_sf,sparse=F)){
      goodPoint=T
    }
  }
  return(c(x,y))
}

# draw a point inside one ellipse and outside a second ellipse
drawPointOne=function(ellipseCoordIn,ellipseCoordOut){
  
  ellipseCoordIn[nrow(ellipseCoordIn),]=ellipseCoordIn[1,]
  polygon1=st_polygon(list(as.matrix(ellipseCoordIn)))
  polygon1_sf <- st_sfc(polygon1)
  
  ellipseCoordOut[nrow(ellipseCoordOut),]=ellipseCoordOut[1,]
  polygon2=st_polygon(list(as.matrix(ellipseCoordOut)))
  polygon2_sf <- st_sfc(polygon2)
  
  bbox=st_bbox(polygon1_sf)
  goodPoint=F
  while(!goodPoint){
    x=runif(1,bbox["xmin"],bbox["xmax"])
    y=runif(1,bbox["ymin"],bbox["ymax"])
    point=st_point(c(x,y))
    
    if(st_within(st_sfc(point),polygon1_sf,sparse=F) & !st_within(st_sfc(point),polygon2_sf,sparse=F)){
      goodPoint=T
    }
  }
  return(c(x,y))
}

# draw a point inside both ellipses
drawPointBoth=function(ellipseCoord1,ellipseCoord2){
  
  ellipseCoord1[nrow(ellipseCoord1),]=ellipseCoord1[1,]
  polygon1=st_polygon(list(as.matrix(ellipseCoord1)))
  polygon1_sf <- st_sfc(polygon1)
  
  ellipseCoord2[nrow(ellipseCoord2),]=ellipseCoord2[1,]
  polygon2=st_polygon(list(as.matrix(ellipseCoord2)))
  polygon2_sf <- st_sfc(polygon2)
  
  bbox=st_bbox(polygon1_sf)
  goodPoint=F
  while(!goodPoint){
    x=runif(1,bbox["xmin"],bbox["xmax"])
    y=runif(1,bbox["ymin"],bbox["ymax"])
    point=st_point(c(x,y))
    
    if(st_within(st_sfc(point),polygon1_sf,sparse=F) & st_within(st_sfc(point),polygon2_sf,sparse=F)){
      goodPoint=T
    }
  }
  return(c(x,y))
}

#' Title Get parameter values from a point on PC axes 1 and 2
#'
#' @param point the point, using PCA axes 1 and 2 
#' @param fullRot the rotation used transform to the PCA axes
#' @param scale the vector of standard deviations used to transform each parameter before determining PCA rotation
#' @param center the vector of means used to transform each parameter before determining PCA rotation
#' @param method the method used in random draws that determine the weighting of the PCA axes beyond 1 and 2
#' @param maxTime cut-off time for finding non-negative parameter values
#'
#' @return vector of fully back-transformed parameter values
#' @export
#'
#' @examples
getFullPCcoords=function(point,fullRot,scale,center,method="uniform",maxTime=60){
  parPC12=fullRot[,1:2]%*%point
  startTime=Sys.time()
  if(method=="uniform"){
    parPC38=fullRot[,3:nrow(fullRot)]%*%runif((nrow(fullRot)-2),-1,1)
    par_untransformed=parPC12+parPC38
    par_transformed=par_untransformed*scale+center
    while((any(par_transformed<0)|par_transformed[3]>1)&(as.numeric(Sys.time()-startTime,units="secs")<maxTime)){
      parPC38=fullRot[,3:nrow(fullRot)]%*%runif((nrow(fullRot)-2),-1,1)
      par_untransformed=parPC12+parPC38
      par_transformed=par_untransformed*scale+center
    }
  }
  if(method=="gaussian"){
    parPC38=fullRot[,3:nrow(fullRot)]%*%rnorm((nrow(fullRot)-2),0,0.5)
    par_untransformed=parPC12+parPC38
    par_transformed=par_untransformed*scale+center
    while((any(parPC12+parPC38<0)|par_transformed[3]>1)&(as.numeric(Sys.time()-startTime,units="secs")<maxTime)){
      parPC38=fullRot[,3:nrow(fullRot)]%*%rnorm((nrow(fullRot)-2),0,0.5)
      par_untransformed=parPC12+parPC38
      par_transformed=par_untransformed*scale+center
    }
  }
  if(any(par_transformed<0)){
    return(NA)
  } else {
    return(par_transformed)
  }
}


NatChangeGrid=function(gridResult,gridResultNoEvo,parameter){
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
        # ok so matrix1 looks like [invasive parameter, native parameter]
      }
    }
  }
  
  # NEEED to check axes still # looking into this Aug 30 2025
  # each gridResult object comes from the sweepSolverModGM3 function
  # when using image(), the columns of the matrix will be plotted as rows of the image
  # this is like x, y
  # so Sept 19 2025- this was verified in the last run through and I think figures have been finalized now
  
  par(xpd=F)
  if(parameter=="sigma") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(sigma["inv"]),ylab=expression(sigma["nat"]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="tau"){
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(tau["inv"]),ylab=expression(tau["nat"]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="b") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(b["inv"]),ylab=expression(b["nat"]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
  } else if(parameter=="m") {
    image(1:gridSize,1:gridSize,z=matrix1,zlim=c(1,4),xlab=expression(m["inv"]),ylab=expression(m["nat"]),col=c("seagreen","thistle1","skyblue","salmon1"),axes="F")
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
  
  # segments(x0=value1-0.5*by1,x1=value1+0.5*by1,y0=value2+0.5*by2,y1=value2+0.5*by2,lwd=3,col="blue")
  # segments(value1-0.5*by1,value1+0.5*by1,value2-0.5*by2,value2-0.5*by2,lwd=3)
  # segments(value1-0.5*by1,value1-0.5*by1,value2-0.5*by2,value2+0.5*by2,lwd=3)
  # segments(value1+0.5*by1,value1+0.5*by1,value2-0.5*by2,value2+0.5*by2,lwd=3)
  # 

  
}

# 0-either
# 1-1
# 2-2
# 3-both
# 4-1 only
# 5-2 only
trimGrid=function(xseq,yseq,ellipseCoords,single=0){
  
  ellipseCoord1=ellipseCoords[[1]]
  ellipseCoord2=ellipseCoords[[2]]
  
  ellipseCoord1[nrow(ellipseCoord1),]=ellipseCoord1[1,]
  polygon1=st_polygon(list(as.matrix(ellipseCoord1)))
  polygon1_sf <- st_sfc(polygon1)
  
  ellipseCoord2[nrow(ellipseCoord2),]=ellipseCoord2[1,]
  polygon2=st_polygon(list(as.matrix(ellipseCoord2)))
  polygon2_sf <- st_sfc(polygon2)
  
  keepPoint=matrix(NA,nrow=length(xseq),ncol=length(yseq))
  if(single==0){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(st_within(st_sfc(point),polygon1_sf,sparse=F)|st_within(st_sfc(point),polygon2_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  } else if(single==1){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(st_within(st_sfc(point),polygon1_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  } else if(single==2){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(st_within(st_sfc(point),polygon2_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  } else if(single==3){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(st_within(st_sfc(point),polygon1_sf,sparse=F)&st_within(st_sfc(point),polygon2_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  } else if(single==4){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(st_within(st_sfc(point),polygon1_sf,sparse=F)&!st_within(st_sfc(point),polygon2_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  } else if(single==5){
    for(i in 1:length(xseq)){ # i is x
      for(j in 1:length(yseq)){ # j is y
        point=st_point(c(xseq[i],yseq[j]))
        if(!st_within(st_sfc(point),polygon1_sf,sparse=F)&st_within(st_sfc(point),polygon2_sf,sparse=F)){
          keepPoint[i,j]=T
        } else {
          keepPoint[i,j]=F
        }
      }
    }
  }
  
  return(keepPoint)
}


extractInvPop=function(modlist,ntrials=100){
  allPops=numeric(ntrials)
  for(i in 1:ntrials){
    allPops[i]=modlist[[i]][25][[1]]
  }
  return(allPops)
}

extractNatPop=function(modlist,ntrials=100){
  allPops=numeric(ntrials)
  for(i in 1:ntrials){
    allPops[i]=modlist[[i]][26][[1]]
  }
  return(allPops)
}

extractInvTrait=function(modlist,ntrials=100){
  allPops=numeric(ntrials)
  for(i in 1:ntrials){
    allPops[i]=modlist[[i]][23][[1]]
  }
  return(allPops)
}

extractNatTrait=function(modlist,ntrials=100){
  allPops=numeric(ntrials)
  for(i in 1:ntrials){
    allPops[i]=modlist[[i]][24][[1]]
  }
  return(allPops)
}

extractA=function(modlist,ntrials=100){
  allPops=numeric(ntrials)
  for(i in 1:ntrials){
    allPops[i]=modlist[[i]][5][[1]]
  }
  return(allPops)
}

# let commMatPointsList be a list of all the points
pcaGridSearch=function(commMatPointsList,commMatAxes,gridLength=16, ntrials=100, trim=T, nsplit=1, splitsec=1){
  
  allx=numeric(0)
  ally=numeric(0)
  # get and plot ellipse coordinates
  for(i in 1:length(commMatPointsList)){
    commMatPoints=commMatPointsList[[i]]
    plotOnAxesDif(commMatPoints,commMatAxes)
    ellipseCoords=pcaEllipseDif(commMatPoints,commMatAxes,col=rainbow(5),verbose=T)
    allx=c(allx,ellipseCoords[[1]]$x,ellipseCoords[[2]]$x)
    ally=c(ally,ellipseCoords[[1]]$y,ellipseCoords[[2]]$y)
  }
  
  # define grid outline by the max and min x and y coordinates of ALL the ellipses
  max_x=max(allx)
  min_x=min(allx)
  max_y=max(ally)
  min_y=min(ally)
  
  xseq=seq(min_x,max_x,length=gridLength)
  yseq=seq(min_y,max_y,length=gridLength)
  
  # trim grid to remove points outside ellipses
  if(trim){
    trimmer=trimGrid(xseq,yseq,ellipseCoords)
  } else {
    trimmer=matrix(TRUE,nrow=length(xseq),ncol=length(yseq))
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
        if(trimmer[k,j]){ # if the point is inside an ellipse, go
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
        } else { # send back just NA for the whole list if we are outside ellipses
          smallList<-NA
          smallListNoEvo<-NA
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

recoverParLine1=function(purplepoint,jperc,focalparam="b",nfig4ALL=nfig4ALL,modellistFig2=modellistFig2,modellistFig2NoEvo=modellistFig2NoEvo,plotout=T){
  i=purplepoint
  j=jperc
  
  resb=readRDS(paste0('/Users/amypatterson/Documents/Grad stuff/A Research etc/GM project/Vasseur paper work/Figures with invasive/big scan Jan 30 2024/PCAgridRes/point_',focalparam,nfig4ALL[i],'_',j,'.rds'))
  resbNoEvo=readRDS(paste0('/Users/amypatterson/Documents/Grad stuff/A Research etc/GM project/Vasseur paper work/Figures with invasive/big scan Jan 30 2024/PCAgridRes/point_',focalparam,nfig4ALL[i],'_',j,'NoEvo.rds'))
  
  focalparam_inital=resb[[16]][1]
  # add in graphical stuff starting here
  if(plotout){
    mat=matrix(c(1,2,4,6,
                 1,3,5,7),nrow=2,byrow=T)
    layout(mat, widths = c(2, 1,1,1), heights = c(1, 1))
    NatChangeGrid(resb,resbNoEvo,"b")
  }

  
  nat150=extractNatPop(modellistFig2[[nfig4ALL[i]]],ntrials=150)
  natnoevo150=extractNatPop(modellistFig2NoEvo[[nfig4ALL[i]]],ntrials=150)
  
  natPopChange1=mapply(function(nat150,natnoevo150){
    if(is.na(natnoevo150)|is.na(nat150)){
      return(NA)
    } else if(natnoevo150==0){
      if(nat150>0){
        return(9999) # just give a high number instead of inf, note that taking the mean will now be pointless
      }
      if(nat150==0){
        return(1) # zero to zero, change should be considered 1
      }
    } else {
      return(nat150/natnoevo150)
    }
  },nat150,natnoevo150)
  
  percentileInd=numeric(3)
  percentileVal=numeric(3)
  
  percentileInd[1]=round(0.05*length(natPopChange1))
  percentileVal[1]=sort(natPopChange1)[percentileInd[1]]
  percentileInd[2]=round(0.50*length(natPopChange1))
  percentileVal[2]=sort(natPopChange1)[percentileInd[2]]
  percentileInd[3]=round(0.95*length(natPopChange1))
  percentileVal[3]=sort(natPopChange1)[percentileInd[3]]
  
  # works only if you don't need to randomly draw- make a new thing for randomly drawing between population change values that are the same\
  
  if(length(which(natPopChange1==percentileVal[j]))==1){
    par_recovered=extractPars(modellistFig2[[nfig4ALL[i]]][which(natPopChange1==percentileVal[j])][[1]])
  } else {
    
    v=which(natPopChange1==percentileVal[j])
    
    results <- lapply(v, function(k) {
      extractPars(modellistFig2[[nfig4ALL[i]]][[k]])
    })
    
    comparepar=which(names(results[[1]])==focalparam)
    
    candidate_pars=sapply(results, function(x) x[comparepar])
    which(candidate_pars==focalparam_inital)
    
    par_recovered=extractPars(modellistFig2[[nfig4ALL[i]]][v[1]][[1]])
    
  }
  
  
  return(par_recovered)
  
}


recoverParLinePlot=function(parameters,focalparam,focalcell,altset1,altset2,tmax=5000){
  
  
  segments(x0=(focalcell[1]+altset1[1])-0.5,x1=(focalcell[1]+altset1[1])+0.5,y0=(focalcell[2]+altset1[2]+1)-0.5,y1=(focalcell[2]+altset1[2]+1)-0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset1[1])-0.5,x1=(focalcell[1]+altset1[1])+0.5,y0=(focalcell[2]+altset1[2]+1)+0.5,y1=(focalcell[2]+altset1[2]+1)+0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset1[1])-0.5,x1=(focalcell[1]+altset1[1])-0.5,y0=(focalcell[2]+altset1[2]+1)-0.5,y1=(focalcell[2]+altset1[2]+1)+0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset1[1])+0.5,x1=(focalcell[1]+altset1[1])+0.5,y0=(focalcell[2]+altset1[2]+1)-0.5,y1=(focalcell[2]+altset1[2]+1)+0.5,lwd=3,col="black")
  
  
  segments(x0=(focalcell[1]+altset2[1])-0.5,x1=(focalcell[1]+altset2[1])+0.5,y0=(focalcell[2]+altset2[2]+1)-0.5,y1=(focalcell[2]+altset2[2]+1)-0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset2[1])-0.5,x1=(focalcell[1]+altset2[1])+0.5,y0=(focalcell[2]+altset2[2]+1)+0.5,y1=(focalcell[2]+altset2[2]+1)+0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset2[1])-0.5,x1=(focalcell[1]+altset2[1])-0.5,y0=(focalcell[2]+altset2[2]+1)-0.5,y1=(focalcell[2]+altset2[2]+1)+0.5,lwd=3,col="black")
  segments(x0=(focalcell[1]+altset2[1])+0.5,x1=(focalcell[1]+altset2[1])+0.5,y0=(focalcell[2]+altset2[2]+1)-0.5,y1=(focalcell[2]+altset2[2]+1)+0.5,lwd=3,col="black")
  
  text(focalcell[1],focalcell[2]+1,1)
  text(focalcell[1]+altset1[1],focalcell[2]+altset1[2]+1,2)
  text(focalcell[1]+altset2[1],focalcell[2]+altset2[2]+1,3)
  
  comparepar=which(names(parameters)==focalparam)
  
  for(i in 1:3){
    if(i==1){
      par_alt=parameters
    }
    if(i==2){
      par_alt=parameters
      par_alt[comparepar]=par_alt[comparepar]+altset1[1]*0.25
      par_alt[comparepar+1]=par_alt[comparepar+1]+altset1[2]*0.25
    }
    if(i==3){
      par_alt=parameters
      par_alt[comparepar]=par_alt[comparepar]+altset2[1]*0.25
      par_alt[comparepar+1]=par_alt[comparepar+1]+altset2[2]*0.25
    }
    # checks 1: set optima for single species to be same as m, n
    parameters[which(names(parameters)=="thetaC")]<-parameters[which(names(parameters)=="m")]
    parameters[which(names(parameters)=="thetaH")]<-parameters[which(names(parameters)=="n")]
    
    #############################################
    # checks 2: solve for native species alone in order to set initial conditions
    # Note: we need to parameters first!
    state<-c(Nf=0, Nv=1, xbar=0,ybar=1) 
    
    stateSolve<-tryCatch({
      namer="withEvo"
      runGMSolver(state,par_alt,tmax=simlengthmin, deltat=0.1, repsiz=100, plotout=F, report=T, save=F,func=GMspecial,extinctlim=extinctlim)
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
        stateSolve<-runGMSolver(state=state,parameters=par_alt,tmax=simlengthstep,deltat=dtcur,report=T,plotout=F,func=func)
        nSampCur<-nSampCur+simlengthstep
      }
      
    }
    
    # use the solution from the native species only for the initial condition
    state<-c(Nf=(stateSolve[5,4]*invP0), Nv=stateSolve[5,4], xbar=(stateSolve[5,2]*invT0),ybar=stateSolve[5,2]) 
    
    
    # then run the final parameter graph
    runGMSolver(state=state,parameters=par_alt,tmax=tmax,deltat=0.1,report=T,plotout=T,func=func,plotoutadd=T)
    
  }
  
}


NatCG_alphas=function(gridResult,gridResultNoEvo,parameter,purplepoint,jperc,showsp="native"){
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
  
  recovered_pars1=recoverParLine1(purplepoint=purplepoint,jperc=jperc,focalparam = parameter,nfig4ALL=nfig4ALL,modellistFig2=modellistFig2,modellistFig2NoEvo=modellistFig2NoEvo,plotout=F) 
  
  altname_param=parameter
  if(altname_param=="tau"){
    altname_param="tauInter1"
  }
  if(altname_param=="sigma"){
    altname_param="sigmax"
  }
  
  incLev=1.1
  decLev=0.9
  matrix1=matrix(NA,nrow=gridSize,ncol=gridSize)
  for(i in 1:gridSize){
    for(j in 1:gridSize){
      param_ind1=which(names(recovered_pars1)==altname_param)
      
      recovered_pars1[param_ind1]=values1[i] #invasive parameter
      recovered_pars1[param_ind1+1]=values2[j] #native parameter
      
      finalInv=gridResult[[i]][[j]][[2]][5,1] # # invasive trait
      finalNat=gridResult[[i]][[j]][[2]][5,2] # # native trait
      print(finalInv)
      if(showsp=="native"){
        matrix1[i,j]=ridgeDif2(finalInv,finalNat,recovered_pars1)
        if(gridResult[[i]][[j]][[2]][5,4]==0){
          matrix1[i,j]=NA
        }
      }
      
      if(showsp=="invasive"){
        matrix1[i,j]=ridgeDif1(finalInv,finalNat,recovered_pars1)
        if(gridResult[[i]][[j]][[2]][5,3]==0){
          matrix1[i,j]=NA
        }
      }
      
      
    }
  }
  
  # NEEED to check axes still # looking into this Aug 30 2025
  # each gridResult object comes from the sweepSolverModGM3 function
  # when using image(), the columns of the matrix will be plotted as rows of the image
  # this is like x, y
  # so Sept 19 2025- this was verified in the last run through and I think figures have been finalized now
  
  par(xpd=F)
  matrix2=matrix1
  matrix2[is.na(matrix2)]=-10
  zlim=c(0,1)
  breaks=seq(0,2,length.out=13)
  breaks=c(-11,breaks)
  cols=c("gray",hcl.colors(12, "YlOrRd", rev = TRUE))
  if(parameter=="sigma") {
    image(1:gridSize,1:gridSize,z=matrix2,xlab=expression(sigma["inv"]),ylab=expression(sigma["nat"]),axes="F",col=cols,breaks=breaks,zlim=c(-11,2))
  } else if(parameter=="tau"){
    image(1:gridSize,1:gridSize,z=matrix2,xlab=expression(tau["inv"]),ylab=expression(tau["nat"]),axes="F",col=cols,breaks=breaks,zlim=c(-11,2))
  } else if(parameter=="b") {
    image(1:gridSize,1:gridSize,z=matrix2,xlab=expression(b["inv"]),ylab=expression(b["nat"]),axes="F",col=cols,breaks=breaks,zlim=c(-11,2))
  } else if(parameter=="m") {
    image(1:gridSize,1:gridSize,z=matrix2,xlab=expression(m["inv"]),ylab=expression(m["nat"]),axes="F",col=cols,breaks=breaks,zlim=c(-11,2))
  } else {
    image(1:gridSize,1:gridSize,z=matrix2,xlab=paste(parameter, "invasive"),ylab=paste(parameter, native),axes="F",col=cols,breaks=breaks)
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
