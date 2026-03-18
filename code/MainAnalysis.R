# New Jun 7 2024 edited November 2024

library(here)
library(ConfidenceEllipse)
library(sf)
library(plotrix)


source(here("code/SolverFunctions.R"))
source(here("code/PcaFunctions.R"))


modellistFeb16par1=readRDS(here("output/modellistFeb16par1.rds"))
modellistFeb16par2=readRDS(here("output/modellistFeb16par2.rds"))
modellistFeb16par3=readRDS(here("output/modellistFeb16par3.rds"))
modellistFeb16par4=readRDS(here("output/modellistFeb16par4.rds"))
modellistFeb16par5=readRDS(here("output/modellistFeb16par5.rds"))
modellistFeb16par6=readRDS(here("output/modellistFeb16par6.rds"))


modellistNoEvoFeb16par1=readRDS(here("output/modellistNoEvoFeb16par1.rds"))
modellistNoEvoFeb16par2=readRDS(here("output/modellistNoEvoFeb16par2.rds"))
modellistNoEvoFeb16par3=readRDS(here("output/modellistNoEvoFeb16par3.rds"))
modellistNoEvoFeb16par4=readRDS(here("output/modellistNoEvoFeb16par4.rds"))
modellistNoEvoFeb16par5=readRDS(here("output/modellistNoEvoFeb16par5.rds"))
modellistNoEvoFeb16par6=readRDS(here("output/modellistNoEvoFeb16par6.rds"))


# Prepare communities by categorizing output
# We initially looked at 
#Axes 4-7
#4
modList=modellistFeb16par1
modListNoEvo=modellistNoEvoFeb16par1
scen2=runScen1(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm4=pcaPlotChange1(pcaMat)
rot4=comm4[[2]]$rotation[,1:2]

#5
modList=modellistFeb16par2
modListNoEvo=modellistNoEvoFeb16par2
scen2=runScen1(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm5=pcaPlotChange1(pcaMat)
rot5=comm5[[2]]$rotation[,1:2]

#6
modList=modellistFeb16par3
modListNoEvo=modellistNoEvoFeb16par3
scen2=runScen1(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm6=pcaPlotChange1(pcaMat)
rot6=comm6[[2]]$rotation[,1:2]

#7
modList=modellistFeb16par4
modListNoEvo=modellistNoEvoFeb16par4
scen2=runScen1(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm7=pcaPlotChange1(pcaMat)
rot7=comm7[[2]]$rotation[,1:2]


# Axes 9-12
#9
modList=modellistFeb16par1
modListNoEvo=modellistNoEvoFeb16par1
scen2=runScen2(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm9=pcaPlotChange2(pcaMat)
rot9=comm9[[2]]$rotation[,1:2]

#10
modList=modellistFeb16par2
modListNoEvo=modellistNoEvoFeb16par2
scen2=runScen2(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm10=pcaPlotChange2(pcaMat)
rot10=comm10[[2]]$rotation[,1:2]

#11
modList=modellistFeb16par3
modListNoEvo=modellistNoEvoFeb16par3
scen2=runScen2(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm11=pcaPlotChange2(pcaMat)
rot11=comm11[[2]]$rotation[,1:2]

#12
modList=modellistFeb16par4
modListNoEvo=modellistNoEvoFeb16par4
scen2=runScen2(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm12=pcaPlotChange2(pcaMat)
rot12=comm12[[2]]$rotation[,1:2]


# Axes 14-17 Scen3
#14
modList=modellistFeb16par1
modListNoEvo=modellistNoEvoFeb16par1
scen2=runScen3(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen2)
commMat=getCommMat(pcaMat)
comm14=pcaPlotChange2(pcaMat)
rot14=comm14[[2]]$rotation[,1:2]

#15
modList=modellistFeb16par2
modListNoEvo=modellistNoEvoFeb16par2
scen1=runScen3(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen1)
commMat=getCommMat(pcaMat)
comm15=pcaPlotChange2(pcaMat)
rot15=comm15[[2]]$rotation[,1:2]

#16
modList=modellistFeb16par3
modListNoEvo=modellistNoEvoFeb16par3
scen1=runScen3(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen1)
commMat=getCommMat(pcaMat)
comm16=pcaPlotChange2(pcaMat)
rot16=comm16[[2]]$rotation[,1:2]

#17
modList=modellistFeb16par4
modListNoEvo=modellistNoEvoFeb16par4
scen1=runScen3(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen1)
commMat=getCommMat(pcaMat)
comm17=pcaPlotChange2(pcaMat)
rot17=comm17[[2]]$rotation[,1:2]


# Axes 21-23 Scen4
#21
modList=modellistFeb16par1
modListNoEvo=modellistNoEvoFeb16par1
scen4=runScen4(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen4)
commMat=getCommMat(pcaMat)
comm21=pcaPlotChange4(pcaMat)
rot21=comm21[[2]]$rotation[,1:2]

#22
modList=modellistFeb16par2
modListNoEvo=modellistNoEvoFeb16par2
scen4=runScen4(modList,modListNoEvo)
pcaMat=makePcaMat(modList,scen4)
commMat=getCommMat(pcaMat)
comm22=pcaPlotChange4(pcaMat)
rot22=comm22[[2]]$rotation[,1:2]

#23
modList=modellistFeb16par3
modListNoEvo=modellistNoEvoFeb16par3
scen4=runScen4(modList,modListNoEvo)
pcaMat4=makePcaMat(modList,scen4)
commMat4=getCommMat(pcaMat4)
comm23=pcaPlotChange4(pcaMat4)
rot23=comm23[[2]]$rotation[,1:2]

#24
modList=modellistFeb16par4
modListNoEvo=modellistNoEvoFeb16par4
scen4=runScen4(modList,modListNoEvo)
pcaMat4=makePcaMat(modList,scen4)
commMat4=getCommMat(pcaMat4)
comm24=pcaPlotChange4(pcaMat4)
rot24=comm24[[2]]$rotation[,1:2]



# add arrows for components of axes

scenList=list(comm4,comm5,comm6,comm7,
              comm9,comm10,comm11,comm12,
              comm14,comm15,comm16,comm17,
              comm21,comm22,comm23,comm24)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))

# Visualize outcomes on axes made with different evaluations (eg native population change vs invasive population change)
tryAllAxes(comm9,scenList,type="n",fancy=T,fancylabels = fancylabels)
tryAllAxes(comm10,scenList,type="n",fancy=T,fancylabels = fancylabels)
tryAllAxes(comm11,scenList)
tryAllAxes(comm12,scenList)
tryAllAxes(comm14,scenList)
tryAllAxes(comm15,scenList)
tryAllAxes(comm16,scenList)
tryAllAxes(comm17,scenList)
tryAllAxes(comm21,scenList)
tryAllAxes(comm22,scenList)
tryAllAxes(comm23,scenList)
tryAllAxes(comm24,scenList)

#################################################
# PCA plots for supplement
# Supplement: Figure S1

# focus on native density increase/decrease
pdf(file=here("plots/PCAforSupplement.pdf"),
    width=14,height=12)
par(mfrow=c(2,2))
par(mar=c(2,3,1,1))
par(xpd=F)
plotOnAxesDif(comm9,comm9,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Default Scenario"),
              ylabopt = expression("PC2 Flatter Scenario"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150
pcaEllipseDif(comm9,comm9,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=1)
pcaEllipseDif(comm10,comm9,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=2)
pcaEllipseDif(comm11,comm9,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=3)
pcaEllipseDif(comm12,comm9,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=4)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))
addArrows(comm9[2][[1]]$rotation,mult=3.5,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)

# custom legend
par(xpd=T)

par(xpd=F)
plotOnAxesDif(comm10,comm10,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Lower m Scenario"),
              ylabopt = expression("PC2 Lower m Scenario"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150
pcaEllipseDif(comm9,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=1)
pcaEllipseDif(comm10,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=2)
pcaEllipseDif(comm11,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=3)
pcaEllipseDif(comm12,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=4)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))
addArrows(comm10[2][[1]]$rotation,mult=3.5,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)

# custom legend
par(xpd=T)

par(xpd=F)
plotOnAxesDif(comm11,comm11,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Flatter Scenario, increased "*tau*", decreased b"),
              ylabopt = expression("PC2 Flatter Scenario, increased "*tau*", decreased b"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150
pcaEllipseDif(comm9,comm11,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=1)
pcaEllipseDif(comm10,comm11,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=2)
pcaEllipseDif(comm11,comm11,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=3)
pcaEllipseDif(comm12,comm11,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=4)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))
addArrows(comm11[2][[1]]$rotation,mult=3.5,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)

# custom legend
par(xpd=T)

par(xpd=F)
plotOnAxesDif(comm12,comm12,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Higher "*sigma*" Scenario"),
              ylabopt = expression("PC2 Higher "*tau*" Scenario"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150
pcaEllipseDif(comm9,comm12,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=1)
pcaEllipseDif(comm10,comm12,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=2)
pcaEllipseDif(comm11,comm12,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=3)
pcaEllipseDif(comm12,comm12,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=4)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))
addArrows(comm12[2][[1]]$rotation,mult=3.5,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)

# custom legend
par(xpd=T)
dev.off()

#################################################

################################################################
# Drawing parameters from grid on PCA axes in a structured way and then solving ODE 
# note, this takes a long time (approximately overnight for me, depending on your machine)
# output from this is saved in the output folder

# I split this into 16 pieces with nsplit, and ran different pieces in different R instances
# Though you could probably parallelize in a more efficient way

# this shows solving for section 13/16, need to manually change for other sections
set.seed(1496)
startTime=Sys.time()
result1=pcaGridSearch(list(comm9,comm10,comm11,comm12),comm10,gridLength = 16,ntrials=150,trim=F, nsplit=16, splitsec=13)
saveRDS(result1,here(paste0("output/scen_all_",13,".rds")))
Sys.time()-startTime

################################################################################

################################################################################
# compile output and prep for next figure

par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1) )
par(mgp=c(3,1,0))

# compile data
split1=readRDS(here("output/scen_all_1.rds"))
split2=readRDS(here("output/scen_all_2.rds"))
split3=readRDS(here("output/scen_all_3.rds"))
split4=readRDS(here("output/scen_all_4.rds"))
split5=readRDS(here("output/scen_all_5.rds"))
split6=readRDS(here("output/scen_all_6.rds"))
split7=readRDS(here("output/scen_all_7.rds"))
split8=readRDS(here("output/scen_all_8.rds"))
split9=readRDS(here("output/scen_all_9.rds"))
split10=readRDS(here("output/scen_all_10.rds"))
split11=readRDS(here("output/scen_all_11.rds"))
split12=readRDS(here("output/scen_all_12.rds"))
split13=readRDS(here("output/scen_all_13.rds"))
split14=readRDS(here("output/scen_all_14.rds"))
split15=readRDS(here("output/scen_all_15.rds"))
split16=readRDS(here("output/scen_all_16.rds"))

gridLength=16


commMatPointsList=list(comm9,comm10,comm11,comm12)
commMatAxes=comm11
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

full=Map(c,split1,split2,split3,split4,split5,split6,split7,split8,split9,split10,split11,split12,split13,split14,split15,split16)
modellistFig2=full[[1]]
modellistFig2NoEvo=full[[2]]

# these will be our chosen points for further exploration (purple points)
nfig4ALL=c(36,41,46,116,121,126,196,201,206)
nfig4=nfig4ALL[1:9] # assign numbers here!

################################################################################

################################################################################
# Compound pie charts figure is a combination of 3 sub figures
# Main Text: Figure 2
########################################################################
########################################################################
########################################################################
########################################################################
# Figure A: Simple Pies Figure
#- 2) Native species new equilibrium w/ invader w/ evolution vs w/o evolution
#- native species equilibrium higher: evo/noevo=1.1
#- native species equilibrium lower: evo/noevo=0.91
#- native species equilibrium same: 0.91<evo/noevo<1.1
#- native species equilibrium totally extinct with and without

modellist_list=list(modellistFeb16par1,modellistFeb16par2,modellistFeb16par3,modellistFeb16par4,modellistFeb16par5,modellistFeb16par6)
modellistNoEvo_list=list(modellistNoEvoFeb16par1,modellistNoEvoFeb16par2,modellistNoEvoFeb16par3,modellistNoEvoFeb16par4,modellistNoEvoFeb16par5,modellistNoEvoFeb16par6)
pname=c(expression(bold("I: default")),
        expression(atop(bold("II: lower native optimal"), bold("trait value (lower "* m[nat]*")"))),
        expression(atop(bold("III: native selective landscape"), bold("flatter  (higher "*tau[nat]*", lower "*b[nat]*")"))),
        expression(atop(bold("IV: higher native trait"), bold("variance  (higher "*sigma[nat]*")"))),
        expression(atop(bold("V: conspecific optimum less"), bold("than heterospecific ("*theta*"<"*m*")"))),
        expression(atop(bold("VI: conspecific optimum greater"),bold("than heterospecific ("*theta*">"*m*")"))))


incLev=1.1
decLev=0.9
################################
pdf(file=here("plots/CompoundPies.pdf"),
    width=14,height=14)

layout(matrix(c(1,1,2,2,3,3,
                4,4,5,5,6,6,
                7,7,7,8,8,8),
              nrow = 3, byrow = TRUE),heights=c(1,1,2))


par(mar=c(0.5, 4.1, 4.1, 2.1))
#par(mfrow=c(2,3))


for(j in 1:6){
  
  modList=modellist_list[[j]]
  modListNoEvo=modellistNoEvo_list[[j]]
  
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
    } else if(modList[[i]][26][[1]]==0&modListNoEvo[[i]][26][[1]]==0) { # Native abundance is same
      listerC=c(listerC,i)
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
  
  
  
  piedat=c(length(listerA),length(listerB),length(listerC))/(length(listerA)+length(listerB)+length(listerC))
  
  pie(x=piedat,col=c("seagreen","thistle1","skyblue"),labels=c("increase","decrease","same"),main=pname[j],cex=2,cex.main=2)
  
  
}

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(xpd=NA)
text(-8.55,3.85,"A",cex=3)


################################################################################
# Figure B: arrows and ellipses


par(xpd=F)
plotOnAxesDif(comm10,comm10,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Flatter Scenario, increased "*tau*", decreased b"),
              ylabopt = expression("PC2 Flatter Scenario, increased "*tau*", decreased b"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150
pcaEllipseDif(comm9,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=1)
pcaEllipseDif(comm10,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=2)
pcaEllipseDif(comm11,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=3)
pcaEllipseDif(comm12,comm10,col=c("seagreen","thistle1"),verbose=F,lwd=3,lty=4)

fancylabels=c(expression(m[inv]),expression(b[inv]),expression(A),
              expression(H),expression("c'"),expression(delta),
              expression(tau[inv]),expression(sigma[inv]))
addArrows(comm10[2][[1]]$rotation,mult=3,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)

# custom legend
par(xpd=T)
# points(-4,5.3,pch=16,cex=2,col="seagreen")
# text(-3.5,5.3,"native increased",adj=0)
# points(-4,4.95,pch=16,cex=2,col="skyblue")
# text(-3.5,4.95,"no change",adj=0)
# points(-4,4.6,pch=16,cex=2,col="thistle1")
# text(-3.5,4.6,"native decreased",adj=0)

segments(-3.8,4.0,-3.3,4.0,lwd=3,lty=1,col="thistle1")
text(-3.2,4.0,"native decreased",adj=0,cex=2)
segments(-3.8,3.7,-3.3,3.7,lwd=3,lty=1,col="seagreen")
text(-3.2,3.7,"native increased",adj=0,cex=2)

segments(-0.3,4.0,-0.8,4.0,lwd=3,lty=1,col="black")
text(-0.2,4.0,"default",adj=0,cex=2)
segments(-0.3,3.7,-0.8,3.7,lwd=3,lty=2,col="black")
text(-0.2,3.8,"lower m",adj=0,cex=2)
segments(1.2,4.0,1.7,4.0,lwd=3,lty=3,col="black")
text(1.8,4.0,expression("flatter (inc "*tau*", dec b)"),adj=0,cex=2)
segments(1.2,3.7,1.7,3.7,lwd=3,lty=4,col="black")
text(1.8,3.7,expression("higher "*sigma),adj=0,cex=2)

text(-3.4,3,"B",cex=3)


################################################################################
# Figure C: arrows and pie charts

par(xpd=F)
plotOnAxesDif(comm10,comm10,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = expression("PC1 Flatter Scenario, increased "*tau*", decreased b"),
              ylabopt = expression("PC2 Flatter Scenario, increased "*tau*", decreased b"),cex=2)
#plotOnAxesDif(comm11,comm11)
ntrials=150

for(i in 1:gridLength){
  for(j in 1:gridLength){
    n=(i-1)*gridLength+j
    if(length(modellistFig2[[n]])==1){
      points(xseq[j],yseq[i],cex=2,pch=16,col="black")
    } else {
      natPop100=extractNatPop(modellistFig2[[n]],ntrials)
      #print(natPop100)
      natPop100[which(natPop100<0)]=NA
      natPopNoEvo100=extractNatPop(modellistFig2NoEvo[[n]],ntrials)
      #print(natPopNoEvo100)
      natPopNoEvo100[which(natPopNoEvo100<0)]=NA
      
      natPopChange=mapply(function(natPop100,natPopNoEvo100){
        if(is.na(natPopNoEvo100)|is.na(natPop100)){
          return(NA)
        } else if(natPopNoEvo100==0){
          if(natPop100>0){
            return("increase")
          }
          if(natPop100==0){
            return("same")
          }
        } else {
          if(natPop100/natPopNoEvo100>1.1){
            return("increase")
          } else if(natPop100/natPopNoEvo100<0.9){
            return("decrease")
          } else {
            return("same")
          }
        }
      },natPop100,natPopNoEvo100)
      
      #print(natPopChange)
      nat_inc=length(which(natPopChange=="increase"))
      nat_dec=length(which(natPopChange=="decrease"))
      nat_sam=length(which(natPopChange=="same"))
      nat_na=length(which(is.na(natPopChange)))
      
      nat_cat=c(nat_inc,nat_sam,nat_dec,nat_na)
      
      # add pie charts
      if(!(i%%gridLength)==1){
        floating.pie(
          xseq[j], 
          yseq[i], 
          nat_cat / sum(nat_cat), # Normalize data for pie chart
          radius = 0.2,
          col=c("seagreen","skyblue","thistle1","black"),
          border=NA
        ) 
      }
      
      
      
      
      
    }
    
  }
}

for(i in 1:gridLength){
  for(j in 1:gridLength){
    n=(i-1)*gridLength+j
    if(n%in%nfig4){
      points(xseq[j],yseq[i],cex=6.5,pch=1,col="darkviolet")
      points(xseq[j],yseq[i],cex=6.7,pch=1,col="darkviolet")
      points(xseq[j],yseq[i],cex=6.1,pch=1,col="darkviolet")
      points(xseq[j],yseq[i],cex=6.3,pch=1,col="darkviolet")
    } else {
      #points(xseq[j],yseq[i],cex=2,pch=1)
    }
    
  }
}
addArrows(comm10[2][[1]]$rotation,mult=3,lwd=3,tsp=1.3,fancy=T,fancylabels=fancylabels,cex=2)
par(xpd=T)
points(-3.3,3.7,pch=16,cex=3,col="seagreen")
text(-3.1,3.7,"native increased",adj=0,cex=2)
points(-0.6,3.7,pch=16,cex=3,col="skyblue")
text(-0.4,3.7,"no change",adj=0,cex=2)
points(1.5,3.7,pch=16,cex=3,col="thistle1")
text(1.7,3.7,"native decreased",adj=0,cex=2)

text(-3.4,3,"C",cex=3)
dev.off()

################################################################################

################################################################################
# heat map of native species both with and without evolution being zero
# Supplement: Figure S2

pdf(file=here("plots/FigPCA_nativezero.pdf"),
    width=7,height=7.5)

plotOnAxesDif(comm11,comm11,type="n",pbuff=0.0,shrinkaxes = 1.1,
              xlabopt = "PC1 Flatter Scenario, increased tau, decreased b",
              ylabopt = "PC2 Flatter Scenario, increased tau, decreased b")
ntrials=150

for(i in 1:gridLength){
  for(j in 1:gridLength){
    n=(i-1)*gridLength+j
    if(length(modellistFig2[[n]])==1){
      points(xseq[j],yseq[i],cex=2,pch=16,col="black")
    } else {
      natPop100=extractNatPop(modellistFig2[[n]],ntrials)
      natPop100[which(natPop100<0)]=NA
      natPopNoEvo100=extractNatPop(modellistFig2NoEvo[[n]],ntrials)
      natPopNoEvo100[which(natPopNoEvo100<0)]=NA
      
      zfreq=length(which(natPop100==0&natPopNoEvo100==0))/ntrials
      print(zfreq)
      # exclude infs and nans
   

        if(zfreq>0.9){
          color1="gray5"
        } else if(zfreq>0.8){
          color1="gray15"
        } else if(zfreq>0.7){
          color1="gray25"
        } else if(zfreq>0.6){
          color1="gray35"
        } else if(zfreq>0.5){
          color1="gray45"
        } else if(zfreq>0.4){
          color1="gray55"
        } else if(zfreq>0.3){
          color1="gray65"
        } else if(zfreq>0.2){
          color1="gray75"
        } else if(zfreq>0.1){
          color1="gray85"
        } else {
          color1="gray95"
        }
        points(xseq[j],yseq[i],cex=2,pch=16,col=color1)
      
      
    }
    
  }
}
dev.off()

################################################################################
# This code runs the structured analysis for chosen points from structured PCA search
# it will run slowly (overnight to a couple days, depending on your machine)
# output from this is saved in the output folder

# needs: draw from several points, across the grid search- utilizing already gathered info!
# m/n, b/u, sigma, tau

# choose 9 points of interest
# choose 3 parameter sets per point (5% percentile, median, 95% percentile)
# choose 4 different parameter plots per parameter set
# 9X9 grids with original point in center BUT need to stop negative values of parameters
# go to zero but if not enough space, extend the parameter range further in the positive direction instead

# identifying purple points, need to get their index number from the big list
# they go from 1 to 265, 

# recalculate the grid
commMatPointsList=c(list(comm9,comm10,comm11,comm12))
commMatAxes=comm10
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

# these will be our chosen points
nfig4ALL=c(36,41,46,116,121,126,196,201,206)
nfig4=nfig4ALL[1:9] # assign numbers here!

for(i in 1:gridLength){
  for(j in 1:gridLength){
    n=(i-1)*gridLength+j
      if(n%in%nfig4){
        points(xseq[j],yseq[i],cex=2.7,pch=1,col="darkviolet")
        points(xseq[j],yseq[i],cex=2.5,pch=1,col="darkviolet")
        points(xseq[j],yseq[i],cex=2.3,pch=1,col="darkviolet")
      } else {
        #points(xseq[j],yseq[i],cex=2,pch=1)
      }
    
  }
}


ntrials=150
# iterate through chosen points
for(i in 1:length(nfig4)){
  # for each point of interest, choose the 5%ile, 50%ile, 95%ile
  # first extract the simulations for the chosen point
  natPop100=extractNatPop(modellistFig2[[nfig4[i]]],ntrials)
  natPop100[which(natPop100<0)]=NA
  natPopNoEvo100=extractNatPop(modellistFig2NoEvo[[nfig4[i]]],ntrials)
  natPopNoEvo100[which(natPopNoEvo100<0)]=NA

  
  natPopChange1=mapply(function(natPop100,natPopNoEvo100){
    if(is.na(natPopNoEvo100)|is.na(natPop100)){
      return(NA)
    } else if(natPopNoEvo100==0){
      if(natPop100>0){
        return(9999) # just give a high number instead of inf, note that taking the mean will be pointless!
      }
      if(natPop100==0){
        return(1) # zero to zero, change should be considered 1 (no change)
      }
    } else {
      return(natPop100/natPopNoEvo100)
    }
  },natPop100,natPopNoEvo100)
  
  print(natPopChange1)
  
  # exclude NAs
  if(any(is.na(natPopChange1))){
    natPopChange=natPopChange1[-which(is.na(natPopChange1))]
  } else {
    natPopChange=natPopChange1
  }
  
  # we have them now, get the percentiles (if percentiles have identical values, choose randomly)
  percentileInd=numeric(3)
  percentileVal=numeric(3)
  
  percentileInd[1]=round(0.05*length(natPopChange))
  percentileVal[1]=sort(natPopChange)[percentileInd[1]]
  percentileInd[2]=round(0.50*length(natPopChange))
  percentileVal[2]=sort(natPopChange)[percentileInd[2]]
  percentileInd[3]=round(0.95*length(natPopChange))
  percentileVal[3]=sort(natPopChange)[percentileInd[3]]
  
  
  for(j in 1:length(percentileVal)){
    
    # with the percentiles, get back the parameters
    if(length(which(natPopChange==percentileVal[j]))>1){
      parCurrent=extractPars(modellistFig2[[nfig4[i]]][sample(which(natPopChange1==percentileVal[j]),1)][[1]])
    } else {
      parCurrent=extractPars(modellistFig2[[nfig4[i]]][which(natPopChange1==percentileVal[j])][[1]])
    }
    
    # set optima for single species to be same as m, n
    parCurrent[which(names(parCurrent)=="thetaC")]<-parCurrent[which(names(parCurrent)=="m")]
    parCurrent[which(names(parCurrent)=="thetaH")]<-parCurrent[which(names(parCurrent)=="n")]
    
    # do a simulation for each grid cell
    # use sweep solver to solve on a gridGMspecialNoEvo
    ranger=1.5
    # b and u
    par1="b"
    par2="u"
    resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial)
    resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo)
    
    saveRDS(resFig,here(paste0('output/point_b',nfig4[i],'_',j,'.rds')))
    saveRDS(resFigNoEvo,here(paste0('output/point_b',nfig4[i],'_',j,'NoEvo.rds')))
    
    # m and n
    par1="m"
    par2="n"
    resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial)
    resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo)
    
    saveRDS(resFig,here(paste0('output/point_m',nfig4[i],'_',j,'.rds')))
    saveRDS(resFigNoEvo,here(paste0('output/point_m',nfig4[i],'_',j,'NoEvo.rds')))
    
    
    # tau1 and tau2
    par1="tauInter1"
    par2="tauInter2"
    resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial)
    resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo)
    
    saveRDS(resFig,here(paste0('output/point_tau',nfig4[i],'_',j,'.rds')))
    saveRDS(resFigNoEvo,here(paste0('output/point_tau',nfig4[i],'_',j,'NoEvo.rds')))
    
    
    # sigmax and sigmay
    par1="sigmax"
    par2="sigmay"
    resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial)
    resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo)
    
    saveRDS(resFig,here(paste0('output/point_sigma',nfig4[i],'_',j,'.rds')))
    saveRDS(resFigNoEvo,here(paste0('output/point_sigma',nfig4[i],'_',j,'NoEvo.rds')))
    
    
  }

}


################################################################################

########################################################
# Produces all 9 parameter space plots
# Supplement: S5

par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
for(i in 1:9){
  for(j in 1:3){
    # b and u
    resb=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'.rds')))
    resbNoEvo=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(resb,resbNoEvo,"b")
    
    # m and n
    resm=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'.rds')))
    resmNoEvo=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(resm,resmNoEvo,"m")
    
    # sigma
    ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'.rds')))
    ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(ress,ressNoEvo,"sigma")
    
    # tau
    rest=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'.rds')))
    restNoEvo=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(rest,restNoEvo,"tau")
    
  }
}
################################################################################################

########################################################
# Detailed structured parameter evaluation
# Main Text: Figure 4

pdf(file="plots/ParameterSpace.pdf",
    width=8.5,height=6)

par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
par(oma=c(1,3,1,1))
for(i in 4){
  for(j in 1:3){
    # sigma
    ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'.rds')))
    ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(ress,ressNoEvo,"sigma")
    
    # b and u
    resb=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'.rds')))
    resbNoEvo=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(resb,resbNoEvo,"b")
    
    
    # m and n
    resm=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'.rds')))
    resmNoEvo=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(resm,resmNoEvo,"m")
    
    
    # tau
    rest=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'.rds')))
    restNoEvo=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'NoEvo.rds')))
    NatChangeGrid(rest,restNoEvo,"tau")
    
  }
}


rowLabels <- c("95th Percentile", "50th Percentile", "5th Percentile")
for(l in 1:3) {
  mtext(rowLabels[l], side = 2, line = 0, outer = TRUE,
        at = 1 - (l - 0.5)/3, las = 0)
}

dev.off()

################################################################################

################################################################################################
# The following code runs more simulations to investigate different initial conditions
# expect longer run times (overnight to a few days depending on your machine)
# output is found in the output folder

par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
for(i in 9){
  for(j in 1:3){
    # b and u
    resb=readRDS(paste0('output/point_b',nfig4ALL[i],'_',j,'.rds'))
    resbNoEvo=readRDS(paste0('output/point_b',nfig4ALL[i],'_',j,'NoEvo.rds'))
    NatChangeGrid(resb,resbNoEvo,"b")
    
    # m and n
    resm=readRDS(paste0('output/point_m',nfig4ALL[i],'_',j,'.rds'))
    resmNoEvo=readRDS(paste0('output/point_m',nfig4ALL[i],'_',j,'NoEvo.rds'))
    NatChangeGrid(resm,resmNoEvo,"m")
    
    # sigma
    ress=readRDS(paste0('output/point_sigma',nfig4ALL[i],'_',j,'.rds'))
    ressNoEvo=readRDS(paste0('output/point_sigma',nfig4ALL[i],'_',j,'NoEvo.rds'))
    NatChangeGrid(ress,ressNoEvo,"sigma")
    
    # tau
    rest=readRDS(paste0('output/point_tau',nfig4ALL[i],'_',j,'.rds'))
    restNoEvo=readRDS(paste0('output/point_tau',nfig4ALL[i],'_',j,'NoEvo.rds'))
    NatChangeGrid(rest,restNoEvo,"tau")
    
  }
}

# Great now grab the code to produce those results

# get our parameters

ntrials=150
nfig4=nfig4ALL[5]
# iterate through chosen points

# for each point of interest, choose the 5%ile, 50%ile, 95%ile
# first extract the simulations for the chosen point
natPop100=extractNatPop(modellistFig2[[nfig4[i]]],ntrials)
natPop100[which(natPop100<0)]=NA
natPopNoEvo100=extractNatPop(modellistFig2NoEvo[[nfig4[i]]],ntrials)
natPopNoEvo100[which(natPopNoEvo100<0)]=NA



# new natPopChange!
natPopChange1=mapply(function(natPop100,natPopNoEvo100){
  if(is.na(natPopNoEvo100)|is.na(natPop100)){
    return(NA)
  } else if(natPopNoEvo100==0){
    if(natPop100>0){
      return(9999) # just give a high number instead of inf, note that taking the mean will now be pointless
    }
    if(natPop100==0){
      return(1) # zero to zero, change should be considered 1
    }
  } else {
    return(natPop100/natPopNoEvo100)
  }
},natPop100,natPopNoEvo100)

print(natPopChange1)

# exclude NAs
if(any(is.na(natPopChange1))){
  natPopChange=natPopChange1[-which(is.na(natPopChange1))]
} else {
  natPopChange=natPopChange1
}

# we have them now, get the percentiles (if percentiles have identical values, choose randomly)
percentileInd=numeric(3)
percentileVal=numeric(3)

percentileInd[1]=round(0.05*length(natPopChange))
percentileVal[1]=sort(natPopChange)[percentileInd[1]]
percentileInd[2]=round(0.50*length(natPopChange))
percentileVal[2]=sort(natPopChange)[percentileInd[2]]
percentileInd[3]=round(0.95*length(natPopChange))
percentileVal[3]=sort(natPopChange)[percentileInd[3]]

invasiveTraitVector=c(1.5,5,10)
invasivePopVector=c(0.005,0.05,0.3)

for(i in 1:3){ # trait loop
  for(k in 1:3){ # population loop
    
    for(j in 1:length(percentileVal)){
      
      # with the percentiles, get back the parameters
      if(length(which(natPopChange==percentileVal[j]))>1){
        parCurrent=extractPars(modellistFig2[[nfig4]][sample(which(natPopChange1==percentileVal[j]),1)][[1]])
      } else {
        parCurrent=extractPars(modellistFig2[[nfig4]][which(natPopChange1==percentileVal[j])][[1]])
      }
      
      # set optima for single species to be same as m, n
      parCurrent[which(names(parCurrent)=="thetaC")]<-parCurrent[which(names(parCurrent)=="m")]
      parCurrent[which(names(parCurrent)=="thetaH")]<-parCurrent[which(names(parCurrent)=="n")]
      
      # do a simulation for each grid cell
      # use sweep solver to solve on a gridGMspecialNoEvo
      ranger=1.5
      # b and u
      par1="b"
      par2="u"
      resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      
      saveRDS(resFig,here(paste0('output/point_b',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'.rds')))
      saveRDS(resFigNoEvo,here(paste0('output/point_b',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'NoEvo.rds')))
      
      # m and n
      par1="m"
      par2="n"
      resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      
      saveRDS(resFig,here(paste0('output/point_m',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'.rds')))
      saveRDS(resFigNoEvo,here(paste0('output/point_m',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'NoEvo.rds')))
      
      
      # tau1 and tau2
      par1="tauInter1"
      par2="tauInter2"
      resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      
      saveRDS(resFig,here(paste0('output/point_tau',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'.rds')))
      saveRDS(resFigNoEvo,here(paste0('output/point_tau',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'NoEvo.rds')))
      
      
      # sigmax and sigmay
      par1="sigmax"
      par2="sigmay"
      resFig=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecial,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      resFigNoEvo=sweepSolverModGM3(parCurrent,param1=par1,param2=par2,start1=(parCurrent[which(names(parCurrent)==par1)]-ranger),end1=(parCurrent[which(names(parCurrent)==par1)]+ranger),by1=0.25,start2=(parCurrent[which(names(parCurrent)==par2)]-ranger),end2=(parCurrent[which(names(parCurrent)==par2)]+ranger),by2=0.25,func=GMspecialNoEvo,invT0=invasiveTraitVector[i], invP0=invasivePopVector[k])
      
      saveRDS(resFig,here(paste0('output/point_sigma',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'.rds')))
      saveRDS(resFigNoEvo,here(paste0('output/point_sigma',nfig4,'_',j,'trait',invasiveTraitVector[i],'pop',invasivePopVector[k],'NoEvo.rds')))
      
      
    }
    
  }
}

################################################################################################

################################################################################################
# The following code produces all the initial conditions plots
# Supplement: Figure S8

par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
for(k in 1:length(invasivePopVector)){
  for(l in 1:length(invasiveTraitVector)){
    for(i in 5){
      for(j in 1:3){
        # b and u
        resb=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'.rds')))
        resbNoEvo=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'NoEvo.rds')))
        NatChangeGrid(resb,resbNoEvo,"b")
        
        # m and n
        resm=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'.rds')))
        resmNoEvo=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'NoEvo.rds')))
        NatChangeGrid(resm,resmNoEvo,"m")
        
        # sigma
        ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'.rds')))
        ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'NoEvo.rds')))
        NatChangeGrid(ress,ressNoEvo,"sigma")
        
        # tau
        rest=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'.rds')))
        restNoEvo=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[l],'pop',invasivePopVector[k],'NoEvo.rds')))
        NatChangeGrid(rest,restNoEvo,"tau")
        
      }
    }
  }
}

################################################################################################

################################################################################################
# Main text initial condition plots
# Figure 5

invasiveTraitVector=c(1.5,5,10)
invasivePopVector=c(0.005,0.05,0.3)

pdf(file=here("plots/InitialConditions.pdf"),
    width=8,height=8)

par(mar=c(2,2,1,1),
    mgp=c(0.5,0.5,0),
    mfrow=c(3,3),
    oma = c(1, 4, 4, 1))

# set second row: 50th percentile
j=2
# set center point out of the 9 explored
i=5

for(k in 1:length(invasiveTraitVector)){ # the 3 initial invasive populations
  for(l in 1:length(invasivePopVector)){ # the 3 initial invasive traits
    
    # sigma
    ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[k],'pop',invasivePopVector[l],'.rds')))
    ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'trait',invasiveTraitVector[k],'pop',invasivePopVector[l],'NoEvo.rds')))
    NatChangeGrid(ress,ressNoEvo,"sigma")
    
  }
}

columnLabels <- c("Invasive population\n 0.5% of native", "Invasive population\n 5% of native", "Invasive population\n 30% of native")
for(k in 1:3) {
  mtext(columnLabels[k], side = 3, line = 0, outer = TRUE,
        at = (k - 0.5)/3)
}


rowLabels <- c("Invasive Trait\n 1.5 times higher", "Invasive Trait\n 5 times higher", "Invasive Trait\n 10 times higher")
for(l in 1:3) {
  mtext(rowLabels[l], side = 2, line = 0, outer = TRUE,
        at = 1 - (l - 0.5)/3, las = 0)
}

dev.off()


################################################################################################

################################################################################################
# Plots of equilibrium Alpha values (competitive ability)
# Supplement: Figures S3 and S4
# Supplement: Figures S6 and S7 are compilations of these figures, with i=1:9

pdf(file=here("output/NativeAlphasSupplement.pdf"),
    width=8.5,height=6)
# native values
par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
for(i in 4){ # select the center image here
  for(j in 1:3){ # the 3 percentiles
    # we have 4 variables (b, m, sigma, tau)
    print(j)
    
    # b and u
    resb=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'.rds')))
    resbNoEvo=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(resb,resbNoEvo,"b",i,j,showsp = "native")
    
    # m and n
    resm=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'.rds')))
    resmNoEvo=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(resm,resmNoEvo,"m",i,j,showsp = "native")
    
    # sigma
    ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'.rds')))
    ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(ress,ressNoEvo,"sigma",i,j,showsp = "native")
    
    # tau
    rest=readRDS(paste0('output/point_tau',nfig4ALL[i],'_',j,'.rds'))
    restNoEvo=readRDS(paste0('output/point_tau',nfig4ALL[i],'_',j,'NoEvo.rds'))
    
    NatCG_alphas(rest,restNoEvo,"tau",i,j,showsp = "native")
    
  }
}

dev.off()


pdf(file=here("plots/InvasiveAlphasSupplement.pdf"),
    width=8.5,height=6)
# invasive values
par(mar=c(2,2,1,1))
par(mgp=c(0.5,0.5,0))
par(mfrow=c(3,4))
for(i in 4){ # select the center image here
  for(j in 1:3){ # the 3 percentiles
    # we have 4 variables (b, m, sigma, tau)
    print(j)
    
    # b and u
    resb=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'.rds')))
    resbNoEvo=readRDS(here(paste0('output/point_b',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(resb,resbNoEvo,"b",i,j,showsp = "invasive")
    
    # m and n
    resm=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'.rds')))
    resmNoEvo=readRDS(here(paste0('output/point_m',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(resm,resmNoEvo,"m",i,j,showsp = "invasive")
    
    # sigma
    ress=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'.rds')))
    ressNoEvo=readRDS(here(paste0('output/point_sigma',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(ress,ressNoEvo,"sigma",i,j,showsp = "invasive")
    
    # tau
    rest=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'.rds')))
    restNoEvo=readRDS(here(paste0('output/point_tau',nfig4ALL[i],'_',j,'NoEvo.rds')))
    
    NatCG_alphas(rest,restNoEvo,"tau",i,j,showsp = "invasive")
    
  }
}

dev.off()