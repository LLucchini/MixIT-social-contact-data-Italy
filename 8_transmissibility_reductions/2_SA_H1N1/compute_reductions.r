# Clear the workspace
rm(list = ls())

# Load required plotting functions
source("../../plotting_functions.R")

CLOSE.RANDOM <- 0
R0 <- 1.44 ## Merler et al, PlOS One
load("input/matr_in_presence.Rdata")
load("input/matr_no_presence.Rdata")

nages <- 15
NBOOT <- 1000

susceptibility <- read.csv("input/susc_age_H1N1.csv")
values <- susceptibility$susc

SUSC.MEAN <- c(rep(values[1],1),rep(values[2],2),rep(values[3],10),rep(values[4],2))

load("input/matr_in_presence.Rdata")
load("input/matr_no_presence.Rdata")
##==============================##
## COMBINED SCENARIOS
##==============================##

scen.school <- c("no_closure","uni","uni+upper_secondary",
                 "uni+secondary","uni+secondary+primary","all_closed")

scen.smart <- c("pre-lockdown","lifting lockdown", "lockdown")

nscen.school <- length(scen.school)
pop.scen.combined <- array(NA,c(length(scen.school),length(scen.smart),nages,2))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    pop.scen.combined[s,t,,] <- as.matrix(read.csv(paste("../0_population_scenarios/pop_scenarios/scenario_smart_",t,"_sclose_",scen.school[s],".csv",sep="")))
  }
}


## ============================== ##
## NEXT GENERATION MATRIX (NGM)
## ============================== ##
# Arrays storing block NGMs
block.inpresence <- block.nopresence <- array(
  NA, c(NBOOT, nscen.school, length(scen.smart), nages, nages)
)

# Store full NGMs in a list
NGM <- vector("list", NBOOT)

# Compute NGM blocks
for(n in 1:NBOOT){
  NGM[[n]] <- list()
  for(s in 1:length(scen.school)){
    NGM[[n]][[s]] <- list()
    for(t in 1:length(scen.smart)){
        for(i in 1:nages){
          for(j in 1:nages){
            if(CLOSE.RANDOM==0){
              block.inpresence[n,s,t,i,j] <-   SUSC.MEAN[i]*matr.inpresence[n,7,i,j]*pop.scen.combined[s,t,i,1]/sum(pop.scen.combined[s,t,j,])
              block.nopresence[n,s,t,i,j] <-   SUSC.MEAN[i]*matr.nopresence[n,7,i,j]*pop.scen.combined[s,t,i,2]/sum(pop.scen.combined[s,t,j,])
            }else{
              tmp1 <- matr.inpresence[n,1,i,j]+ matr.inpresence[n,2,i,j]+ matr.inpresence[n,3,i,j] ##solo casa scuola e lavoro
            
              block.inpresence[n,s,t,i,j] <-   SUSC.MEAN[i]*tmp1*pop.scen.combined[s,t,i,1]/sum(pop.scen.combined[s,t,j,])
              block.nopresence[n,s,t,i,j] <-   SUSC.MEAN[i]*matr.nopresence[n,1,i,j]*pop.scen.combined[s,t,i,2]/sum(pop.scen.combined[s,t,j,])
            }
          } 
        }
    
      ROW.inpresence <- cbind(block.inpresence[n,s,t,,],block.inpresence[n,s,t,,])
      ROW.nopresence <- cbind(block.nopresence[n,s,t,,],block.nopresence[n,s,t,,])
      NGM[[n]][[s]][[t]] <- rbind(ROW.inpresence,ROW.nopresence)
    }
  }
}
## ============================== ##
## DOMINANT EIGENVALUES (Rho)
## ============================== ##
rho <- array(NA,c(NBOOT,length(scen.school),length(scen.smart)))
for(n in 1:NBOOT){
  for(s in 1:length(scen.school)){
    for(t in 1:length(scen.smart)){
      rho[n,s,t] <-max(Mod(eigen(NGM[[n]][[s]][[t]], only.values = TRUE)$values))
    }
  }
}

## ============================== ##
## PERCENTAGE REDUCTION VS BASE CASE
## ============================== ##
reductions.combined <- array(NA,c(NBOOT,length(scen.school),length(scen.smart)))
for(n in 1:NBOOT){
  for(s in 1:length(scen.school)){
    for(t in 1:length(scen.smart)){
      reductions.combined[n,s,t] <- (1-rho[n,s,t]/rho[n,1,1])*100
    }
  }
}

# Store mean and 95% CI of reductions

stat.reductions.boot <- array(NA,c(length(scen.school),length(scen.smart),3))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    stat.reductions.boot[s,t,1] <- mean(reductions.combined[,s,t])
    stat.reductions.boot[s,t,2] <- quantile(reductions.combined[,s,t],0.025)
    stat.reductions.boot[s,t,3] <- quantile(reductions.combined[,s,t],0.975)
    
  }
}

## ============================== ##
## PLOT VISUALIZATION
## ============================== ##

# Load the libraries
library(sysfonts)
library(showtext)

# Add the Barlow font
font_add_google("Barlow Semi Condensed")
showtext_auto()
par(family = "Barlow Semi Condensed")


verde <- "#1a9850"
cols <- c(lighten(verde, 0.6),lighten(verde,0.3),verde)



lab.scen.smart <- c("baseline (pre-pandemic)",
                    "sustainable (no closure of economic activities)",
                    "minimum (only essential activities)")


lab.scen.school <- c("no level\nclosed",
                     "ISCED 5-8",
                     "ISCED 5-8\nISCED 3",
                     "ISCED 5-8\nISCED 2-3",
                     "ISCED 5-8\nISCED 1-3",
                     "all levels\nclosed")
makefigs=T
if(makefigs==T){
  pdf(paste("Figure_SI14A.pdf",sep=""),
      width=16,height=8)
  par(mfrow=c(1,1),mar=c(8,7,1,1),mgp=c(6,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 5
  BY=1
  xcoord.labx<- seq(1,nscen.school,by=STEP)
  MIN <- 0
  MAX <- 65
  plot(c(0.5,nscen.school+0.5),c(MIN,MAX), type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=2)
  min.RW <- stat.reductions.boot[,1,]
  sus.RW <- stat.reductions.boot[,2,]
  max.RW <- stat.reductions.boot[,3,]
  
  for(s in 1:nscen.school){
    if(min.RW[s,1]>0){
      makerect(s-0.2,min.RW[s,1],0.2,col=cols[1])
      segments(s-0.2,min.RW[s,2],s-0.2,min.RW[s,3],col=darken(cols[1],0.2),lwd=3)
    }
    makerect(s,sus.RW[s,1],0.2,col=cols[2])
    makerect(s+0.2,max.RW[s,1],0.2,col=cols[3])
    segments(s,sus.RW[s,2],s,sus.RW[s,3],col=darken(cols[2],0.2),lwd=3)
    segments(s+0.2,max.RW[s,2],s+0.2,max.RW[s,3],col=darken(cols[3],0.2),lwd=3)
    
  }
  
  text(x=-nscen.school*0.04,y=MAX/2,
       "Transmissibility reduction of H1N1pdm09 (%)",
       cex=2,
       srt=90)
  # 
  ycoord <- c(3,2.,3.2,3.2,3)*1.8
  text(x=(1:nscen.school),y=-ycoord,
       lab.scen.school,
       cex=2,
       srt=0,adj=0.5)
  
  text(x=nscen.school/1.8,y=-17,
       paste("Level of school closure"),
       cex=2,
       srt=0)
  text(x=1.1,y=65,"Level of in-person work:",cex=2)
  legend(x=0.2,y=64,lab.scen.smart,
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2)
  
 
  
  dev.off()
}

##==============================================================================##
##WITHOUT PHI - accounting for transmissibility reduction due to other infection-precautions (Marziano et al, PNAS 2021)
##==============================================================================##

PHI <- rep(0,NBOOT)
R0 <- (1-PHI)*R0

Reff <- reductions.combined*0
for(n in 1:NBOOT){
  for(s in 1:length(scen.school)){
    for(t in 1:length(scen.smart)){
      Reff[n,s,t] <- (1 - reductions.combined[n,s,t]/100)*R0[n]
    }
  }
}


stat.Reff.boot <- array(NA,c(length(scen.school),length(scen.smart),3))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    
    stat.Reff.boot[s,t,1] <- mean(Reff[,s,t])
    stat.Reff.boot[s,t,2] <- quantile(Reff[,s,t],0.025)
    stat.Reff.boot[s,t,3] <- quantile(Reff[,s,t],0.975)
  }
}

R0.baseline <- stat.Reff.boot[1,1,1]

verde <- "#88419d"

cols <- c(lighten(verde, 0.6),lighten(verde,0.3),verde)

makefigs=T
if(makefigs==T){
  pdf(paste("Figure_SI15A.pdf",sep=""),
      width=8,height=8)
  par(mfrow=c(1,1),mar=c(8,7,1,1),mgp=c(6.5,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 0.2
  BY=1
  xcoord.labx<- seq(1,nscen.school,by=STEP)
  MIN <- 0
  MAX <- 2
  plot(c(0.5,nscen.school+0.5),c(MIN,MAX), type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=2)
  min.RW <- stat.Reff.boot[,1,]
  sus.RW <- stat.Reff.boot[,2,]
  max.RW <- stat.Reff.boot[,3,]
  
  for(s in 1:nscen.school){
    if(min.RW[s,1]>0){
      makerect(s-0.2,min.RW[s,1],0.2,col=cols[1])
      segments(s-0.2,min.RW[s,2],s-0.2,min.RW[s,3],col=darken(cols[1],0.2),lwd=3)
    }
    makerect(s,sus.RW[s,1],0.2,col=cols[2])
    makerect(s+0.2,max.RW[s,1],0.2,col=cols[3])
    segments(s,sus.RW[s,2],s,sus.RW[s,3],col=darken(cols[2],0.2),lwd=3)
    segments(s+0.2,max.RW[s,2],s+0.2,max.RW[s,3],col=darken(cols[3],0.2),lwd=3)
    
  }
  segments(0.2,1,nscen.school+0.7,1,col=darken("#bababa",0.6),lty=2,lwd=2)
  segments(0.2,R0.baseline,nscen.school+0.7,R0.baseline,col="dark red",lty=2,lwd=2)
  text(x=-nscen.school*0.15,y=MAX/2,
       "Effective reproduction number of H1N1pdm09",
       cex=2,
       srt=90)
  # 
  ycoord <- c(3,2.,3.2,3.2,3)*0.05
  ycoord <- 2*0.03
  
  text(x=(1:nscen.school),y=-ycoord,
       lab.scen.school,
       cex=2,
       srt=90,adj=1)
  
  text(x=nscen.school/1.8,y=-0.59,
       paste("Level of school closure"),
       cex=2,
       srt=0)
  text(x=2.1,y=2.08,"Level of in-person work:",cex=2)
  legend(x=0.1,y=2.05,lab.scen.smart,
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2, x.intersp = 0.1)
  
  
  dev.off()
}


##==============================================================================##
##WITH PHI - accounting for transmissibility reduction due to other infection-precautions (Marziano et al, PNAS 2021)
##==============================================================================##
load("input/phi_PNAS.Rdata") ## riduzione trasmissibilita' PNAS
phi <- phi2[,1]
idsel <- sample(10000, NBOOT)
PHI <- phi[idsel]  

R0 <- (1-PHI)*R0


Reff <- reductions.combined*0
for(n in 1:NBOOT){
  for(s in 1:length(scen.school)){
    for(t in 1:length(scen.smart)){
      Reff[n,s,t] <- (1 - reductions.combined[n,s,t]/100)*R0[n]
    }
  }
}


stat.Reff.boot <- array(NA,c(length(scen.school),length(scen.smart),3))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    
    stat.Reff.boot[s,t,1] <- mean(Reff[,s,t])
    stat.Reff.boot[s,t,2] <- quantile(Reff[,s,t],0.025)
    stat.Reff.boot[s,t,3] <- quantile(Reff[,s,t],0.975)
  }
}

R0.baseline <- 1.44

verde <- "#88419d"

cols <- c(lighten(verde, 0.6),lighten(verde,0.3),verde)

makefigs=T
if(makefigs==T){
  pdf(paste("Figure_SI15B.pdf",sep=""),
      width=8,height=8)
  par(mfrow=c(1,1),mar=c(8,7,1,1),mgp=c(6.5,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 0.2
  BY=1
  xcoord.labx<- seq(1,nscen.school,by=STEP)
  MIN <- 0
  MAX <- 2
  plot(c(0.5,nscen.school+0.5),c(MIN,MAX), type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=2)
  min.RW <- stat.Reff.boot[,1,]
  sus.RW <- stat.Reff.boot[,2,]
  max.RW <- stat.Reff.boot[,3,]
  
  for(s in 1:nscen.school){
    if(min.RW[s,1]>0){
      makerect(s-0.2,min.RW[s,1],0.2,col=cols[1])
      segments(s-0.2,min.RW[s,2],s-0.2,min.RW[s,3],col=darken(cols[1],0.2),lwd=3)
    }
    makerect(s,sus.RW[s,1],0.2,col=cols[2])
    makerect(s+0.2,max.RW[s,1],0.2,col=cols[3])
    segments(s,sus.RW[s,2],s,sus.RW[s,3],col=darken(cols[2],0.2),lwd=3)
    segments(s+0.2,max.RW[s,2],s+0.2,max.RW[s,3],col=darken(cols[3],0.2),lwd=3)
    
  }
  segments(0.2,1,nscen.school+0.7,1,col=darken("#bababa",0.6),lty=2,lwd=2)
  segments(0.2,R0.baseline,nscen.school+0.7,R0.baseline,col="dark red",lty=2,lwd=2)
  text(x=-nscen.school*0.15,y=MAX/2,
       "Effective reproduction number of H1N1pdm09",
       cex=2,
       srt=90)
  # 
  ycoord <- c(3,2.,3.2,3.2,3)*0.05
  ycoord <- 2*0.03
  
  text(x=(1:nscen.school),y=-ycoord,
       lab.scen.school,
       cex=2,
       srt=90,adj=1)
  
  text(x=nscen.school/1.8,y=-0.59,
       paste("Level of school closure"),
       cex=2,
       srt=0)
  text(x=2.1,y=2.08,"Level of in-person work:",cex=2)
  legend(x=0.1,y=2.05,lab.scen.smart,
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2, x.intersp = 0.1)
  
  
  dev.off()
}

