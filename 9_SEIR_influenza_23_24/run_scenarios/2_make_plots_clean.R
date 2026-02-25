rm(list=ls())
library(data.table)
source("../../plotting_functions.R")

# Load scenario table
scenarios <- read.csv("readme_scenarios.csv")
idscen <- scenarios[,1]
nscen <- nrow(scenarios)

# MCMC input
mcmcrerun <- read.table("common_input/mcmcrerun", sep="\t")
nsim <- nrow(mcmcrerun)

# General settings
nages <- 15
nweeks <- 24
idsim <- 0:(nsim-1)

# Population by age
popbyage <- as.matrix(read.table("scenario_1/output_files/populationByAge_0", sep="\t")[,1:nweeks])
totp <- colSums(popbyage)[1]

# Load total infections
load("Rdata/totinf.Rdata")

# Reporting ratios from MCMC
reporting <- mcmcrerun$V2

# Compute cumulative infections and attack rates per scenario
cuminf <- artot <- matrix(NA, nscen, nsim)
for(s in 1:nscen){
  cuminf[s,] <- rowSums(totinf[s,,])
  artot[s,] <- cuminf[s,] / totp
}

# Scenario groups: school closure vs smart working
scen.school <- unique(scenarios[,2])
scen.smart <- unique(scenarios[,3])
nscen.school <- length(scen.school)
nscen.smart <- length(scen.smart)

lab.scen.school <- c("no level\nclosed",
                     "ISCED 5-8",
                     "ISCED 5-8\nISCED 3",
                     "ISCED 5-8\nISCED 2-3",
                     "ISCED 5-8\nISCED 1-3",
                     "all levels\nclosed")

lab.scen.smart <- c("baseline (pre-pandemic)",
                    "sustainable (no closure of economic activities)",
                    "minimum (only essential activities)")

# Compute statistics (mean, 95% CrI)
stat.artot.boot <- array(NA,c(nscen.school,nscen.smart,3))
for(s in 1:nscen.school){
  for(t in 1:nscen.smart){
    sel <- scenarios[scenarios[,2]==scen.school[s] &
                       scenarios[,3]==scen.smart[t], 1]
    stat.artot.boot[s,t,] <- c(
      mean(artot[sel,]),
      quantile(artot[sel,],0.025),
      quantile(artot[sel,],0.975)
    )
  }
}

# Load Google font
library(sysfonts)
library(showtext)
font_add_google("Barlow Semi Condensed")
showtext_auto()
par(family = "Barlow Semi Condensed")

save(stat.artot.boot, file="Rdata/stat_artot.Rdata")
blue <- "#66B2D9"
cols <- c(lighten(blue, 0.6),lighten(blue,0.3),blue)
col1 <- "#C0C0C0"

# Create figs folder if not existing
dir.create("figs", showWarnings = FALSE)
makefigs <- TRUE

### ================= FIGURE: attack rate total ================= ###

makefigs=T
if(makefigs==T){
  pdf(paste("figs/artot.pdf",sep=""),
      width=8,height=8)
  par(mfrow=c(1,1),mar=c(8,7,1,1),mgp=c(6.5,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 5
  BY=1
  xcoord.labx<- seq(1,nscen.school,by=STEP)
  MIN <- 0
  MAX <- 23
  plot(c(0.5,nscen.school+0.5),c(MIN,MAX), type="l",col=NA,xlab="",
       ylab="",axes=F)
  ##axis(1,cex.axis = 1.5,at=seq(1,nscen.school,by=BY),labels = lab.agegr, lwd=1,las=2)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=2)
  min.RW <- stat.artot.boot[,1,]*100
  sus.RW <- stat.artot.boot[,2,]*100
  max.RW <- stat.artot.boot[,3,]*100
  
  for(s in 1:nscen.school){
    if(min.RW[s,1]>0){
      if(s==1){
        rect(s-0.2-0.2/2,0,s-0.2+0.2/2,max(0,min.RW[s,1]),col=lighten(col1[1],0.5),
             border=FALSE,density=NA)
        segments(s-0.2,min.RW[s,2],s-0.2,min.RW[s,3],col=darken(col1,0.3),lwd=3)
        
      }else{
        makerect(s-0.2,min.RW[s,1],0.2,col=cols[1])
        segments(s-0.2,min.RW[s,2],s-0.2,min.RW[s,3],col=darken(cols[1],0.2),lwd=3)
        
      }
    }
    makerect(s,sus.RW[s,1],0.2,col=cols[2])
    makerect(s+0.2,max.RW[s,1],0.2,col=cols[3])
    segments(s,sus.RW[s,2],s,sus.RW[s,3],col=darken(cols[2],0.2),lwd=3)
    segments(s+0.2,max.RW[s,2],s+0.2,max.RW[s,3],col=darken(cols[3],0.2),lwd=3)
    
  }
  text(x=-nscen.school*0.13,y=MAX/2,
       "Infection attack rate of influenza A\nin the 2023-2024 season (%)",
       cex=2,
       srt=90)
  # 
  ycoord <- c(3,2.,3.2,3.2,3)*0.05
  ycoord <- 2*0.4
  
  text(x=(1:nscen.school),y=-ycoord,
       lab.scen.school,
       cex=2,
       srt=90,adj=1)
  
  text(x=nscen.school/1.8,y=-6.8,
       paste("Level of school closure"),
       cex=2,
       srt=0)
  text(x=2.1,y=23.8,"Level of in-person work:",cex=2)
  legend(x=0.1,y=23.4,lab.scen.smart,
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2, x.intersp = 0.1)
  
  
  dev.off()
}

### =============== FIT TO DATA =============== ###
inc.ILI <- read.table("scenario_1/input_files/ILI_A_2023_2024", sep="\t")

# Weekly reported cases from model
reported.cases <- totinf[1,,] * reporting

# Summary statistics
inc.fit.summary <- cbind(
  colMeans(reported.cases)/totp*1000,
  apply(reported.cases,2,quantile,0.025)/totp*1000,
  apply(reported.cases,2,quantile,0.975)/totp*1000
)

labweeks <- c(46:52,1:17)
if(makefigs==T){
  pdf("figs/fit.pdf",width=10,height=8)
  par(mfrow=c(1,1),mar=c(8,7,1,1),mgp=c(3.5,1.,0.),xpd=T,family="Barlow Semi Condensed")
  BY=1
  xcoord.labx<- seq(1,nweeks,by=1)
  MIN <- 0
  MAX <- 12
  plot(c(0.5,nweeks+0.5),c(MIN,MAX), type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(1,cex.axis = 2,at=seq(1,nweeks,by=BY),labels = labweeks,
       lwd=1,las=2)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=2)
  
  text(-3.5,MAX/2,
       "Weekly incidence of ILI due to influenza A
       in the season 2023-2024 (per 1,000 individuals)",
       cex=2,srt=90)
  text(nweeks/2,-3,
       "Week of the year",
       cex=2,srt=0)
  text(3.5,-2,
       "2023",
       cex=2,srt=0)
  text(15.5,-2,
       "2024",
       cex=2,srt=0)
  polygon(c(1:24,24:1),c(inc.fit.summary[,3],rev(inc.fit.summary[,2])),
          border=NA,col=lighten(col1,0.5))
  lines(inc.fit.summary[,1],col=darken(col1,0.3),lwd=2)
  points(inc.ILI[,2],cex=2,pch=16,col="#595959")
  
  dev.off()
}

### ========== ATTACK RATE BY AGE ========== ###

##======================================================##
## ATTACK RATE BY AGE                                   ##
##======================================================##
load(file="Rdata/cuminfage.Rdata")
dim(cuminfage)
lab.age.groups <- c("0-14 years", "15-64 years", "65 years and older")

groups <- list(
  g1 = 1:3,
  g2 = 4:13,
  g3 = 14:15
)
tmp <- popbyage[,1]
popagegroups <- c(sum(tmp[1:3]),sum(tmp[4:13]),sum(tmp[14:15]))

cuminf.agegroups <- array(0, dim = c(18, 80000, 3))

for (i in seq_along(groups)) {
  cuminf.agegroups[,,i] <- apply(cuminfage[,,groups[[i]]], c(1,2), sum)
}

nagegroups <- 3
mean.values <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  mean.values[,a] <- rowMeans(cuminf.agegroups[,,a])
}

mean.arage <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  mean.arage[,a] <- mean.values[,a]/popagegroups[a]
}

upper.ci.values <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  upper.ci.values[,a] <- apply(cuminf.agegroups[,,a],1,quantile,0.975)
}

upper.ci.arage <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  upper.ci.arage[,a] <- upper.ci.values[,a]/popagegroups[a]
}

lower.ci.values <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  lower.ci.values[,a] <- apply(cuminf.agegroups[,,a],1,quantile,0.025)
}

lower.ci.arage <- matrix(NA,nscen,nagegroups)
for(a in 1:nagegroups){
  lower.ci.arage[,a] <- lower.ci.values[,a]/popagegroups[a]
}

library(colorspace)
breaks.ar <- c(0,5,10,15,20,25,30,35)

labels <- as.character(breaks.ar)
nb <- length(breaks.ar)
breaks.ar.sub <- breaks.ar[seq(1,nb,by=1)]
label.sub <- as.character(breaks.ar.sub)

mean.ar.age.matrix <- array(NA,c(nscen.school,nscen.smart,nagegroups))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    sel <- scenarios[which(scenarios[,2]==scen.school[s] & scenarios[,3]==scen.smart[t]),1]
    mean.ar.age.matrix[s,t,1] <- mean.arage[sel,1]
    mean.ar.age.matrix[s,t,2] <- mean.arage[sel,2]
    mean.ar.age.matrix[s,t,3] <- mean.arage[sel,3]
  }
}

paletta <- diverge_hcl(nb,c=100,l=c(50,90),power=1)

lab.x <- c("baseline","sustainable","minimum")
lab.y <- lab.scen.school
nx <- nscen.smart
ny <- nscen.school

cols_mat=hcl.colors((length(breaks.ar))+3,
                    "BluGrn", rev = FALSE)[(length(breaks.ar)-1+2):(1+2)]

for(h in 1:nagegroups){
  data <- matrix(data = mean.ar.age.matrix[,,h]*100,
                 nrow = nrow(mean.ar.age.matrix[,,h]),
                 ncol = ncol(mean.ar.age.matrix[,,h]))
  data <- t(data)
  system("mkdir figs")
  makefigs=T
  if(makefigs==T){
    pdf(paste("figs/image_AR_age_group_",h,".pdf",sep=""),
        width=6,height=8,
        paper="special")#,horizontal=F)
    par(mfrow=c(1,1),mar=c(5,9,9,1),mgp=c(7,1.3,0.5),xpd=T,family="Barlow Semi Condensed")  # margine sup
    image(x=c(1:nx),y=c(1:ny),
          z=data,breaks=breaks.ar,
          col=cols_mat,axes=F,xlab="",
          ylab="",
          main="",cex.lab=2,cex.main=2)
    par(mgp=c(7,4.3,0.5))
    par(mar=c(1, 1, 1, 0))
    for(i in 1:nx){
      for(j in 1:ny){
        text(x=i,y=j,round(data[i,j],digits=1),cex=1.8)
      }
    }
    # Loop per aggiungere bordini
    for(i in 1:nx){
      for(j in 1:ny){
        rect(xleft = i - 0.5, ybottom = j - 0.5,
             xright = i + 0.5, ytop = j + 0.5,
             border = "white", lwd=1)
      }
    }
    text(x=0.35,y=1:ny,lab.y,cex=1.8,srt=0,adj=1)
    text(x=-0.7,y=ny/2+0.5,"Level of school closure",cex=2,srt=90)
    text(x=nx/2+0.5,y=-0.4,"Level of in-person work",cex=2,srt=0)
    
    text(x=1:nx,y=0.2,lab.x,cex=1.8,srt=0)
    text(x=nx/2+0.5,y=ny*1.15,paste("Infection attack rate (%)",sep=""),cex=2,srt=0)
    
    xcoord.legend <- seq(0.75,nx+0.25,length.out=length(breaks.ar))
    for(i in 1:(length(breaks.ar)-1)){
      rect(xcoord.legend[i],ny+1.2,xcoord.legend[i+1],ny+1.7,col=cols_mat[i],border=NA)
    }
    
    text(xcoord.legend,ny+1.9,breaks.ar,cex=1.8)
    title(paste(lab.age.groups[h], sep=""),cex.main=2.5,line=-0.5)
    dev.off()
  }
}

## REDUCTIONS IN AGE SPECIFIC- ATTACK RATES


reduction.sim <- array(NA,c(nscen,nsim,nagegroups))
for (i in 1:nsim){
  for(s in 1:nscen){
    reduction.sim[s,i,1] <- (cuminf.agegroups[1,i,1]-cuminf.agegroups[s,i,1])/cuminf.agegroups[1,i,1]
    reduction.sim[s,i,2] <- (cuminf.agegroups[1,i,2]-cuminf.agegroups[s,i,2])/cuminf.agegroups[1,i,2]
    reduction.sim[s,i,3] <- (cuminf.agegroups[1,i,3]-cuminf.agegroups[s,i,3])/cuminf.agegroups[1,i,3]
    
  }
}

reduction.age.matrix <- array(NA,c(nscen.school,nscen.smart,nsim,nagegroups))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    for(i in 1:nsim){
      sel <- scenarios[which(scenarios[,2]==scen.school[s] & scenarios[,3]==scen.smart[t]),1]
      reduction.age.matrix[s,t,i,1] <- reduction.sim[sel,i,1] 
      reduction.age.matrix[s,t,i,2] <- reduction.sim[sel,i,2] 
      reduction.age.matrix[s,t,i,3] <- reduction.sim[sel,i,3] 
    }
  }
}



bis.reduction.matrix <- array(NA,c(nscen.school,nscen.smart,nagegroups))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    bis.reduction.matrix[s,t,1] <- mean(reduction.age.matrix[s,t,,1])
    bis.reduction.matrix[s,t,2] <- mean(reduction.age.matrix[s,t,,2])
    bis.reduction.matrix[s,t,3] <- mean(reduction.age.matrix[s,t,,3])
  }
}



breaks.ar <- seq(0,100,by=10)
labels <- as.character(breaks.ar)
nb <- length(breaks.ar)

breaks.ar.sub <- breaks.ar[seq(1,nb,by=1)]
label.sub <- as.character(breaks.ar.sub)

paletta <- diverge_hcl(nb,c=100,l=c(50,90),power=1)

lab.x <- c("baseline","sustainable","minimum")
lab.y <- lab.scen.school
nx <- nscen.smart
ny <- nscen.school

cols_mat=hcl.colors((length(breaks.ar))+3,
                    "BluGrn", rev = FALSE)[(length(breaks.ar)-1+2):(1+2)]

for(h in 1:nagegroups){
  data <- matrix(data = bis.reduction.matrix[,,h]*100,
                 nrow = nrow(bis.reduction.matrix[,,h]),
                 ncol = ncol(bis.reduction.matrix[,,h]))
  data <- t(data)
  system("mkdir figs")
  makefigs=T
  if(makefigs==T){
    pdf(paste("figs/image_AR_reduction_age_group_",h,".pdf",sep=""),
        width=6,height=8,
        paper="special")#,horizontal=F)
    par(mfrow=c(1,1),mar=c(5,9,9,1.5),mgp=c(7,1.3,0.5),xpd=T,family="Barlow Semi Condensed")  # margine sup
    image(x=c(1:nx),y=c(1:ny),
          z=data,breaks=breaks.ar,
          col=cols_mat,axes=F,xlab="",
          ylab="",
          main="",cex.lab=2,cex.main=2)
    par(mgp=c(7,4.3,0.5))
    par(mar=c(1, 1, 1, 0))
    for(i in 1:nx){
      for(j in 1:ny){
        text(x=i,y=j,round(data[i,j],digits=1),cex=1.8)
      }
    }
    for(i in 1:nx){
      for(j in 1:ny){
        rect(xleft = i - 0.5, ybottom = j - 0.5,
             xright = i + 0.5, ytop = j + 0.5,
             border = "white", lwd=1)
      }
    }
    text(x=0.35,y=1:ny,lab.y,cex=1.8,srt=0,adj=1)
    text(x=-0.7,y=ny/2+0.5,"Level of school closure",cex=2,srt=90)
    text(x=nx/2+0.5,y=-0.4,"Level of in-person work",cex=2,srt=0)
    
    text(x=1:nx,y=0.2,lab.x,cex=1.8,srt=0)
    text(x=nx/2+0.5,y=ny*1.15,paste("Infection attack rate reduction (%)",sep=""),cex=2,srt=0)
    
    xcoord.legend <- seq(0.75,nx+0.25,length.out=length(breaks.ar))
    for(i in 1:(length(breaks.ar)-1)){
      rect(xcoord.legend[i],ny+1.2,xcoord.legend[i+1],ny+1.7,col=cols_mat[i],border=NA)
    }
    
    text(xcoord.legend,ny+1.9,breaks.ar,cex=1.8,srt=90)
    title(paste(lab.age.groups[h], sep=""),cex.main=2.5,line=-0.5)
    dev.off()
  }
}

