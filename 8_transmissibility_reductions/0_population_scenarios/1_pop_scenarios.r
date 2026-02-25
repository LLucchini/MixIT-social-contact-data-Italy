rm(list=ls())
# Load required plotting functions
source("../../plotting_functions.R")

# Load the libraries
library(sysfonts)
library(showtext)

# Add the Barlow font
font_add_google("Barlow Semi Condensed")
##"Barlow Condensed"   "Barlow Semi Condensed"   
# Enable showtext
showtext_auto()

# Set the base font to Barlow
par(family = "Barlow Semi Condensed")


## FIGURE EDUCATION
cols <- c("#fbdb86","#f8af7b","#ed938e","#d1e6f6","#a8ceeb","#8caeca")


ns.levels <- 5
lab.levels <- c("Early childhood\n(ISCED 0)","Primary\n(ISCED 1)",
                "Lower secondary\n(ISCED 2)","Upper secondary\n(ISCED 3)",
                "Tertiary\n(ISCED 5-8)")
matr <- read.table("input/pop_age_education_level.txt",sep="\t")


cumsum.matr <- apply(matr,2,cumsum)
cumsum.matr
agegr <- c("0-4","5-9","10-14","15-19","20-24","25+")

tot <- colSums(matr)
system("mkdir figs/")
makefigs=T
if(makefigs==T){
  pdf(paste("figs/Figure_3A.pdf",sep=""),
      width=8,height=8)
  par(mfrow=c(1,1),mar=c(8,12,8,1),mgp=c(6,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 0.5
  xcoord.labx<- seq(1,ns.levels,by=1)
  MIN <- 0
  MAX <- 3
  plot(c(MIN,MAX), c(0.5,ns.levels+0.5),type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(3,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=1)
  
  text(x=MAX/2,y=ns.levels*1.4,
       "Population (mln) enrolled to different\nlevels of education (ISCED 2011)",
       cex=2,
       srt=0)
  
  text(x=0,y=(ns.levels:1)-0.05,
       rev(lab.levels),
       cex=2,
       srt=0,srt=0,pos=2)
  
  for(s in 1:ns.levels){
    for(a in 6:1){
      if(cumsum.matr[a,s]>0){
        ##rect(0,ns.levels+1-s-0.3,cumsum.matr[a,s]/10^6,ns.levels+1-s+0.3,col=cols[a],border=(cols[a]))
        rect(0,s-0.3,cumsum.matr[a,s]/10^6,s+0.3,col=cols[a],border=(cols[a]))
      }
    }
  }
  text(x=1,y=0.25,"Age group (years):",cex=2)
  legend(x=-0.5,y=0.2,agegr[1:3],
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2)
  legend(x=1,y=0.2,agegr[4:6],
         pch=15,col=cols[4:6],bty='n',lwd=NA,cex=2)
  dev.off()
}


## POPULATION WORKING IN PERSON

nw.levels <- 3
lab.legend <- c("baseline (pre−pandemic)",
                "sustainable (no closure of economic activities)",
                "minimum (only essential activities)")

pop.workers <- read.table("input/pop_in_person_workers.txt",sep="\t")
baseline.v <- pop.workers[,1]
sustainable.v <- pop.workers[,2]
minimum.v <- pop.workers[,3]

## aggrego in age groups 
baseline.v.aggr <- c(sum(baseline.v[4:6]),sum(baseline.v[7:8]),sum(baseline.v[9:10]),
                     sum(baseline.v[11:12]),sum(baseline.v[13:15]))
sustainable.v.aggr <- c(sum(sustainable.v[4:6]),sum(sustainable.v[7:8]),sum(sustainable.v[9:10]),
                        sum(sustainable.v[11:12]),sum(sustainable.v[13:15]))
minimum.v.aggr <- c(sum(minimum.v[4:6]),sum(minimum.v[7:8]),sum(minimum.v[9:10]),
                    sum(minimum.v[11:12]),sum(minimum.v[13:15]))


viola <- "#88419d"
cols <- c(lighten(viola, 0.6),lighten(viola,0.3),viola)

agegr <- c("<30","30-39","40-49","50-59","60+")
nagegr <- length(agegr)
makefigs=T
if(makefigs==T){
  pdf(paste("figs/Figure_3B.pdf",sep=""),
      width=8,height=8)
  par(mfrow=c(1,1),mar=c(12,7,1,1),mgp=c(6,1.5,0.5),xpd=T,family="Barlow Semi Condensed")
  STEP <- 0.5
  xcoord.labx<- seq(1,nagegr,by=1)
  MIN <- 0
  MAX <- 6.5
  plot(c(0.5,nagegr+0.5),c(MIN,MAX),type="l",col=NA,xlab="",
       ylab="",axes=F)
  axis(2,cex.axis = 2,at=seq(MIN,MAX,by=STEP),las=1)
  
  for(a in 1:nagegr){
    makerect(a-0.2,baseline.v.aggr[a]/10^6,0.2,col=cols[1])
    makerect(a,sustainable.v.aggr[a]/10^6,0.2,col=cols[2])
    makerect(a+0.2,minimum.v.aggr[a]/10^6,0.2,col=cols[3])
    
  }
  text(x=-0.5,y=MAX/2,
       "Population working in person (mln)",
       cex=2,
       srt=90)
  
  text(x=(1:nagegr),y=-0.1,
       agegr,
       cex=2,
       srt=0,pos=1)
  text(x=nagegr/2+0.5,y=-0.8,"Age group (years)",cex=2)
  text(x=1,y=-1.4,"Level of in-person work:",cex=2)
  # legend(x=-0.5,y=0.2,agegr[1:3],
  #        pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2)
  legend(x=-0.5,y=-1.4,lab.legend,
         pch=15,col=cols[1:3],bty='n',lwd=NA,cex=2)
  
  dev.off()
}


##============================================##
## MAKE POP FILES COMBINED SCENARIOS      ##
##============================================##
## define scenarios with different levels of smart working and different levels of school closure
breaks <- c(seq(0,70,by=5),101)
agegr <- unique(cut(0:100,breaks,include.lowest = T,right=F))
tmp <- cut(0:100,breaks,include.lowest = T,right=F,labels = F)
nagegr <- max(tmp)
scen.school <- c("no_closure","uni","uni+upper_secondary",
                 "uni+secondary","uni+secondary+primary","all_closed")
scen.smart <- c("pre-lockdown","lifting lockdown", "lockdown")
nscen.smart <- length(scen.smart)
nscen.school <- length(scen.school)
pop.scenario.max.presence <- read.csv("input/pop_scenario_max_presence.csv")

# ##a=cbind(pop.scenario.max.presence[1:15,1],pop.scenario.medium.presence[1:15,1],pop.scenario.min.presence[1:15,1])
# a=pop.workers
# active.workers <- matrix(NA,nrow(pop.workers),nscen.smart)
# active.workers[,1] <- a[,1]
# active.workers[,2] <- a[,2]
# active.workers[,nscen.smart] <- a[,3]
# 
active.workers <- pop.workers

pop.scenarios.combined <- array(0,c(length(scen.school),length(scen.smart),nagegr,2))


for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    ##prima colonna: popolazione attiva.
    pop.scenarios.combined[s,t,,1] <- active.workers[,t]
    if(scen.school[s]=="no_closure"){
      students.at.home <- 0
      pop.scenarios.combined[s,t,,2] <- 0
      students.at.school <- rowSums(pop.scenario.max.presence[,2:6])
      pop.scenarios.combined[s,t,,1] <- (pop.scenarios.combined[s,t,,1]+
                                           students.at.school)                                         
    }else if(scen.school[s]=="uni"){
      students.at.home <- pop.scenario.max.presence[,6]
      pop.scenarios.combined[s,t,,2] <- students.at.home
      students.at.school <- rowSums(pop.scenario.max.presence[,2:5])
      pop.scenarios.combined[s,t,,1] <- (pop.scenarios.combined[s,t,,1]+
                                           students.at.school)    
    }else if(scen.school[s]=="uni+upper_secondary"){
      students.at.home <- rowSums(pop.scenario.max.presence[,5:6])
      pop.scenarios.combined[s,t,,2] <- students.at.home
      students.at.school <- rowSums(pop.scenario.max.presence[,2:4])
      pop.scenarios.combined[s,t,,1] <- (pop.scenarios.combined[s,t,,1]+
                                           students.at.school) 
    }else if(scen.school[s]=="uni+secondary"){
      students.at.home <- rowSums(pop.scenario.max.presence[,4:6])
      pop.scenarios.combined[s,t,,2] <- students.at.home
      students.at.school <- rowSums(pop.scenario.max.presence[,2:3])
      pop.scenarios.combined[s,t,,1] <- (pop.scenarios.combined[s,t,,1]+
                                           students.at.school) 
    }else if(scen.school[s]=="uni+secondary+primary"){
      students.at.home <- rowSums(pop.scenario.max.presence[,3:6])
      pop.scenarios.combined[s,t,,2] <- students.at.home
      students.at.school <- pop.scenario.max.presence[,2]
      pop.scenarios.combined[s,t,,1] <- (pop.scenarios.combined[s,t,,1]+
                                           students.at.school) 
    }else if(scen.school[s]=="all_closed"){
      students.at.home <- rowSums(pop.scenario.max.presence[,2:6])
      pop.scenarios.combined[s,t,,2] <- students.at.home
    }
    inactive <- (pop.scenario.max.presence[,7]+pop.scenario.max.presence[,1]-active.workers[,t])
    pop.scenarios.combined[s,t,,2] <- pop.scenarios.combined[s,t,,2]+inactive
  }
}
    
system("mkdir pop_scenarios")
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    write.csv(pop.scenarios.combined[s,t,,],
            paste("pop_scenarios/scenario_smart_",t,"_sclose_",scen.school[s],".csv",sep=""),
            row.names = F)
  }
}

