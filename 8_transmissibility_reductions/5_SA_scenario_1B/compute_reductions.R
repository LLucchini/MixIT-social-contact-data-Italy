# Clear the environment
rm(list = ls())

# Load required plotting functions
source("../../plotting_functions.R")

# Number of age groups
nages <- 15

# Number of bootstrap samples
NBOOT <- 1000

# Load contact matrices
load("input/matr_in_presence.Rdata")
load("input/matr_no_presence.Rdata")

## ============================== ##
## COMBINED SCENARIOS
## ============================== ##

# School closure scenarios
scen.school <- c("no_closure", "uni", "uni+upper_secondary",
                 "uni+secondary", "uni+secondary+primary", "all_closed")

# Smart-working (in presence / remote) scenarios
scen.smart <- c("pre-lockdown", "lifting lockdown", "lockdown")

nscen.school <- length(scen.school)

# Population distributions per scenario:
# dim = school scenario × smart-working scenario × age × presence (1=in presence, 2=remote)
pop.scen.combined <- array(NA, c(nscen.school, length(scen.smart), nages, 2))

# Load population data for each combined scenario
for(s in 1:nscen.school){
  for(t in 1:length(scen.smart)){
    pop.scen.combined[s, t, , ] <- as.matrix(
      read.csv(paste("../0_population_scenarios/pop_scenarios/scenario_smart_",
                     t, "_sclose_", scen.school[s], ".csv", sep = ""))
    )
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
  NGM[[n]] <- vector("list", nscen.school)
  
  for(s in 1:nscen.school){
    NGM[[n]][[s]] <- vector("list", length(scen.smart))
    
    for(t in 1:length(scen.smart)){
      # Compute contact blocks weighted by population
      for(i in 1:nages){
        for(j in 1:nages){
          block.inpresence[n,s,t,i,j] <-
            matr.inpresence[n,1,i,j] *
            pop.scen.combined[s,t,i,1] / sum(pop.scen.combined[s,t,j,])
          block.nopresence[n,s,t,i,j] <-
            matr.nopresence[n,1,i,j] *
            pop.scen.combined[s,t,i,2] / sum(pop.scen.combined[s,t,j,])
        }
      }
      
      # Build full 2x2 block NGM matrix
      NGM[[n]][[s]][[t]] <- rbind(
        cbind(block.inpresence[n,s,t,,], block.inpresence[n,s,t,,]),
        cbind(block.nopresence[n,s,t,,], block.nopresence[n,s,t,,])
      )
    }
  }
}

## ============================== ##
## DOMINANT EIGENVALUES (Rho)
## ============================== ##

rho <- array(NA, c(NBOOT, nscen.school, length(scen.smart)))

for(n in 1:NBOOT){
  for(s in 1:nscen.school){
    for(t in 1:length(scen.smart)){
      # Dominant eigenvalue
      rho[n,s,t] <- max(Mod(eigen(NGM[[n]][[s]][[t]], only.values = TRUE)$values))
    }
  }
}

## ============================== ##
## PERCENTAGE REDUCTION VS BASE CASE
## ============================== ##

reductions.combined <- array(NA, c(NBOOT, nscen.school, length(scen.smart)))

for(n in 1:NBOOT){
  for(s in 1:nscen.school){
    for(t in 1:length(scen.smart)){
      reductions.combined[n,s,t] <-
        (1 - rho[n,s,t] / rho[n,1,1]) * 100
    }
  }
}

# Store mean and 95% CI of reductions
stat.reductions.boot <- array(NA, c(nscen.school, length(scen.smart), 3))

for(s in 1:nscen.school){
  for(t in 1:length(scen.smart)){
    stat.reductions.boot[s,t,1] <- mean(reductions.combined[,s,t])
    stat.reductions.boot[s,t,2] <- quantile(reductions.combined[,s,t], 0.025)
    stat.reductions.boot[s,t,3] <- quantile(reductions.combined[,s,t], 0.975)
  }
}

## ============================== ##
## PLOT VISUALIZATION
## ============================== ##

library(sysfonts)
library(showtext)
font_add_google("Barlow Semi Condensed")
showtext_auto()

par(family = "Barlow Semi Condensed")

verde <- "#1a9850"
cols <- c(lighten(verde, 0.6), lighten(verde,0.3), verde)

# Labels for plots
lab.scen.smart <- c("baseline (pre-pandemic)",
                    "sustainable (no closure of economic activities)",
                    "minimum (only essential activities)")

lab.scen.school <- c("no level\nclosed",
                     "ISCED 5-8",
                     "ISCED 5-8\nISCED 3",
                     "ISCED 5-8\nISCED 2-3",
                     "ISCED 5-8\nISCED 1-3",
                     "all levels\nclosed")

makefigs <- TRUE

if(makefigs){
  pdf(paste("Figure_SI13.pdf",sep=""),
      width=16,height=8)
  par(mfrow = c(1,1), mar = c(8,7,1,1), mgp = c(6,1.5,0.5),
      xpd = TRUE, family = "Barlow Semi Condensed")
  
  STEP <- 5
  MIN <- 0
  MAX <- 40
  
  plot(c(0.5,nscen.school+0.5), c(MIN,MAX), type = "n",
       xlab = "", ylab = "", axes = FALSE)
  axis(2, cex.axis = 2, at = seq(MIN, MAX, by = STEP), las = 2)
  
  # Extract stats per smart-working level
  min.RW <- stat.reductions.boot[,1,]
  sus.RW <- stat.reductions.boot[,2,]
  max.RW <- stat.reductions.boot[,3,]
  
  # Draw bars per scenario
  for(s in 1:nscen.school){
    if(min.RW[s,1] > 0){
      makerect(s-0.2, min.RW[s,1], 0.2, col = cols[1])
      segments(s-0.2, min.RW[s,2], s-0.2, min.RW[s,3],
               col = darken(cols[1],0.2), lwd = 3)
    }
    
    makerect(s, sus.RW[s,1], 0.2, col = cols[2])
    makerect(s+0.2, max.RW[s,1], 0.2, col = cols[3])
    
    segments(s, sus.RW[s,2], s, sus.RW[s,3],
             col = darken(cols[2],0.2), lwd = 3)
    segments(s+0.2, max.RW[s,2], s+0.2, max.RW[s,3],
             col = darken(cols[3],0.2), lwd = 3)
  }
  
  text(x = -nscen.school*0.04, y = MAX/2,
       "Transmissibility reduction (%)", cex = 2, srt = 90)
  
  text(x = 1:nscen.school, y = -3,
       lab.scen.school, cex = 2, adj = 0.5)
  
  text(x = nscen.school/1.8, y = -9,
       "Level of school closure", cex = 2)
  text(x = 1.1, y = 35,
       "Level of in-person work:", cex = 2)
  
  legend(x = 0.2, y = 34.5, lab.scen.smart,
         pch = 15, col = cols, bty = "n", cex = 2)
  
  dev.off()
}

