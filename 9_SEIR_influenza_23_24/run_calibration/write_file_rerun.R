## RSCRIPT TO rerun simulations
rm(list=ls())
mcmcprogress <- read.table("output_files/mcmcprogress",sep="\t")
burn.in <- 20000
idsel <- (burn.in+1):nrow(mcmcprogress)
nsim <- length(idsel)

npar <- 3 ## beta; reporting; initial infections;
mcmcrerun <- matrix(NA,nsim,npar)

mcmcrerun[,1] <- mcmcprogress$V2[idsel]
mcmcrerun[,2] <- mcmcprogress$V3[idsel]
mcmcrerun[,3] <- mcmcprogress$V4[idsel]

write.table(mcmcrerun, "mcmcrerun", sep="\t", col.names = F, row.names = F)

write.table( mcmcprogress$V5[idsel], "overdisp", sep="\t", col.names = F, row.names = F)
