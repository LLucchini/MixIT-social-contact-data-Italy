rm(list=ls())
library(data.table)

scenarios <- read.csv("readme_scenarios.csv")
idscen <- scenarios[,1]
nscen <- nrow(scenarios)
mcmcrerun <- read.table("common_input/mcmcrerun",sep="\t")
folders <- paste("scenario_",idscen,sep="")
nsim <- nrow(mcmcrerun)
nages <- 15
nweeks <- 24
idsim <- 0:(nsim-1)

popbyage <- array(NA,c(nages,nweeks))
popbyage <- as.matrix(read.table(paste(folders[1],"/output_files/populationByAge_0",sep=""),sep="\t")[,1:nweeks])
totp <- colSums(popbyage)[1]

cuminfbyageP <- array(NA,c(nsim,nages,nweeks))
cuminfbyageNP <- array(NA,c(nsim,nages,nweeks))



totinf <- array(NA, c(nscen, nsim, nweeks))
cuminfage <- array(NA, c(nscen, nsim, nages))

pb <- txtProgressBar(min = 0, max = nscen * nsim, style = 3)
counter <- 0

for (s in 1:nscen) {
  for (i in 1:nsim) {
    outputP <- fread(file.path(folders[s], "output_files", paste0("CumInfByAgeP_", idsim[i])))[, 1:nweeks]
    outputNP <- fread(file.path(folders[s], "output_files", paste0("CumInfByAgeNP_", idsim[i])))[, 1:nweeks]
    
    cuminfage[s, i, ] <- rowSums(outputP) + rowSums(outputNP)
    totinf[s, i, ]    <- colSums(outputP) + colSums(outputNP)
    
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
}
close(pb)
system("mkdir Rdata")
save(totinf,file="Rdata/totinf.Rdata")
save(cuminfage,file="Rdata/cuminfage.Rdata")

