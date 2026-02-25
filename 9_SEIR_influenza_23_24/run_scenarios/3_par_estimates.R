rm(list=ls())
library(data.table)

mcmcrerun <- read.table("../calibration/sw/mcmcrerun",sep="\t")
overdisp <- read.table("../calibration/sw/overdisp",sep="\t")

# Esempio: tre vettori di 80.000 elementi ciascuno
set.seed(123)
v1 <- mcmcrerun[,1]
v2 <- mcmcrerun[,2]
v3 <- mcmcrerun[,3]
v4 <- overdisp[,1]
library(sysfonts)
library(showtext)
# Add the Barlow font
font_add_google("Barlow Semi Condensed")
##"Barlow Condensed"   "Barlow Semi Condensed"   
# Enable showtext
showtext_auto()

matr <- cbind(v1,v2,v3,v4)

ESTIMATES <- cbind(colMeans(matr),t(apply(matr,2,quantile,c(0.025,0.975))))
row.names(ESTIMATES) <- c("beta","reporting","i0","overdisp")
write.table(ESTIMATES,file="PAR_ESTIMATES",sep="\t",row.names=T)
# Imposta layout: 2 righe x 2 colonne
pdf("posterior_distributions.pdf",   width=8,height=8,
    paper="special")

par(mfrow = c(2, 2),mgp=c(7,1.3,0.5),
    mar = c(4, 4, 2, 1),xpd=T,family="Barlow Semi Condensed")  # margini (bottom, left, top, right)

# Istogramma 1
hist(v1,
     main = "Scaling factor for transmissibility",
     xlab = "Values",
     col = "lightblue",
     border = "white")

# Istogramma 2
hist(v2,
     main = "Reporting ratio",
     xlab = "Values",
     col = "lightgreen",
     border = "white")

# Istogramma 3
hist(v3,
     main = "Infectious at simulation start",
     xlab = "Values",
     col = "lightcoral",
     border = "white")

# Istogramma 4
hist(v4,
     main = "Overdispersion",
     xlab = "Values",
     col = "lightgoldenrod",
     border = "white")

dev.off()

