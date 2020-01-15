#!/usr/bin/env Rscript

the.data <- read.table("summary_vary_autocorrelation_new.csv",header=T,sep=";")

init_feedback <- 0.0
ad <- 0.5
aP <- 0.5
r <- 1
u <- 0.5
sNP2P_1 <- 0.01
sP2NP_1 <- 0.01


the.subset <- the.data[
    the.data$init_feedback == init_feedback
        &  the.data$ad == ad
        &  the.data$aP == aP
        &  the.data$r == r
        &  the.data$u == u
        &  the.data$sNP2P_1 == sNP2P_1
        &  the.data$sP2NP_1 == sP2NP_1,]

print("how many graphs we are printing?")
print(nrow(the.subset))


command = "../../src/ibm/plot_simulation.py"

for (i in 1:nrow(the.subset))
{
    system(paste(command,as.character(the.subset[i,"file"])))
}
