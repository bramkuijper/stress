library("lattice")

the.data <- read.table("../../data/summary_stress.csv"
                       ,sep=";"
                       ,header=T)


the.names <- names(the.data)
the.names[the.names == "r"] <- "damage_decay"
the.names[the.names == "u"] <- "hormone_damage"
the.names[the.names == "mean_feedback"] <- "hormone_decay"

names(the.data) <- the.names

print(the.names)

str(the.data)

pdf("feedback_vary_autocorr.pdf")
print(
      cloud(
             hormone_decay ~ sNP2P_1 + sP2NP_1 | cue_P * cue_NP * damage_decay * hormone_damage
             ,data=the.data
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,default.scales=list(arrows=F)
             )
      )
dev.off()


pdf("influx_vary_autocorr.pdf")
print(
      cloud(
             (mean_stress_influx - mean_influx) ~ sNP2P_1 + sP2NP_1 | cue_P * cue_NP * damage_decay * hormone_damage
             ,data=the.data
             ,auto.key=T
             ,pch=21
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,default.scales=list(arrows=F)
             )
      )
dev.off()

# create dataset that has 30 timepoints
time.data <- NULL

nrep <- 1

for (row_i in seq(1,nrow(the.data),1))
{
    # get hormone level
    mean_hormone <- the.data[row_i,"mean_hormone"]
    var_hormone <- the.data[row_i,"var_hormone"]

    influx <- the.data[row_i,"mean_influx"]
    var_influx <- the.data[row_i,"var_influx"]

    stress_influx <- the.data[row_i,"mean_stress_influx"]
    var_stress_influx <- the.data[row_i,"var_stress_influx"]

    hormone_decay <- the.data[row_i,"hormone_decay"]
    var_hormone_decay <- the.data[row_i,"var_feedback"]

    sNP2P_1 <- the.data[row_i,"sNP2P_1"]
    sP2NP_1 <- the.data[row_i,"sP2NP_1"]
    cue_P <- the.data[row_i,"cue_P"] 
    cue_NP <-the.data[row_i,"cue_NP"]
    damage_decay <- the.data[row_i,"damage_decay"]
    hormone_damage <- the.data[row_i,"hormone_damage"]

    zt_vals <- rnorm(n=nrep
                     ,mean=mean_hormone
                     ,sd=sqrt(var_hormone))

    decay_vals <- rnorm(n=nrep
                        ,mean=hormone_decay
                        ,sd=sqrt(var_hormone_decay))

    influx_vals <- rnorm(n=nrep
                         ,mean=influx
                         ,sd=sqrt(var_influx))

    stress_influx_vals <- rnorm(n=nrep
                                ,mean=stress_influx
                                ,sd=sqrt(var_stress_influx))

    for (i in seq(1,nrep,1))
    {
        stress_tplus1 = 0

        stress_t = zt_vals[[i]]

        for (t in seq(0,30,1))
        {
            stress_tplus1 = stress_t * (1.0 - decay_vals[[i]]) + influx_vals[[i]]

            if (t == 10)
            {
                stress_tplus1 = stress_tplus1 + stress_influx_vals[[i]]
            }

            the_list <- list(
                             "sNP2P_1"=sNP2P_1
                             ,"sP2NP_1"=sP2NP_1
                             ,"cue_P"=cue_P
                             ,"cue_NP"=cue_NP
                             ,"damage_decay"=damage_decay
                             ,"hormone_damage"=hormone_damage
                             ,"replicate"=i
                             ,"time"=t
                             ,"stress"=stress_tplus1)

            stress_t <- stress_tplus1

            if (is.null(time.data))
            {
                time.data <- as.data.frame(the_list)
            }
            else
            {
                time.data <- rbind(time.data,the_list)
            }

        }
    }
}

str(time.data)


pdf("lines_stress.pdf")
print(xyplot(
             stress ~ time | sNP2P_1 * sP2NP_1 * cue_P * cue_NP * damage_decay * hormone_damage
             ,data=time.data
             ,groups=replicate
             ,ylim=c(0,10)
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ))
dev.off()


