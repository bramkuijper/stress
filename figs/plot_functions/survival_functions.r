#!/usr/bin/env Rscript 


library("lattice")

z <- seq(0,1,0.01)

baseline_mort <- c(0,0.1,0.2,0.3)

zmax <- c(1)

baseline_survival <- c(0,0.1,0.2,0.3)

zpower <- c(0.5,1,2)

hormone.grid <- as.data.frame(
        expand.grid(
                z=z
                ,baseline_mort=baseline_mort
                ,baseline_survival=baseline_survival
                ,zmax=zmax
                ,zpower=zpower))

mortality <- function(row)
{
    zpower <- row["zpower"]
    baseline_mort <- row["baseline_mort"]
    baseline_survival <- row["baseline_survival"]
    z <- row["z"]
    zmax <- row["zmax"]

    val <- baseline_mort + (1.0 - baseline_mort) * (1.0 - (z / zmax)^zpower)

    if (val <= 0.0)
    {
        val <- 0.0
    }

    survival <- baseline_survival + (1.0 - baseline_survival) * (1.0 - val)

    return(1.0 - survival)
}

hormone.grid[,"mortality"] <- apply(X=hormone.grid,FUN=mortality,MARGIN=1)

pdf("overview_survival.pdf",width=10,height=10)
print(
        xyplot(mortality ~ z | baseline_mort * baseline_survival
                ,type="l"
                ,groups=zpower
                ,auto.key=T
                ,data=hormone.grid
                ,ylim=c(-0.1,1.1)
                ,xlab="Hormone level"
                ,ylab="Probability of mortality in presence of predator"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                )
        )
dev.off()

################# damage affecting fecundity #################

damage <- seq(0,5,0.01)

dmax <- c(1,5,10,100)

ad <- c(0.5,1,2.0)

damage.grid <- as.data.frame(
        expand.grid(
                damage=damage
                ,dmax=dmax
                ,ad=ad))


fecundity <- function(row)
{
    ad <- row["ad"]
    damage <- row["damage"]
    dmax <- row["dmax"]

    val <- 1.0 - (damage/dmax)^ad

    if (val < 0.0)
    {
        val <- 0.0
    }

    return(val)
}

damage.grid[,"fecundity"] <- apply(
        X=damage.grid
        ,FUN=fecundity
        ,MARGIN=1)


pdf("overview_damage.pdf")
print(
        xyplot(fecundity ~ damage | dmax
                ,type="l"
                ,groups=ad
                ,auto.key=T
                ,data=damage.grid
                ,xlab="Amount of damage"
                ,ylab="Share of fecundity"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                )
)
dev.off()




