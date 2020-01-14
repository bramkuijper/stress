library("lattice")

z <- seq(0,5,0.01)

baseline <- c(0.1,0.2,0.3)

zmax <- c(1,5, 10,100)

zpower <- c(0.5,1,2)

hormone.grid <- as.data.frame(
        expand.grid(
                z=z
                ,baseline=baseline
                ,zmax=zmax
                ,zpower=zpower))

survival <- function(row)
{
    zpower <- row["zpower"]
    baseline <- row["baseline"]
    z <- row["z"]
    zmax <- row["zmax"]

    val <- baseline + (1.0 - baseline) * (1.0 - (z / zmax)^zpower)

    if (val <= 0.0)
    {
        val <- 0.0
    }

    return(val)
}

hormone.grid[,"survival"] <- apply(X=hormone.grid,FUN=survival,MARGIN=1)

pdf("overview_survival.pdf")
print(
        xyplot(survival ~ z | baseline * zmax
                ,type="l"
                ,groups=zpower
                ,auto.key=T
                ,data=hormone.grid
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

    return(1.0 - (damage/dmax)^ad)
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




