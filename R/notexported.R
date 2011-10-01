## 2011-04-07 mash(), cluster.count(), and tile.popn  not exported yet

## 2011-04-07, 2011-04-08
## not exported

tile.popn <- function (popn, method = "reflect") {
    bbox <- attr(popn, 'boundingbox')
    if (method== "reflect") {
        p2 <- rbind.popn(popn, flip(popn,lr=min(bbox$x)), flip(popn,lr=max(bbox$x)))
        rbind.popn(p2, flip(p2,tb=min(bbox$y)), flip(p2,tb=max(bbox$y)))
    }
    else if (method == "copy") {
        ht <- max(bbox$y) - min(bbox$y)
        wd <- max(bbox$x) - min(bbox$x)
        p2 <- rbind.popn(popn, shift(popn,c(-wd,0)), shift(popn, c(wd,0)))
        rbind.popn(p2, shift(p2,c(0,-ht)), shift(p2, c(0,ht)))
    }
    else
        stop ("unrecognised method")
}

## demonstration
# pop <- sim.popn(D=100, core=make.grid(),model2D='coastal')
# par(mfrow=c(1,2), mar=c(2,2,2,2))
# plot(tile.popn(pop, 'copy'))
# polygon(cbind(-100,200,200,-100),c(-100,-100,200,200), col='red', density=0)
# plot(tile.popn(pop, 'reflect'))
# polygon(cbind(-100,200,200,-100),c(-100,-100,200,200), col='red', density=0)
