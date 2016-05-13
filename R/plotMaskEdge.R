plotMaskEdge <- function (mask, add = FALSE, ...) {
    ## plots outer border of a mask, however irregular
    ## assumes integer x,y OK
    draw <- function (x,y, vertical) {
        if (vertical)
            segments (x, y+spc/2, x, y-spc/2, ...)
        else
            segments (x+spc/2, y, x-spc/2, y, ...)
    }
    inmask <- function (x,y) {
        !is.na(match (paste(round(x),round(y),sep='.'), maskvector))
    }
    onecell <- function (xy) {
        x <- xy[1]
        y <- xy[2]
        if (!inmask(x-spc, y)) draw (x-spc/2, y, TRUE)
        if (!inmask(x+spc, y)) draw (x+spc/2, y, TRUE)
        if (!inmask(x, y-spc)) draw (x, y-spc/2, FALSE)
        if (!inmask(x, y+spc)) draw (x, y+spc/2, FALSE)
    }
    if (!add)
        plot(mask, dots = FALSE)
    spc <- spacing(mask)
    maskvector <- apply(round(mask),1,paste, collapse='.')
    apply(mask, 1, onecell)
    invisible()
}

#########################################
