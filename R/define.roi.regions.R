#' Function to assign samples as low intermediate or high based on ROI
#' 
#' This is a function to samples as low intermediate or high based on ROI
#' @param bresat Full return data returned from the sig.ranksum function.
#' @param random.rank.counts ROI region returned from random.ranks function.
#' @param middle.range The ROI range).
.define.roi.regions <- function (bresat, random.rank.counts, middle.range = 0.95) 
{
    random.ranks.cdf <- cumsum(random.rank.counts)/sum(random.rank.counts)
    left <- max(c(0, which(random.ranks.cdf < (1 - middle.range)/2)))
    right <- min(which(random.ranks.cdf > 1 - ((1 - middle.range)/2))) - 
        1
    ret <- rep(2, length(bresat$rank))
    ret[bresat$rank < left] <- 1
    ret[bresat$rank > right] <- 3
    ret
}
