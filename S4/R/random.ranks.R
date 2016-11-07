#' Function to assign the random region for a expression matrix
#' 
#' This is a function to assign the random region for a dataset
#' @param bs Full return data returned from the sig.ranksum function.
#' @param n Number of random samples.
#' @param seed Seed to generate random data (for reproducibility).
#' @param mc.cores An integer specifying how many cores to use.
#' @param renice A boolean to specify if the child processes should be reniced to the lowest nice priority.
#' @export
#' @import parallel
random.ranks <- function (bs, n = 1000, seed = 123456, mc.cores = 2, renice = T) 
{
    datrank.up <- datrank.dn <- NULL
    if (any(bs$up.dn > 0)) 
        datrank.up <- t(apply(bs$dat[bs$up.dn > 0, , drop = FALSE], 
            1, function(x) {
                rank(x, "average", na.last = "keep")
            }))
    if (any(bs$up.dn < 0)) 
        datrank.dn <- ncol(bs$dat) - t(apply(bs$dat[bs$up.dn < 
            0, , drop = FALSE], 1, function(x) {
            rank(x, "average", na.last = "keep")
        })) + 1
    datrank <- rbind(datrank.up, datrank.dn)
    nvals.cols <- nrow(datrank) - c(colSums(is.na(datrank)), 
        0)
    nvals.rows <- ncol(datrank) - rowSums(is.na(datrank))
    set.seed(seed)
    random.cols <- t(sapply(1:nrow(datrank), function(i) {
        runif(n, 0, nvals.rows[i] + 1)
    }))
    cl.rr <- makeCluster(mc.cores)
    if (renice == T) 
        temp <- clusterCall(cl.rr, function(x) system(paste("renice", 
            19, Sys.getpid()), ignore.stdout = T))
    clusterExport(cl.rr, c("datrank", "nvals.cols"), envir = environment())
    rand.dist <- unlist(parCapply(cl.rr, random.cols, function(i) {
        datrank.rand.col <- cbind(datrank, i)
        datrank.rand.col <- t(apply(datrank.rand.col, 1, function(x) {
            rank(x, "average", na.last = "keep")
        }))
        ranksums <- colSums(datrank.rand.col, na.rm = TRUE)/nvals.cols
        tail(rank(ranksums, "average", na.last = TRUE), 1)
    })) - 1
    stopCluster(cl.rr)
    ret <- hist(rand.dist, breaks = 0:(ncol(datrank) + 1), plot = FALSE)$counts
    return(ret)
}
