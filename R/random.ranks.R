#' Function to assign the random region for a expression matrix
#'
#' This is a function to assign the random region for a dataset
#' @param bs Full return data returned from the sig.ranksum function.
#' @param n Number of random samples.
#' @param seed Seed to generate random data (for reproducibility).
#' @param mc.cores An integer specifying how many cores to use.
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom stats runif
#' @importFrom utils tail
#' @importFrom graphics hist
.random.ranks <- function (bs, n = 1000, seed = 123456, workers = snowWorkers())
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
    nvals.cols <- nrow(datrank) - c(colSums(is.na(datrank)), 0)
    nvals.rows <- ncol(datrank) - rowSums(is.na(datrank))

    set.seed(seed)
    random.cols <- t(sapply(seq_len(nrow(datrank)), function(i) {
        runif(n, 0, nvals.rows[i] + 1)
    }))

    snow <- SnowParam(workers = workers, type = "SOCK")

    rand.dist <- unlist(bplapply(seq_len(ncol(random.cols)), datrank=datrank, nvals.cols=nvals.cols, BPPARAM = snow, FUN=function(i, datrank, nvals.cols){

      datrank.rand.col <- cbind(datrank, random.cols[, i])
      matrixStats::rowRanks(datrank.rand.col)
      ranksums <- colSums(datrank.rand.col, na.rm = TRUE)/nvals.cols
      tail(rank(ranksums, "average", na.last = TRUE), 1)

    })) - 1

    ret <- hist(rand.dist, breaks = seq(0,(ncol(datrank) + 1)), plot = FALSE)$counts
    return(ret)
}
