#' Function to assign ranksums to an expression matrix
#'
#' This is a function to assign ranksums to a dataset
#' @param exprdata An expression matrix. Rows correspont to genes, columns to samples. Genes with multiple probes should be collapsed to a single row e.g. using most variable probe.
#' @param up A vector specifying rows of genes in the up direction.
#' @param dn A vector specifying rows of genes in the down direction.
#' @param ns A vector specifying rows of genes with an unspecified direction. These will be split using the pam function from the cluster package.
#' @param full.return A boolean specifying if the full result should be returned (T) or only the ranksum values (F).
#' @import cluster
#' @importFrom matrixStats rowSds
#' @importFrom stats cor
#' @importFrom matrixStats rowSds
.sig.ranksum <- function (exprdata, up = NULL, dn = NULL, ns = NULL, full.return = FALSE)
{
    if (length(dim(exprdata)) != 2)
        stop("'exprdata' must be a matrix")
    if (length(up) == 0 && length(dn) == 0 && length(ns) == 0)
        stop("no indices were specified")
    if (length(ns) > 0 && (length(up) > 0 || length(dn) > 0))
        stop("both directional and non-directional indices were specified")
    if (is.logical(up))
        up <- which(up)
    if (is.logical(dn))
        dn <- which(dn)
    if (is.logical(ns))
        ns <- which(ns)
    if (ncol(exprdata) < 2) {
        if (isTRUE(full.return)) {
            ret <- list()
            up <- c(up, ns)
            ret$rank <- 1
            ret$pat.order <- 1
            ret$gene.order <- c(up, dn)
            ret$dat <- exprdata[ret$gene.order, , drop = FALSE]
            ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
            return(ret)
        }
        else {
            return(1)
        }
    }
    if (length(ns) > 0) {
        row.sds <- rowSds(exprdata[ns, , drop = FALSE], na.rm = TRUE)
        tmp <- row.sds == 0 | is.na(row.sds)
        zero.sd.idx <- ns[tmp]
        ns <- ns[!tmp]
        if (length(ns) == 1) {
            up <- ns
        }
        else if (length(ns) == 2) {
            if (cor(exprdata[ns[1], ], exprdata[ns[2], ], use = "pairwise") <
                0) {
                up <- ns[1]
                dn <- ns[2]
            }
            else {
                up <- ns
            }
        }
        else if (length(ns) > 2) {
            diss <- 1 - cor(t(exprdata[ns, , drop = FALSE]),
                use = "pairwise")
            diss[which(is.na(diss))] <- 1
            diss[diss >= 1] <- diss[diss >= 1] + 1
            clustering <- pam(diss, k = 2, diss = TRUE, cluster.only = TRUE)
            up.cluster <- which.max(table(clustering))
            up <- ns[which(clustering == up.cluster)]
            dn <- ns[which(clustering != up.cluster)]
            up.dn.cor <- cor(t(exprdata[up, , drop = FALSE]),
                t(exprdata[dn, , drop = FALSE]), use = "pairwise")
            if (sum(up.dn.cor < 0, na.rm = TRUE) < length(up) *
                length(dn)/2) {
                up <- ns
                dn <- NULL
            }
            rm(diss, clustering, up.cluster, up.dn.cor)
        }
        if (length(zero.sd.idx) > 0) {
            up <- c(up, zero.sd.idx)
        }
    }
    ranksum <- double(ncol(exprdata))
    col.counts <- rep(0, ncol(exprdata))
    if (length(up) != 0) {
        dat <- exprdata[up, , drop = FALSE]
        ranksum <- rowSums(apply(dat, 1, function(x) {
            rank(x, "average", na.last = "keep")
        }), na.rm = TRUE)
        col.counts <- colSums(!is.na(dat))
    }
    if (length(dn) != 0) {
        dat <- exprdata[dn, , drop = FALSE]
        ranksum <- ranksum + rowSums(ncol(exprdata) - apply(dat,
            1, function(x) {
                rank(x, "average", na.last = "keep")
            }) + 1, na.rm = TRUE)
        col.counts <- col.counts + colSums(!is.na(dat))
    }
    ranksum <- ranksum/col.counts
    rank <- rank(ranksum, "average", na.last = TRUE)
    if (full.return == FALSE)
        return(rank)
    if (length(up) == 0)
        up <- NULL
    if (length(dn) == 0)
        dn <- NULL
    gene.cor <- function(gene.idx, is.up, exprdata, pat.order) {
        gene.expr <- exprdata[gene.idx, pat.order]
        if (is.up == TRUE) 
            gene.expr <- -gene.expr
        cor(gene.expr, seq_len(ncol(exprdata)), method = "spearman",
            use = "pairwise")
    }
    ret <- list()
    ret$rank <- rank
    ret$pat.order <- order(-rank)
    up.cor <- sapply(up, gene.cor, TRUE, exprdata, ret$pat.order)
    dn.cor <- sapply(dn, gene.cor, FALSE, exprdata, ret$pat.order)
    ret$gene.order <- c(up[rev(order(up.cor))], dn[order(dn.cor)])
    ret$up <- up
    ret$dn <- dn
    ret$dat <- exprdata[ret$gene.order, ret$pat.order, drop = FALSE]
    ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
    ret$ranksum <- ranksum
    return(ret)
}
