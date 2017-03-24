#' Function to assign properties to an expression matrix
#'
#' This is a function to assign stromal property and TNBCType generative property levels to a TNBC dataset
#' @param ESet An ExpressionSet object. Rows correspond to genes, columns to samples. If there are genes with multiple probes, they will be collapsed to a single row using the function defined in var.method.
#' @param geneID.column Integer or column name corresponding to column in featureData table corresponding to HGNC ID.
#' @param genelists A vector with either Stroma4, TNBCType, or both to specify which genelists to use.
#' @param n An integer value specifying the number of random samples to generate.
#' @param seed An integer value specifying a seed for the random ranks function. Default value is 123456.
#' @param mc.cores An integer specifying how many cores to use. Defaults to use the function snowWorkers().
#' @param var.method Function for assessing variance to collapse probes. Default is IQR
#' @export
#' @import Biobase
#' @importFrom BiocParallel snowWorkers
#' @importFrom matrixStats rowIQRs
#' @importFrom utils data
#' @return The function returns a list of assignments for each property. Assignments are from 1-3, with 1 being low, 2 being intermediate, and 3 being high assignment respectively.
#' @examples
#' library(breastCancerMAINZ)
#' data(mainz, package='breastCancerMAINZ')
#' all.properties <- assign.properties(ESet=mainz, geneID.column='Gene.symbol',
#'		genelists=c('Stroma4', 'TNBCType'), n=10, mc.cores=1)
assign.properties <- function (ESet, geneID.column = 1, genelists = c("Stroma4", "TNBCType"),
    n = 1000, seed = 123456, mc.cores = snowWorkers(), var.method = function(x) rowIQRs(x,
        na.rm = TRUE))
{
    message("--Assigning properties to expression data--")
    if (!class(ESet) == "ExpressionSet")
        stop("Error in assigning property: exprs is not of class \"ExpressionSet\"")
    exprs <- exprs(ESet)
    genes <- fData(ESet)[, geneID.column]
    if (any(duplicated(genes))) {
        message("--There are duplicated genes. Using most variable to collapse--")
        var.estimate <- order(var.method(exprs), decreasing = TRUE)
        exprs <- exprs[var.estimate, ]
        genes <- genes[var.estimate]
        to.keep <- !duplicated(genes)
        genes <- genes[to.keep]
        exprs <- exprs[to.keep, ]
    }

    temp.envir <- new.env()
    if ("Stroma4" %in% genelists) {
        data("B.stroma.property", package = "STROMA4", envir = temp.envir)
        data("E.stroma.property", package = "STROMA4", envir = temp.envir)
        data("D.stroma.property", package = "STROMA4", envir = temp.envir)
        data("T.stroma.property", package = "STROMA4", envir = temp.envir)
    }

    if ("TNBCType" %in% genelists) {
        data("BL1.property", package = "STROMA4", envir = temp.envir)
        data("BL2.property", package = "STROMA4", envir = temp.envir)
        data("IM.property", package = "STROMA4", envir = temp.envir)
        data("LAR.property", package = "STROMA4", envir = temp.envir)
        data("M.property", package = "STROMA4", envir = temp.envir)
        data("MSL.property", package = "STROMA4", envir = temp.envir)
    }

    if (!("TNBCType" %in% genelists) & !("Stroma4" %in% genelists))
        stop("Need to specify either Stroma4 or TNBCType as genelists")

    ret <- vector("list", length(temp.envir))
    names(ret) <- names(temp.envir)

    for (i in names(temp.envir)) {
        match.genes <- temp.envir[[i]]
        match.genes <- match.genes[which(match.genes[, 1] %in%
            genes), , drop = FALSE]
        if (nrow(match.genes) == 0) {
            ret[[i]] <- "Error: No matching genes in expression matrix"
        }
        else {
            message("----", nrow(match.genes), " out of ", nrow(temp.envir[[i]]),
                " total genes matching for ", i, "----")
            up.genes <- which(genes %in% match.genes[which(match.genes[,
                2] == "up"), 1])
            dn.genes <- which(genes %in% match.genes[which(match.genes[,
                2] == "down"), 1])
            match.exprs <- exprs[c(up.genes, dn.genes), ]
            directions <- rep(c("up", "down"), c(length(up.genes),
                length(dn.genes)))
            ranksum.object <- .sig.ranksum(exprdata = match.exprs,
                up = which(directions == "up"), dn = which(directions ==
                  "down"), full.return = TRUE)
            roi <- .random.ranks(ranksum.object, n = n, seed = seed,
                workers = mc.cores)
            ret[[i]] <- .define.roi.regions(ranksum.object, roi)
        }
    }
    rm(temp.envir)

    df <- as.data.frame(Map(factor, ret, MoreArgs=list(labels=c("low", "intermediate", "high"))))

    pData(ESet) <- cbind(pData(ESet), df)

    return(ESet)
}
