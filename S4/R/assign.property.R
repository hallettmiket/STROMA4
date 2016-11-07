#' Function to assign properties to an expression matrix
#' 
#' This is a function to assign stromal property and TNBCType generative property levels to a TNBC dataset
#' @param exprs An expression matrix. Rows correspont to genes, columns to samples. Genes with multiple probes should be collapsed to a single row e.g. using most variable probe.
#' @param genes A vector with gene symbols for exprs. Should be the same length as the number of rows in exprs.
#' @param genelists A vector with either S4, TNBCType, or both to specify which genelists to use.
#' @param n An integer value specifying the number of random samples to generate.
#' @param seed An integer value specifying a seed for the random ranks function. Defauly value is 123456.
#' @param mc.cores An integer specifying how many cores to use. Default value is 1.
#' @export
assign.properties <- function (exprs, genes, genelists = c("S4", "TNBCType"), n = 1000, 
    seed = 123456, mc.cores = 1) 
{
    cat("--Assigning properties to expression data--", "\n")
    if (!class(exprs) %in% c("data.frame", "matrix")) 
        stop("Error in assigning property: exprs is not a matrix or dataframe")
    if (length(genes) != nrow(exprs)) 
        stop("Error in assigning property: Length of genes vector not equal to number of rows of expression data")
    if (any(duplicated(genes))) {
        cat("--There are duplicated genes. Using most variable (IQR) to collapse--", 
            "\n")
        cat()
        iqrs <- order(apply(exprs, 1, IQR), decreasing = T)
        exprs <- exprs[iqrs, ]
        genes <- genes[iqrs]
        to.keep <- !duplicated(genes)
        genes <- genes[to.keep]
        exprs <- exprs[to.keep, ]
    }
    temp.envir <- new.env()
    if ("S4" %in% genelists) {
        data(B.stroma.property, package = "S4", envir = temp.envir)
        data(E.stroma.property, package = "S4", envir = temp.envir)
        data(F.stroma.property, package = "S4", envir = temp.envir)
        data(T.stroma.property, package = "S4", envir = temp.envir)
    }
    if ("TNBCType" %in% genelists) {
        data(BL1.property, package = "S4", envir = temp.envir)
        data(BL2.property, package = "S4", envir = temp.envir)
        data(IM.property, package = "S4", envir = temp.envir)
        data(LAR.property, package = "S4", envir = temp.envir)
        data(M.property, package = "S4", envir = temp.envir)
        data(MSL.property, package = "S4", envir = temp.envir)
    }
    if (!("TNBCType" %in% genelists) & !("S4" %in% genelists)) 
        stop("Need to specify either S4 or TNBCType as genelists")
    ret <- list()
    for (i in names(temp.envir)) {
        match.genes <- temp.envir[[i]]
        match.genes <- match.genes[which(match.genes[, 1] %in% 
            genes), , drop = F]
        if (nrow(match.genes) == 0) {
            ret[[i]] <- "Error: No matching genes in expression matrix"
        }
        else {
            cat("----", nrow(match.genes), "out of", nrow(temp.envir[[i]]), 
                "total genes matching for", i, "----", "\n")
            up.genes <- which(genes %in% match.genes[which(match.genes[, 
                2] == "up"), 1])
            dn.genes <- which(genes %in% match.genes[which(match.genes[, 
                2] == "down"), 1])
            match.exprs <- exprs[c(up.genes, dn.genes), ]
            directions <- rep(c("up", "down"), c(length(up.genes), 
                length(dn.genes)))
            ranksum.object <- sig.ranksum(exprdata = match.exprs, 
                up = which(directions == "up"), dn = which(directions == 
                  "down"), full.return = T)
            roi <- random.ranks(ranksum.object, n = n, seed = seed, 
                mc.cores = mc.cores)
            ret[[i]] <- define.roi.regions(ranksum.object, roi)
        }
    }
    rm(temp.envir)
    return(ret)
}
