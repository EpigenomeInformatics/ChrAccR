# ==============================================================================
# CLASS DEFINITION
# ==============================================================================
#' @include DsAcc-class.R
#' 
#' DsASC Class
#'
#' A class for storing Allele-Specific Chromatin (ASC) data
#' inherits from \code{\linkS4class{DsAcc}}
#' @author Irem B. GUNDUZ
#' @section Slots:
#' \describe{
#'   \item{\code{counts}}{
#'      List containing two matrices: \code{'ref'} (Reference allele counts) and \code{'alt'} (Alternative allele counts).
#'   }
#'   \item{\code{diskDump}}{
#'      Flag indicating whether count data is stored on disk (HDF5) rather than in main memory.
#'   }
#' }
#' @name DsASC-class
#' @rdname DsASC-class
#' @exportClass DsASC
setClass("DsASC",
    slots = list(
        counts = "list",
        diskDump = "logical"
    ),
    contains = "DsAcc",
    package = "ChrAccR"
)

# ==============================================================================
# Initialize
# ==============================================================================

setMethod("initialize", "DsASC",
    function(
        .Object,
        counts,
        coord,
        sampleAnnot,
        genome,
        diskDump
    ) {
        .Object@counts      <- counts
        .Object@coord       <- coord
        .Object@sampleAnnot <- sampleAnnot
        .Object@genome      <- genome
        .Object@diskDump    <- diskDump
        .Object@pkgVersion  <- packageVersion("ChrAccR")
        .Object
    }
)

# ==============================================================================
# Construcutor
# ==============================================================================

#' Create a new DsASC object
#' @author Irem B. GUNDUZ
#' @param sampleAnnot data.frame with sample annotation.
#' @param genome      Character string containing genome assembly (e.g. "hg38").
#' @param diskDump    Logical. If TRUE, matrices are realized as HDF5 arrays.
#' 
#' @return A \code{\linkS4class{DsASC}} object.
#' @export
DsASC <- function(sampleAnnot, genome, diskDump=FALSE) {
    obj <- new("DsASC",
        list(),       # counts
        list(),       # coord
        sampleAnnot,
        genome,
        diskDump
    )
    return(obj)
}

# ==============================================================================
# Show 
# ==============================================================================

setMethod("show", "DsASC", function(object) {
    cat("DsASC allele-specific chromatin accessibility dataset\n")
    cat("contains:\n")
    cat(sprintf(" * %d samples\n", nrow(object@sampleAnnot)))
    
    if (length(object@coord) > 0) {
        nSites <- length(object@coord[[1]])
        cat(sprintf(" * %d heterozygous SNP sites\n", nSites))
    } else {
        cat(" * no sites defined (empty)\n")
    }
    
    if (length(object@counts) > 0) {
        cat(" * assays: ref, alt\n")
        if (object@diskDump) cat(" * data stored on disk (HDF5)\n")
    }
})

# ==============================================================================
# Getters
# ==============================================================================

if (!isGeneric("getCounts")) {
    setGeneric("getCounts", function(.object, ...) standardGeneric("getCounts"), signature=c(".object"))
}

#' @export
setMethod("getCounts", signature(.object="DsASC"),
    function(.object, type, i=NULL, j=NULL, asMatrix=TRUE, naIsZero=TRUE) {
        if (!type %in% c("ref", "alt")) stop(paste("Unsupported allele type:", type))
        
        res <- .object@counts[[type]]

        if (!is.null(i)) res <- res[i,,drop=FALSE]
        if (!is.null(j)) res <- res[,j,drop=FALSE]
        
        if (asMatrix && !is.matrix(res)){
            if (.object@diskDump) {
                res <- as.matrix(res)
            }
        }

        if (naIsZero){
            res[is.na(res)] <- 0
        }
        return(res)
    }
)
# ==============================================================================
# Ref Counts Getter
# ==============================================================================
#' @export
setGeneric("getRefCounts", function(object, ...) standardGeneric("getRefCounts"))

#' @export
setMethod("getRefCounts", "DsASC", 
    function(object, naIsZero = TRUE, ...) {
        # Explicitly pass naIsZero to the internal getCounts
        getCounts(object, type = "ref", naIsZero = naIsZero, ...)
    }
)

# ==============================================================================
# Alt Counts Getter
# ==============================================================================
#' @export
setGeneric("getAltCounts", function(object, ...) standardGeneric("getAltCounts"))

#' @export
setMethod("getAltCounts", "DsASC", 
    function(object, naIsZero = TRUE, ...) {
        # Explicitly pass naIsZero to the internal getCounts
        getCounts(object, type = "alt", naIsZero = naIsZero, ...)
    }
)
#' @export
setGeneric("getAllelicBalance", function(object, minCoverage=0, ...) standardGeneric("getAllelicBalance"))
setMethod("getAllelicBalance", "DsASC", function(object, minCoverage=0, ...) {
    alt <- getAltCounts(object, ...)
    total <- getRefCounts(object, ...) + getAltCounts(object, ...)
    ratio <- alt / total
    
    if (minCoverage > 0) {
        ratio[total < minCoverage] <- NA
    }
    ratio[total == 0] <- NA
    
    return(ratio)
})

#' @export
setMethod("[", signature(x = "DsASC", i = "ANY", j = "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        if (missing(i)) i <- TRUE
        if (missing(j)) j <- TRUE

        # Subset counts
        x@counts <- lapply(x@counts, function(mat) {
             return(mat[i, j, drop = FALSE])
        })

        # Subset sample annotation
        x@sampleAnnot <- x@sampleAnnot[j, , drop = FALSE]

        # Subset coordinates (coord is a list of GRanges)
        if (.hasSlot(x, "coord")) {
            if (is.list(x@coord) || inherits(x@coord, "List")) {
                x@coord <- lapply(x@coord, function(gr) gr[i])
            } else {
                x@coord <- x@coord[i]
            }
        }

        return(x)
    }
)

# ==============================================================================
# Filtering
# ==============================================================================
#' @export
setGeneric("filterLowCovg", function(.object, ...) standardGeneric("filterLowCovg"))

#' filterLowCovg-methods (DsASC)
#' @param .object   A \code{\linkS4class{DsASC}} object.
#' @param thresh    Integer. Minimum total coverage required per site per sample.
#' @param reqSamples Numeric. Minimum number or fraction of samples that must meet the coverage threshold.
#' @author Irem B. GUNDUZ
#' @export
setMethod("filterLowCovg", signature(.object = "DsASC"),
  function(.object, thresh = 1L, reqSamples = 0.75) {
    
    # Determin e allowed failing samples
    N <- length(getSamples(.object))
    numAllowed <- reqSamples
    if (numAllowed < 1 && numAllowed >= 0) numAllowed <- as.integer(ceiling(numAllowed * N))
    percAllowed <- round(numAllowed / N, 2)
    
    logger.status(c("Removing sites with total coverage <", thresh, "in >", N - numAllowed, "samples"))

    # Get Coverage 
    ref <- getRefCounts(.object, naIsZero = TRUE)
    alt <- getAltCounts(.object, naIsZero = TRUE)
    totalCov <- ref + alt

    # Identify sites to remove
    passing_samples <- rowSums(totalCov >= thresh)
    rem <- passing_samples < numAllowed
    nRem <- sum(rem)
    nSites <- nrow(ref)

    # Filter the object
    if (nRem > 0) {
      .object <- .object[!rem, ]
    }
    
    logger.info(c("Removed", nRem, "sites", paste0("(", round(nRem / nSites, 4) * 100, "%)")))
    return(.object)
  }
)

#' @export
setGeneric("filterChroms", function(.object, ...) standardGeneric("filterChroms"))

#' filterChroms-methods (DsASC)
#'
#' Filter the dataset to retain or exclude specific chromosomes.
#'
#' @param .object   \code{\linkS4class{DsASC}} object
#' @param keep      Character vector of chromosomes to keep. If NULL (default), all are kept unless 'remove' is specified.
#'                  Special value "standard" keeps chr1-22, X, Y (human) or chr1-19, X, Y.
#' @param remove    Character vector of chromosomes to remove (e.g., "chrM").
#' @return          A new \code{\linkS4class{DsASC}} object with the specified chromosomes filtered.
#'
#' @rdname filterChroms-DsASC-method
#' @docType methods
#' @aliases filterChroms
#' @export
setMethod("filterChroms", signature(.object = "DsASC"),
    function(.object, keep = NULL, remove = NULL) {
        
        if (.hasSlot(.object, "coord")) {
            gr <- unlist(.object@coord)$snps
        } else {
            stop("Could not find coordinate slot (loci/gr/coord) in DsASC object.")
        }

        current_chroms <- as.character(GenomeInfoDb::seqnames(gr))
        unique_chroms <- unique(current_chroms)
        to_keep <- unique_chroms

        if (!is.null(keep)) {
            if (length(keep) == 1 && keep == "standard") {
                std_chrs <- paste0("chr", c(1:22, "X", "Y")) 
                std_nums <- c(1:22, "X", "Y")
                to_keep <- intersect(unique_chroms, c(std_chrs, std_nums))
            } else {
                to_keep <- intersect(unique_chroms, keep)
            }
        }

        if (!is.null(remove)) {
            to_keep <- setdiff(to_keep, remove)
        }

        keep_idx <- which(current_chroms %in% to_keep)
        
        n_before <- length(current_chroms)
        n_after <- length(keep_idx)
        n_removed <- n_before - n_after

        if (n_after == 0) {
            logger.warning("Filtering would remove ALL sites. Returning empty object.")
        }

        .object <- .object[keep_idx, ]

        logger.info(c("Kept", length(to_keep), "chromosomes. Removed", n_removed, "sites", 
                      paste0("(", round(n_removed/n_before*100, 2), "%)")))

        return(.object)
    }
)
