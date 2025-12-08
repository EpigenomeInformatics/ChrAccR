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

setMethod(
  "initialize", "DsASC",
  function(.Object,
           counts,
           coord,
           sampleAnnot,
           genome,
           diskDump) {
    .Object@counts <- counts
    .Object@coord <- coord
    .Object@sampleAnnot <- sampleAnnot
    .Object@genome <- genome
    .Object@diskDump <- diskDump
    .Object@pkgVersion <- packageVersion("ChrAccR")
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
DsASC <- function(sampleAnnot, genome, diskDump = FALSE) {
  obj <- new(
    "DsASC",
    list(), # counts
    list(), # coord
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
  setGeneric("getCounts", function(.object, ...) standardGeneric("getCounts"), signature = c(".object"))
}

#' @export
setMethod(
  "getCounts", signature(.object = "DsASC"),
  function(.object, type, i = NULL, j = NULL, asMatrix = TRUE, naIsZero = TRUE) {
    if (!type %in% c("ref", "alt")) stop(paste("Unsupported allele type:", type))

    res <- .object@counts[[type]]

    if (!is.null(i)) res <- res[i, , drop = FALSE]
    if (!is.null(j)) res <- res[, j, drop = FALSE]

    if (asMatrix && !is.matrix(res)) {
      if (.object@diskDump) {
        res <- as.matrix(res)
      }
    }

    if (naIsZero) {
      res[is.na(res)] <- 0
    }
    return(res)
  }
)

#' @export
setGeneric("getRefCounts", function(object, ...) standardGeneric("getRefCounts"))

#' @export
setMethod(
  "getRefCounts", "DsASC",
  function(object, naIsZero = TRUE, ...) {
    # Explicitly pass naIsZero to the internal getCounts
    getCounts(object, type = "ref", naIsZero = naIsZero, ...)
  }
)

#' @export
setGeneric("getAltCounts", function(object, ...) standardGeneric("getAltCounts"))

#' @export
setMethod(
  "getAltCounts", "DsASC",
  function(object, naIsZero = TRUE, ...) {
    # Explicitly pass naIsZero to the internal getCounts
    getCounts(object, type = "alt", naIsZero = naIsZero, ...)
  }
)
#' @export
setGeneric("getAllelicBalance", function(object, minCoverage = 0, ...) standardGeneric("getAllelicBalance"))
setMethod("getAllelicBalance", "DsASC", function(object, minCoverage = 0, ...) {
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
setMethod(
  "[", signature(x = "DsASC", i = "ANY", j = "ANY"),
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
setMethod(
  "filterLowCovg", signature(.object = "DsASC"),
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
setMethod(
  "filterChroms", signature(.object = "DsASC"),
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

    logger.info(c(
      "Kept", length(to_keep), "chromosomes. Removed", n_removed, "sites",
      paste0("(", round(n_removed / n_before * 100, 2), "%)")
    ))

    return(.object)
  }
)

#' Annotate DsASC SNPs with Consensus Peaks
#'
#' Overlaps the SNPs in a DsASC object with a set of Consensus Peaks (GRanges).
#' Assigns the Peak ID to the SNP metadata.
#'
#' @param dsObj A DsASC object.
#' @param peakGr A GRanges object containing consensus peaks. Must have unique names().
#' @return A DsASC object with updated @coord metadata.
#' @export
annotateASCSites <- function(dsObj, peakGr) {
  logger.start("Annotating SNPs with Consensus Peaks")

  snpGr <- dsObj@coord$snps

  # Find Overlaps
  ov <- GenomicRanges::findOverlaps(snpGr, peakGr)

  # Initialize peakId column as NA
  GenomicRanges::mcols(snpGr)$peakId <- NA_character_

  # If a SNP overlaps multiple peaks, this takes the first one.
  GenomicRanges::mcols(snpGr)$peakId[S4Vectors::queryHits(ov)] <- names(peakGr)[S4Vectors::subjectHits(ov)]

  # Update Object
  dsObj@coord$snps <- snpGr

  n_overlap <- sum(!is.na(GenomicRanges::mcols(snpGr)$peakId))
  pct <- round(n_overlap / length(snpGr) * 100, 2)

  logger.info(paste0(n_overlap, " SNPs (", pct, "%) overlap with consensus peaks."))
  logger.completed()

  return(dsObj)
}

#' Filter DsASC to Peak-Overlapping SNPs
#'
#' Removes SNPs that do not overlap with any consensus peak.
#' This reduces dataset size and focuses analysis on regulatory regions.
#'
#' @param dsObj A DsASC object (must be annotated first).
#' @return A subsetted DsASC object.
#' @export
filterASCByPeaks <- function(dsObj) {
  logger.start("Filtering SNPs outside peaks")

  snpGr <- dsObj@coord$snps

  if (!"peakId" %in% colnames(GenomicRanges::mcols(snpGr))) {
    logger.error("DsASC object not annotated. Run annotateASCSites() first.")
    stop("Missing peakId annotation")
  }

  # Identify indices to keep
  keepIdx <- which(!is.na(GenomicRanges::mcols(snpGr)$peakId))

  if (length(keepIdx) == 0) {
    logger.warning("No SNPs inside peaks found!")
    return(dsObj)
  }

  # Subset Coordinates
  dsObj@coord$snps <- snpGr[keepIdx]

  # Subset Count Matrices
  dsObj@counts$ref <- dsObj@counts$ref[keepIdx, , drop = FALSE]
  dsObj@counts$alt <- dsObj@counts$alt[keepIdx, , drop = FALSE]

  logger.info(paste("Retained", length(keepIdx), "SNPs inside peaks."))
  logger.completed()

  return(dsObj)
}

# ==============================================================================
# Downstream Analysis
# ==============================================================================

#' Calculate ASC Statistics (Binomial Test)
#'
#' Performs a two-sided Binomial Test against p=0.5 for every SNP.
#' Adds P-values and FDR (Benjamini-Hochberg) to the result.
#'
#' @param dsObj A DsASC object (filtered).
#' @param minCoverage Minimum coverage to perform test (default 10).
#' @return A data.table containing statistics for every SNP-Sample pair.
#' @export
calcASCStatistics <- function(dsObj, minCoverage = 10) {
  logger.start("Calculating ASC Statistics (Binomial Test)")

  #  Fetch Data
  refMat <- as.matrix(getRefCounts(dsObj))
  altMat <- as.matrix(getAltCounts(dsObj))
  totalMat <- refMat + altMat

  # Initialize Results Table
  dt <- data.table::as.data.table(as.table(totalMat))
  setnames(dt, c("snpId", "sampleId", "total"))

  # Add Ref/Alt counts
  dt$ref <- as.vector(refMat)
  dt$alt <- as.vector(altMat)

  # Filter for coverage
  dt <- dt[total >= minCoverage]

  if (nrow(dt) == 0) {
    logger.warning("No sites passed coverage threshold.")
    return(NULL)
  }

  # Calculate Ratios
  dt$ratio <- dt$alt / dt$total
  dt$log2FC <- log2((dt$alt + 1) / (dt$ref + 1))

  # Vectorized Binomial Test (Two-sided)
  # The probability of observing k successes in n trials with p=0.5
  # P-value = 2 * pbinom(min(k, n-k), n, 0.5)
  logger.info("Running Binomial Tests...")

  min_k <- pmin(dt$ref, dt$alt)
  dt$pVal <- 2 * pbinom(min_k, size = dt$total, prob = 0.5)

  # 6. FDR Correction (Per Sample)
  logger.info("Applying FDR correction...")
  dt[, fdr := p.adjust(pVal, method = "BH"), by = sampleId]

  # 7. Add Consensus Peak Info if available
  if (!is.null(dsObj@coord$snps$peakId)) {
    peak_map <- setNames(dsObj@coord$snps$peakId, names(dsObj@coord$snps))
    dt$peakId <- peak_map[dt$snpId]
  }

  logger.info(paste("Calculated stats for", nrow(dt), "site-sample pairs."))
  logger.completed()

  return(dt)
}
