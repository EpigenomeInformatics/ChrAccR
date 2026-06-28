# ==============================================================================
# CLASS DEFINITION
# ==============================================================================
#' @include DsAcc-class.R
#'
#' DsASC Class
#'
#' A class for storing Allele-Specific Chromatin (ASC) data together with the
#' underlying peak-level chromatin accessibility signal.
#' Inherits from \code{\linkS4class{DsAcc}}.
#'
#' @author Irem B. GUNDUZ
#' @section Slots:
#' \describe{
#'   \item{\code{counts}}{
#'      List of two matrices: \code{'ref'} (reference-allele read counts) and
#'      \code{'alt'} (alternative-allele read counts), each [SNP x sample].
#'   }
#'   \item{\code{accessibility}}{
#'      Matrix of peak-level ATAC-seq counts [peak x sample]. This is the
#'      "total" chromatin accessibility signal, independent of allele.
#'      May be a DelayedArray/HDF5Matrix when \code{diskDump = TRUE}.
#'   }
#'   \item{\code{diskDump}}{
#'      Flag indicating whether count data is stored on disk (HDF5) rather than
#'      in main memory.
#'   }
#' }
#' The inherited \code{coord} slot is a list holding two GRanges:
#'   \code{coord$snps}  (heterozygous SNP positions, annotated with peakId) and
#'   \code{coord$peaks} (consensus peak coordinates).
#'
#' @name DsASC-class
#' @rdname DsASC-class
#' @exportClass DsASC
setClass("DsASC",
  slots = list(
    counts        = "list",
    accessibility = "ANY",     # peak x sample matrix (or HDF5 when diskDump)
    diskDump      = "logical"
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
           accessibility,
           coord,
           sampleAnnot,
           genome,
           diskDump) {
    .Object@counts        <- counts
    .Object@accessibility <- accessibility
    .Object@coord         <- coord
    .Object@sampleAnnot   <- sampleAnnot
    .Object@genome        <- genome
    .Object@diskDump      <- diskDump
    .Object@pkgVersion    <- packageVersion("ChrAccR")
    .Object
  }
)

# ==============================================================================
# Constructor
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
    list(),                                   # counts (ref / alt)
    matrix(numeric(0), nrow = 0, ncol = 0),   # accessibility (empty peak matrix)
    list(),                                   # coord (snps / peaks)
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

  # SNP sites
  if (length(object@coord) > 0 && !is.null(object@coord$snps)) {
    cat(sprintf(" * %d heterozygous SNP sites\n", length(object@coord$snps)))
  } else {
    cat(" * no SNP sites defined (empty)\n")
  }

  # Peaks / accessibility
  if (!is.null(object@coord$peaks)) {
    cat(sprintf(" * %d consensus peaks\n", length(object@coord$peaks)))
  }
  if (length(object@accessibility) > 0 && nrow(object@accessibility) > 0) {
    cat(sprintf(" * accessibility matrix: %d peaks x %d samples\n",
                nrow(object@accessibility), ncol(object@accessibility)))
  } else {
    cat(" * accessibility: not loaded\n")
  }

  # Allele assays
  if (length(object@counts) > 0) {
    cat(" * allele assays: ref, alt\n")
    if (object@diskDump) cat(" * data stored on disk (HDF5)\n")
  }
})

# ==============================================================================
# Getters  (allele counts)
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
      if (.object@diskDump) res <- as.matrix(res)
    }

    if (naIsZero) res[is.na(res)] <- 0
    return(res)
  }
)

#' @export
setGeneric("getRefCounts", function(object, ...) standardGeneric("getRefCounts"))
#' @export
setMethod(
  "getRefCounts", "DsASC",
  function(object, naIsZero = TRUE, ...) {
    getCounts(object, type = "ref", naIsZero = naIsZero, ...)
  }
)

#' @export
setGeneric("getAltCounts", function(object, ...) standardGeneric("getAltCounts"))
#' @export
setMethod(
  "getAltCounts", "DsASC",
  function(object, naIsZero = TRUE, ...) {
    getCounts(object, type = "alt", naIsZero = naIsZero, ...)
  }
)

#' @export
setGeneric("getAllelicBalance", function(object, minCoverage = 0, ...) standardGeneric("getAllelicBalance"))
setMethod("getAllelicBalance", "DsASC", function(object, minCoverage = 0, ...) {
  alt   <- getAltCounts(object, ...)
  total <- getRefCounts(object, ...) + getAltCounts(object, ...)
  ratio <- alt / total

  if (minCoverage > 0) ratio[total < minCoverage] <- NA
  ratio[total == 0] <- NA

  return(ratio)
})

# ==============================================================================
# Getters / setter  (accessibility)
# ==============================================================================

#' Get the peak-level accessibility matrix
#'
#' Returns the total (allele-agnostic) ATAC-seq count matrix [peak x sample].
#' This is the chromatin accessibility signal used for differential
#' accessibility and for filtering SNPs to open regions.
#'
#' @param object   A \code{\linkS4class{DsASC}} object.
#' @param asMatrix Logical. Realize HDF5-backed arrays into a dense matrix.
#' @return A matrix of peak accessibility counts.
#' @author Irem B. GUNDUZ
#' @export
setGeneric("getAccessibility", function(object, ...) standardGeneric("getAccessibility"))
setMethod("getAccessibility", "DsASC", function(object, asMatrix = TRUE, ...) {
  res <- object@accessibility
  if (asMatrix && !is.matrix(res) && object@diskDump) res <- as.matrix(res)
  return(res)
})

#' Get peak coordinates
#' @export
setGeneric("getPeaks", function(object, ...) standardGeneric("getPeaks"))
setMethod("getPeaks", "DsASC", function(object, ...) object@coord$peaks)

#' Attach a peak accessibility matrix (+ peak coordinates) to a DsASC object
#'
#' Populates the \code{accessibility} slot and \code{coord$peaks}. The peak
#' matrix is typically produced by ChrAccR / your 02_ChrAccR.R step.
#'
#' @param dsObj    A \code{\linkS4class{DsASC}} object.
#' @param peakMat  Matrix of peak counts [peak x sample]. rownames = peak IDs,
#'                 colnames must contain every \code{getSamples(dsObj)}.
#' @param peakGr   GRanges of consensus peaks. names(peakGr) must match
#'                 rownames(peakMat).
#' @return The updated DsASC object.
#' @author Irem B. GUNDUZ
#' @export
setAccessibility <- function(dsObj, peakMat, peakGr) {
  logger.start("Attaching accessibility matrix to DsASC")

  # --- identifier checks ---
  if (is.null(rownames(peakMat)) || is.null(names(peakGr))) {
    logger.error("peakMat must have rownames and peakGr must have names().")
    stop("Missing peak identifiers.")
  }
  if (!identical(rownames(peakMat), names(peakGr))) {
    logger.warning("Reordering peakGr to match peakMat rownames.")
    peakGr <- peakGr[rownames(peakMat)]
  }

  # --- sample alignment ---
  sampleIds <- getSamples(dsObj)
  missing   <- setdiff(sampleIds, colnames(peakMat))
  if (length(missing) > 0) {
    logger.error(paste("peakMat is missing columns for samples:",
                       paste(missing, collapse = ", ")))
    stop("Sample mismatch between DsASC and peakMat.")
  }
  peakMat <- peakMat[, sampleIds, drop = FALSE]   # align column order

  dsObj@accessibility <- peakMat
  dsObj@coord$peaks    <- peakGr

  logger.info(paste0("Stored ", nrow(peakMat), " peaks x ", ncol(peakMat), " samples."))
  logger.completed()
  return(dsObj)
}

# ==============================================================================
# Subsetting   (FIXED: peaks are not indexed by SNP rows; accessibility cols
#               follow the sample index j)
# ==============================================================================

#' @export
setMethod(
  "[", signature(x = "DsASC", i = "ANY", j = "ANY"),
  function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    # Allele counts: SNP rows (i) x sample cols (j)
    x@counts <- lapply(x@counts, function(mat) mat[i, j, drop = FALSE])

    # Accessibility: peak rows x sample cols -> only sample cols move with j
    if (length(x@accessibility) > 0 && nrow(x@accessibility) > 0) {
      x@accessibility <- x@accessibility[, j, drop = FALSE]
    }

    # Sample annotation follows j
    x@sampleAnnot <- x@sampleAnnot[j, , drop = FALSE]

    # Coordinates: ONLY snps are indexed by i; peaks are independent and stay put
    if (!is.null(x@coord$snps)) {
      x@coord$snps <- x@coord$snps[i]
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
#' @param .object    A \code{\linkS4class{DsASC}} object.
#' @param thresh     Integer. Minimum total coverage required per site per sample.
#' @param reqSamples Numeric. Minimum number (>=1) or fraction (<1) of samples
#'                   that must meet the coverage threshold.
#' @author Irem B. GUNDUZ
#' @export
setMethod(
  "filterLowCovg", signature(.object = "DsASC"),
  function(.object, thresh = 1L, reqSamples = 0.75) {
    N <- length(getSamples(.object))
    numAllowed <- reqSamples
    if (numAllowed < 1 && numAllowed >= 0) numAllowed <- as.integer(ceiling(numAllowed * N))
    percAllowed <- round(numAllowed / N, 2)

    logger.status(c("Removing sites with total coverage <", thresh, "in >", N - numAllowed, "samples"))

    ref <- getRefCounts(.object, naIsZero = TRUE)
    alt <- getAltCounts(.object, naIsZero = TRUE)
    totalCov <- ref + alt

    passing_samples <- rowSums(totalCov >= thresh)
    rem <- passing_samples < numAllowed
    nRem <- sum(rem)
    nSites <- nrow(ref)

    if (nRem > 0) .object <- .object[!rem, ]

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
#'                  Special value "standard" keeps chr1-22, X, Y.
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
    # FIXED: read coord$snps directly (unlist() merges snps+peaks and breaks $snps)
    gr <- .object@coord$snps
    if (is.null(gr)) stop("Could not find coord$snps in DsASC object.")

    current_chroms <- as.character(GenomeInfoDb::seqnames(gr))
    unique_chroms  <- unique(current_chroms)
    to_keep        <- unique_chroms

    if (!is.null(keep)) {
      if (length(keep) == 1 && keep == "standard") {
        std_chrs <- paste0("chr", c(1:22, "X", "Y"))
        std_nums <- c(1:22, "X", "Y")
        to_keep  <- intersect(unique_chroms, c(std_chrs, std_nums))
      } else {
        to_keep <- intersect(unique_chroms, keep)
      }
    }

    if (!is.null(remove)) to_keep <- setdiff(to_keep, remove)

    keep_idx <- which(current_chroms %in% to_keep)

    n_before  <- length(current_chroms)
    n_after   <- length(keep_idx)
    n_removed <- n_before - n_after

    if (n_after == 0) logger.warning("Filtering would remove ALL sites. Returning empty object.")

    .object <- .object[keep_idx, ]

    logger.info(c(
      "Kept", length(to_keep), "chromosomes. Removed", n_removed, "sites",
      paste0("(", round(n_removed / n_before * 100, 2), "%)")
    ))

    return(.object)
  }
)

# ==============================================================================
# Peak annotation & overlap filtering
# ==============================================================================

#' Annotate DsASC SNPs with Consensus Peaks
#'
#' Overlaps the SNPs in a DsASC object with a set of Consensus Peaks (GRanges).
#' Assigns the Peak ID to the SNP metadata.
#'
#' @param dsObj  A DsASC object.
#' @param peakGr A GRanges object containing consensus peaks. Must have unique names().
#'               If NULL, uses the peaks already stored in \code{coord$peaks}.
#' @return A DsASC object with updated @coord metadata.
#' @author Irem B. GUNDUZ
#' @export
annotateASCSites <- function(dsObj, peakGr = NULL) {
  logger.start("Annotating SNPs with Consensus Peaks")

  if (is.null(peakGr)) peakGr <- dsObj@coord$peaks
  if (is.null(peakGr)) {
    logger.error("No peaks supplied and coord$peaks is empty.")
    stop("Missing peaks for annotation.")
  }
  if (is.null(names(peakGr))) {
    logger.error("peakGr has no names(); peakId assignment would be all NA.")
    stop("peakGr must have names().")
  }

  snpGr <- dsObj@coord$snps

  ov <- GenomicRanges::findOverlaps(snpGr, peakGr)

  GenomicRanges::mcols(snpGr)$peakId <- NA_character_
  # If a SNP overlaps multiple peaks, this takes the first one.
  GenomicRanges::mcols(snpGr)$peakId[S4Vectors::queryHits(ov)] <-
    names(peakGr)[S4Vectors::subjectHits(ov)]

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
#' @author Irem B. GUNDUZ
#' @export
filterASCByPeaks <- function(dsObj) {
  logger.start("Filtering SNPs outside peaks")

  snpGr <- dsObj@coord$snps

  if (!"peakId" %in% colnames(GenomicRanges::mcols(snpGr))) {
    logger.error("DsASC object not annotated. Run annotateASCSites() first.")
    stop("Missing peakId annotation")
  }

  keepIdx <- which(!is.na(GenomicRanges::mcols(snpGr)$peakId))

  if (length(keepIdx) == 0) {
    logger.warning("No SNPs inside peaks found!")
    return(dsObj)
  }

  # Subset via the [ method so counts + coord$snps stay in sync
  dsObj <- dsObj[keepIdx, ]

  logger.info(paste("Retained", length(keepIdx), "SNPs inside peaks."))
  logger.completed()

  return(dsObj)
}

# ==============================================================================
# Downstream Analysis
# ==============================================================================

#' Calculate ASC Statistics (Conditional Bias Correction)
#'
#' Performs a two-sided Binomial Test. If global residual bias is detected,
#' the null hypothesis is automatically adjusted to correct for the bias.
#'
#' @param dsObj A DsASC object (filtered).
#' @param minCoverage Minimum coverage to perform test (default 10).
#' @param globalProbTolerance Tolerance for deviation from p=0.5 (default 0.02).
#' @return A data.table containing statistics, including log2FC_norm which is bias-corrected only if needed.
#' @export
calcASCStatistics <- function(dsObj, minCoverage = 10, globalProbTolerance = 0.02) {

  logger.start("Calculating ASC Statistics (Conditional Correction)")

  refMat   <- as.matrix(getRefCounts(dsObj))
  altMat   <- as.matrix(getAltCounts(dsObj))
  totalMat <- refMat + altMat

  dt <- data.table::as.data.table(as.table(totalMat))
  data.table::setnames(dt, c("snpId", "sampleId", "total"))
  dt$ref <- as.vector(refMat)
  dt$alt <- as.vector(altMat)

  dt <- dt[dt$total >= minCoverage]

  if (nrow(dt) == 0) {
    logger.warning("No sites passed coverage threshold.")
    return(NULL)
  }

  # Conditional Bias Check
  global_total    <- sum(as.numeric(dt$total))
  global_alt_prob <- sum(as.numeric(dt$alt)) / global_total

  if (abs(global_alt_prob - 0.5) > globalProbTolerance) {
    test_prob     <- global_alt_prob
    is_normalized <- TRUE
    logger.warning(paste0("Bias Detected! Global Alt Prob is ", round(global_alt_prob, 3),
                          ". Test probability set to observed bias."))
  } else {
    test_prob     <- 0.5
    is_normalized <- FALSE
    logger.info("Bias not significant. Test probability set to 0.5.")
  }

  # LFC Metrics
  dt$log2FC_raw <- log2((dt$alt + 1) / (dt$ref + 1))

  if (is_normalized) {
    global_odds      <- (1 - global_alt_prob) / global_alt_prob
    dt$log2FC_report <- log2((dt$alt / dt$ref) / global_odds)
  } else {
    dt$log2FC_report <- dt$log2FC_raw
  }

  # Binomial Test
  logger.info(paste("Running Binomial Test against p =", round(test_prob, 3)))
  p_lower <- pbinom(dt$alt,     size = dt$total, prob = test_prob)
  p_upper <- pbinom(dt$alt - 1, size = dt$total, prob = test_prob, lower.tail = FALSE)

  dt$pVal <- 2 * pmin(p_lower, p_upper)
  dt$pVal[dt$pVal > 1] <- 1.0

  # FDR per sample
  dt[, fdr := p.adjust(pVal, method = "BH"), by = sampleId]

  # Peak metadata
  if (!is.null(dsObj@coord$snps$peakId)) {
    peak_map   <- setNames(dsObj@coord$snps$peakId, names(dsObj@coord$snps))
    dt$peakId  <- peak_map[dt$snpId]
  }

  data.table::setnames(dt, "log2FC_report", "log2FC_norm")

  logger.completed()
  return(dt)
}


#' Estimate the Proportion of Shared ASC Imbalance Effects
#'
#' Estimates the proportion of ASC sites (significant in Cell Type A) that still
#' show a non-zero effect in Cell Type B, using the Beta-Binomial posterior and ashR.
#'
#' @param dsObjA DsASC object for Cell Type A.
#' @param dsObjB DsASC object for Cell Type B.
#' @param fdrCutoff Numeric. FDR threshold for defining significance in Cell Type A (Default 0.01).
#' @return A list with the ashR model output and the estimated proportion of shared effects.
#' @export
estimateSharedImbalance <- function(dsObjA, dsObjB, fdrCutoff = 0.01) {

  logger.start("Estimating Shared Imbalance (ashR Method)")

  if (!requireNamespace("ashr", quietly = TRUE)) {
    logger.error("The 'ashr' package is required for this analysis.")
    stop("Missing package: ashr")
  }

  calculate_posterior <- function(dsObj) {
    ref   <- rowSums(as.matrix(getRefCounts(dsObj)), na.rm = TRUE)
    alt   <- rowSums(as.matrix(getAltCounts(dsObj)), na.rm = TRUE)
    total <- ref + alt

    valid_sites <- total > 0
    ref <- ref[valid_sites]; alt <- alt[valid_sites]; total <- total[valid_sites]

    # Posterior Beta(r+1, a+1)
    mu     <- (ref + 1) / (total + 2)
    sigma2 <- ((ref + 1) * (alt + 1)) / ((total + 2)^2 * (total + 3))

    pVal <- 2 * pmin(pbinom(ref, total, 0.5), pbinom(alt, total, 0.5))
    fdr  <- p.adjust(pVal, method = "BH")

    return(data.table(snpId = names(ref), mu = mu, sigma2 = sigma2, fdr = fdr))
  }

  dtA <- calculate_posterior(dsObjA)
  dtB <- calculate_posterior(dsObjB)

  sigA <- dtA[fdr < fdrCutoff]
  logger.info(paste("Identified", nrow(sigA), "significant ASC sites in Cell Type A (FDR <", fdrCutoff, ")."))

  if (nrow(sigA) < 100) logger.warning("Fewer than 100 significant sites found. ashR power may be low.")

  merged <- merge(sigA[, .(snpId)], dtB, by = "snpId", suffixes = c("_A", "_B"))

  merged[, effect_size := mu - 0.5]
  merged[, sem := sqrt(sigma2)]
  merged <- merged[sem > 0]

  logger.info("Running ashR to estimate shared proportion...")
  ash_res <- ashr::ash(
    Bhat       = merged$effect_size,
    SEbetahat  = merged$sem,
    mixcompdist = "normal"
  )

  pi_nonzero <- 1 - ash_res$pi[1]

  logger.info(paste0("Estimated Proportion of Shared Effects (Non-Zero): ", round(pi_nonzero, 3)))
  logger.completed()

  return(list(
    ash_model        = ash_res,
    sharing_estimate = pi_nonzero,
    sites_tested     = nrow(merged)
  ))
}


#' Filter for Recurrent Allele-Specific Chromatin Events
#'
#' Identifies SNPs that are significantly imbalanced (ASC) in a minimum number
#' of independent donors.
#'
#' @param stats_dt A data.table containing statistics (output of calcASCStatistics).
#' @param minDonors Integer. Minimum number of independent donors required (default 3).
#' @param fdrCutoff Numeric. FDR threshold to define significance (default 0.05).
#' @return A data.table containing only the recurrent ASC events.
#' @export
filterForRecurrence <- function(stats_dt, minDonors = 3, fdrCutoff = 0.05) {

  logger.start("Identifying Recurrent ASC Sites")

  if (!"donor" %in% colnames(stats_dt)) {
    stats_dt[, donor := stringr::str_extract(sampleId, "\\d+")]
    if (any(is.na(stats_dt$donor))) {
      logger.error("Could not parse donor ID from sampleId.")
      stop("Donor ID parsing failed.")
    }
  }

  sig_dt <- stats_dt[fdr < fdrCutoff]

  if (nrow(sig_dt) == 0) {
    logger.warning("No significant ASC events found below FDR cutoff.")
    return(data.table(snpId = character(0)))
  }

  recurrence_counts <- sig_dt[,
    list(
      n_donors_sig = length(unique(donor)),
      mean_LFC     = mean(log2FC_norm)
    ),
    by = snpId
  ]

  recurrent_dt <- recurrence_counts[n_donors_sig >= minDonors]

  logger.info(paste("Initial significant events:", nrow(sig_dt)))
  logger.info(paste("Recurrent events (N>=", minDonors, "donors):", nrow(recurrent_dt)))
  logger.completed()

  return(recurrent_dt)
}