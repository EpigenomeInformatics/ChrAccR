#' @include DsAcc-class.R
NULL

# ==============================================================================
# CLASS DEFINITION
# ==============================================================================
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

if (!exists("getCounts", envir = topenv(environment()), inherits = FALSE)) {
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

  # Guard: empty SNP set (e.g. a low-coverage cell type where filterLowCovg
  # removed everything). Return unchanged rather than erroring on a 0-length GRanges.
  if (is.null(snpGr) || length(snpGr) == 0) {
    logger.warning("No SNPs present (empty object) - skipping peak annotation.")
    logger.completed()
    return(dsObj)
  }

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

#' Calculate ASC Statistics (Calderon et al. 2019)
#'
#' Tests for allele-specific chromatin accessibility at each heterozygous site,
#' following Calderon et al., Nat. Genet. 2019. After WASP filtering (which
#' removes reference mapping bias), each site is tested with a two-sided
#' binomial test against the null hypothesis of a 50:50 ref:alt split.
#' Multiple-testing correction (Benjamini-Hochberg) is applied PER SAMPLE.
#'
#' Note: the null is a fixed p = 0.5. WASP is responsible for removing mapping
#' bias upstream; re-centering the null on the observed (pooled) allele fraction
#' is NOT done, because (a) it conflates allele identity across donors, and
#' (b) it would calibrate the null against the very ASC signal being detected.
#'
#' @param dsObj A DsASC object (WASP-filtered, peak-annotated).
#' @param minCoverage Minimum total coverage (ref+alt) required to test a site
#'                    in a sample (default 10).
#' @param minAllele Minimum reads required on EACH allele (min(ref, alt)) for a
#'                    site-sample to be tested (default 2). At low coverage the
#'                    only calls that reach significance against 0.5 are the most
#'                    extreme imbalances, with minor allele = 0; those are
#'                    low-confidence "pseudo-ASC" that inflate the ref skew and
#'                    the significant count. Requiring >= 2 reads on both alleles
#'                    is the standard ASE/ASC guard and keeps calls two-sided.
#' @return A data.table with one row per tested (snp x sample): ref, alt, total,
#'         log2FC (alt vs ref effect size), pVal, fdr (per-sample BH), and peakId.
#' @export
calcASCStatistics <- function(dsObj, minCoverage = 10, minAllele = 2) {

  logger.start("Calculating ASC Statistics (binomial test vs p = 0.5)")

  refMat   <- as.matrix(getRefCounts(dsObj))
  altMat   <- as.matrix(getAltCounts(dsObj))
  totalMat <- refMat + altMat

  dt <- data.table::as.data.table(as.table(totalMat))
  data.table::setnames(dt, c("snpId", "sampleId", "total"))
  dt$ref <- as.vector(refMat)
  dt$alt <- as.vector(altMat)

  # Only test sites with sufficient coverage AND >= minAllele reads on BOTH
  # alleles in that sample (drops minor-allele = 0 low-confidence pseudo-ASC).
  dt <- dt[dt$total >= minCoverage & pmin(dt$ref, dt$alt) >= minAllele]

  if (nrow(dt) == 0) {
    logger.warning("No sites passed coverage / min-allele threshold.")
    return(NULL)
  }

  # Effect size: log2 ratio of alt vs ref (pseudocount-stabilised).
  # Reported for interpretation only; not used in the test.
  dt$log2FC <- log2((dt$alt + 1) / (dt$ref + 1))

  # Two-sided binomial test against a fixed null of p = 0.5 (paper's approach).
  logger.info("Running two-sided binomial test against p = 0.5")
  p_lower <- pbinom(dt$alt,     size = dt$total, prob = 0.5)
  p_upper <- pbinom(dt$alt - 1, size = dt$total, prob = 0.5, lower.tail = FALSE)

  dt$pVal <- 2 * pmin(p_lower, p_upper)
  dt$pVal[dt$pVal > 1] <- 1.0

  # Benjamini-Hochberg FDR computed PER SAMPLE (as in the paper).
  dt[, fdr := p.adjust(pVal, method = "BH"), by = sampleId]

  # Attach peak annotation if present
  if (!is.null(dsObj@coord$snps$peakId)) {
    peak_map  <- setNames(dsObj@coord$snps$peakId, names(dsObj@coord$snps))
    dt$peakId <- peak_map[dt$snpId]
  }

  n_sig <- nrow(dt[fdr < 0.1])
  logger.info(paste0("Tested ", nrow(dt), " site-sample observations; ",
                     n_sig, " significant at FDR < 0.1."))
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

  pi_nonzero <- 1 - ashr::get_pi0(ash_res)   # FIXED: was ash_res$pi[1] (NULL in current ashr)

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

#' @include DsASC-class.R
NULL

# ==============================================================================
# DsASC-analysis.R
# Reusable analysis helpers for allele-specific chromatin (DsASC) objects.
# These were previously duplicated across validation scripts; centralising them
# keeps method definitions in one place.
#
# Style: plain exported functions (matches ChrAccR plotting/analysis helpers).
# ==============================================================================

#' @import data.table
#' @importFrom GenomicRanges seqnames start
NULL

# ------------------------------------------------------------------------------
# Loading / merging per-cell-type array outputs
# ------------------------------------------------------------------------------

#' Merge per-cell-type DsASC objects into unified count matrices
#'
#' Reads ds_asc_*.rds produced by the per-cell-type SLURM array and rebuilds
#' aligned ref/alt matrices (union of SNPs x all samples), plus combined
#' sample annotation. Empties (_EMPTY.rds) are skipped.
#'
#' @param dir          Directory containing ds_asc_*.rds.
#' @param pattern      File pattern (default "^ds_asc_.*\\.rds$").
#' @return list(ref, alt, annot, snpIds) with ref/alt integer matrices.
#' @export
mergeDsASCArray <- function(dir, pattern = "^ds_asc_.*\\.rds$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("_EMPTY\\.rds$", files)]
  if (length(files) == 0) stop(paste("No DsASC array files in", dir))

  annot_list <- list(); ref_list <- list(); alt_list <- list(); all_snps <- character()
  for (f in files) {
    ds <- readRDS(f)
    annot_list[[length(annot_list)+1]] <- data.table::as.data.table(getSampleAnnot(ds))
    r <- as.matrix(getRefCounts(ds)); a <- as.matrix(getAltCounts(ds))
    ref_list[[length(ref_list)+1]] <- r; alt_list[[length(alt_list)+1]] <- a
    all_snps <- c(all_snps, rownames(r))
  }
  all_snps <- unique(all_snps)
  annot <- unique(data.table::rbindlist(annot_list, fill = TRUE))
  annot[, donor := as.character(donor)]

  ref <- matrix(0L, length(all_snps), nrow(annot), dimnames = list(all_snps, annot$sampleId))
  alt <- matrix(0L, length(all_snps), nrow(annot), dimnames = list(all_snps, annot$sampleId))
  for (i in seq_along(ref_list)) {
    r <- ref_list[[i]]; a <- alt_list[[i]]
    ref[rownames(r), colnames(r)] <- r
    alt[rownames(a), colnames(a)] <- a
  }
  list(ref = ref, alt = alt, annot = annot, snpIds = all_snps)
}

# ------------------------------------------------------------------------------
# Core statistics
# ------------------------------------------------------------------------------

#' Per-site binomial posterior + FDR for a ref/alt vector
#'
#' @param refV,altV  integer vectors of reference / alternative counts.
#' @return list(total, fdr) with BH-adjusted two-sided binomial p.
#' @export
ascSitePosterior <- function(refV, altV) {
  total <- refV + altV
  pV <- 2 * pmin(stats::pbinom(refV, total, 0.5), stats::pbinom(altV, total, 0.5))
  pV[total == 0] <- 1
  list(total = total, fdr = stats::p.adjust(pV, method = "BH"))
}

#' Aggregate ref/alt across a set of samples (columns)
#'
#' @param ref,alt   count matrices [snp x sample].
#' @param sampleIds samples to sum over.
#' @return list(ref, alt) row sums.
#' @export
ascAggregate <- function(ref, alt, sampleIds) {
  vs <- intersect(sampleIds, colnames(ref))
  if (length(vs) == 0) return(list(ref = numeric(nrow(ref)), alt = numeric(nrow(alt))))
  list(ref = rowSums(ref[, vs, drop = FALSE], na.rm = TRUE),
       alt = rowSums(alt[, vs, drop = FALSE], na.rm = TRUE))
}

#' Estimate shared-imbalance proportion between two count sets (ashR)
#'
#' Given A-significant sites, estimates the proportion with a non-zero effect in
#' B using ashR (uses get_pi0; robust across ashr versions).
#'
#' @param refA,altA,refB,altB  count vectors aligned by SNP.
#' @param fdrA       FDR cutoff defining significance in A (default 0.01).
#' @param minReadsB  reads in B below which a site is "inaccessible" (default 4).
#' @return list(n_total, Inaccessible, NotShared, Shared) proportions.
#' @export
ascSharing <- function(refA, altA, refB, altB, fdrA = 0.01, minReadsB = 4) {
  postA <- ascSitePosterior(refA, altA)
  sigA  <- which(postA$fdr < fdrA & postA$total > 0)
  if (length(sigA) == 0) return(NULL)

  totalB <- (refB + altB)[sigA]
  inaccessible <- totalB < minReadsB
  accIdx <- sigA[!inaccessible]
  shared_prop <- NA_real_

  if (length(accIdx) >= 10) {
    tB <- refB[accIdx] + altB[accIdx]
    effect <- ((refB[accIdx] + 1) / (tB + 2)) - 0.5
    sem <- sqrt(((refB[accIdx] + 1) * (altB[accIdx] + 1)) / ((tB + 2)^2 * (tB + 3)))
    keep <- sem > 0
    if (sum(keep) >= 10) {
      sp <- tryCatch(1 - ashr::get_pi0(ashr::ash(effect[keep], sem[keep], mixcompdist = "normal")),
                     error = function(e) NA_real_)
      if (length(sp) == 1 && is.finite(sp)) shared_prop <- sp
    }
  }
  n_total <- length(sigA); n_inacc <- sum(inaccessible); n_acc <- n_total - n_inacc
  n_shared <- if (is.na(shared_prop)) 0L else round(shared_prop * n_acc)
  list(n_total = n_total,
       Inaccessible = n_inacc / n_total,
       NotShared = (n_acc - n_shared) / n_total,
       Shared = n_shared / n_total)
}

# ------------------------------------------------------------------------------
# Genotype-quality cleanup
# ------------------------------------------------------------------------------

#' Drop homozygous / mis-genotyped sites using donor-pooled evidence
#'
#' A genuinely heterozygous site must show BOTH alleles in a BALANCED way once
#' reads are pooled across all of a donor's samples. Two kinds of genotyping
#' error are removed (they generate the homozygous "leak" that
#' 04_asc_benchmarking.R quantifies):
#' \enumerate{
#'   \item pooled minor-allele count below \code{minMinor} (strictly monoallelic);
#'   \item pooled minor-allele FRACTION below \code{minMinorFrac} -- a hom-ref (or
#'         hom-alt) site whose few minor reads are just sequencing error scattered
#'         across many samples. An absolute count alone misses these, because
#'         ~0.3\% error across ~50 deep samples can sum to >=2 minor reads while
#'         the pooled fraction stays near zero; the fraction test catches them.
#' }
#' For each donor, offending sites have their counts ZEROED in that donor's sample
#' columns (ref = alt = 0), giving total coverage 0 so they are dropped by the
#' downstream coverage filters and never tested. Zeroing (rather than NA) keeps
#' plain \code{sum()} / \code{rowSums()} calls valid. Genuine ASC survives: even
#' strong allele-specific sites pool to a minor fraction well above
#' \code{minMinorFrac} across a donor's many samples.
#'
#' @param refMat,altMat count matrices [snp x sample].
#' @param annot         sample annotation with \code{sampleId}, \code{donor}.
#' @param minMinor      minimum donor-pooled minor-allele reads (default 2).
#' @param minMinorFrac  minimum donor-pooled minor-allele fraction (default 0.05);
#'                      below this a covered site is treated as homozygous+error.
#' @param minTotal      donor-pooled coverage above which the tests apply rather
#'                      than treating the site as merely low-coverage (default 10).
#' @param verbose       log how many donor x site cells were masked (default TRUE).
#' @return list(ref, alt) with donor-specific zeroing applied.
#' @author Irem B. GUNDUZ
#' @export
ascDropHomozygous <- function(refMat, altMat, annot,
                              minMinor = 2L, minMinorFrac = 0.05,
                              minTotal = 10L, verbose = TRUE) {
  annot <- data.table::as.data.table(annot)
  masked <- 0L; badSites <- 0L
  for (d in unique(annot$donor)) {
    s <- intersect(annot[donor == d, sampleId], colnames(refMat))
    if (length(s) == 0) next
    pr    <- rowSums(refMat[, s, drop = FALSE], na.rm = TRUE)
    pa    <- rowSums(altMat[, s, drop = FALSE], na.rm = TRUE)
    minor <- pmin(pr, pa); tot <- pr + pa
    frac  <- ifelse(tot > 0, minor / tot, 0)
    # covered, but monoallelic by count OR by fraction (error-only minor allele)
    bad   <- which(tot >= minTotal & (minor < minMinor | frac < minMinorFrac))
    if (length(bad) > 0) {
      refMat[bad, s] <- 0L; altMat[bad, s] <- 0L         # zero -> total 0 -> not tested
      masked <- masked + length(bad) * length(s); badSites <- badSites + length(bad)
    }
  }
  if (verbose) logger.info(paste0("ascDropHomozygous: masked ", badSites,
                 " donor x site loci (", masked, " cells) as homozygous/mis-genotyped."))
  list(ref = refMat, alt = altMat)
}

#' Drop ASC sites in blacklist / low-mappability regions
#'
#' Repetitive and low-mappability regions are where WASP fails to remove reference
#' mapping bias AND where GATK makes false-heterozygous genotype calls (reads from
#' a paralog map in). Both push significant ASC toward the reference allele. Since
#' we lack the paper's imputed-genotype QC, removing these regions is the standard,
#' faithful mitigation (Calderon et al. also excluded blacklist regions in QC).
#'
#' @param snpIds   character vector of "chr:pos_REF_ALT" ids (rownames of the
#'                 count matrices).
#' @param regionsGr GRanges of regions to remove (ENCODE blacklist and/or
#'                 low-mappability). Seqnames must be UCSC-style ("chr1", ...).
#' @param verbose  log how many sites were removed (default TRUE).
#' @return the subset of \code{snpIds} NOT overlapping \code{regionsGr}.
#' @author Irem B. GUNDUZ
#' @export
ascFilterRegions <- function(snpIds, regionsGr, verbose = TRUE) {
  chrpos <- sub("_.*$", "", snpIds)                 # "chr:pos"
  chrom  <- sub(":.*$", "", chrpos)
  pos    <- suppressWarnings(as.integer(sub("^.*:", "", chrpos)))
  ok     <- !is.na(pos)
  gr     <- GenomicRanges::GRanges(chrom[ok], IRanges::IRanges(pos[ok], pos[ok]))
  ov     <- GenomicRanges::findOverlaps(gr, regionsGr)
  badLocal <- unique(S4Vectors::queryHits(ov))
  bad    <- which(ok)[badLocal]
  keep   <- setdiff(seq_along(snpIds), bad)
  if (verbose) logger.info(paste0("ascFilterRegions: removed ", length(bad), " / ",
                 length(snpIds), " sites in blacklist/low-mappability regions."))
  snpIds[keep]
}

#' Pool technical-replicate libraries into biological samples
#'
#' Calderon et al. count and test allele-specific chromatin at the level of a
#' biological sample (a donor x cell-type x condition), merging technical
#' replicate libraries. Our per-library ATAC objects instead carry each replicate
#' as its own column, which halves the reads per binomial test and sharply
#' reduces power (median per-test coverage ~17 vs a paper-comparable ~34). This
#' function sums ref/alt counts across libraries that share the same
#' \code{groupCols}, reproducing the paper's per-sample unit before testing.
#'
#' @param refMat,altMat count matrices [snp x library].
#' @param annot         sample annotation with \code{sampleId} plus \code{groupCols}.
#' @param groupCols     columns defining a biological sample
#'                      (default c("cellType","stimulus","donor")).
#' @param verbose       log the collapse (default TRUE).
#' @return list(ref, alt, annot): pooled matrices [snp x biological-sample] and a
#'         one-row-per-sample annotation whose \code{sampleId} is the group key.
#' @author Irem B. GUNDUZ
#' @export
ascPoolReplicates <- function(refMat, altMat, annot,
                              groupCols = c("cellType", "stimulus", "donor"),
                              verbose = TRUE) {
  annot <- data.table::as.data.table(annot)
  if (!all(groupCols %in% names(annot)))
    stop(paste("annot lacks grouping columns:",
               paste(setdiff(groupCols, names(annot)), collapse = ", ")))

  a   <- annot[match(colnames(refMat), sampleId)]
  grp <- do.call(paste, c(a[, ..groupCols], sep = "_"))
  idxByGroup <- split(seq_along(grp), grp)
  ug  <- names(idxByGroup)

  poolRef <- vapply(ug, function(g)
    rowSums(refMat[, idxByGroup[[g]], drop = FALSE], na.rm = TRUE), numeric(nrow(refMat)))
  poolAlt <- vapply(ug, function(g)
    rowSums(altMat[, idxByGroup[[g]], drop = FALSE], na.rm = TRUE), numeric(nrow(altMat)))
  rownames(poolRef) <- rownames(refMat); rownames(poolAlt) <- rownames(altMat)

  poolAnnot <- unique(a[, ..groupCols])
  poolAnnot[, sampleId := do.call(paste, c(.SD, sep = "_")), .SDcols = groupCols]
  poolAnnot <- poolAnnot[match(ug, sampleId)]              # align annot to matrix cols
  data.table::setcolorder(poolAnnot, c("sampleId", groupCols))

  if (verbose) logger.info(paste0("ascPoolReplicates: pooled ", ncol(refMat),
                 " libraries into ", length(ug), " biological samples (",
                 paste(groupCols, collapse = " x "), ")."))
  list(ref = poolRef, alt = poolAlt, annot = poolAnnot)
}

# ------------------------------------------------------------------------------
# Donor utilities
# ------------------------------------------------------------------------------

#' Rank donors by sequencing depth at tested sites
#' @export
ascDonorDepth <- function(ref, alt, annot) {
  data.table::rbindlist(lapply(unique(annot$donor), function(d) {
    s <- intersect(annot[donor == d, sampleId], colnames(ref))
    data.table::data.table(donor = d,
      total_reads = sum(ref[, s, drop = FALSE]) + sum(alt[, s, drop = FALSE]))
  }))[order(-total_reads)]
}

#' Best donor for a SNP: het there (per master membership) and most reads
#'
#' @param sn          snpId.
#' @param ref,alt     count matrices.
#' @param annot       sample annotation (with donor).
#' @param membership  optional named list snpId -> donors het there.
#' @return donor id, or NA.
#' @export
ascBestDonor <- function(sn, ref, alt, annot, membership = NULL) {
  donor_ids <- unique(annot$donor)
  cand <- if (!is.null(membership) && !is.null(membership[[sn]]))
            as.character(membership[[sn]]) else donor_ids
  cand <- intersect(cand, donor_ids)
  if (length(cand) == 0) return(NA_character_)
  depths <- sapply(cand, function(d) {
    s <- intersect(annot[donor == d, sampleId], colnames(ref))
    if (length(s) == 0) return(0)
    sum(ref[sn, s]) + sum(alt[sn, s])
  })
  if (max(depths) <= 0) return(NA_character_)
  cand[which.max(depths)]
}

#' Stratify a pair of (cellType, stimulus) groups into the 4 sharing strata
#' @export
ascStratum <- function(cellA, condA, cellB, condB) {
  same_lin <- cellA == cellB; same_cond <- condA == condB
  data.table::fifelse(same_lin & same_cond, "Same lineage & condition",
   data.table::fifelse(same_lin & !same_cond, "Same lineage, diff condition",
    data.table::fifelse(!same_lin & same_cond, "Diff lineage, same condition",
                                               "Diff lineage & condition")))
}

#' @include DsASC-class.R
NULL

# ==============================================================================
# DsASC-plots.R
# Plotting functions for allele-specific chromatin (DsASC) objects.
#
# Style: plain exported functions returning ggplot objects (matches ChrAccR's
# existing plotting helpers, which are functions rather than S4 methods).
# Each takes a DsASC object as the first argument.
#
# Depends on calcASCStatistics() from DsASC-class.R.
# ==============================================================================

#' @import ggplot2
#' @importFrom data.table as.data.table setnames := data.table
NULL

# ------------------------------------------------------------------------------
# Internal: build a per-sample allelic-balance table for a set of sites
# ------------------------------------------------------------------------------
.ascSiteTable <- function(dsObj, snpIds, minReads = 4) {
  ref <- as.matrix(getRefCounts(dsObj))
  alt <- as.matrix(getAltCounts(dsObj))

  miss <- setdiff(snpIds, rownames(ref))
  if (length(miss) > 0) {
    logger.warning(paste(length(miss), "requested SNP(s) not in object; dropping."))
    snpIds <- intersect(snpIds, rownames(ref))
  }
  if (length(snpIds) == 0) stop("None of the requested SNPs are present.")

  annot <- as.data.table(getSampleAnnot(dsObj))
  rows <- list()
  for (sn in snpIds) {
    r <- ref[sn, ]; a <- alt[sn, ]
    total <- r + a
    refFrac <- r / total
    # Wilson-ish 95% CI from read depth (binomial)
    se <- sqrt(refFrac * (1 - refFrac) / total)
    dt <- data.table(
      snpId    = sn,
      sampleId = colnames(ref),
      ref = r, alt = a, total = total,
      refFrac = refFrac,
      lo = pmax(0, refFrac - 1.96 * se),
      hi = pmin(1, refFrac + 1.96 * se)
    )
    rows[[sn]] <- dt
  }
  tab <- data.table::rbindlist(rows)
  tab <- merge(tab, annot, by = "sampleId", all.x = TRUE)
  tab <- tab[!is.na(total) & total >= minReads]
  tab
}

# ------------------------------------------------------------------------------
# 1. Figure 4a — allelic-balance forest plot
# ------------------------------------------------------------------------------
#' Plot allele-specific chromatin balance across samples (Fig 4a style)
#'
#' For one or more heterozygous SNPs, plots the proportion of reads mapping to
#' the reference allele in each sample, with binomial confidence intervals.
#' Samples can be coloured by any annotation column (e.g. stimulus).
#'
#' @param dsObj      A \code{\linkS4class{DsASC}} object.
#' @param snpIds     Character vector of snpIds (rownames of the count matrix).
#' @param colorBy    Sample-annotation column to colour points by (default "stimulus").
#' @param orderBy    Sample-annotation column to order the y axis by (default "cellType").
#' @param minReads   Minimum reads to display a sample at a site (default 4).
#' @param sigStats   Optional output of \code{calcASCStatistics}; significant
#'                   (fdr < sigCut) sample-site points are drawn solid, others faded.
#' @param sigCut     FDR threshold for "significant" shading (default 0.1).
#' @return A ggplot object.
#' @export
plotASCBalance <- function(dsObj, snpIds, colorBy = "stimulus",
                           orderBy = "cellType", minReads = 4,
                           sigStats = NULL, sigCut = 0.1) {
  tab <- .ascSiteTable(dsObj, snpIds, minReads = minReads)

  # significance shading
  tab$sig <- TRUE
  if (!is.null(sigStats)) {
    key <- paste(sigStats$snpId, sigStats$sampleId)
    sigset <- key[sigStats$fdr < sigCut]
    tab$sig <- paste(tab$snpId, tab$sampleId) %in% sigset
  }

  if (!is.null(orderBy) && orderBy %in% colnames(tab)) {
    tab$sampleId <- factor(tab$sampleId,
                           levels = unique(tab$sampleId[order(tab[[orderBy]])]))
  }

  aes_color <- if (!is.null(colorBy) && colorBy %in% colnames(tab)) colorBy else NULL

  p <- ggplot(tab, aes(x = refFrac, y = sampleId)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi, alpha = sig), height = 0) +
    geom_point(aes_string(colour = aes_color, alpha = "sig"), size = 2) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    facet_wrap(~ snpId, nrow = 1) +
    labs(x = "Proportion of reads mapping to reference allele",
         y = NULL, colour = colorBy,
         title = "Allele-specific chromatin accessibility") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6))
  p
}

# ------------------------------------------------------------------------------
# 2. Per-sample ASC volcano
# ------------------------------------------------------------------------------
#' Volcano plot of ASC effect size vs significance
#'
#' @param dsObj      A \code{\linkS4class{DsASC}} object.
#' @param stats      Optional output of \code{calcASCStatistics}. If NULL it is
#'                   computed with the given minCoverage.
#' @param sample     Optional single sampleId to restrict the plot to.
#' @param minCoverage Passed to calcASCStatistics if stats is NULL (default 10).
#' @param sigCut     FDR threshold for highlighting (default 0.1).
#' @return A ggplot object.
#' @export
plotASCVolcano <- function(dsObj, stats = NULL, sample = NULL,
                           minCoverage = 10, sigCut = 0.1) {
  if (is.null(stats)) stats <- calcASCStatistics(dsObj, minCoverage = minCoverage)
  if (is.null(stats)) stop("No statistics available to plot.")

  dt <- data.table::as.data.table(stats)
  if (!is.null(sample)) dt <- dt[sampleId == sample]
  if (nrow(dt) == 0) stop("No rows to plot (check 'sample').")

  dt[, negLogP := -log10(pmax(pVal, .Machine$double.xmin))]
  dt[, sig := fdr < sigCut]

  p <- ggplot(dt, aes(x = log2FC, y = negLogP, colour = sig)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c("TRUE" = "#C0392B", "FALSE" = "grey70"),
                        name = paste0("FDR < ", sigCut)) +
    labs(x = expression(log[2]~"(alt / ref)"),
         y = expression(-log[10]~"(p)"),
         title = if (is.null(sample)) "ASC volcano (all samples)" else paste("ASC volcano:", sample)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  p
}

# ------------------------------------------------------------------------------
# 3. Manhattan / genome-position view of ASC
# ------------------------------------------------------------------------------
#' Manhattan-style plot of ASC significance along the genome
#'
#' @param dsObj      A \code{\linkS4class{DsASC}} object (provides snp coordinates).
#' @param stats      Optional output of \code{calcASCStatistics}. If NULL it is
#'                   computed with the given minCoverage.
#' @param sample     Optional single sampleId to restrict the plot to. If NULL,
#'                   the most significant observation per SNP is shown.
#' @param minCoverage Passed to calcASCStatistics if stats is NULL (default 10).
#' @param sigCut     FDR threshold for the significance line (default 0.1).
#' @return A ggplot object.
#' @export
plotASCManhattan <- function(dsObj, stats = NULL, sample = NULL,
                             minCoverage = 10, sigCut = 0.1) {
  if (is.null(stats)) stats <- calcASCStatistics(dsObj, minCoverage = minCoverage)
  if (is.null(stats)) stop("No statistics available to plot.")

  dt <- data.table::as.data.table(stats)
  if (!is.null(sample)) dt <- dt[sampleId == sample]

  # coordinates from the object
  gr <- dsObj@coord$snps
  coord <- data.table(
    snpId = names(gr),
    chrom = as.character(GenomicRanges::seqnames(gr)),
    pos   = GenomicRanges::start(gr)
  )
  dt <- merge(dt, coord, by = "snpId")
  if (nrow(dt) == 0) stop("No rows with coordinates to plot.")

  # collapse to most significant per SNP if multiple samples
  if (is.null(sample)) {
    dt <- dt[dt[, .I[which.min(fdr)], by = snpId]$V1]
  }

  # order chromosomes naturally and build cumulative x position
  chrom_order <- paste0("chr", c(1:22, "X", "Y"))
  dt <- dt[chrom %in% chrom_order]
  dt[, chrom := factor(chrom, levels = chrom_order)]
  setkey(dt, chrom, pos)

  chrom_len <- dt[, .(maxpos = max(pos)), by = chrom]
  chrom_len[, offset := cumsum(as.numeric(maxpos)) - maxpos]
  dt <- merge(dt, chrom_len[, .(chrom, offset)], by = "chrom")
  dt[, xpos := pos + offset]
  dt[, negLogFDR := -log10(pmax(fdr, .Machine$double.xmin))]

  axis_df <- dt[, .(center = mean(range(xpos))), by = chrom]

  p <- ggplot(dt, aes(x = xpos, y = negLogFDR, colour = chrom)) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_hline(yintercept = -log10(sigCut), linetype = "dashed", colour = "#C0392B") +
    scale_colour_manual(values = rep(c("#34495E", "#95A5A6"),
                                     length.out = nlevels(dt$chrom)), guide = "none") +
    scale_x_continuous(breaks = axis_df$center, labels = sub("chr", "", axis_df$chrom)) +
    labs(x = "Chromosome", y = expression(-log[10]~"(FDR)"),
         title = if (is.null(sample)) "ASC across the genome (best per SNP)" else paste("ASC:", sample)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  p
}