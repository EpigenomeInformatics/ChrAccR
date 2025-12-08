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

#' Calculate ASC Statistics (Conditional Bias Correction)
#'
#' Performs a two-sided Binomial Test. If global residual bias is detected, 
#' the null hypothesis is automatically adjusted to correct for the bias.
#'
#' @param dsObj A DsASC object (filtered).
#' @param minCoverage Minimum coverage to perform test (default 10).
#' @param globalProbTolerance Tolerance for deviation from p=0.5 (default 0.02).
#' @return A data.table containing statistics, including log2FC_report which is bias-corrected only if needed.
#' @export
calcASCStatistics <- function(dsObj, minCoverage=10, globalProbTolerance=0.02) {
  
  logger.start("Calculating ASC Statistics (Conditional Correction)")
  
  # Fetch Data
  refMat <- as.matrix(getRefCounts(dsObj))
  altMat <- as.matrix(getAltCounts(dsObj))
  totalMat <- refMat + altMat
  
  dt <- data.table::as.data.table(as.table(totalMat))
  setnames(dt, c("snpId", "sampleId", "total"))
  dt$ref <- as.vector(refMat)
  dt$alt <- as.vector(altMat)
  
  dt <- dt[dt$total >= minCoverage]
  
  if (nrow(dt) == 0) {
      logger.warning("No sites passed coverage threshold.")
      return(NULL)
  }
  
  # Conditional Bias Check
  global_ref <- sum(as.numeric(dt$ref))
  global_total <- sum(as.numeric(dt$total))
  
  # Global Alt Proportion (our observed p)
  global_alt_prob <- sum(as.numeric(dt$alt)) / global_total
  
  # Decide whether to normalize
  if (abs(global_alt_prob - 0.5) > globalProbTolerance) {
    test_prob <- global_alt_prob
    is_normalized <- TRUE
    logger.warning(paste0("Bias Detected! Global Alt Prob is ", round(global_alt_prob, 3), ". Test probability set to observed bias."))
  } else {
    test_prob <- 0.5
    is_normalized <- FALSE
    logger.info("Bias not significant. Test probability set to 0.5.")
  }
  
  # Calculate LFC Metrics
  dt$log2FC_raw <- log2((dt$alt + 1) / (dt$ref + 1))
  
  if (is_normalized) {
    # Calculate normalized LFC to center plots at 0
    global_odds <- (1 - global_alt_prob) / global_alt_prob
    dt$log2FC_report <- log2((dt$alt / dt$ref) / global_odds)
  } else {
    # Use raw LFC if no bias was detected
    dt$log2FC_report <- dt$log2FC_raw
  }
  
  # Binomial Test using the determined test_prob
  logger.info(paste("Running Binomial Test against p =", round(test_prob, 3)))
  
  # We test Alt count (successes) vs Total count (trials)
  p_lower <- pbinom(dt$alt, size=dt$total, prob=test_prob)
  p_upper <- pbinom(dt$alt - 1, size=dt$total, prob=test_prob, lower.tail=FALSE)
  
  dt$pVal <- 2 * pmin(p_lower, p_upper)
  dt$pVal[dt$pVal > 1] <- 1.0
  
  # FDR Correction
  dt[, fdr := p.adjust(pVal, method = "BH"), by = sampleId]
  
  # Add Metadata
  if (!is.null(dsObj@coord$snps$peakId)) {
    peak_map <- setNames(dsObj@coord$snps$peakId, names(dsObj@coord$snps))
    dt$peakId <- peak_map[dt$snpId]
  }
  
  # Rename the reporting LFC column
  data.table::setnames(dt, "log2FC_report", "log2FC_norm")
  
  logger.completed()
  return(dt)
}


#' Estimate the Proportion of Shared ASC Imbalance Effects 
#'
#' This function implements the Corces 2016 method to estimate the proportion of ASC 
#' sites (significant in Cell Type A) that still show a non-zero effect in Cell Type B.
#' It uses the Beta-Binomial posterior to calculate effect size and variance.
#'
#' @param dsObjA DsASC object for Cell Type A.
#' @param dsObjB DsASC object for Cell Type B.
#' @param fdrCutoff Numeric. FDR threshold for defining significance in Cell Type A (Default is 0.01).
#' @return A list containing the ashR model output and the estimated proportion of shared effects.
#' @export
estimateSharedImbalance <- function(dsObjA, dsObjB, fdrCutoff=0.01) {
  
  logger.start("Estimating Shared Imbalance (ashR Method)")
  
  if (!requireNamespace("ashr", quietly = TRUE)) {
    logger.error("The 'ashr' package is required for this analysis.")
    stop("Missing package: ashr")
  }
  
  # Helper function to calculate Bayesian Posterior Stats
  calculate_posterior <- function(dsObj) {
    # Sum counts across all samples/donors within the cell type
    ref <- rowSums(as.matrix(getRefCounts(dsObj)), na.rm=TRUE)
    alt <- rowSums(as.matrix(getAltCounts(dsObj)), na.rm=TRUE)
    total <- ref + alt
    
    # Filter sites with zero coverage
    valid_sites <- total > 0
    ref <- ref[valid_sites]; alt <- alt[valid_sites]; total <- total[valid_sites]
    
    # Posterior Mean (mu) and Variance (sigma2) of the proportion of reference reads
    # Posterior is Beta(r+1, a+1) 
    mu <- (ref + 1) / (total + 2)
    sigma2 <- ((ref + 1) * (alt + 1)) / ((total + 2)^2 * (total + 3))
    
    # Simple Binomial Test to identify significant sites for filtering A
    pVal <- 2 * pmin(pbinom(ref, total, 0.5), pbinom(alt, total, 0.5))
    fdr <- p.adjust(pVal, method="BH")
    
    return(data.table(
      snpId = names(ref),
      mu = mu,
      sigma2 = sigma2,
      fdr = fdr
    ))
  }
  
  # Calculate Full Bayesian Statistics for Both Cell Types
  dtA <- calculate_posterior(dsObjA)
  dtB <- calculate_posterior(dsObjB)
  
  # Define Set of ASC Sites in Cell Type A (Discovery Set)
  sigA <- dtA[fdr < fdrCutoff]
  logger.info(paste("Identified", nrow(sigA), "significant ASC sites in Cell Type A (FDR <", fdrCutoff, ")."))
  
  if (nrow(sigA) < 100) {
    logger.warning("Fewer than 100 significant sites found. ashR power may be low.")
  }
  
  # Collect Effects in Cell Type B (Test Set)
  merged <- merge(sigA[, .(snpId)], dtB, by="snpId", suffixes = c("_A", "_B"))
  
  # Define the Effect Size (Bhat) and Standard Error (SE) 
  # Effect Size (E) = mu - 0.5 (measured in Cell Type B)
  merged[, effect_size := mu - 0.5]
  merged[, sem := sqrt(sigma2)]
  
  # Filter out sites where variance is zero
  merged <- merged[sem > 0]
  
  # Estimate Proportion of Nonzero Effects using ashR
  logger.info("Running ashR to estimate shared proportion...")
  
  ash_res <- ashr::ash(
    Bhat = merged$effect_size, # Effect size vector (mu - 0.5)
    SEbetahat = merged$sem,    # Standard Error of the Effect Size
    mixcompdist = "normal"     
  )
  
  # We want the proportion that is NON-ZERO (shared/true effect).
  pi_nonzero <- 1 - ash_res$pi[1]
  
  logger.info(paste0("Estimated Proportion of Shared Effects (Non-Zero): ", round(pi_nonzero, 3)))
  logger.completed()
  
  return(list(
    ash_model = ash_res,
    sharing_estimate = pi_nonzero,
    sites_tested = nrow(merged)
  ))
}


#' Filter for Recurrent Allele-Specific Chromatin Events
#'
#' Identifies SNPs that are significantly imbalanced (ASC) in a minimum number 
#' of independent donors.
#'
#' @param stats_dt A data.table containing statistics (output of calcASCStatistics).
#' @param minDonors Integer. Minimum number of independent donors required for recurrence (default 3).
#' @param fdrCutoff Numeric. FDR threshold to define significance (default 0.05).
#' @return A data.table containing only the recurrent ASC events.
#' @export
filterForRecurrence <- function(stats_dt, minDonors = 3, fdrCutoff = 0.05) {
  
  logger.start("Identifying Recurrent ASC Sites")
  
  # 1. Derive Donor ID from sampleId (If not already present)
  if (!"donor" %in% colnames(stats_dt)) {
    # Assuming donor ID is numeric and embedded in sampleId (e.g., 1008 in Mono_S_1008_a)
    stats_dt[, donor := stringr::str_extract(sampleId, "\\d+")]
    if (any(is.na(stats_dt$donor))) {
      logger.error("Could not parse donor ID from sampleId.")
      stop("Donor ID parsing failed.")
    }
  }

  # 2. Identify Significant Events Per Donor
  sig_dt <- stats_dt[fdr < fdrCutoff]
  
  if (nrow(sig_dt) == 0) {
    logger.warning("No significant ASC events found below FDR cutoff.")
    return(data.table(snpId=character(0)))
  }

  # 3. Count Recurrence Across Independent Donors
  recurrence_counts <- sig_dt[, 
    list(
      n_donors_sig = length(unique(donor)),
      mean_LFC = mean(log2FC_norm)
    ), 
    by = snpId
  ]
  
  # 4. Apply Recurrence Filter
  recurrent_dt <- recurrence_counts[n_donors_sig >= minDonors]
  
  logger.info(paste("Initial significant events:", nrow(sig_dt)))
  logger.info(paste("Recurrent events (N>=", minDonors, "donors):", nrow(recurrent_dt)))

  logger.completed()
  
  # Return the table of recurrent SNP IDs and their mean LFC
  return(recurrent_dt)
}