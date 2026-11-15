#' @include DsASC-class.R
NULL

################################################################################
# Analysis of allele counts after cleanup
################################################################################
# The cleanup functions (ascDropHomozygous, ascFilterRegions, ascPoolReplicates)
# return count matrices. The functions in this file continue from those matrices:
# the allele frequency filter, the test for allelic imbalance, the motif preference
# of each SNP, and the test of binding preference against allelic accessibility.

# "chr:pos_REF_ALT" ids as a table of single-nucleotide sites.
.ascParseSnpIds <- function(snpIds) {
  m  <- regmatches(snpIds, regexec("^([^:]+):([0-9]+)_([^_]+)_([^_]+)$", snpIds))
  ok <- lengths(m) == 5
  dt <- data.table::data.table(
    snpId = snpIds[ok],
    chrom = vapply(m[ok], `[`, character(1), 2),
    pos   = as.integer(vapply(m[ok], `[`, character(1), 3)),
    REF   = vapply(m[ok], `[`, character(1), 4),
    ALT   = vapply(m[ok], `[`, character(1), 5))
  dt[nchar(REF) == 1 & nchar(ALT) == 1]
}

#' Keep SNPs that are common in a reference population
#'
#' When genotypes are called from the ATAC-seq reads, a homozygous reference site
#' that is called heterozygous carries almost no alternative reads and looks like
#' a strong allelic effect. Such errors are frequent at rare and unknown variants,
#' so the analysis can be restricted to variants that are common in a reference
#' population, for example the 1000 Genomes Project.
#'
#' @param snpIds character vector of "chr:pos_REF_ALT" ids (rownames of the
#'               count matrices).
#' @param mafDb  a \code{MafDb} object of the \code{GenomicScores} package, for
#'               example \code{MafDb.1Kgenomes.phase3.GRCh38}.
#' @param minMaf SNPs with an allele frequency above this value are kept
#'               (default 0.01). SNPs without an entry in \code{mafDb} are removed.
#' @param pop    name of the population in \code{mafDb}. By default the frequency
#'               over all populations ("AF" or "AF_ALL").
#' @param verbose log how many sites were removed (default TRUE).
#' @return the subset of \code{snpIds} with an allele frequency above \code{minMaf}.
#' @author Irem B. GUNDUZ
#' @export
ascFilterCommonVariants <- function(snpIds, mafDb, minMaf = 0.01, pop = NULL, verbose = TRUE) {
  if (!requireNamespace("GenomicScores", quietly = TRUE)) {
    stop("The 'GenomicScores' package is required for the allele frequency filter.")
  }
  if (is.null(pop)) pop <- intersect(c("AF", "AF_ALL"), GenomicScores::populations(mafDb))[1]
  if (is.na(pop)) stop("No population given and neither 'AF' nor 'AF_ALL' in mafDb.")

  key <- sub("_.*$", "", snpIds)                    # "chr:pos"
  gr  <- GenomicRanges::GRanges(sub(":.*$", "", key),
                                IRanges::IRanges(as.integer(sub("^.*:", "", key)), width = 1L))
  try(GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(mafDb)[1], silent = TRUE)
  maf <- as.numeric(S4Vectors::mcols(GenomicScores::gscores(mafDb, gr, pop = pop))[[pop]])

  keep <- !is.na(maf) & maf > minMaf
  if (verbose) logger.info(paste0("ascFilterCommonVariants: removed ", sum(!keep), " of ",
                 length(snpIds), " sites with an allele frequency of ", minMaf,
                 " or less, or without an entry."))
  snpIds[keep]
}

#' Test allele count matrices for allelic imbalance
#'
#' The test of \code{\link{calcASCStatistics}} for a pair of count matrices, for
#' example the output of \code{\link{ascPoolReplicates}}: a two-sided binomial
#' test against 0.5 for every site in every sample, with Benjamini-Hochberg
#' correction within each sample.
#'
#' @param refMat,altMat reference and alternative count matrices (sites x samples)
#'                    with identical dimnames.
#' @param minCoverage minimum number of reads (ref + alt) for a site to be tested
#'                    in a sample (default 10).
#' @param minAllele   minimum number of reads on each allele (default 2).
#' @return A data.table with one row per tested site and sample: snpId, sampleId,
#'         ref, alt, total, log2FC (alt vs ref), pVal and fdr. \code{NULL} if no
#'         site passes the thresholds.
#' @author Irem B. GUNDUZ
#' @export
ascTestCounts <- function(refMat, altMat, minCoverage = 10, minAllele = 2) {
  refMat <- as.matrix(refMat); altMat <- as.matrix(altMat)
  if (!identical(dimnames(refMat), dimnames(altMat))) {
    stop("refMat and altMat must have identical dimnames.")
  }
  dt <- data.table::data.table(
    snpId    = rep(rownames(refMat), times = ncol(refMat)),
    sampleId = rep(colnames(refMat), each  = nrow(refMat)),
    ref      = as.vector(refMat),
    alt      = as.vector(altMat))
  dt <- dt[!is.na(ref) & !is.na(alt)]
  dt[, total := ref + alt]
  dt <- dt[total >= minCoverage & pmin(ref, alt) >= minAllele]
  if (nrow(dt) == 0) {
    logger.warning("No sites passed coverage / min-allele threshold.")
    return(NULL)
  }

  dt[, log2FC := log2((alt + 1) / (ref + 1))]
  pLower <- stats::pbinom(dt$alt,     size = dt$total, prob = 0.5)
  pUpper <- stats::pbinom(dt$alt - 1, size = dt$total, prob = 0.5, lower.tail = FALSE)
  dt[, pVal := pmin(2 * pmin(pLower, pUpper), 1)]
  dt[, fdr  := stats::p.adjust(pVal, method = "BH"), by = sampleId]
  dt[]
}

#' Motif preference of each SNP
#'
#' Scores the sequence around each SNP with a TF motif, once with the reference
#' and once with the alternative allele. The window reaches from (w - 1) bp
#' upstream to (w - 1) bp downstream of the SNP for a motif of width w, so that
#' every placement of the motif overlaps the SNP. \code{delta} is the best motif
#' score with the alternative allele minus the best score with the reference
#' allele. A SNP is assigned a preference only if a motif match overlaps it on at
#' least one allele.
#'
#' @param snpIds    character vector of "chr:pos_REF_ALT" ids. Only
#'                  single-nucleotide variants are scored.
#' @param pwms      a \code{PWMatrixList} or \code{PFMatrixList} of the
#'                  \code{TFBSTools} package, or a list of such matrices.
#' @param genome    a \code{BSgenome} object, for example
#'                  \code{BSgenome.Hsapiens.UCSC.hg38}.
#' @param prefDelta a SNP prefers an allele if the score differs by more than
#'                  this value (default 1).
#' @param matchPval p-value threshold of the motif match (default 2e-3).
#' @return A data.table with the columns snpId, tf (motif name and id), delta and
#'         preference ("Reference", "Alternative" or "No preference").
#' @author Irem B. GUNDUZ
#' @export
ascMotifPreference <- function(snpIds, pwms, genome, prefDelta = 1, matchPval = 2e-3) {
  for (pkg in c("TFBSTools", "motifmatchr", "Biostrings")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The '", pkg, "' package is required for motif scoring."))
    }
  }
  sites <- .ascParseSnpIds(snpIds)
  if (nrow(sites) == 0) stop("No single-nucleotide sites to score.")
  pwms <- lapply(seq_along(pwms), function(i) {
    p <- pwms[[i]]
    if (is(p, "PFMatrix")) TFBSTools::toPWM(p) else p
  })

  scan <- function(pwm, seqs) {
    mm <- motifmatchr::matchMotifs(TFBSTools::PWMatrixList(pwm), seqs, genome = genome,
                                   out = "scores", bg = "even", p.cutoff = matchPval)
    list(score = as.vector(motifmatchr::motifScores(mm)),
         match = as.vector(motifmatchr::motifMatches(mm)))
  }

  data.table::rbindlist(lapply(pwms, function(pwm) {
    w  <- ncol(TFBSTools::Matrix(pwm))
    gr <- GenomicRanges::GRanges(sites$chrom,
                                 IRanges::IRanges(sites$pos - (w - 1), sites$pos + (w - 1)))
    seqRef <- Biostrings::getSeq(genome, gr)
    at <- matrix(FALSE, nrow = length(seqRef), ncol = 2 * w - 1)
    at[, w] <- TRUE
    seqAlt <- Biostrings::replaceLetterAt(seqRef, at, letter = sites$ALT)

    r <- scan(pwm, seqRef); a <- scan(pwm, seqAlt)
    delta <- a$score - r$score
    hit   <- r$match | a$match
    data.table::data.table(
      snpId = sites$snpId,
      tf    = paste0(TFBSTools::name(pwm), "_", TFBSTools::ID(pwm)),
      delta = delta,
      preference = data.table::fifelse(hit & delta >  prefDelta, "Alternative",
                   data.table::fifelse(hit & delta < -prefDelta, "Reference", "No preference")))
  }))
}

#' Reference allele fraction of each SNP in a set of samples
#'
#' Sums the reads over the samples and returns the fraction of reads on the
#' reference allele, joined with the motif preference of the SNPs. This is the
#' table behind a plot of binding preference against allelic accessibility.
#'
#' @param refMat,altMat reference and alternative count matrices (sites x samples).
#' @param sampleIds     samples to sum over.
#' @param prefDt        the output of \code{\link{ascMotifPreference}}, or any
#'                      table with the columns snpId, tf and preference.
#' @param minReads      SNPs with fewer reads in the summed samples are dropped
#'                      (default 4).
#' @return A data.table with one row per SNP and motif: snpId, refFrac, total and
#'         the columns of \code{prefDt}.
#' @author Irem B. GUNDUZ
#' @export
ascPreferenceTable <- function(refMat, altMat, sampleIds, prefDt, minReads = 4) {
  agg  <- ascAggregate(refMat, altMat, sampleIds)
  frac <- data.table::data.table(snpId = rownames(refMat), total = agg$ref + agg$alt)
  frac[, refFrac := agg$ref / total]
  frac <- frac[total >= minReads]
  merge(frac[, .(snpId, refFrac, total)], data.table::as.data.table(prefDt),
        by = "snpId", allow.cartesian = TRUE)
}

#' Test binding preference against allelic accessibility
#'
#' For each motif, compares the reference allele fraction of the SNPs that prefer
#' the reference allele with that of the SNPs that prefer the alternative allele.
#' If the TF opens chromatin, the first group has the higher reference allele
#' fraction.
#'
#' @param prefTable the output of \code{\link{ascPreferenceTable}}.
#' @param method    "wilcoxon" for a two-sided Wilcoxon rank-sum test, in which
#'                  every SNP has the same weight, or "quasibinomial" for a
#'                  logistic regression of the read counts, in which SNPs are
#'                  weighted by their reads.
#' @param pAdjust   method of \code{\link[stats]{p.adjust}} applied over the
#'                  motifs (default "bonferroni").
#' @param minSites  motifs with fewer SNPs in either group are not tested
#'                  (default 2).
#' @return A data.table with one row per motif: tf, the number of SNPs and the
#'         mean reference allele fraction of both groups, their difference
#'         (reference minus alternative preference), p, pAdj, and z, the standard
#'         normal quantile of the two-sided p-value with the sign of the
#'         difference.
#' @author Irem B. GUNDUZ
#' @export
ascPreferenceTest <- function(prefTable, method = c("wilcoxon", "quasibinomial"),
                              pAdjust = "bonferroni", minSites = 2) {
  method <- match.arg(method)
  dt  <- data.table::as.data.table(prefTable)
  dt  <- dt[preference %in% c("Reference", "Alternative") & !is.na(refFrac)]
  res <- dt[, {
    isRef <- preference == "Reference"
    nRef  <- sum(isRef); nAlt <- sum(!isRef)
    p <- NA_real_
    if (min(nRef, nAlt) >= minSites) {
      g <- factor(preference, levels = c("Alternative", "Reference"))
      p <- tryCatch(
        if (method == "wilcoxon") {
          stats::wilcox.test(refFrac ~ g)$p.value
        } else {
          k   <- round(refFrac * total)
          fit <- stats::glm(cbind(k, total - k) ~ g, family = stats::quasibinomial())
          stats::coef(summary(fit))[2, 4]
        }, error = function(e) NA_real_)
    }
    list(sitesReference     = nRef,
         sitesAlternative   = nAlt,
         refFracReference   = mean(refFrac[isRef]),
         refFracAlternative = mean(refFrac[!isRef]),
         p                  = p)
  }, by = tf]
  res[, diff := refFracReference - refFracAlternative]
  res[, pAdj := stats::p.adjust(p, method = pAdjust)]
  res[, z    := sign(diff) * stats::qnorm(p / 2, lower.tail = FALSE)]
  res[order(p)]
}
