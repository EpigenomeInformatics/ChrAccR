#' @include DsASC-class.R
NULL

# ==============================================================================
# DsASC-tfactivity.R
# Allele-specific transcription-factor ACTIVITY for DsASC objects.
#
# Motivation
# ----------
# calcASCStatistics() and the 03_tf_binding.R workflow operate at the level of a
# single heterozygous SNP: they test allele-specific chromatin (ASC) per site and
# correlate a per-site PWM binding delta with the per-site reference fraction.
#
# This file aggregates ONE level up. For a transcription factor T, it pools every
# ASC SNP that disrupts a T motif and asks a single question:
#
#     Across all of T's motif-disrupting variants, is the allele PREDICTED to
#     bind T better also the allele that is MORE ACCESSIBLE?
#
# The answer is a per-(TF x context) "allele-specific TF activity" score. It is a
# chromVAR-analogue, but instead of aggregating accessibility deviations across
# background-matched peaks (with the sample as the contrasted axis), it aggregates
# allelic imbalance across a motif's het SNPs, with the two ALLELES as the
# contrasted axis. The null is obtained by permuting each site's predicted-binder
# orientation, so every site is its own internal control (ref vs alt at the same
# locus in the same individuals). This exploits the allelic design directly and
# avoids chromVAR's GC/background matching.
#
# A positive score means: the allele predicted to bind T better is preferentially
# accessible -> evidence of allele-specific TF activity in that context. Computing
# the score separately in resting vs stimulated samples reveals TFs whose
# allele-specific activity switches with cell state ("dynamic TF handoff").
#
# Style: plain exported functions returning data.tables / ggplots, matching the
# rest of the DsASC analysis helpers. The core statistic is pure-numeric so it can
# be unit-tested without motif or genome dependencies.
#
# @author Irem B. GUNDUZ
# ==============================================================================

#' @import data.table
NULL

# ------------------------------------------------------------------------------
# Core statistic (pure numeric; no motif / genome dependency)
# ------------------------------------------------------------------------------

#' Orient allele counts toward the predicted stronger-binding allele
#'
#' Given reference / alternative counts and a per-site PWM binding delta
#' (\code{delta = score(ALT) - score(REF)}), returns the read counts on the
#' predicted "high-binding" (\code{H}) and "low-binding" (\code{L}) allele.
#' Sites with \code{delta > 0} bind better on the ALT allele, so \code{H = alt};
#' sites with \code{delta < 0} bind better on REF, so \code{H = ref}.
#'
#' @param refV,altV integer vectors of reference / alternative counts (aligned).
#' @param delta     numeric vector, PWM score difference \code{alt - ref}.
#' @return list(H, L) oriented count vectors.
#' @author Irem B. GUNDUZ
#' @export
ascOrientReads <- function(refV, altV, delta) {
  altStronger <- delta > 0
  H <- ifelse(altStronger, altV, refV)
  L <- ifelse(altStronger, refV, altV)
  list(H = as.numeric(H), L = as.numeric(L))
}

#' Allele-specific TF-activity statistic for one motif's sites
#'
#' The observed activity is the net read fraction toward the predicted
#' stronger-binding allele:
#' \deqn{obs = \frac{\sum_i (H_i - L_i)}{\sum_i (H_i + L_i)}}
#' bounded in \code{[-1, 1]}; \code{obs > 0} means the predicted-binder allele is
#' preferentially accessible. Significance is assessed with an
#' orientation-permutation null: each site's H/L assignment is flipped with
#' probability 0.5 (a Rademacher sign flip), which is the exact null of "predicted
#' binding is unrelated to which allele is accessible". This is vectorised as a
#' single matrix multiply over \code{nPerm} sign draws.
#'
#' @param H,L    oriented count vectors (from \code{ascOrientReads}).
#' @param nPerm  number of orientation permutations (default 2000).
#' @param seed   optional RNG seed for reproducibility.
#' @return list(obs, z, p, meanNull, sdNull, nSites).
#' @author Irem B. GUNDUZ
#' @export
ascActivityStat <- function(H, L, nPerm = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tot   <- H + L
  keep  <- tot > 0
  H <- H[keep]; L <- L[keep]; tot <- tot[keep]
  n <- length(H)
  if (n == 0 || sum(tot) == 0) {
    return(list(obs = NA_real_, z = NA_real_, p = NA_real_,
                meanNull = NA_real_, sdNull = NA_real_, nSites = 0L))
  }
  d       <- H - L
  sumTot  <- sum(tot)
  obs     <- sum(d) / sumTot

  # Vectorised orientation-permutation null: flips is n x nPerm of +/-1.
  flips   <- matrix(sample(c(-1, 1), n * nPerm, replace = TRUE), nrow = n, ncol = nPerm)
  nullV   <- as.numeric(crossprod(d, flips)) / sumTot     # length nPerm

  meanNull <- mean(nullV); sdNull <- stats::sd(nullV)
  z <- if (isTRUE(sdNull > 0)) (obs - meanNull) / sdNull else NA_real_
  # two-sided empirical p (with +1 smoothing)
  p <- (sum(abs(nullV - meanNull) >= abs(obs - meanNull)) + 1) / (nPerm + 1)

  list(obs = obs, z = z, p = p, meanNull = meanNull, sdNull = sdNull, nSites = n)
}

# ------------------------------------------------------------------------------
# Driver over a full delta table + count matrices
# ------------------------------------------------------------------------------

#' Compute allele-specific TF activity for a set of samples
#'
#' Pools reads across the requested samples, orients each ASC SNP toward its
#' predicted stronger-binding allele using \code{deltaDt}, and computes the
#' \code{\link{ascActivityStat}} per transcription factor. Only sites whose
#' absolute PWM delta exceeds \code{prefDelta} (i.e. the motif is actually
#' disrupted) and whose pooled coverage reaches \code{minReads} are used.
#'
#' @param refMat,altMat  count matrices [snp x sample] (e.g. from
#'                        \code{mergeDsASCArray}).
#' @param deltaDt         data.table with columns \code{snpId}, \code{tf},
#'                        \code{delta} (PWM score alt - ref), as produced by the
#'                        motif-scoring step of the driver script.
#' @param sampleIds       samples to pool over (a condition, lineage, etc.).
#' @param prefDelta       |delta| threshold to count a site as motif-disrupting
#'                        (default 1.0, matching 03_tf_binding.R).
#' @param minReads        minimum pooled coverage (H+L) to keep a site (default 4).
#' @param minSites        minimum number of usable sites to report a TF (default 10).
#' @param nPerm           permutations for the null (default 2000).
#' @param seed            RNG seed (default 42).
#' @param label           optional context label copied into the output column
#'                        \code{context}.
#' @return data.table: tf, context, nSites, obs, z, p, fdr (BH over TFs),
#'         meanNull, sdNull.
#' @author Irem B. GUNDUZ
#' @export
ascTFActivity <- function(refMat, altMat, deltaDt, sampleIds,
                          prefDelta = 1.0, minReads = 4L, minSites = 10L,
                          nPerm = 2000L, seed = 42L, label = NA_character_) {
  logger.start(paste0("Allele-specific TF activity",
                      if (!is.na(label)) paste0(" [", label, "]") else ""))

  stopifnot(all(c("snpId", "tf", "delta") %in% colnames(deltaDt)))
  set.seed(seed)

  sids <- intersect(sampleIds, colnames(refMat))
  if (length(sids) == 0) {
    logger.warning("No requested samples present in count matrices.")
    logger.completed(); return(NULL)
  }

  # Pool reads across the chosen samples once.
  refPool <- rowSums(refMat[, sids, drop = FALSE], na.rm = TRUE)
  altPool <- rowSums(altMat[, sids, drop = FALSE], na.rm = TRUE)

  dd <- data.table::as.data.table(deltaDt)[abs(delta) > prefDelta]
  dd <- dd[snpId %in% names(refPool)]
  if (nrow(dd) == 0) {
    logger.warning("No motif-disrupting sites present after filtering.")
    logger.completed(); return(NULL)
  }

  dd[, refN := refPool[snpId]]
  dd[, altN := altPool[snpId]]
  dd <- dd[(refN + altN) >= minReads]

  tfs <- sort(unique(dd$tf))
  rows <- vector("list", length(tfs))
  for (k in seq_along(tfs)) {
    sub <- dd[tf == tfs[k]]
    if (nrow(sub) < minSites) next
    or  <- ascOrientReads(sub$refN, sub$altN, sub$delta)
    st  <- ascActivityStat(or$H, or$L, nPerm = nPerm)
    rows[[k]] <- data.table::data.table(
      tf = tfs[k], context = label, nSites = st$nSites,
      obs = st$obs, z = st$z, p = st$p,
      meanNull = st$meanNull, sdNull = st$sdNull)
  }
  res <- data.table::rbindlist(rows)
  if (nrow(res) == 0) {
    logger.warning(paste0("No TF reached minSites=", minSites, "."))
    logger.completed(); return(NULL)
  }
  res[, fdr := stats::p.adjust(p, method = "BH")]
  data.table::setorder(res, -z)
  logger.info(paste0("Scored ", nrow(res), " TFs; ",
                     nrow(res[fdr < 0.1]), " significant at FDR<0.1."))
  logger.completed()
  res[]
}

#' Allele-specific TF activity across conditions, with a switch (handoff) score
#'
#' Runs \code{\link{ascTFActivity}} separately in two contexts (e.g. resting vs
#' stimulated) and joins them, adding the change in activity
#' \code{deltaObs = obs_B - obs_A} and a \code{flip} flag for TFs whose activity
#' changes sign between contexts. These are candidate "dynamic TF handoffs".
#'
#' @param refMat,altMat count matrices [snp x sample].
#' @param deltaDt       motif delta table (snpId, tf, delta).
#' @param samplesA,samplesB  sample vectors for the two contexts.
#' @param labelA,labelB      context labels (default "A"/"B").
#' @param ...           passed to \code{ascTFActivity} (prefDelta, minReads, ...).
#' @return data.table with per-TF activity in both contexts plus deltaObs, flip.
#' @author Irem B. GUNDUZ
#' @export
ascTFActivitySwitch <- function(refMat, altMat, deltaDt,
                                samplesA, samplesB,
                                labelA = "A", labelB = "B", ...) {
  a <- ascTFActivity(refMat, altMat, deltaDt, samplesA, label = labelA, ...)
  b <- ascTFActivity(refMat, altMat, deltaDt, samplesB, label = labelB, ...)
  if (is.null(a) || is.null(b)) return(NULL)
  m <- merge(a[, .(tf, obs_A = obs, z_A = z, p_A = p, fdr_A = fdr, nSites_A = nSites)],
             b[, .(tf, obs_B = obs, z_B = z, p_B = p, fdr_B = fdr, nSites_B = nSites)],
             by = "tf")
  m[, deltaObs := obs_B - obs_A]
  m[, flip := sign(obs_A) != sign(obs_B) &
        ((fdr_A < 0.1) | (fdr_B < 0.1))]
  data.table::setorder(m, -deltaObs)
  m[]
}

# ------------------------------------------------------------------------------
# TF family classification (for colouring)
# ------------------------------------------------------------------------------

#' Assign a broad structural family to a (JASPAR) TF name
#'
#' Coarse buckets used for colouring: "AP-1 / bZIP" (FOS/JUN/BATF/BACH/MAF/JDP/
#' ATF/NFE2, incl. dimers), "ZNF / GC" (ZNF/VEZF/RREB/GLIS/PLAGL/ZBTB/SP/KLF/EGR/
#' ZIC/WT1/CTCF), else "Other". Matching is case-insensitive and matches either
#' side of a dimer name (e.g. "FOSL1::JUND").
#'
#' @param tf character vector of TF names.
#' @return character vector of family labels.
#' @author Irem B. GUNDUZ
#' @export
ascTFFamily <- function(tf) {
  fam <- rep("Other", length(tf))
  ap1 <- grepl("FOS|JUN|BATF|BACH|MAF[FGK]?|JDP|ATF3|NFE2", tf, ignore.case = TRUE)
  znf <- grepl("ZNF|VEZF|RREB|GLIS|PLAGL|ZBTB|(^|:)SP[0-9]|KLF|EGR|ZIC|WT1|CTCF",
               tf, ignore.case = TRUE)
  fam[znf] <- "ZNF / GC"
  fam[ap1] <- "AP-1 / bZIP"   # AP-1 wins if a name somehow matches both
  fam
}

# ------------------------------------------------------------------------------
# Per-donor robustness
# ------------------------------------------------------------------------------

#' Allele-specific TF activity computed within each donor
#'
#' Repeats \code{\link{ascTFActivity}} separately for each donor's samples in a
#' given condition, so a pooled result can be checked for donor robustness (i.e.
#' that it is not driven by one deep donor). Per-donor coverage is lower, so some
#' TFs will drop below \code{minSites} in some donors -- that is expected and is
#' part of the robustness read-out.
#'
#' @param refMat,altMat count matrices [snp x sample].
#' @param deltaDt       motif delta table (snpId, tf, delta); subset to the TFs
#'                      of interest beforehand to keep this fast.
#' @param annot         sample annotation with \code{sampleId}, \code{donor} and
#'                      the condition column.
#' @param condition     value of \code{stimCol} to select (e.g. "U" or "S").
#' @param stimCol       name of the condition column (default "stimulus").
#' @param minDonorSamples minimum samples a donor must have in the condition
#'                      (default 1).
#' @param ...           passed to \code{ascTFActivity} (prefDelta, minReads,
#'                      minSites, nPerm, ...).
#' @return data.table with per-(TF, donor) activity plus \code{donor},
#'         \code{condition}; NULL if nothing computable.
#' @author Irem B. GUNDUZ
#' @export
ascTFActivityByDonor <- function(refMat, altMat, deltaDt, annot, condition,
                                 stimCol = "stimulus", minDonorSamples = 1L, ...) {
  logger.start(paste0("Per-donor TF activity [", condition, "]"))
  annot <- data.table::as.data.table(annot)
  donors <- unique(annot$donor)
  out <- list()
  for (d in donors) {
    sids <- annot[donor == d & get(stimCol) == condition, sampleId]
    if (length(sids) < minDonorSamples) next
    a <- ascTFActivity(refMat, altMat, deltaDt, sids,
                       label = paste0(condition, ":", d), ...)
    if (!is.null(a)) {
      a[, donor := d]; a[, condition := condition]
      out[[as.character(d)]] <- a
    }
  }
  logger.completed()
  if (length(out) == 0) return(NULL)
  data.table::rbindlist(out, fill = TRUE)
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

#' Lollipop / dot plot of allele-specific TF activity across contexts
#'
#' @param actDt   rbind of \code{ascTFActivity} outputs (needs tf, context, obs, z, fdr).
#' @param topN    keep the top-N TFs by |z| (default 25).
#' @param sigCut  FDR threshold used for the significance aesthetic (default 0.1).
#' @return a ggplot.
#' @author Irem B. GUNDUZ
#' @export
plotASCTFActivity <- function(actDt, topN = 25L, sigCut = 0.1) {
  requireNamespace("ggplot2")
  dt <- data.table::as.data.table(actDt)
  keep <- dt[, .(mz = max(abs(z), na.rm = TRUE)), by = tf][order(-mz)][seq_len(min(topN, .N))]$tf
  dt <- dt[tf %in% keep]
  dt[, tf := factor(tf, levels = rev(keep))]
  dt[, sig := fdr < sigCut]
  ggplot2::ggplot(dt, ggplot2::aes(x = obs, y = tf, colour = context)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = obs, yend = tf),
                          position = ggplot2::position_dodge(width = 0.6),
                          linewidth = 0.4, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = nSites, alpha = sig),
                        position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    ggplot2::labs(x = "Allele-specific TF activity  (net fraction toward predicted binder)",
                  y = NULL, colour = "Context", size = "ASC sites",
                  title = "Allele-specific transcription-factor activity") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Resting-vs-stimulated allele-specific TF activity scatter
#'
#' Each point is a TF: x = activity in resting, y = activity in stimulated.
#' Points are coloured by TF family (\code{\link{ascTFFamily}}) and faded when not
#' significant in either context. Points above the diagonal gain allele-specific
#' activity on stimulation (amplification); the two families typically separate
#' into opposite quadrants. (This replaces the earlier "handoff/flip" framing,
#' which the data did not support.)
#'
#' @param switchDt output of \code{ascTFActivitySwitch}.
#' @param topN     number of TFs (by |deltaObs|) to label (default 18).
#' @param sigCut   FDR threshold for the significance fade (default 0.1).
#' @return a ggplot.
#' @author Irem B. GUNDUZ
#' @export
plotASCTFActivityScatter <- function(switchDt, topN = 18L, sigCut = 0.1) {
  requireNamespace("ggplot2")
  dt <- data.table::as.data.table(switchDt)
  dt[, family := ascTFFamily(tf)]
  dt[, sig := (fdr_A < sigCut) | (fdr_B < sigCut)]
  dt[, lab := ""]
  ord <- order(-abs(dt$deltaObs))
  dt$lab[ord[seq_len(min(topN, nrow(dt)))]] <- dt$tf[ord[seq_len(min(topN, nrow(dt)))]]
  fam_cols <- c("AP-1 / bZIP" = "#C0392B", "ZNF / GC" = "#2C7FB8", "Other" = "grey70")
  ggplot2::ggplot(dt, ggplot2::aes(x = obs_A, y = obs_B, colour = family)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(ggplot2::aes(size = pmin(nSites_A, nSites_B), alpha = sig)) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.6, vjust = -0.8, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = fam_cols, name = "TF family") +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3),
                                name = paste0("FDR < ", sigCut)) +
    ggplot2::labs(x = "Allele-specific activity (resting)",
                  y = "Allele-specific activity (stimulated)",
                  size = "min ASC sites",
                  title = "Allele-specific TF activity: resting vs stimulated",
                  subtitle = "Above diagonal = amplified by stimulation; AP-1/bZIP (up) and ZNF/GC (down) separate") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Per-donor robustness of the stimulation effect
#'
#' Boxplot of per-donor Delta activity (stimulated - resting) for each TF, so you
#' can see the pooled amplification is consistent across donors rather than driven
#' by one deep donor. Each point is one donor.
#'
#' @param donorDeltaDt data.table with columns tf, donor, deltaObs (per-donor
#'                     stimulated - resting activity), e.g. built by joining the
#'                     two \code{\link{ascTFActivityByDonor}} outputs.
#' @param topN         number of TFs (by donor count) to show (default 20).
#' @return a ggplot.
#' @author Irem B. GUNDUZ
#' @export
plotASCTFDonorDelta <- function(donorDeltaDt, topN = 20L) {
  requireNamespace("ggplot2")
  dt <- data.table::as.data.table(donorDeltaDt)
  dt[, family := ascTFFamily(tf)]
  keep <- dt[, .N, by = tf][order(-N)][seq_len(min(topN, .N))]$tf
  dt <- dt[tf %in% keep]
  ord <- dt[, .(m = stats::median(deltaObs, na.rm = TRUE)), by = tf][order(m)]$tf
  dt[, tf := factor(tf, levels = ord)]
  fam_cols <- c("AP-1 / bZIP" = "#C0392B", "ZNF / GC" = "#2C7FB8", "Other" = "grey70")
  ggplot2::ggplot(dt, ggplot2::aes(x = deltaObs, y = tf, colour = family)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    ggplot2::geom_point(size = 0.9, alpha = 0.6,
                        position = ggplot2::position_jitter(height = 0.12, width = 0)) +
    ggplot2::scale_colour_manual(values = fam_cols, name = "TF family") +
    ggplot2::labs(x = "Per-donor change in activity (stimulated minus resting)", y = NULL,
                  title = "Per-donor robustness of stimulation amplification",
                  subtitle = "Each point = one donor; box = across-donor distribution") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
