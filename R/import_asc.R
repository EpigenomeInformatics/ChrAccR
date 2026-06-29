#' @include DsASC-class.R
NULL

#' Build a shared master heterozygous-SNP list across donors
#'
#' Reads donor VCFs once and returns a sorted, de-duplicated GRanges of
#' biallelic heterozygous SNPs spanning all requested donors. The returned
#' object carries, per site, the set of donors heterozygous there
#' (\code{$donorMembership}, a CharacterList), so the donor mask can be applied
#' later without re-reading VCFs. Save this once and pass it to
#' \code{DsASC.gatk(..., masterGr = )} in each per-cell-type job to keep all
#' objects row-aligned.
#'
#' @param donors  Character vector of donor IDs.
#' @param vcfDir  Directory with \code{\{donor\}_heterozygous.vcf.gz} files.
#' @param genome  Genome assembly string (e.g. "hg38").
#' @return A GRanges with REF, ALT, snpId, donorMembership metadata; names() = snpId.
#' @export
buildMasterSNPs <- function(donors, vcfDir, genome) {
  logger.start("Constructing Master SNP List from VCFs")

  grList      <- list()
  donorSnpIds <- list()

  for (d in donors) {
    vcfFile <- file.path(vcfDir, paste0(d, "_heterozygous.vcf.gz"))
    if (file.exists(vcfFile)) {
      vcfParam <- VariantAnnotation::ScanVcfParam(fixed = c("ALT"), info = NA, geno = NA)
      vcf <- suppressWarnings(VariantAnnotation::readVcf(vcfFile, genome, param = vcfParam))
      vcf <- vcf[VariantAnnotation::isSNV(vcf) & elementNROWS(VariantAnnotation::alt(vcf)) == 1]

      if (length(vcf) > 0) {
        gr <- SummarizedExperiment::rowRanges(vcf)
        GenomicRanges::mcols(gr)$REF <- as.character(VariantAnnotation::ref(vcf))
        GenomicRanges::mcols(gr)$ALT <- as.character(unlist(VariantAnnotation::alt(vcf)))
        gr$snpId <- paste0(GenomicRanges::seqnames(gr), ":", GenomicRanges::start(gr),
                           "_", gr$REF, "_", gr$ALT)
        donorSnpIds[[d]] <- gr$snpId
        grList[[d]]      <- gr
      }
    } else {
      logger.warning(c("VCF missing for donor:", d))
    }
  }

  grList <- grList[!sapply(grList, is.null)]
  if (length(grList) == 0) {
    logger.error("No valid VCFs found / all filtered to zero SNPs.")
    stop("Missing VCFs")
  }

  masterGr <- unique(unlist(GenomicRanges::GRangesList(grList)))
  masterGr <- sort(masterGr)
  masterGr$snpId <- paste0(GenomicRanges::seqnames(masterGr), ":",
                           GenomicRanges::start(masterGr), "_",
                           masterGr$REF, "_", masterGr$ALT)
  names(masterGr) <- masterGr$snpId

  # Per-site donor membership: which donors are het at each master site.
  # Match on snpId STRING via match() -- this is the same mapping the counting
  # loop uses successfully, and both masterGr and each donor gr build snpId
  # identically (chr:pos_REF_ALT), so they align exactly. findOverlaps(type=
  # "equal") was unreliable here after unlist/sort changed seqlevels ordering.
  masterIds <- names(masterGr)
  qhAll <- integer(0)   # master row indices
  dvAll <- character(0) # donor for that match

  for (d in names(grList)) {
    ids <- grList[[d]]$snpId
    rows <- match(ids, masterIds)        # master row for each donor het site
    rows <- rows[!is.na(rows)]
    qhAll <- c(qhAll, rows)
    dvAll <- c(dvAll, rep(d, length(rows)))
  }

  if (length(qhAll) == 0) {
    logger.error("Donor membership is EMPTY: no donor snpId matched the master.")
    logger.error("This means snpId construction differs between master and donor VCFs.")
    stop("Empty donor membership; aborting before producing all-NA counts.")
  }

  memb <- IRanges::CharacterList(rep(list(character(0)), length(masterGr)))
  filled <- S4Vectors::splitAsList(dvAll, factor(qhAll, levels = seq_along(masterGr)))
  memb <- IRanges::CharacterList(lapply(filled, unique))
  masterGr$donorMembership <- memb

  nOrphan <- sum(S4Vectors::elementNROWS(masterGr$donorMembership) == 0)
  if (nOrphan > 0) {
    logger.warning(c(nOrphan, "of", length(masterGr),
                     "master sites matched no donor."))
  }
  logger.info(c("Membership built: ",
                length(masterGr) - nOrphan, "of", length(masterGr),
                "sites have >=1 donor."))

  logger.info(c("Identified", length(masterGr),
                "unique heterozygous SNPs across", length(grList), "donors."))
  logger.completed()
  return(masterGr)
}

#' Import Allele-Specific Chromatin Data (GATK VCFs + BAM Pileup)
#'
#' This function imports ASC data by counting reads at heterozygous SNP sites
#' defined in donor-specific VCF files. It uses \code{Rsamtools::pileup} to
#' quantify reference and alternative alleles.
#'
#' Counts are stored against a union ("master") SNP list spanning all donors,
#' but each sample is counted ONLY at sites heterozygous in its own donor.
#' Sites outside a sample's donor het set are set to NA so they are never tested
#' (prevents homozygous genotypes from being mis-called as allele-specific).
#'
#' @importFrom DelayedArray close
#' @import HDF5Array
#' @param sampleAnnot  Data frame with sample annotation. Must contain 'sampleId', 'bamFilename', 'donor'.
#' @param vcfDir       Directory containing donor VCFs named \code{\{donor\}_heterozygous.vcf.gz}.
#' @param genome       Character string containing genome assembly (e.g. "hg38").
#' @param diskDump     Logical. If TRUE, matrices are realized as HDF5 arrays.
#'
#' @return A \code{\linkS4class{DsASC}} object.
#' @param masterGr   Optional pre-built master SNP GRanges (from
#'                   \code{buildMasterSNPs}). If supplied, the per-VCF master
#'                   construction is skipped and counts are stored against this
#'                   shared site set. Use this so per-cell-type objects built in
#'                   separate jobs remain row-aligned for cross-cell-type
#'                   comparison (Fig 4d). Must carry REF, ALT, snpId metadata and
#'                   names() set to snpId (as buildMasterSNPs produces).
#' @export
DsASC.gatk <- function(sampleAnnot, vcfDir, genome, diskDump = FALSE, masterGr = NULL) {
  # Validation & Setup
  reqCols <- c("sampleId", "bamFilename", "donor")
  if (!all(reqCols %in% colnames(sampleAnnot))) {
    logger.error(c("sampleAnnot must contain columns:", paste(reqCols, collapse = ", ")))
    stop("Invalid annotation columns")
  }

  if (any(duplicated(sampleAnnot$sampleId))) {
    logger.error("Duplicate sampleIds found in annotation.")
    stop("Duplicate sampleIds")
  }
  rownames(sampleAnnot) <- sampleAnnot$sampleId

  # Config Defaults
  minMapq <- getConfigElement("minMapq")
  minBaseq <- getConfigElement("minBaseq")
  maxDepth <- getConfigElement("maxDepth")

  # Check Backend
  be <- DelayedArray::getAutoRealizationBackend()
  if (isTRUE(diskDump) || (!is.null(be) && be == "HDF5Array")) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
      logger.error("HDF5Array package required for HDF5 backend")
      stop("Missing package: HDF5Array")
    }
    logger.info("Using HDF5Array backend for disk dumping")
  }

  # 2. Master SNP List  (build, or reuse a supplied one)
  if (is.null(masterGr)) {
    masterGr <- buildMasterSNPs(unique(sampleAnnot$donor), vcfDir, genome)
  } else {
    logger.info("Using supplied master SNP list (shared across jobs).")
    if (is.null(masterGr$donorMembership)) {
      logger.error("Supplied masterGr lacks $donorMembership; rebuild with buildMasterSNPs().")
      stop("Invalid masterGr: missing donorMembership metadata.")
    }
  }

  nSites   <- as.integer(length(masterGr))
  nSamples <- as.integer(nrow(sampleAnnot))

  # [DONOR-MASK] per-donor row indices, derived from membership carried in masterGr.
  # Only donors present in THIS job's annotation are needed.
  # NOTE: force character keys throughout. If donor is integer, list[[donor]]
  # would be interpreted as a POSITIONAL index -> "subscript out of bounds".
  membership  <- masterGr$donorMembership   # CharacterList: donors het at each site
  jobDonors   <- as.character(unique(sampleAnnot$donor))
  memb_unl    <- as.character(unlist(membership, use.names = FALSE))
  memb_rowmap <- rep(seq_along(membership), times = S4Vectors::elementNROWS(membership))
  donorRowIdx <- lapply(jobDonors, function(d) unique(memb_rowmap[memb_unl == d]))
  names(donorRowIdx) <- jobDonors

  # diagnostic: warn if any job donor is absent from the master membership
  missingDonors <- setdiff(jobDonors, as.character(memb_unl))
  if (length(missingDonors) > 0) {
    logger.warning(c("Donors with NO het sites in master (will get empty mask):",
                     paste(missingDonors, collapse = ", ")))
  }

  logger.info("Per-donor het-site membership resolved from master (donor mask).")
  logger.info(c("Master SNP list:", nSites, "sites; job covers", length(jobDonors), "donor(s)."))

  # 3. Initialize Data Structures
  if (diskDump) {
    logger.info("Initializing HDF5 sinks...")

    sinkRef <- HDF5Array::HDF5RealizationSink(
      dim = c(nSites, nSamples),
      dimnames = list(names(masterGr), sampleAnnot$sampleId),
      type = "integer",
      name = "refCounts",
      level = 6
    )

    sinkAlt <- HDF5Array::HDF5RealizationSink(
      dim = c(nSites, nSamples),
      dimnames = list(names(masterGr), sampleAnnot$sampleId),
      type = "integer",
      name = "altCounts",
      level = 6
    )

    # Define Grid: Spacings = (All Rows, 1 Column)
    chunkGrid <- DelayedArray::RegularArrayGrid(
      refdim = c(nSites, nSamples),
      spacings = c(nSites, 1L)
    )
  } else {
    matRef <- matrix(0L, nrow = nSites, ncol = nSamples, dimnames = list(names(masterGr), sampleAnnot$sampleId))
    matAlt <- matrix(0L, nrow = nSites, ncol = nSamples, dimnames = list(names(masterGr), sampleAnnot$sampleId))
  }

  # 4. Process Samples
  logger.start(c("Counting alleles for", nSamples, "samples"))

  pParam <- Rsamtools::PileupParam(
    max_depth = maxDepth,
    min_base_quality = minBaseq,
    min_mapq = minMapq,
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    include_deletions = FALSE,
    include_insertions = FALSE
  )

  for (i in seq_len(nSamples)) {
    sid <- sampleAnnot$sampleId[i]
    bam <- sampleAnnot$bamFilename[i]
    donor_i <- as.character(sampleAnnot$donor[i])   # [DONOR-MASK] character key, not positional

    logger.status(c("Processing sample", paste0(i, "/", nSamples), ":", sid))

    # [DONOR-MASK] Initialize as NA, not 0. Only this donor's het sites will be
    # filled with real counts (incl. genuine zeros); everything else stays NA
    # and is excluded from testing downstream.
    vecRef <- rep(NA_integer_, nSites)
    vecAlt <- rep(NA_integer_, nSites)

    keepRows <- donorRowIdx[[donor_i]]
    if (is.null(keepRows)) keepRows <- integer(0)

    # Real (possibly zero) counts only at this donor's het sites
    if (length(keepRows) > 0) {
      vecRef[keepRows] <- 0L
      vecAlt[keepRows] <- 0L
    }

    if (file.exists(bam) && length(keepRows) > 0) {
      # [DONOR-MASK] restrict pileup targets to THIS donor's het sites
      donorGr <- masterGr[keepRows]

      bamHeader <- Rsamtools::scanBamHeader(bam)
      commonSeq <- intersect(GenomeInfoDb::seqlevels(donorGr), names(bamHeader[[1]]$targets))

      if (length(commonSeq) > 0) {
        targetGr <- GenomeInfoDb::keepSeqlevels(donorGr, commonSeq, pruning.mode = "coarse")
        sbP <- Rsamtools::ScanBamParam(which = targetGr, flag = Rsamtools::scanBamFlag(isDuplicate = FALSE))

        pileupRes <- Rsamtools::pileup(bam, scanBamParam = sbP, pileupParam = pParam)

        if (nrow(pileupRes) > 0) {
          dt <- data.table::as.data.table(pileupRes)
          counts <- data.table::dcast(dt, seqnames + pos ~ nucleotide, value.var = "count", fill = 0)

          targetDt <- data.table::as.data.table(data.frame(
            seqnames = as.character(GenomicRanges::seqnames(targetGr)),
            pos = GenomicRanges::start(targetGr),
            REF = targetGr$REF,
            ALT = targetGr$ALT,
            id = targetGr$snpId,
            stringsAsFactors = FALSE
          ))

          merged <- merge(targetDt, counts, by = c("seqnames", "pos"), all.x = TRUE)

          for (nuc in c("A", "C", "G", "T")) if (!nuc %in% names(merged)) merged[[nuc]] <- 0L
          merged[is.na(merged)] <- 0L

          merged[, r_c := 0L]
          merged[, a_c := 0L]
          for (nuc in c("A", "C", "G", "T")) {
            merged[REF == nuc, r_c := get(nuc)]
            merged[ALT == nuc, a_c := get(nuc)]
          }

          idx <- match(merged$id, names(masterGr))
          validIdx <- !is.na(idx)
          if (any(validIdx)) {
            vecRef[idx[validIdx]] <- merged$r_c[validIdx]
            vecAlt[idx[validIdx]] <- merged$a_c[validIdx]
          }
        }
      }
    } else if (!file.exists(bam)) {
      logger.warning(c("BAM file not found:", bam))
    }

    # Write Block (Column)
    if (diskDump) {
      # Format block as matrix (nSites x 1)
      blockRef <- matrix(vecRef, nrow = nSites, ncol = 1L)
      blockAlt <- matrix(vecAlt, nrow = nSites, ncol = 1L)

      # Use linear indexing [[i]] for RegularArrayGrid
      DelayedArray::write_block(sinkRef, chunkGrid[[i]], blockRef)
      DelayedArray::write_block(sinkAlt, chunkGrid[[i]], blockAlt)

      rm(blockRef, blockAlt, vecRef, vecAlt)
      gc()
    } else {
      matRef[, i] <- vecRef
      matAlt[, i] <- vecAlt
    }
  }

  logger.completed()

  # 5. Populate DsASC Object
  logger.start("Creating DsASC object")

  if (diskDump) {
    DelayedArray::close(sinkRef)
    DelayedArray::close(sinkAlt)

    countsList <- list(
      ref = as(sinkRef, "HDF5Array"),
      alt = as(sinkAlt, "HDF5Array")
    )
  } else {
    countsList <- list(
      ref = matRef,
      alt = matAlt
    )
  }
  # Create DsASC Object
  dsObj <- DsASC(sampleAnnot, genome, diskDump)
  dsObj@coord <- list(snps = masterGr)
  dsObj@counts <- countsList

  logger.completed()

  return(dsObj)
}