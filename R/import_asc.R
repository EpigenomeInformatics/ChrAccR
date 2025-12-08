#' Import Allele-Specific Chromatin Data (GATK VCFs + BAM Pileup)
#'
#' This function imports ASC data by counting reads at heterozygous SNP sites
#' defined in donor-specific VCF files. It uses \code{Rsamtools::pileup} to
#' quantify reference and alternative alleles.
#' 
#' @importFrom DelayedArray close
#' @import HDF5Array
#' @param sampleAnnot  Data frame with sample annotation. Must contain 'sampleId', 'bamFilename', 'donor'.
#' @param vcfDir       Directory containing donor VCFs named \code{\{donor\}_heterozygous.vcf.gz}.
#' @param genome       Character string containing genome assembly (e.g. "hg38").
#' @param diskDump     Logical. If TRUE, matrices are realized as HDF5 arrays.
#'
#' @return A \code{\linkS4class{DsASC}} object.
#' @export
DsASC.gatk <- function(sampleAnnot, vcfDir, genome, diskDump = FALSE) {
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

  # 2. Master SNP List
  logger.start("Constructing Master SNP List from VCFs")

  donors <- unique(sampleAnnot$donor)
  grList <- list()

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
        grList[[d]] <- gr
      }
    } else {
      logger.warning(c("VCF missing for donor:", d))
    }
  }
  # Remove NULLs
  if (length(grList) == 0) {
    logger.error("No valid VCFs found.")
    stop("Missing VCFs")
  }
  grList <- grList[!sapply(grList, is.null)]
  
  if (length(grList) == 0) {
    muLogR::logger.error("VCFs filtered to zero SNPs for all donors in this cell type. Cannot proceed.")
    stop("VCFs filtered to zero SNPs.")
  }

  masterGr <- unique(unlist(GenomicRanges::GRangesList(grList)))
  masterGr <- sort(masterGr)
  masterGr$snpId <- paste0(GenomicRanges::seqnames(masterGr), ":", GenomicRanges::start(masterGr), "_", masterGr$REF, "_", masterGr$ALT)
  names(masterGr) <- masterGr$snpId

  nSites <- as.integer(length(masterGr))
  nSamples <- as.integer(nrow(sampleAnnot))

  logger.info(c("Identified", nSites, "unique heterozygous SNPs across", length(grList), "donors."))
  logger.completed()

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

    logger.status(c("Processing sample", paste0(i, "/", nSamples), ":", sid))

    vecRef <- integer(nSites)
    vecAlt <- integer(nSites)

    if (file.exists(bam)) {
      bamHeader <- Rsamtools::scanBamHeader(bam)
      commonSeq <- intersect(GenomeInfoDb::seqlevels(masterGr), names(bamHeader[[1]]$targets))

      if (length(commonSeq) > 0) {
        targetGr <- GenomeInfoDb::keepSeqlevels(masterGr, commonSeq, pruning.mode = "coarse")
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
    } else {
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
