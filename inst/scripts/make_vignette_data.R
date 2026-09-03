#!/usr/bin/env Rscript

#####################################################################
# make_vignette_data.R
#
# Builds the small example object used by vignettes/allelespec.Rmd from
# a completed AlleleSpeC run. Run it once, on the machine that holds the
# per-cell-type ds_asc_*.rds files; the result is committed to the
# package so the vignette knits during R CMD check without any download.
#
#   Rscript inst/scripts/make_vignette_data.R
#   ASC_DIR=/path/to/asc_counts Rscript inst/scripts/make_vignette_data.R
#
# Writes a single list to inst/extdata/dsASC_example.rds:
#   $ds       DsASC restricted to a few cell types, one chromosome, a
#             few thousand sites, resting and stimulated samples only
#   $peakGr   consensus peaks overlapping those sites
#   $peakMat  peak accessibility for those peaks and samples
#   $delta    per-site PWM delta table (snpId, tf, delta)
#
# Thymic populations are excluded throughout: they are the shallowest
# libraries in the dataset, several of them produce _EMPTY.rds at the
# counting stage, and none of them survive the coverage filter in
# numbers that would make a useful example.
#
# The point of the subsetting is size. The full merged object is
# hundreds of megabytes and has no business in a git repository; this
# script exists so that the example in the package is reproducible from
# the full run rather than being an unexplained binary.
#####################################################################

suppressPackageStartupMessages({
  library(ChrAccR)
  library(muLogR)
  library(GenomicRanges)
  library(data.table)
})

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
ASC_DIR  <- Sys.getenv("ASC_DIR",  "/icbb_triton/scratch/igunduz/allelSpec/asc_counts")
OUT_FILE <- Sys.getenv("OUT_FILE", file.path("inst", "extdata", "dsASC_example.rds"))

EXCLUDE_CELLTYPES <- c("^Thy", "^mTEC$", "^TpreDN$")
KEEP_CELLTYPES    <- c("TeffNaive", "Bmem", "Mono")
KEEP_STIMULI      <- c("U", "S")
KEEP_CHROM        <- "chr20"

MAX_SITES    <- 4000L
MIN_COVERAGE <- 10L
MIN_SAMPLES  <- 5L

TF_NAMES <- c("JUNB", "BACH1", "FOSL1", "CTCF", "SPI1",
              "RUNX1", "NFKB1", "IRF4", "GATA3", "TCF7")

SIZE_BUDGET_MB <- 1
SEED <- 42L

set.seed(SEED)
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

if (!"DsASC" %in% getNamespaceExports("ChrAccR")) {
  stop("The installed ChrAccR does not export DsASC. Install this working copy first: ",
       "devtools::install('.')")
}

# The stored objects were written by scripts that source() the class definitions
# on top of the package, so their class can be stamped .GlobalEnv rather than
# ChrAccR and S4 dispatch will not find the DsAcc methods. Slots are read
# directly here so the script works regardless of how the objects were created.
ascSlots <- function(ds) {
  list(annot = data.table::as.data.table(ds@sampleAnnot),
       snps  = ds@coord$snps,
       peaks = ds@coord$peaks,
       acc   = ds@accessibility)
}
ascCounts <- function(ds, type, i, j) {
  m <- as.matrix(ds@counts[[type]][i, j, drop = FALSE])
  m[is.na(m)] <- 0L
  m
}

# ---------------------------------------------------------------------
# 1. Pick the input objects
# ---------------------------------------------------------------------
logger.start("Selecting per-cell-type ASC objects")

files <- list.files(ASC_DIR, pattern = "^ds_asc_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_EMPTY\\.rds$", files)]
if (length(files) == 0) stop("No ds_asc_*.rds found in ", ASC_DIR)

cellTypes <- sub("^ds_asc_(.*)\\.rds$", "\\1", basename(files))
logger.info(paste("Available:", paste(sort(cellTypes), collapse = ", ")))

isThymic <- Reduce(`|`, lapply(EXCLUDE_CELLTYPES, function(p) grepl(p, cellTypes)))
if (any(isThymic)) {
  logger.info(paste("Excluding thymic populations:",
                    paste(cellTypes[isThymic], collapse = ", ")))
  files <- files[!isThymic]; cellTypes <- cellTypes[!isThymic]
}

if (!is.null(KEEP_CELLTYPES)) {
  sel <- cellTypes %in% KEEP_CELLTYPES
  missing <- setdiff(KEEP_CELLTYPES, cellTypes)
  if (length(missing) > 0) logger.warning(paste("Requested but absent:",
                                                paste(missing, collapse = ", ")))
  files <- files[sel]; cellTypes <- cellTypes[sel]
}
if (length(files) == 0) stop("No cell types left after filtering.")
logger.info(paste("Using:", paste(cellTypes, collapse = ", ")))
logger.completed()

# ---------------------------------------------------------------------
# 2. Read, restricting rows and columns before anything is stacked
# ---------------------------------------------------------------------
logger.start(paste("Reading objects, restricted to", KEEP_CHROM))

refL <- list(); altL <- list(); annL <- list()
coordAll <- NULL; peaksAll <- NULL; accAll <- NULL

for (i in seq_along(files)) {
  ds <- readRDS(files[i])
  sl <- ascSlots(ds)

  keepCols <- which(sl$annot$stimulus %in% KEEP_STIMULI)
  if (length(keepCols) == 0) { rm(ds, sl); next }

  keepRows <- which(as.character(seqnames(sl$snps)) %in% KEEP_CHROM)
  if (length(keepRows) == 0) { rm(ds, sl); next }

  r <- ascCounts(ds, "ref", keepRows, keepCols)
  a <- ascCounts(ds, "alt", keepRows, keepCols)

  refL[[cellTypes[i]]] <- r
  altL[[cellTypes[i]]] <- a
  annL[[cellTypes[i]]] <- sl$annot[keepCols]

  newSnps <- sl$snps[keepRows]
  coordAll <- if (is.null(coordAll)) newSnps else
    c(coordAll, newSnps[!names(newSnps) %in% names(coordAll)])

  if (is.null(peaksAll) && !is.null(sl$peaks)) peaksAll <- sl$peaks
  if (length(sl$acc) > 0 && nrow(sl$acc) > 0) {
    acc <- as.matrix(sl$acc)
    accAll <- if (is.null(accAll)) acc else
      cbind(accAll, acc[rownames(accAll), setdiff(colnames(acc), colnames(accAll)), drop = FALSE])
  }

  logger.info(paste0(cellTypes[i], ": ", nrow(r), " sites x ", ncol(r), " samples"))
  rm(ds, sl, r, a); gc(verbose = FALSE)
}
if (length(refL) == 0) stop("Nothing survived the row/column restriction.")
logger.completed()

# ---------------------------------------------------------------------
# 3. Align into one pair of matrices
# ---------------------------------------------------------------------
logger.start("Aligning cell types onto a common site set")

allSites   <- Reduce(union, lapply(refL, rownames))
allSamples <- unlist(lapply(refL, colnames), use.names = FALSE)

REF <- matrix(0L, length(allSites), length(allSamples),
              dimnames = list(allSites, allSamples))
ALT <- REF
for (nm in names(refL)) {
  REF[rownames(refL[[nm]]), colnames(refL[[nm]])] <- refL[[nm]]
  ALT[rownames(altL[[nm]]), colnames(altL[[nm]])] <- altL[[nm]]
}
rm(refL, altL); gc(verbose = FALSE)
logger.info(paste(nrow(REF), "sites x", ncol(REF), "samples before site selection"))
logger.completed()

# ---------------------------------------------------------------------
# 4. Choose the sites to keep
# ---------------------------------------------------------------------
logger.start("Selecting example sites")

covered <- rowSums((REF + ALT) >= MIN_COVERAGE) >= MIN_SAMPLES
candidates <- rownames(REF)[covered]
logger.info(paste(length(candidates), "sites covered >=", MIN_COVERAGE,
                  "in >=", MIN_SAMPLES, "samples"))
if (length(candidates) == 0) stop("No sites meet the coverage requirement.")

chosen <- if (length(candidates) > MAX_SITES) sample(candidates, MAX_SITES) else candidates

siteGr <- coordAll[names(coordAll) %in% chosen]
siteGr <- sort(siteGr)
chosen <- names(siteGr)

REF <- REF[chosen, , drop = FALSE]
ALT <- ALT[chosen, , drop = FALSE]
logger.info(paste("Keeping", length(chosen), "sites"))
logger.completed()

# ---------------------------------------------------------------------
# 5. Assemble the DsASC
# ---------------------------------------------------------------------
logger.start("Assembling the example DsASC")

annotCols <- intersect(c("sampleId", "donor", "cellType", "stimulus", "replicate"),
                       names(rbindlist(annL, fill = TRUE)))
annot <- as.data.frame(rbindlist(annL, fill = TRUE)[, ..annotCols])
annot$donor <- as.character(annot$donor)
rownames(annot) <- annot$sampleId
annot <- annot[colnames(REF), , drop = FALSE]

dsSub <- DsASC(annot, genome = "hg38")
dsSub@coord  <- list(snps = siteGr)
dsSub@counts <- list(ref = REF, alt = ALT)
dsSub@diskDump <- FALSE
print(dsSub)
logger.completed()

# ---------------------------------------------------------------------
# 6. Peaks for the sites that were kept
# ---------------------------------------------------------------------
logger.start("Subsetting peaks and accessibility")

peakGr <- NULL; peakMat <- NULL
if (!is.null(peaksAll) && !is.null(siteGr$peakId)) {
  peakIds <- unique(stats::na.omit(siteGr$peakId))
  peakIds <- intersect(peakIds, names(peaksAll))
  peakGr  <- peaksAll[peakIds]
  if (!is.null(accAll)) {
    shared  <- intersect(colnames(REF), colnames(accAll))
    peakMat <- accAll[peakIds, shared, drop = FALSE]
    peakMat <- peakMat[, colnames(REF)[colnames(REF) %in% shared], drop = FALSE]
  }
  logger.info(paste("Kept", length(peakGr), "peaks"))
} else {
  logger.warning("No peak annotation on the source objects; skipping peaks.")
}
logger.completed()

# ---------------------------------------------------------------------
# 7. Motif deltas for the same sites
# ---------------------------------------------------------------------
logger.start("Scoring allelic motif disruption")

motifPkgs <- c("JASPAR2020", "TFBSTools", "motifmatchr", "BSgenome.Hsapiens.UCSC.hg38")
missingPkgs <- motifPkgs[!vapply(motifPkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missingPkgs) > 0) {
  stop("The vignette needs the motif delta table, so these packages are required: ",
       paste(missingPkgs, collapse = ", "))
}

{
  suppressPackageStartupMessages({
    library(JASPAR2020); library(TFBSTools); library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
  })

  pfms <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
  nms  <- vapply(pfms, name, "")
  pfms <- pfms[toupper(nms) %in% toupper(TF_NAMES)]
  if (length(pfms) == 0) stop("None of TF_NAMES found in JASPAR2020 CORE.")

  flank  <- max(vapply(pfms, function(p) ncol(as.matrix(p)), 1L))
  win    <- GRanges(as.character(seqnames(siteGr)),
                    IRanges(start(siteGr) - flank, start(siteGr) + flank))
  seqRef <- getSeq(BSgenome.Hsapiens.UCSC.hg38, win)
  seqAlt <- replaceLetterAt(
    seqRef,
    at = matrix(seq_len(width(seqRef)[1]) == (flank + 1L),
                nrow = length(seqRef), ncol = width(seqRef)[1], byrow = TRUE),
    letter = siteGr$ALT)

  scoreSet <- function(s) as.matrix(motifScores(
    matchMotifs(pfms, s, out = "scores", bg = "even", p.cutoff = 1)))

  delta <- scoreSet(seqAlt) - scoreSet(seqRef)
  colnames(delta) <- vapply(pfms, name, "")

  deltaDt <- rbindlist(lapply(colnames(delta), function(tfName)
    data.table(snpId = names(siteGr), tf = tfName,
               delta = round(as.numeric(delta[, tfName]), 4))))
  logger.info(paste(nrow(deltaDt), "site-motif deltas for",
                    length(unique(deltaDt$tf)), "factors"))
}
logger.completed()

# ---------------------------------------------------------------------
# 8. Save and check the size
# ---------------------------------------------------------------------
logger.start("Writing example object")

saveRDS(list(ds = dsSub, peakGr = peakGr, peakMat = peakMat, delta = deltaDt),
        OUT_FILE, compress = "xz")

sizeMb <- file.size(OUT_FILE) / 1024^2
logger.info(paste0(OUT_FILE, ": ", round(sizeMb, 2), " MB"))
if (sizeMb > SIZE_BUDGET_MB) {
  logger.warning(paste0(
    "Over the ", SIZE_BUDGET_MB, " MB budget. Lower MAX_SITES, drop a cell type ",
    "from KEEP_CELLTYPES, or restrict KEEP_STIMULI, and rerun."))
} else {
  logger.info("Within budget; safe to commit.")
}
logger.completed()
