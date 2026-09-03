# R/AllGenerics.R
# Centralized generic definitions to guarantee correct load order and prevent overwriting

################################################################################
# Base DsAcc Generics
################################################################################
setGeneric("getSamples", function(.object) standardGeneric("getSamples"), signature = c(".object"))
setGeneric("getSampleAnnot", function(.object) standardGeneric("getSampleAnnot"), signature=c(".object"))
setGeneric("getGenome", function(.object) standardGeneric("getGenome"), signature=c(".object"))
setGeneric("getRegionTypes", function(.object, ...) standardGeneric("getRegionTypes"), signature=c(".object"))
setGeneric("getCoord", function(.object, ...) standardGeneric("getCoord"), signature=c(".object"))
setGeneric("getNRegions", function(.object, ...) standardGeneric("getNRegions"), signature=c(".object"))
setGeneric("addSampleAnnotCol", function(.object, ...) standardGeneric("addSampleAnnotCol"), signature=c(".object"))
setGeneric("removeRegions", function(.object, ...) standardGeneric("removeRegions"), signature=c(".object"))
setGeneric("getComparisonTable", function(.object, ...) standardGeneric("getComparisonTable"), signature=c(".object"))

################################################################################
# DsATAC Generics
################################################################################
setGeneric("getCounts", function(.object, ...) standardGeneric("getCounts"), signature=c(".object"))
setGeneric("getCountsSE", function(.object, ...) standardGeneric("getCountsSE"), signature=c(".object"))
setGeneric("getFragmentGr", function(.object, ...) standardGeneric("getFragmentGr"), signature=c(".object"))
setGeneric("getFragmentGrl", function(.object, ...) standardGeneric("getFragmentGrl"), signature=c(".object"))
setGeneric("getFragmentNum", function(.object, ...) standardGeneric("getFragmentNum"), signature=c(".object"))
setGeneric("getInsertionSites", function(.object, ...) standardGeneric("getInsertionSites"), signature=c(".object"))
setGeneric("getCoverage", function(.object, ...) standardGeneric("getCoverage"), signature=c(".object"))
setGeneric("regionAggregation", function(.object, ...) standardGeneric("regionAggregation"), signature=c(".object"))
setGeneric("mergeSamples", function(.object, ...) standardGeneric("mergeSamples"), signature=c(".object"))
setGeneric("join", function(.object, ...) standardGeneric("join"), signature=c(".object"))
setGeneric("removeFragmentData", function(object) standardGeneric("removeFragmentData"), signature=c("object"))
setGeneric("undiskFragmentData", function(object) standardGeneric("undiskFragmentData"), signature=c("object"))
setGeneric("subsampleFragmentData", function(object, ...) standardGeneric("subsampleFragmentData"), signature=c("object"))
setGeneric("addCountDataFromBam", function(.object, ...) standardGeneric("addCountDataFromBam"), signature=c(".object"))
setGeneric("addCountDataFromGRL", function(.object, ...) standardGeneric("addCountDataFromGRL"), signature=c(".object"))
setGeneric("addSignalDataFromGRL", function(.object, ...) standardGeneric("addSignalDataFromGRL"), signature=c(".object"))
setGeneric("addInsertionDataFromBam", function(.object, ...) standardGeneric("addInsertionDataFromBam"), signature=c(".object"))
setGeneric("removeRegionType", function(.object, ...) standardGeneric("removeRegionType"), signature=c(".object"))
setGeneric("removeRegionData", function(.object) standardGeneric("removeRegionData"), signature=c(".object"))
setGeneric("removeSamples", function(.object, ...) standardGeneric("removeSamples"), signature=c(".object"))
setGeneric("transformCounts", function(.object, ...) standardGeneric("transformCounts"), signature=c(".object"))
setGeneric("filterLowCovg", function(.object, ...) standardGeneric("filterLowCovg"), signature=c(".object"))
setGeneric("filterChroms", function(.object, ...) standardGeneric("filterChroms"), signature=c(".object"))
setGeneric("filterByGRanges", function(.object, ...) standardGeneric("filterByGRanges"), signature=c(".object"))
setGeneric("regionSetCounts", function(.object, ...) standardGeneric("regionSetCounts"), signature=c(".object"))
setGeneric("getInsertionKmerFreq", function(.object, ...) standardGeneric("getInsertionKmerFreq"), signature=c(".object"))
setGeneric("aggregateRegionCounts", function(.object, ...) standardGeneric("aggregateRegionCounts"), signature=c(".object"))
setGeneric("getMotifEnrichment", function(.object, ...) standardGeneric("getMotifEnrichment"), signature=c(".object"))
setGeneric("getChromVarDev", function(.object, ...) standardGeneric("getChromVarDev"), signature=c(".object"))
setGeneric("getMotifFootprints", function(.object, ...) standardGeneric("getMotifFootprints"), signature=c(".object"))
setGeneric("getDESeq2Dataset", function(.object, ...) standardGeneric("getDESeq2Dataset"), signature=c(".object"))
setGeneric("getDiffAcc", function(.object, ...) standardGeneric("getDiffAcc"), signature=c(".object"))
setGeneric("exportCountTracks", function(.object, ...) standardGeneric("exportCountTracks"), signature=c(".object"))
setGeneric("callPeaks", function(.object, ...) standardGeneric("callPeaks"), signature=c(".object"))
setGeneric("plotInsertSizeDistribution", function(.object, ...) standardGeneric("plotInsertSizeDistribution"), signature=c(".object"))
setGeneric("getTssEnrichment", function(.object, ...) standardGeneric("getTssEnrichment"), signature=c(".object"))
setGeneric("getTssEnrichmentBatch", function(.object, ...) standardGeneric("getTssEnrichmentBatch"), signature=c(".object"))
setGeneric("getQuickTssEnrichment", function(.object, ...) standardGeneric("getQuickTssEnrichment"), signature=c(".object"))
setGeneric("getMonocleCellDataSet", function(.object, ...) standardGeneric("getMonocleCellDataSet"), signature=c(".object"))
setGeneric("getCiceroGeneActivities", function(.object, ...) standardGeneric("getCiceroGeneActivities"), signature=c(".object"))
setGeneric("getRBFGeneActivities", function(.object, ...) standardGeneric("getRBFGeneActivities"), signature=c(".object"))

################################################################################
# DsATACsc (Single-Cell) Generics
################################################################################
setGeneric("simulateDoublets", function(.object, ...) standardGeneric("simulateDoublets"), signature=c(".object"))
setGeneric("getScQcStatsTab", function(.object, ...) standardGeneric("getScQcStatsTab"), signature=c(".object"))
setGeneric("filterCellsTssEnrichment", function(.object, ...) standardGeneric("filterCellsTssEnrichment"), signature=c(".object"))
setGeneric("unsupervisedAnalysisSc", function(.object, ...) standardGeneric("unsupervisedAnalysisSc"), signature=c(".object"))
setGeneric("dimRed_UMAP", function(.object, ...) standardGeneric("dimRed_UMAP"), signature=c(".object"))
setGeneric("iterativeLSI", function(.object, ...) standardGeneric("iterativeLSI"), signature=c(".object"))
setGeneric("mergePseudoBulk", function(.object, ...) standardGeneric("mergePseudoBulk"), signature=c(".object"))
setGeneric("samplePseudoBulk", function(.object, ...) standardGeneric("samplePseudoBulk"), signature=c(".object"))

################################################################################
# DsNOMe Generics
################################################################################
setGeneric("getMeth", function(.object, ...) standardGeneric("getMeth"), signature=c(".object"))
setGeneric("getCovg", function(.object, ...) standardGeneric("getCovg"), signature=c(".object"))
setGeneric("mergeStrands", function(.object, ...) standardGeneric("mergeStrands"), signature=c(".object"))
setGeneric("getRegionMapping", function(.object, ...) standardGeneric("getRegionMapping"), signature=c(".object"))
setGeneric("maskMethNA", function(.object, ...) standardGeneric("maskMethNA"), signature=c(".object"))
setGeneric("normalizeMeth", function(.object, ...) standardGeneric("normalizeMeth"), signature=c(".object"))

################################################################################
# Report Generics
################################################################################
setGeneric("createReport", function(.object, ...) standardGeneric("createReport"), signature=c(".object"))
setGeneric("createReport_summary", function(.object, ...) standardGeneric("createReport_summary"), signature=c(".object"))
setGeneric("createReport_exploratory", function(.object, ...) standardGeneric("createReport_exploratory"), signature=c(".object"))
setGeneric("createReport_differential", function(.object, ...) standardGeneric("createReport_differential"), signature=c(".object"))
setGeneric("createReport_filtering", function(.object, ...) standardGeneric("createReport_filtering"), signature=c(".object"))
setGeneric("createReport_normalization", function(.object, ...) standardGeneric("createReport_normalization"), signature=c(".object"))
