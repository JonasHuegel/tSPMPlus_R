// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// createTransitiveSequences
size_t createTransitiveSequences(DataFrame& df_dbMart, size_t numOfPatients, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads);
RcppExport SEXP _tSPMPlus_createTransitiveSequences(SEXP df_dbMartSEXP, SEXP numOfPatientsSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< size_t >::type numOfPatients(numOfPatientsSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(createTransitiveSequences(df_dbMart, numOfPatients, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads));
    return rcpp_result_gen;
END_RCPP
}
// tSPMPlus
DataFrame tSPMPlus(DataFrame& df_dbMart, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_tSPMPlus(SEXP df_dbMartSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseSequences(removeSparseSequencesSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type createTemporalBuckets(createTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< bool >::type durationSparsity(durationSparsitySEXP);
    Rcpp::traits::input_parameter< double >::type durationSparsityValue(durationSparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseTemporalBuckets(removeSparseTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< int >::type patIdLength(patIdLengthSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(tSPMPlus(df_dbMart, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// extractNonSparseSequences
DataFrame extractNonSparseSequences(DataFrame& df_dbMart, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, double sparsityValue, int numOfThreads, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_extractNonSparseSequences(SEXP df_dbMartSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP sparsityValueSEXP, SEXP numOfThreadsSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(extractNonSparseSequences(df_dbMart, storeSeqDuringCreation, outputDir, outputFilePrefix, sparsityValue, numOfThreads, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// extractAllTransiviteSequences
DataFrame extractAllTransiviteSequences(DataFrame& df_dbMart, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_extractAllTransiviteSequences(SEXP df_dbMartSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(extractAllTransiviteSequences(df_dbMart, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// getSequencesWithCandidateEnd
DataFrame getSequencesWithCandidateEnd(DataFrame& df_dbMart, std::uint64_t minDuration, unsigned int bitShift, unsigned int lengthOfPhenx, IntegerVector& lowerBucketThresholds, IntegerVector& startPhenxOfInterrest, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_getSequencesWithCandidateEnd(SEXP df_dbMartSEXP, SEXP minDurationSEXP, SEXP bitShiftSEXP, SEXP lengthOfPhenxSEXP, SEXP lowerBucketThresholdsSEXP, SEXP startPhenxOfInterrestSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::uint64_t >::type minDuration(minDurationSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bitShift(bitShiftSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lengthOfPhenx(lengthOfPhenxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThresholds(lowerBucketThresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type startPhenxOfInterrest(startPhenxOfInterrestSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseSequences(removeSparseSequencesSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type createTemporalBuckets(createTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< bool >::type durationSparsity(durationSparsitySEXP);
    Rcpp::traits::input_parameter< double >::type durationSparsityValue(durationSparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseTemporalBuckets(removeSparseTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< int >::type patIdLength(patIdLengthSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequencesWithCandidateEnd(df_dbMart, minDuration, bitShift, lengthOfPhenx, lowerBucketThresholds, startPhenxOfInterrest, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tSPMPlus_createTransitiveSequences", (DL_FUNC) &_tSPMPlus_createTransitiveSequences, 6},
    {"_tSPMPlus_tSPMPlus", (DL_FUNC) &_tSPMPlus_tSPMPlus, 15},
    {"_tSPMPlus_extractNonSparseSequences", (DL_FUNC) &_tSPMPlus_extractNonSparseSequences, 9},
    {"_tSPMPlus_extractAllTransiviteSequences", (DL_FUNC) &_tSPMPlus_extractAllTransiviteSequences, 8},
    {"_tSPMPlus_getSequencesWithCandidateEnd", (DL_FUNC) &_tSPMPlus_getSequencesWithCandidateEnd, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_tSPMPlus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
