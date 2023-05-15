// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// createTransitiveSequences
size_t createTransitiveSequences(DataFrame& df_dbMart, size_t numOfPatients, std::string& outputDir, std::string& outputFilePrefix, int numOfThreads);
RcppExport SEXP _tSPMPlus_createTransitiveSequences(SEXP df_dbMartSEXP, SEXP numOfPatientsSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< size_t >::type numOfPatients(numOfPatientsSEXP);
    Rcpp::traits::input_parameter< std::string& >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string& >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(createTransitiveSequences(df_dbMart, numOfPatients, outputDir, outputFilePrefix, numOfThreads));
    return rcpp_result_gen;
END_RCPP
}
// tSPMPlus
DataFrame tSPMPlus(DataFrame& df_dbMart, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_tSPMPlus(SEXP df_dbMartSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
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
    rcpp_result_gen = Rcpp::wrap(tSPMPlus(df_dbMart, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// extractNonSparseSequences
DataFrame extractNonSparseSequences(DataFrame& df_dbMart, std::string outputDir, std::string outputFilePrefix, double sparsityValue, int numOfThreads, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_extractNonSparseSequences(SEXP df_dbMartSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP sparsityValueSEXP, SEXP numOfThreadsSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(extractNonSparseSequences(df_dbMart, outputDir, outputFilePrefix, sparsityValue, numOfThreads, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// extractAllTransiviteSequences
DataFrame extractAllTransiviteSequences(DataFrame& df_dbMart, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_extractAllTransiviteSequences(SEXP df_dbMartSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    rcpp_result_gen = Rcpp::wrap(extractAllTransiviteSequences(df_dbMart, outputDir, outputFilePrefix, numOfThreads, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// getSequencesWithCandidateEnd
DataFrame getSequencesWithCandidateEnd(DataFrame& df_dbMart, std::string outputDir, std::string outputFilePrefix, unsigned long minDuration, unsigned int bitShift, unsigned int lengthOfPhenx, IntegerVector& lowerBucketThresholds, IntegerVector& startPhenxOfInterrest, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, bool durationInWeeks, bool durationInMonths);
RcppExport SEXP _tSPMPlus_getSequencesWithCandidateEnd(SEXP df_dbMartSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP minDurationSEXP, SEXP bitShiftSEXP, SEXP lengthOfPhenxSEXP, SEXP lowerBucketThresholdsSEXP, SEXP startPhenxOfInterrestSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationInWeeksSEXP, SEXP durationInMonthsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type minDuration(minDurationSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bitShift(bitShiftSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lengthOfPhenx(lengthOfPhenxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThresholds(lowerBucketThresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type startPhenxOfInterrest(startPhenxOfInterrestSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseSequences(removeSparseSequencesSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type createTemporalBuckets(createTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< bool >::type durationSparsity(durationSparsitySEXP);
    Rcpp::traits::input_parameter< double >::type durationSparsityValue(durationSparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseTemporalBuckets(removeSparseTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< int >::type patIdLength(patIdLengthSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDuration(returnDurationSEXP);
    Rcpp::traits::input_parameter< bool >::type durationInWeeks(durationInWeeksSEXP);
    Rcpp::traits::input_parameter< bool >::type durationInMonths(durationInMonthsSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequencesWithCandidateEnd(df_dbMart, outputDir, outputFilePrefix, minDuration, bitShift, lengthOfPhenx, lowerBucketThresholds, startPhenxOfInterrest, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationInWeeks, durationInMonths));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tSPMPlus_createTransitiveSequences", (DL_FUNC) &_tSPMPlus_createTransitiveSequences, 5},
    {"_tSPMPlus_tSPMPlus", (DL_FUNC) &_tSPMPlus_tSPMPlus, 14},
    {"_tSPMPlus_extractNonSparseSequences", (DL_FUNC) &_tSPMPlus_extractNonSparseSequences, 8},
    {"_tSPMPlus_extractAllTransiviteSequences", (DL_FUNC) &_tSPMPlus_extractAllTransiviteSequences, 7},
    {"_tSPMPlus_getSequencesWithCandidateEnd", (DL_FUNC) &_tSPMPlus_getSequencesWithCandidateEnd, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_tSPMPlus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
