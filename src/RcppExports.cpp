// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// createTransitiveSequences
size_t createTransitiveSequences(DataFrame& df_dbMart, std::string outputDir, std::string outputFilePrefix, int numOfThreads);
RcppExport SEXP _tSPMPlus_createTransitiveSequences(SEXP df_dbMartSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(createTransitiveSequences(df_dbMart, outputDir, outputFilePrefix, numOfThreads));
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
// getSequencesWithEndPhenx
DataFrame getSequencesWithEndPhenx(DataFrame& df_dbMart, unsigned int bitShift, unsigned int lengthOfPhenx, IntegerVector& lowerBucketThresholds, IntegerVector& endPhenx, bool includeCorBuckets, std::uint64_t minDuration, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence, bool returnSummary, bool summaryOnPatientLevel);
RcppExport SEXP _tSPMPlus_getSequencesWithEndPhenx(SEXP df_dbMartSEXP, SEXP bitShiftSEXP, SEXP lengthOfPhenxSEXP, SEXP lowerBucketThresholdsSEXP, SEXP endPhenxSEXP, SEXP includeCorBucketsSEXP, SEXP minDurationSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP, SEXP returnSummarySEXP, SEXP summaryOnPatientLevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bitShift(bitShiftSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lengthOfPhenx(lengthOfPhenxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThresholds(lowerBucketThresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type endPhenx(endPhenxSEXP);
    Rcpp::traits::input_parameter< bool >::type includeCorBuckets(includeCorBucketsSEXP);
    Rcpp::traits::input_parameter< std::uint64_t >::type minDuration(minDurationSEXP);
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
    Rcpp::traits::input_parameter< bool >::type returnSummary(returnSummarySEXP);
    Rcpp::traits::input_parameter< bool >::type summaryOnPatientLevel(summaryOnPatientLevelSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequencesWithEndPhenx(df_dbMart, bitShift, lengthOfPhenx, lowerBucketThresholds, endPhenx, includeCorBuckets, minDuration, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence, returnSummary, summaryOnPatientLevel));
    return rcpp_result_gen;
END_RCPP
}
// getSequencesWithStartPhenx
DataFrame getSequencesWithStartPhenx(DataFrame& df_dbMart, unsigned int bitShift, unsigned int lengthOfPhenx, IntegerVector& lowerBucketThresholds, IntegerVector& startPhenx, bool includeCorBuckets, std::uint64_t minDuration, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence, bool returnSummary, bool summaryOnPatientLevel);
RcppExport SEXP _tSPMPlus_getSequencesWithStartPhenx(SEXP df_dbMartSEXP, SEXP bitShiftSEXP, SEXP lengthOfPhenxSEXP, SEXP lowerBucketThresholdsSEXP, SEXP startPhenxSEXP, SEXP includeCorBucketsSEXP, SEXP minDurationSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP, SEXP returnSummarySEXP, SEXP summaryOnPatientLevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bitShift(bitShiftSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lengthOfPhenx(lengthOfPhenxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThresholds(lowerBucketThresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type startPhenx(startPhenxSEXP);
    Rcpp::traits::input_parameter< bool >::type includeCorBuckets(includeCorBucketsSEXP);
    Rcpp::traits::input_parameter< std::uint64_t >::type minDuration(minDurationSEXP);
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
    Rcpp::traits::input_parameter< bool >::type returnSummary(returnSummarySEXP);
    Rcpp::traits::input_parameter< bool >::type summaryOnPatientLevel(summaryOnPatientLevelSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequencesWithStartPhenx(df_dbMart, bitShift, lengthOfPhenx, lowerBucketThresholds, startPhenx, includeCorBuckets, minDuration, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence, returnSummary, summaryOnPatientLevel));
    return rcpp_result_gen;
END_RCPP
}
// getCandidateSequencesForPOI
DataFrame getCandidateSequencesForPOI(DataFrame& df_dbMart, std::uint64_t minDuration, unsigned int bitShift, unsigned int lengthOfPhenx, IntegerVector& lowerBucketThresholds, IntegerVector& startPhenxOfInterest, bool storeSeqDuringCreation, std::string outputDir, std::string outputFilePrefix, int numOfThreads, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, bool returnDuration, double durationPeriods, unsigned int daysForCoOoccurence);
RcppExport SEXP _tSPMPlus_getCandidateSequencesForPOI(SEXP df_dbMartSEXP, SEXP minDurationSEXP, SEXP bitShiftSEXP, SEXP lengthOfPhenxSEXP, SEXP lowerBucketThresholdsSEXP, SEXP startPhenxOfInterestSEXP, SEXP storeSeqDuringCreationSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP numOfThreadsSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP returnDurationSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< std::uint64_t >::type minDuration(minDurationSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bitShift(bitShiftSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lengthOfPhenx(lengthOfPhenxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThresholds(lowerBucketThresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type startPhenxOfInterest(startPhenxOfInterestSEXP);
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
    rcpp_result_gen = Rcpp::wrap(getCandidateSequencesForPOI(df_dbMart, minDuration, bitShift, lengthOfPhenx, lowerBucketThresholds, startPhenxOfInterest, storeSeqDuringCreation, outputDir, outputFilePrefix, numOfThreads, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, returnDuration, durationPeriods, daysForCoOoccurence));
    return rcpp_result_gen;
END_RCPP
}
// getStartPhenxFromSequence
unsigned int getStartPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength);
RcppExport SEXP _tSPMPlus_getStartPhenxFromSequence(SEXP sequenceSEXP, SEXP phenxLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::uint64_t >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type phenxLength(phenxLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(getStartPhenxFromSequence(sequence, phenxLength));
    return rcpp_result_gen;
END_RCPP
}
// getEndPhenxFromSequence
unsigned int getEndPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength);
RcppExport SEXP _tSPMPlus_getEndPhenxFromSequence(SEXP sequenceSEXP, SEXP phenxLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::uint64_t >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type phenxLength(phenxLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(getEndPhenxFromSequence(sequence, phenxLength));
    return rcpp_result_gen;
END_RCPP
}
// createSequence
std::uint64_t createSequence(unsigned int firstPhenx, unsigned int secondPhenx, unsigned int phenxLength);
RcppExport SEXP _tSPMPlus_createSequence(SEXP firstPhenxSEXP, SEXP secondPhenxSEXP, SEXP phenxLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type firstPhenx(firstPhenxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type secondPhenx(secondPhenxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type phenxLength(phenxLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(createSequence(firstPhenx, secondPhenx, phenxLength));
    return rcpp_result_gen;
END_RCPP
}
// sequenceAndSummarize
DataFrame sequenceAndSummarize(DataFrame df_dbMart, IntegerVector& lowerBucketThreshold, bool storeSeqDuringCreation, bool includeDurations, int numOfThreads, std::string outputDir, std::string outputFilePrefix, bool removeSparseSequences, double sparsityValue, bool createTemporalBuckets, bool durationSparsity, double durationSparsityValue, bool removeSparseTemporalBuckets, int patIdLength, double durationPeriods, unsigned int daysForCoOoccurence, bool summaryOnPatientLevel);
RcppExport SEXP _tSPMPlus_sequenceAndSummarize(SEXP df_dbMartSEXP, SEXP lowerBucketThresholdSEXP, SEXP storeSeqDuringCreationSEXP, SEXP includeDurationsSEXP, SEXP numOfThreadsSEXP, SEXP outputDirSEXP, SEXP outputFilePrefixSEXP, SEXP removeSparseSequencesSEXP, SEXP sparsityValueSEXP, SEXP createTemporalBucketsSEXP, SEXP durationSparsitySEXP, SEXP durationSparsityValueSEXP, SEXP removeSparseTemporalBucketsSEXP, SEXP patIdLengthSEXP, SEXP durationPeriodsSEXP, SEXP daysForCoOoccurenceSEXP, SEXP summaryOnPatientLevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df_dbMart(df_dbMartSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type lowerBucketThreshold(lowerBucketThresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type storeSeqDuringCreation(storeSeqDuringCreationSEXP);
    Rcpp::traits::input_parameter< bool >::type includeDurations(includeDurationsSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputDir(outputDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilePrefix(outputFilePrefixSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseSequences(removeSparseSequencesSEXP);
    Rcpp::traits::input_parameter< double >::type sparsityValue(sparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type createTemporalBuckets(createTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< bool >::type durationSparsity(durationSparsitySEXP);
    Rcpp::traits::input_parameter< double >::type durationSparsityValue(durationSparsityValueSEXP);
    Rcpp::traits::input_parameter< bool >::type removeSparseTemporalBuckets(removeSparseTemporalBucketsSEXP);
    Rcpp::traits::input_parameter< int >::type patIdLength(patIdLengthSEXP);
    Rcpp::traits::input_parameter< double >::type durationPeriods(durationPeriodsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type daysForCoOoccurence(daysForCoOoccurenceSEXP);
    Rcpp::traits::input_parameter< bool >::type summaryOnPatientLevel(summaryOnPatientLevelSEXP);
    rcpp_result_gen = Rcpp::wrap(sequenceAndSummarize(df_dbMart, lowerBucketThreshold, storeSeqDuringCreation, includeDurations, numOfThreads, outputDir, outputFilePrefix, removeSparseSequences, sparsityValue, createTemporalBuckets, durationSparsity, durationSparsityValue, removeSparseTemporalBuckets, patIdLength, durationPeriods, daysForCoOoccurence, summaryOnPatientLevel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tSPMPlus_createTransitiveSequences", (DL_FUNC) &_tSPMPlus_createTransitiveSequences, 4},
    {"_tSPMPlus_tSPMPlus", (DL_FUNC) &_tSPMPlus_tSPMPlus, 15},
    {"_tSPMPlus_extractNonSparseSequences", (DL_FUNC) &_tSPMPlus_extractNonSparseSequences, 9},
    {"_tSPMPlus_extractAllTransiviteSequences", (DL_FUNC) &_tSPMPlus_extractAllTransiviteSequences, 8},
    {"_tSPMPlus_getSequencesWithEndPhenx", (DL_FUNC) &_tSPMPlus_getSequencesWithEndPhenx, 23},
    {"_tSPMPlus_getSequencesWithStartPhenx", (DL_FUNC) &_tSPMPlus_getSequencesWithStartPhenx, 23},
    {"_tSPMPlus_getCandidateSequencesForPOI", (DL_FUNC) &_tSPMPlus_getCandidateSequencesForPOI, 20},
    {"_tSPMPlus_getStartPhenxFromSequence", (DL_FUNC) &_tSPMPlus_getStartPhenxFromSequence, 2},
    {"_tSPMPlus_getEndPhenxFromSequence", (DL_FUNC) &_tSPMPlus_getEndPhenxFromSequence, 2},
    {"_tSPMPlus_createSequence", (DL_FUNC) &_tSPMPlus_createSequence, 3},
    {"_tSPMPlus_sequenceAndSummarize", (DL_FUNC) &_tSPMPlus_sequenceAndSummarize, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_tSPMPlus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
