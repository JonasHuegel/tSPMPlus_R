#include <Rcpp.h>
#include "tspm_cpp_backend/utils/sequencing.cpp"
#include "tspm_cpp_backend/utils/utils.cpp"
#include "tspm_cpp_backend/utils/sorter.cpp"
#include "tspm_cpp_backend/utils/dbMartEntry.h"
#include "tspm_cpp_backend/utils/workflows.cpp"


using namespace Rcpp;


//' transform the dbmart data frame to c++ structure
//' Function to create all transitive sequences
//' an store in binary format them to files.
//' @name transformDataFrameToStruct
//' @returns dbMart as c++ struct
//' @param  df_dbMart The dataframe that stores the numeric data mart.
std::vector<tspm::dbMartEntry> transformDataFrameToStruct(DataFrame &dfDbMart){
  IntegerVector patientIds = dfDbMart[0];
  IntegerVector phenxIds = dfDbMart[1];
  DateVector startDates = dfDbMart[2];
  size_t numOfEntries = dfDbMart.nrows();
  Rcout<<"numOfEntries: " << numOfEntries <<"\n";
  Rcout.flush();
  std::vector<tspm::dbMartEntry> dbMart;
  dbMart.reserve(numOfEntries);
  
  for(size_t i = 0; i < numOfEntries; ++i){
    tspm::dbMartEntry entry;
    entry.patID = patientIds[i];
    entry.phenID = phenxIds[i];
    entry.date = tspm::getTimeFromString(((Date)startDates[i]).format().c_str());
    dbMart.emplace_back(entry);
  }
  return dbMart;
}


//' Transforms the C++ representation used to store the sequences into R
//' 
//' Function to transform a vector of the C++ representation used to store the sequences into an R dataframe.
//' @name transformToDefaultDataFrame
//' @param sequences std::vector of the sequences that should be transformed into a R dataframe
//' @param returnDuration Boolean, if the duration should be included into the returned dataframe
DataFrame transformToDefaultDataFrame(std::vector<tspm::temporalSequence> sequences, bool returnDuration){
 Rcout << "transform sequences from c++ structure in R DataFrame!\n";                                                                      
 std::vector<std::uint64_t> seqIDs;
 seqIDs.reserve(sequences.size());
 std::vector<int> patIDs;
 patIDs.reserve(sequences.size());
 std::vector<unsigned int> durations;
 if(returnDuration){
   durations.reserve(sequences.size());
 }
 
 for(tspm::temporalSequence seq: sequences){
   int patId = seq.patientID;
   std::int64_t seqId = seq.seqID;
   unsigned int duration = seq.duration;
   patIDs.emplace_back(patId);
   seqIDs.emplace_back(seqId);
   if(returnDuration){
     durations.emplace_back(duration);
   }
 }
 DataFrame sequenceDataFrame;
 if(returnDuration){
   sequenceDataFrame= DataFrame::create(Named("patient_num") = patIDs, Named("sequence") = seqIDs, Named("duration") = durations);
 }else{
   sequenceDataFrame = DataFrame::create(Named("patient_num") = patIDs, Named("sequence") = seqIDs);
 }
 return sequenceDataFrame;
}

//' create all transitive Sequences
//' 
//' Function to create all transitive sequences
//' an store in binary format them to files.
//' 
//' @returns the overall number of sequences stored.
//' @param  df_dbMart The dataframe that stores the numeric data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @export
// [[Rcpp::export]]
size_t createTransitiveSequences(DataFrame &df_dbMart,
                                 std::string outputDir, 
                                 std::string outputFilePrefix,
                                 int numOfThreads = 1){
  
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  
  Rcout <<"determine start positions!\n";
  Rcout.flush();
  std::vector<size_t> startPositions = tspm::extractStartPositions(dbMart);
  Rcout<< "Creating sequences!\n";
  Rcout.flush();

  size_t numOfCreatedSequences = tspm::writeSequencesFromArrayToFile(dbMart,
                                                           startPositions,
                                                           outputDir,
                                                           outputFilePrefix,
                                                           7,
                                                           numOfThreads);
  Rcout<< "Created " << numOfCreatedSequences << " sequences.\n";
  return numOfCreatedSequences;
  
}




//' summarize
//'
//' Function to summarize a set of sequences and  returing a R data frame that contains the summary.
//' that calls this functions with predefined parameters.
//' @name summarize
//' @returns The summary as data frame.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param lowerBucketThresholds IntegerVector, lower duration Thresholds for the duration buckets of the candidate sequences
//' @param includeDurations include duration buckets in the summary
//' @param sequences the vector containing the sequences to summarize
DataFrame summarize(std::vector<tspm::temporalSequence> &sequences,
                    std::vector<unsigned int> lowerBucketThreshold,
                    bool includeDurations = false,
                    bool summaryOnPatientLevel = false,
                    unsigned int numOfThreads = 1){
  std::vector<std::pair<tspm::temporalSequence, size_t>> summary =
    tspm::summarizeSequencesAsVector(sequences,
                                     includeDurations,
                                     lowerBucketThreshold,
                                     summaryOnPatientLevel,
                                     numOfThreads);
  sequences.clear();
  sequences.shrink_to_fit();
  std::vector<std::uint64_t> seqIDs;
  seqIDs.reserve(summary.size());
  std::vector<std::uint64_t> counts;
  counts.reserve(summary.size());
  std::vector<std::uint32_t> durationBuckets;
  durationBuckets.reserve(summary.size());
  std::vector<std::uint32_t> patientIDs;
  patientIDs.reserve(summary.size());
  
  
  for(std::pair<tspm::temporalSequence, size_t> entry: summary){
    std::uint64_t seq = entry.first.seqID;
    seqIDs.emplace_back(seq);
    
    if(includeDurations){
      std::uint32_t dur = entry.first.duration;
      durationBuckets.emplace_back(dur);
    }
    if(summaryOnPatientLevel){
      std::uint32_t patID = entry.first.patientID;
      patientIDs.emplace_back(patID);
    }
    counts.emplace_back(entry.second);
  }
  DataFrame summaryDF = NULL; 
  if(includeDurations){
    if(summaryOnPatientLevel){
      summaryDF = DataFrame::create(Named("patientID") = patientIDs, Named("sequence") = seqIDs, Named("durBucket") = durationBuckets, Named("count") = counts);
    }else{
      summaryDF = DataFrame::create(Named("sequence") = seqIDs, Named("durBucket") = durationBuckets, Named("count") = counts);
    }
  }else{
    if(summaryOnPatientLevel){
      summaryDF = DataFrame::create(Named("patientID") = patientIDs, Named("sequence") = seqIDs,Named("count") = counts);
    }else{
      summaryDF = DataFrame::create(Named("sequence") = seqIDs,Named("count") = counts);
    }
  }

  return summaryDF;
}


//' tSPMPlus
//'
//' Function to call the tSPM+ workflow, most of the other function provide by this package are wrapper functions
//' that calls this functions with predefined parameters.
//' 
//' @returns The sequences as data frame.
//' @param  df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param removeSparseSequences  Boolean parameter to control if the sparsity should be applied
//' @param sparsityValue          The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets  Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param durationSparsity  Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param patIdLength describes the number of digits that are used for the patient number.
//' @param returnDuration Boolean, controls if the data frame that is returns contains. 
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @export
// [[Rcpp::export]]
DataFrame tSPMPlus(DataFrame &df_dbMart,
             bool storeSeqDuringCreation = false,
             std::string outputDir = "",
             std::string outputFilePrefix = "",
             int numOfThreads = 1,
             bool removeSparseSequences = true,
             double sparsityValue = 0.05,
             bool createTemporalBuckets = false,
             bool durationSparsity = false,
             double durationSparsityValue = 0,
             bool removeSparseTemporalBuckets = false,
             int patIdLength= 7,
             bool returnDuration = true,
             double durationPeriods = 30.437,
             unsigned int daysForCoOoccurence = 14 ){
  
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                              storeSeqDuringCreation,
                                                              outputDir,
                                                              outputFilePrefix,
                                                              removeSparseSequences,
                                                              sparsityValue,
                                                              createTemporalBuckets,
                                                              durationPeriods,
                                                              daysForCoOoccurence,
                                                              durationSparsity,
                                                              durationSparsityValue,
                                                              removeSparseTemporalBuckets,
                                                              patIdLength,
                                                              numOfThreads);
  
  Rcout << "created " << sequences.size() << " transitive sequences!\n";
  
  return transformToDefaultDataFrame(sequences, returnDuration);
} 



//' extract Non Sparse Sequences
//'
//' Extracts all non sparse transitive sequences from the data mart, after converting it to numeric.
//'
//' @returns The sequences as data frame.
//' @param df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param sparsityValue          The numeric value for the sparsity. DEFAULT = 0.05.
//' @param returnDuration Boolean, controls if the data frame that is returns contains.
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, eg. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @export
// [[Rcpp::export]]
DataFrame extractNonSparseSequences(DataFrame &df_dbMart,
                                      bool storeSeqDuringCreation = false,
                                      std::string outputDir = "",
                                      std::string outputFilePrefix = "",
                                      double sparsityValue = 0.05,
                                      int numOfThreads = 1,
                                      bool returnDuration = true,
                                      double durationPeriods = 30.437,
                                      unsigned int daysForCoOoccurence = 14 ){
    bool removeSparseSequences = true;
    bool createTemporalBuckets = false;
    bool removeSparseTemporalBuckets = false;
    int patIdLength= 7;
    bool durationSparsity = false;
    double durationSparsityValue = 0;
    return tSPMPlus(df_dbMart,
                    storeSeqDuringCreation,
                    outputDir,
                    outputFilePrefix,
                    numOfThreads,
                    removeSparseSequences,
                    sparsityValue,
                    createTemporalBuckets,
                    durationSparsity,
                    durationSparsityValue,
                    removeSparseTemporalBuckets,
                    patIdLength,
                    returnDuration,
                    durationPeriods,
                    daysForCoOoccurence);

}

//' Extract All Transivite Sequences
//' 
//' Function to extract all transitive sequences.
//' 
//' @returns The sequences as data frame.
//' @param  df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param returnDuration Boolean, controls if the data frame that is returns contains.
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @export
// [[Rcpp::export]]
DataFrame extractAllTransiviteSequences(DataFrame &df_dbMart,
                                    bool storeSeqDuringCreation = false,
                                    std::string outputDir = "",
                                    std::string outputFilePrefix = "",
                                    int numOfThreads = 1,
                                    bool returnDuration = true,
                                    double durationPeriods = 30.437,
                                    unsigned int daysForCoOoccurence = 14 ){
  double sparsityValue = 0;
  bool removeSparseSequences = false;
  bool createTemporalBuckets = false;
  bool removeSparseTemporalBuckets = false;
  int patIdLength= 7;
  bool durationSparsity = false;
  double durationSparsityValue = 0;
  return tSPMPlus(df_dbMart,
                  storeSeqDuringCreation,
                  outputDir,
                  outputFilePrefix,
                  numOfThreads,
                  removeSparseSequences,
                  sparsityValue,
                  createTemporalBuckets,
                  durationSparsity,
                  durationSparsityValue,
                  removeSparseTemporalBuckets,
                  patIdLength,
                  returnDuration,
                  durationPeriods,
                  daysForCoOoccurence);
  
}



std::vector<tspm::temporalSequence> extractCandidatesSequences(std::vector<tspm::temporalSequence> &originalSequences,
                                    std::uint64_t minDuration, unsigned int bitShift,
                                    unsigned int lengthOfPhenx, unsigned int numOfBuckets,
                                    std::vector<unsigned int> lowerBucketThresholds,
                                    std::vector<unsigned int> startPhenxOfInterrest,
                                    int &numOfThreads){
  
  std::set<unsigned int> candidateEndPhenx;
  candidateEndPhenx= extractEndPhenxWithGivenStartPhenx(originalSequences,
                                                        minDuration,
                                                        bitShift,
                                                        lengthOfPhenx,
                                                        startPhenxOfInterrest,
                                                        numOfThreads);
  
  Rcout<<"Extracted CandidatePhenx: " << candidateEndPhenx.size() << std::endl;
  Rcout.flush();
  
  std::vector<tspm::temporalSequence> sequencesOfInterest;
  sequencesOfInterest = extractSequencesWithEnd(originalSequences,
                                                bitShift,
                                                lengthOfPhenx,
                                                candidateEndPhenx,
                                                numOfThreads);
  
  Rcout<<"Extracted Candidate Sequences: " << sequencesOfInterest.size() << std::endl;
  Rcout.flush();
  return sequencesOfInterest;
}



DataFrame transformToCandidateDataFrame(std::vector<tspm::temporalSequence> &sequencesOfInterest,
                                        std::vector<unsigned int> lowerBucketThresholds,
                                        unsigned int lengthOfPhenx){
  
  std::vector<std::uint64_t> seqIDs;
  seqIDs.reserve(sequencesOfInterest.size());
  std::vector<int> patIDs;
  patIDs.reserve(sequencesOfInterest.size());
  std::vector<unsigned int> durations;
  durations.reserve(sequencesOfInterest.size());
  
  for(tspm::temporalSequence seq : sequencesOfInterest){
    int patId = seq.patientID;
    std::uint64_t seqId = seq.seqID;
    int duration = seq.duration;
    patIDs.emplace_back(patId);
    seqIDs.emplace_back(seqId);
    durations.emplace_back(duration);
  }
  
  sequencesOfInterest.clear();
  sequencesOfInterest.shrink_to_fit();
  
  std::vector<unsigned int> endPhenxVector;
  endPhenxVector.reserve(patIDs.size());
  std::vector<unsigned int> durationBuckets;
  durationBuckets.reserve(patIDs.size());
  for(size_t i = 0; i < patIDs.size();++i){
    tspm::temporalSequence seq;
    seq.seqID = seqIDs[i];

    unsigned int duration = durations[i];
    unsigned int endPhenx = getEndPhenx(seq, lengthOfPhenx);
    unsigned int durationBucket = tspm::getCandidateBucket(duration, lowerBucketThresholds);

    endPhenxVector.emplace_back(endPhenx);
    durationBuckets.emplace_back(durationBucket);
  }
  
  
  DataFrame  candidates = DataFrame::create(Named("patient_num") = patIDs, 
                                            Named("sequence") = seqIDs,
                                            Named("duration") = durations,
                                            Named("endPhenx") = endPhenxVector,
                                            Named("durationBucket") = durationBuckets);
  return candidates;
}



//' Get Sequences With Specific End Phenx
//' 
//' Function to extract all transitive sequences that end with given endPhenxs.
//' 
//' @returns The sequences as data frame.
//' @param df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param endPhenx IntegerVector, contains the phenx that sequences should end with.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param returnDuration Boolean, controls if the data frame that is returns contains. 
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @param minDuration the minimum duration a sequence must have, for j to be considered a candidate. Not Implemented at the moment!
//' @param bitShift  Integer, the number of bits used to shift the duration into sequnceID.
//' @param lengthOfPhenx describes the number of digits that represents a phenx in the sequence.
//' @param lowerBucketThresholds the lower thresholds for the temporal buckets, that are stored when the includeCorBuckets flag is set. 
//' @param includeCorBuckets Boolean, flag to control if the R data frame that is returned should contain columns for for the endPhenx and the buckets set in lowerBucketThresholds 
//' @param removeSparseSequences Boolean parameter to control if the sparsity should be applied.
//' @param sparsityValue The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param durationSparsity Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param patIdLength Integer, describes the number of digits that are used for the patient number.
//' @param returnSummary Boolean, if return a summary of the sequences instead of the sequences
//' @param summaryOnPatientLevel bool, that defines if the summary should be on the patient level (counting occurrences for each patient) or on the dbMart level
//' @param returnCandidateDataFrame Boolean to controll is a candidate dataframe should be returned, if returnSummary and returnCandidateDataFrame are both false the default sequence dataframe is returned
//' @export
// [[Rcpp::export]]
DataFrame getSequencesWithEndPhenx(DataFrame &df_dbMart,
                                   unsigned int bitShift,
                                   unsigned int lengthOfPhenx,
                                   IntegerVector &lowerBucketThresholds,
                                   IntegerVector &endPhenx,
                                   bool includeCorBuckets = false,
                                   std::uint64_t minDuration = 0,
                                   bool storeSeqDuringCreation = false,
                                   std::string outputDir = "",
                                   std::string outputFilePrefix = "",
                                   int numOfThreads = 1,
                                   bool removeSparseSequences = true,
                                   double sparsityValue = 0.05,
                                   bool createTemporalBuckets = false,
                                   bool durationSparsity = false,
                                   double durationSparsityValue = 0,
                                   bool removeSparseTemporalBuckets = false,
                                   int patIdLength= 7,
                                   bool returnDuration = true,
                                   double durationPeriods = 30.437,
                                   unsigned int daysForCoOoccurence = 14,
                                   bool returnSummary = false,
                                   bool summaryOnPatientLevel = false,
                                   bool returnCandidateDataFrame = true){
  
  
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                              storeSeqDuringCreation,
                                                              outputDir,
                                                              outputFilePrefix,
                                                              removeSparseSequences,
                                                              sparsityValue,
                                                              createTemporalBuckets,
                                                              durationPeriods,
                                                              daysForCoOoccurence,
                                                              durationSparsity,
                                                              durationSparsityValue,
                                                              removeSparseTemporalBuckets,
                                                              patIdLength,
                                                              numOfThreads);
  

  ips4o::parallel::sort(sequences.begin(), sequences.end(),tspm::timedSequencesSorter, numOfThreads);
  std::set<unsigned int > endPhenxSet;
  endPhenxSet.insert(endPhenx.begin(), endPhenx.end());
  sequences = tspm::extractSequencesWithEnd(sequences, bitShift, lengthOfPhenx, endPhenxSet, numOfThreads);
  Rcout<< sequences.size() << std::endl;
  Rcout.flush();
  if(returnSummary){
    return summarize(sequences,
                     as<std::vector<unsigned int>>(lowerBucketThresholds),
                     returnDuration,
                     summaryOnPatientLevel,
                     (unsigned int) numOfThreads);
  }else if(returnCandidateDataFrame){
   return transformToCandidateDataFrame(sequences, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);
  } else {
    return transformToDefaultDataFrame(sequences, returnDuration);
  }

}



//' Get Sequences With Specific Start Phenx
//' 
//' Function to extract all transitive sequences that end with given endPhenxs.
//' @name getSequencesWithStartPhenx
//' @returns The sequences as data frame.
//' @param df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param startPhenx IntegerVector, contains the phenx that sequences should end with.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param returnDuration Boolean, controls if the data frame that is returns contains. 
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @param minDuration the minimum duration a sequence must have, for j to be considered a candidate. Not Implemented at the moment!
//' @param bitShift  Integer, the number of bits used to shift the duration into sequnceID.
//' @param lengthOfPhenx describes the number of digits that represents a phenx in the sequence.
//' @param lowerBucketThresholds the lower thresholds for the temporal buckets, that are stored when the includeCorBuckets flag is set. 
//' @param includeCorBuckets Boolean, flag to control if the R data frame that is returned should contain columns for for the endPhenx and the buckets set in lowerBucketThresholds 
//' @param removeSparseSequences Boolean parameter to control if the sparsity should be applied.
//' @param sparsityValue The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param durationSparsity Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param patIdLength Integer, describes the number of digits that are used for the patient number.
//' @param summarize Boolean, if return a summary of the sequences instead of the sequences
//' @param returnCandidateDataFrame Boolean to controll is a candidate dataframe should be returned, if returnSummary and returnCandidateDataFrame are both false the default sequence dataframe is returned
//' @export
// [[Rcpp::export]]
DataFrame getSequencesWithStartPhenx(DataFrame &df_dbMart,
                                  unsigned int bitShift,
                                  unsigned int lengthOfPhenx,
                                  IntegerVector &lowerBucketThresholds,
                                  IntegerVector &startPhenx,
                                  bool includeCorBuckets = false,
                                  std::uint64_t minDuration = 0,
                                  bool storeSeqDuringCreation = false,
                                  std::string outputDir = "",
                                  std::string outputFilePrefix = "",
                                  int numOfThreads = 1,
                                  bool removeSparseSequences = true,
                                  double sparsityValue = 0.05,
                                  bool createTemporalBuckets = false,
                                  bool durationSparsity = false,
                                  double durationSparsityValue = 0,
                                  bool removeSparseTemporalBuckets = false,
                                  int patIdLength= 7,
                                  bool returnDuration = true,
                                  double durationPeriods = 30.437,
                                  unsigned int daysForCoOoccurence = 14,
                                  bool returnSummary = false,
                                  bool summaryOnPatientLevel = false,
                                  bool returnCandidateDataFrame = true){
 
 
 if(numOfThreads <= 0){
   numOfThreads = 1;
 }
 Rcout <<"Preparing data!\n";
 Rcout.flush();
 std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
 Rcout <<"Data prepared!\n";
 Rcout.flush();
 std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                                         storeSeqDuringCreation,
                                                                         outputDir,
                                                                         outputFilePrefix,
                                                                         removeSparseSequences,
                                                                         sparsityValue,
                                                                         createTemporalBuckets,
                                                                         durationPeriods,
                                                                         daysForCoOoccurence,
                                                                         durationSparsity,
                                                                         durationSparsityValue,
                                                                         removeSparseTemporalBuckets,
                                                                         patIdLength,
                                                                         numOfThreads);
 
 
 ips4o::parallel::sort(sequences.begin(), sequences.end(),tspm::timedSequencesSorter, numOfThreads);
 std::set<unsigned int > startPhenxSet;
 startPhenxSet.insert(startPhenx.begin(), startPhenx.end());
 sequences = tspm::extractSequencesWithSpecificStart(sequences, minDuration, bitShift, lengthOfPhenx, startPhenxSet, numOfThreads);
 Rcout<< sequences.size() << std::endl;
 Rcout.flush();
 if(returnSummary){
   return summarize(sequences,
                    as<std::vector<unsigned int>>(lowerBucketThresholds),
                    returnDuration,
                    summaryOnPatientLevel,
                    (unsigned int) numOfThreads);
 }else if(returnCandidateDataFrame){
   return transformToCandidateDataFrame(sequences, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);
 }else{
   return transformToDefaultDataFrame(sequences, returnDuration);
 }
 
}

//' Get Sequences that start or end with a vector of given phenx
//' 
//' Function to extract all transitive sequences that end with given endPhenxs.
//' @name getSequencesContainingPhenx
//' @returns The sequences as data frame.
//' @param df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param phenxOfInterest IntegerVector, contains the phenx that sequences should end with.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param returnDuration Boolean, controls if the data frame that is returns contains. 
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @param minDuration the minimum duration a sequence must have, for j to be considered a candidate. Not Implemented at the moment!
//' @param bitShift  Integer, the number of bits used to shift the duration into sequnceID.
//' @param lengthOfPhenx describes the number of digits that represents a phenx in the sequence.
//' @param lowerBucketThresholds the lower thresholds for the temporal buckets, that are stored when the includeCorBuckets flag is set. 
//' @param includeCorBuckets Boolean, flag to control if the R data frame that is returned should contain columns for for the endPhenx and the buckets set in lowerBucketThresholds 
//' @param removeSparseSequences Boolean parameter to control if the sparsity should be applied.
//' @param sparsityValue The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param durationSparsity Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param patIdLength Integer, describes the number of digits that are used for the patient number.
//' @param returnSummary Boolean, if return a summary of the sequences instead of the sequences
//' @param summaryOnPatientLevel bool, that defines if the summary should be on the patient level (counting occurrences for each patient) or on the dbMart level
//' @param returnCandidateDataFrame Boolean to controll is a candidate dataframe should be returned, if returnSummary and returnCandidateDataFrame are both false the default sequence dataframe is returned
//' @export
// [[Rcpp::export]]
DataFrame getSequencesContainingPhenx(DataFrame &df_dbMart,
                                    unsigned int bitShift,
                                    unsigned int lengthOfPhenx,
                                    IntegerVector &lowerBucketThresholds,
                                    IntegerVector &phenxOfInterest,
                                    bool includeCorBuckets = false,
                                    std::uint64_t minDuration = 0,
                                    bool storeSeqDuringCreation = false,
                                    std::string outputDir = "",
                                    std::string outputFilePrefix = "",
                                    int numOfThreads = 1,
                                    bool removeSparseSequences = true,
                                    double sparsityValue = 0.05,
                                    bool createTemporalBuckets = false,
                                    bool durationSparsity = false,
                                    double durationSparsityValue = 0,
                                    bool removeSparseTemporalBuckets = false,
                                    int patIdLength= 7,
                                    bool returnDuration = true,
                                    double durationPeriods = 30.437,
                                    unsigned int daysForCoOoccurence = 14,
                                    bool returnSummary = false,
                                    bool summaryOnPatientLevel = false,
                                    bool returnCandidateDataFrame = false){
 
 
 if(numOfThreads <= 0){
   numOfThreads = 1;
 }
 Rcout <<"Preparing data!\n";
 Rcout.flush();
 std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
 Rcout <<"Data prepared!\n";
 Rcout.flush();
 std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                                         storeSeqDuringCreation,
                                                                         outputDir,
                                                                         outputFilePrefix,
                                                                         removeSparseSequences,
                                                                         sparsityValue,
                                                                         createTemporalBuckets,
                                                                         durationPeriods,
                                                                         daysForCoOoccurence,
                                                                         durationSparsity,
                                                                         durationSparsityValue,
                                                                         removeSparseTemporalBuckets,
                                                                         patIdLength,
                                                                         numOfThreads);
 
 
 ips4o::parallel::sort(sequences.begin(), sequences.end(),tspm::timedSequencesSorter, numOfThreads);
 std::set<unsigned int > phenxSet;
 phenxSet.insert(phenxOfInterest.begin(), phenxOfInterest.end());
 sequences = tspm::extractSequencesWithPhenx(sequences, minDuration, bitShift, lengthOfPhenx, phenxSet, numOfThreads);
 Rcout<< "Reduced to" << sequences.size() << "containing a phenx of interest" << std::endl;
 Rcout.flush();
 if(returnSummary){
   return summarize(sequences,
                    as<std::vector<unsigned int>>(lowerBucketThresholds),
                    returnDuration,
                    summaryOnPatientLevel,
                    (unsigned int) numOfThreads);
 }else if(returnCandidateDataFrame){
   return transformToCandidateDataFrame(sequences, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);
 }else{
   return transformToDefaultDataFrame(sequences, returnDuration);
 }
 
}



//' Get Candidate Sequences For Phenx Of Interest (POI)
//' 
//' Function to extract all sequences that end with phenx j, that is a end phenx for each phenx in startPhenxOfInterest.
//' @name getCandidateSequencesForPOI
//' @returns The sequences as data frame.
//' @param  df_dbMart The data frame that stores the data mart.
//' @param minDuration the minimum duration at least one sequence from type poi->j must have, for j to be considered a candidate.
//' @param bitShift  Integer, the number of bits used to shift the duration into sequnceID.
//' @param lengthOfPhenx Integer, number of digits used to represent the second phenx in the sequence.
//' @param lowerBucketThresholds IntegerVector, lower duration Thresholds for the duration buckets of the candidate sequences
//' @param startPhenxOfInterest IntegerVector, numeric representation of the POI.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param removeSparseSequences  Boolean parameter to control if the sparsity should be applied.
//' @param sparsityValue          The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets  Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param durationSparsity  Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param patIdLength Integer, describes the number of digits that are used for the patient number.
//' @param returnDuration Boolean, controls if the data frame that is returns contains.
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @export
// [[Rcpp::export]]
DataFrame getCandidateSequencesForPOI(DataFrame &df_dbMart,
                                       std::uint64_t minDuration,
                                       unsigned int bitShift,
                                       unsigned int lengthOfPhenx,
                                       IntegerVector &lowerBucketThresholds,
                                       IntegerVector &startPhenxOfInterest,
                                       bool storeSeqDuringCreation = false,
                                       std::string outputDir = "",
                                       std::string outputFilePrefix = "",
                                       int numOfThreads = 1,
                                       bool removeSparseSequences = true,
                                       double sparsityValue = 0.05,
                                       bool createTemporalBuckets = false,
                                       bool durationSparsity = false,
                                       double durationSparsityValue = 0,
                                       bool removeSparseTemporalBuckets = false,
                                       int patIdLength= 7,
                                       bool returnDuration = true,
                                       double durationPeriods = 30.437,
                                       unsigned int daysForCoOoccurence = 14 ){
  
  
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                              storeSeqDuringCreation,
                                                              outputDir,
                                                              outputFilePrefix,
                                                              removeSparseSequences,
                                                              sparsityValue,
                                                              createTemporalBuckets,
                                                              durationPeriods,
                                                              daysForCoOoccurence,
                                                              durationSparsity,
                                                              durationSparsityValue,
                                                              removeSparseTemporalBuckets,
                                                              patIdLength,
                                                              numOfThreads);
  
  
  Rcout << "Reached sequenceofInterest" << std::endl;
  Rcout.flush();
  std::vector<tspm::temporalSequence> sequencesOfInterest = extractCandidatesSequences(sequences,
                                                                                 minDuration,
                                                                                 bitShift,
                                                                                 lengthOfPhenx,
                                                                                 as<std::vector<unsigned int>>(lowerBucketThresholds).size(),
                                                                                 as<std::vector<unsigned int>>(lowerBucketThresholds),
                                                                                 as<std::vector<unsigned int>>(startPhenxOfInterest),
                                                                                 numOfThreads);
  Rcout<< sequencesOfInterest.size() << std::endl;
  Rcout.flush();
  sequences.clear();
  sequences.shrink_to_fit();
  return transformToCandidateDataFrame(sequencesOfInterest, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);
}



//' Get Start Phenx From Sequence
//' 
//' Function to extract the start Phenx from a sequence.
//' @name getStartPhenxFromSequence
//' @returns The startPhenx of a sequence
//' @param  sequence  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
//' @export
// [[Rcpp::export]]
unsigned int getStartPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength = 7){
  return tspm::getStartPhenx(sequence, phenxLength);
}

//' Get End Phenx From Sequence
//' 
//' Function to extract the end Phenx from a sequence.
//' @name getEndPhenxFromSequence
//' @returns The endPhenx of a sequence.
//' @param  sequence  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
//' @export
// [[Rcpp::export]]
unsigned int getEndPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength = 7){
  return tspm::getEndPhenx(sequence, phenxLength); 
}


//' Create a sequence from two phenx
//' 
//' Function create a sequence from two phenx.
//' @name createSequence
//' @returns Integer.
//' @param  firstPhenx  Integer.
//' @param  secondPhenx  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
//' @export
// [[Rcpp::export]]
std::uint64_t createSequence(unsigned int firstPhenx, unsigned int secondPhenx, unsigned int phenxLength = 7){
  return tspm::createSequence(firstPhenx, secondPhenx, phenxLength);
}


//' Sequences and summarizes all sequences in a dbmart
//' 
//' Summarizes the sequences in a dbmart
//' @name sequenceAndSummarize
//' @returns The summary as data frame.
//' @param  df_dbMart The data frame that stores the data mart.
//' @param outputDir The path as string to the directory where the sequences should be stored.
//' @param outputFilePrefix The string file prefix for the patient files storing the sequences.
//' @param numOfThreads The number of threads that should be used during sequencing.
//' @param storeSeqDuringCreation  Boolean parameter to control if the duration should be included in the sequence ID during creation, DEFAULT = FALSE.
//' @param removeSparseSequences  Boolean parameter to control if the sparsity should be applied.
//' @param sparsityValue          The numeric value for the sparsity. DEFAULT = 0.05.
//' @param createTemporalBuckets  Boolean flag if the the the sequences should be split up in dynamic buckets. Number of buckets min(4, max_duration(sequence)).
//' @param removeSparseTemporalBuckets Boolean, to control if the sparsity should be applied on the dynamic temporal buckets.
//' @param durationSparsity  Boolean flag to control if sparse sequences should be removed considering the duration periods of a sequence.
//' @param durationSparsityValue Numeric value.
//' @param patIdLength Integer, describes the number of digits that are used for the patient number.
//' @param durationPeriods Numeric, Upper threshold, stores the number of day in the time period, e.g. 30.471 for months, 364.25 for years. 
//' @param daysForCoOoccurence Integer, sets the upper threshold for the sequence duration so that they are counted as co-occurrence (meaning a duration of 0).
//' @param lowerBucketThresholds IntegerVector, lower duration Thresholds for the duration buckets of the candidate sequences
//' @param includeDurations include duration buckets in the summary
//' @param summaryOnPatientLevel bool, that defines if the summary should be on the patient level (counting occurrences for each patient) or on the dbMart level
//' @export
// [[Rcpp::export]]
DataFrame sequenceAndSummarize(DataFrame df_dbMart,
                               IntegerVector &lowerBucketThresholds,
                               bool storeSeqDuringCreation = false,
                               bool includeDurations = false,
                               int numOfThreads = 1,
                               std::string outputDir = "",
                               std::string outputFilePrefix = "",
                               bool removeSparseSequences = true,
                               double sparsityValue = 0.05,
                               bool createTemporalBuckets = false,
                               bool durationSparsity = false,
                               double durationSparsityValue = 0,
                               bool removeSparseTemporalBuckets = false,
                               int patIdLength= 7,
                               double durationPeriods = 30.437,
                               unsigned int daysForCoOoccurence = 14,
                               bool summaryOnPatientLevel = false){
  
  std::vector<tspm::dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<tspm::temporalSequence> sequences =  tspm::sequenceWorkflow(dbMart,
                                                                          storeSeqDuringCreation,
                                                                          outputDir,
                                                                          outputFilePrefix,
                                                                          removeSparseSequences,
                                                                          sparsityValue,
                                                                          createTemporalBuckets,
                                                                          durationPeriods,
                                                                          daysForCoOoccurence,
                                                                          durationSparsity,
                                                                          durationSparsityValue,
                                                                          removeSparseTemporalBuckets,
                                                                          patIdLength,
                                                                          numOfThreads);
  
  return summarize(sequences,
                   as<std::vector<unsigned int>>(lowerBucketThresholds),
                   includeDurations,
                   summaryOnPatientLevel,
                   (unsigned int) numOfThreads);
}

  


