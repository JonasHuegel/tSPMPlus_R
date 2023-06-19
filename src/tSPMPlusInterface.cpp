#include <Rcpp.h>
#include "tspm_cpp_backend/utils/sequencing.cpp"
#include "tspm_cpp_backend/utils/utils.cpp"
#include "tspm_cpp_backend/utils/sorter.cpp"
#include "tspm_cpp_backend/utils/dbMartEntry.h"
#include "tspm_cpp_backend/utils/workflows.cpp"


using namespace Rcpp;


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
// [[Rcpp::export]]
size_t createTransitiveSequences(DataFrame &df_dbMart,
                                 std::string outputDir, 
                                 std::string outputFilePrefix,
                                 int numOfThreads = 1){
  
  bool storeSeqDuringCreation = false;
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
  

  
  std::set<unsigned int > endPhenxSet;
  endPhenxSet.insert(endPhenx.begin(), endPhenx.end());
  sequences = tspm::extractSequencesWithEnd(sequences, bitShift, lengthOfPhenx, endPhenxSet, numOfThreads);
  Rcout<< sequences.size() << std::endl;
  Rcout.flush();
  return transformToCandidateDataFrame(sequences, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);

}


//' Get Candidate Sequences For Phenx Of Interest (POI)
//' 
//' Function to extract all sequences that end with phenx j, that is a end phenx for each phenx in startPhenxOfInterest.
//' 
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
//' 
//' @returns The startPhenx of a sequence
//' @param  sequence  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
// [[Rcpp::export]]
unsigned int getStartPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength = 7){
  return tspm::getStartPhenx(sequence, phenxLength);
}

//' Get End Phenx From Sequence
//' 
//' Function to extract the end Phenx from a sequence.
//' 
//' @returns The endPhenx of a sequence.
//' @param  sequence  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
// [[Rcpp::export]]
unsigned int getEndPhenxFromSequence(std::uint64_t sequence, unsigned int phenxLength = 7){
  return tspm::getEndPhenx(sequence, phenxLength); 
}


//' Create a sequence from toi phenx
//' 
//' Function create a sequence from toi phenx.
//' 
//' @returns Integer.
//' @param  firstPhenx  Integer.
//' @param  secondPhenx  Integer.
//' @param phenxLength Integer, number of digits used to represent the second phenx in the sequence.
// [[Rcpp::export]]
std::uint64_t createSequence(unsigned int firstPhenx, unsigned int secondPhenx, unsigned int phenxLength = 7){
  return tspm::createSequence(firstPhenx, secondPhenx, phenxLength);
}
