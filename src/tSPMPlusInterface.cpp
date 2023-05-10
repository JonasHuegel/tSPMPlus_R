
#include <Rcpp.h>
#include "tspm_cpp_backend/utils/sequencing.cpp"
#include "tspm_cpp_backend/utils/utils.cpp"
#include "tspm_cpp_backend/utils/sorter.cpp"
#include "tspm_cpp_backend/utils/dbMartEntry.h"
#include "tspm_cpp_backend/utils/workflows.cpp"


using namespace Rcpp;


std::vector<dbMartEntry> transformDataFrameToStruct(DataFrame &dfDbMart){
  IntegerVector patientIds = dfDbMart[0];
  IntegerVector phenxIds = dfDbMart[1];
  DateVector startDates = dfDbMart[2];
  size_t numOfEntries = dfDbMart.nrows();
  Rcout<<"numOfEntries: " << numOfEntries <<"\n";
  Rcout.flush();
  std::vector<dbMartEntry> dbMart;
  dbMart.reserve(numOfEntries);
  
  for(size_t i = 0; i < numOfEntries; ++i){
    dbMartEntry entry;
    entry.patID = patientIds[i];
    entry.phenID = phenxIds[i];
    entry.date = getTimeFromString(((Date)startDates[i]).format().c_str());
    dbMart.emplace_back(entry);
  }
  return dbMart;
}




// [[Rcpp::export]]
size_t createTransitiveSequences(DataFrame df_dbMart,size_t numOfPatients, std::string &outputDir, 
                                 std::string &outputFilePrefix, int numOfThreads){
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  
  Rcout <<"determine start positions!\n";
  Rcout.flush();
  std::vector<size_t> startPositions = extractStartPositions(dbMart);
  Rcout<< "Creating sequences!\n";
  Rcout.flush();

  size_t numOfCreatedSequences = extractSequencesFromArray(dbMart,
                                                           numOfPatients,
                                                           startPositions.data(),
                                                           outputDir,
                                                           outputFilePrefix,
                                                           7,
                                                           numOfThreads);
  Rcout<< "Created " << numOfCreatedSequences << " sequences.\n";
  return numOfCreatedSequences;
  
}

// [[Rcpp::export]]
DataFrame tSPMPlus(DataFrame &df_dbMart,
             std::string outputDir,
             std::string outputFilePrefix,
             int numOfThreads,
             bool removeSparseSequences = true,
             double sparsityValue = 0.05,
             bool createTemporalBuckets = false,
             bool durationSparsity = false,
             double durationSparsityValue = 0,
             bool removeSparseTemporalBuckets = false,
             int patIdLength= 7,
             bool addDuration = true,
             bool durationInWeeks = false,
             bool durationInMonths = false ){
  
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<temporalSequence> sequences =  sequenceWorkflow(dbMart,
                                                              outputDir,
                                                              outputFilePrefix,
                                                              removeSparseSequences,
                                                              sparsityValue,
                                                              createTemporalBuckets,
                                                              durationInWeeks,
                                                              durationInMonths,
                                                              durationSparsity,
                                                              durationSparsityValue,
                                                              removeSparseTemporalBuckets,
                                                              patIdLength,
                                                              numOfThreads);
  
  Rcout << "created " << sequences.size() << " transitive sequences!\n";
  
  Rcout << "transform sequences from c++ structure in R DataFrame!\n";                                                                      
  std::vector<unsigned long> seqIDs;
  seqIDs.reserve(sequences.size());
  std::vector<int> patIDs;
  patIDs.reserve(sequences.size());
  std::vector<unsigned int> durations;
  if(addDuration){
    durations.reserve(sequences.size());
  }
  
  for(temporalSequence seq: sequences){
    int patId = seq.patientID;
    long seqId = seq.seqID;
    long duration = seq.duration;
    patIDs.emplace_back(patId);
    seqIDs.emplace_back(seqId);
    if(addDuration){
      durations.emplace_back(duration);
    }
  }
  DataFrame sequenceDataFrame;
  if(addDuration){
    sequenceDataFrame= DataFrame::create(Named("patient_num") = patIDs, Named("sequence") = seqIDs, Named("duration") = durations);
  }else{
    sequenceDataFrame = DataFrame::create(Named("patient_num") = patIDs, Named("sequence") = seqIDs);
  }
  return sequenceDataFrame;
} 
  
// [[Rcpp::export]]
DataFrame extractNonSparseSequences(DataFrame df_dbMart,
                                      std::string outputDir,
                                      std::string outputFilePrefix,
                                      double sparsityValue,
                                      int numOfThreads,
                                      bool addDuration = true,
                                      bool durationInWeeks = false,
                                      bool durationInMonths = false ){
    bool removeSparseSequences = true;
    bool createTemporalBuckets = false;
    bool removeSparseTemporalBuckets = false;
    int patIdLength= 7;
    bool durationSparsity = false;
    double durationSparsityValue = 0;
    return tSPMPlus(df_dbMart,
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
                    addDuration,
                    durationInWeeks,
                    durationInMonths);

}

// [[Rcpp::export]]
DataFrame extractAllTransiviteSequences(DataFrame df_dbMart,
                                    std::string outputDir,
                                    std::string outputFilePrefix,
                                    int numOfThreads,
                                    bool addDuration = true,
                                    bool durationInWeeks = false,
                                    bool durationInMonths = false ){
  double sparsityValue = 0;
  bool removeSparseSequences = false;
  bool createTemporalBuckets = false;
  bool removeSparseTemporalBuckets = false;
  int patIdLength= 7;
  bool durationSparsity = false;
  double durationSparsityValue = 0;
  return tSPMPlus(df_dbMart,
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
                  addDuration,
                  durationInWeeks,
                  durationInMonths);
  
}

std::vector<temporalSequence> extractCandidatesSequences(std::vector<temporalSequence> &originalSequences,
                                    unsigned long minDuration, unsigned int bitShift,
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
  
  Rcout<<"ExtractedCandidatePhenx: " << candidateEndPhenx.size() << std::endl;
  
  std::vector<temporalSequence> sequencesOfInterest;
  sequencesOfInterest = extractSequencesWithEnd(originalSequences,
                                                bitShift,
                                                lengthOfPhenx,
                                                candidateEndPhenx,
                                                numOfThreads);

  return sequencesOfInterest;
}



DataFrame transformToCandidateDataFrame(std::vector<temporalSequence> sequencesOfInterest,
                                        std::vector<unsigned int> lowerBucketThresholds,
                                        unsigned int lengthOfPhenx){
  
  std::vector<unsigned int> endPhenxVector;
  endPhenxVector.reserve(sequencesOfInterest.size());
  std::vector<unsigned int> durationBuckets;
  durationBuckets.reserve(sequencesOfInterest.size());
  std::vector<unsigned long> seqIDs;
  seqIDs.reserve(sequencesOfInterest.size());
  std::vector<int> patIDs;
  patIDs.reserve(sequencesOfInterest.size());
  std::vector<unsigned int> durations;
  durations.reserve(sequencesOfInterest.size());
  
  for(temporalSequence seq : sequencesOfInterest){
    int patId = seq.patientID;
    long seqId = seq.seqID;
    long duration = seq.duration;
    unsigned int endPhenx = getEndPhenx(seq, lengthOfPhenx);
    unsigned int durationBucket = getCandidateBucket(duration, lowerBucketThresholds);
    patIDs.emplace_back(patId);
    seqIDs.emplace_back(seqId);
    durations.emplace_back(duration);
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



// [[Rcpp::export]]
DataFrame getSequencesWithCandidateEnd(DataFrame &df_dbMart,
                                       std::string outputDir,
                                       std::string outputFilePrefix,
                                       unsigned long minDuration,
                                       unsigned int bitShift,
                                       unsigned int lengthOfPhenx,
                                       unsigned int numOfBuckets,
                                       IntegerVector &lowerBucketThresholds,
                                       IntegerVector &startPhenxOfInterrest,
                                       int numOfThreads = 1,
                                       bool removeSparseSequences = true,
                                       double sparsityValue = 0.05,
                                       bool createTemporalBuckets = false,
                                       bool durationSparsity = true,
                                       double durationSparsityValue = 0,
                                       bool removeSparseTemporalBuckets = false,
                                       int patIdLength= 7,
                                       bool addDuration = true,
                                       bool durationInWeeks = false,
                                       bool durationInMonths = false ){
  
  
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  Rcout <<"Preparing data!\n";
  Rcout.flush();
  std::vector<dbMartEntry> dbMart = transformDataFrameToStruct(df_dbMart);
  Rcout <<"Data prepared!\n";
  Rcout.flush();
  std::vector<temporalSequence> sequences =  sequenceWorkflow(dbMart,
                                                              outputDir,
                                                              outputFilePrefix,
                                                              removeSparseSequences,
                                                              sparsityValue,
                                                              createTemporalBuckets,
                                                              durationInWeeks,
                                                              durationInMonths,
                                                              durationSparsity,
                                                              durationSparsityValue,
                                                              removeSparseTemporalBuckets,
                                                              patIdLength,
                                                              numOfThreads);
  
  
  std::vector<temporalSequence> sequencesOfInterest = extractCandidatesSequences(sequences,
                                                                                 minDuration,
                                                                                 bitShift,
                                                                                 lengthOfPhenx,
                                                                                 numOfBuckets,
                                                                                 as<std::vector<unsigned int>>(lowerBucketThresholds),
                                                                                 as<std::vector<unsigned int>>(startPhenxOfInterrest),
                                                                                 numOfThreads);
  Rcout<< sequencesOfInterest.size() << std::endl;
    
  return transformToCandidateDataFrame(sequencesOfInterest, as< std::vector<unsigned int> >(lowerBucketThresholds), lengthOfPhenx);
}
