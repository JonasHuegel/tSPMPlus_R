#include <Rcpp.h>
#include "tspm_cpp_backend/utils/sequencing.cpp"
#include "tspm_cpp_backend/utils/utils.cpp"
#include "tspm_cpp_backend/utils/sorter.cpp"



using namespace Rcpp;
// [[Rcpp::export]]
size_t createTransitiveSequences(DataFrame df_dbMart,size_t numOfPatients, std::string outputDir, 
                              std::string outputFilesDescription, int numOfThreads){
  Rcout <<"Preparing data!\n";
  size_t numOfEntries = df_dbMart.nrows();
  Rcout<<"numOfEntries: " << numOfEntries <<"\n";
  Rcout<<"numOfPatients: " << numOfPatients <<"\n";
  size_t startPositions[numOfPatients];
  IntegerVector patientIds = df_dbMart[0];
  IntegerVector phenxIds = df_dbMart[1];
  DateVector startDates = df_dbMart[2];
  std::cout << "try to alloc " <<(sizeof(dbMartEntry) * numOfEntries)/ (1024 *1024) << " MB\n";
  dbMartEntry* dbMart = (dbMartEntry *) malloc(sizeof(dbMartEntry) * numOfEntries);
  if(dbMart == nullptr){
    Rcout << "Error! could not allocate memory for dbmart \n";
    Rcout.flush();
    return 21;
  }
  Rcout<< "Allocated memory for DBMart. Filling it from DF and determine start position for each Patient!\n";
  Rcout.flush();
  startPositions[0] = 0;
  for(size_t i = 0; i < numOfEntries; ++i){
    dbMart[i].patID = patientIds[i];
    dbMart[i].phenID = phenxIds[i];
    dbMart[i].date = getTimeFromString(((Date)startDates[i]).format().c_str());;

    if(i>0 && dbMart[i].patID != dbMart[i-1].patID){
      startPositions[dbMart[i].patID] = i; 
    }
  }

  Rcout<< "Creating sequences.\n";
  size_t numOfCreatedSequences = extractSequencesFromArray(dbMart,
                                                           numOfPatients,
                                                           startPositions,
                                                           numOfEntries,
                                                           outputDir,
                                                           outputFilesDescription,
                                                           7,
                                                           numOfThreads);
  free(dbMart);
  Rcout<< "Created " << numOfCreatedSequences << " sequences.\n";
  return numOfCreatedSequences;
  //return 0;
}



// [[Rcpp::export]]
int tSPMPlus(int numOfThreads, std::string dbMartCsv,char inputFileDelimiter, std::string outputDir, 
             std::string outputFilesDescription, int patIDColumn, int phenxColumn,
             int dateColumn, int patientCount, bool createDuration,
             bool removeSparseBuckets, double sparsity ){
  if(numOfThreads <= 0){
    numOfThreads = 1;
  }
  omp_set_num_threads(numOfThreads);
  std::vector<std::string> inputFilePaths;
  inputFilePaths.emplace_back(dbMartCsv);
  
  int patIdColumns[1] = {patIDColumn};
  int phenxIDColumns[1] = {phenxColumn};
  int dateColumns[1] = {dateColumn};
  return sequenceWorkflow(createDuration, removeSparseBuckets, inputFilePaths, inputFileDelimiter,
                          outputDir, outputFilesDescription, patIdColumns, phenxIDColumns,
                          dateColumns, patientCount, sparsity);
  
  
}
// [[Rcpp::export]]
int testSequencingworkflow(){
  omp_set_num_threads(16);
  std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_fourtimes_processed.csv";
  //    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_debug.csv";
  std::string description = "tspm_test_";
  std::string outputDir = "/home/jonas/CLionProjects/tspm_cpp_backend/out/data/";
  int patIdColumn = 0;
  int phenotypeIdColumn = 1;
  int dateColumn = 3;
  int patientCount = 1168; //TODO add function to determine (extract from last line in file)
  bool createDuration = true;
  bool removeSparseBuckets = true;
  double sparsity = 0.005;
  
  std::vector<std::string> inputFilePaths;
  inputFilePaths.push_back(fileName);
  
  char inputFileDelimiter = ',';
  int patIdColumns[1] = {0};
  int phenxIDColumns[1] = {1};
  int dateColumns[1] = {3};
  return sequenceWorkflow(createDuration, removeSparseBuckets, inputFilePaths, inputFileDelimiter,
                   outputDir, description, patIdColumns, phenxIDColumns,
                   dateColumns, patientCount, sparsity);
}
