utils::globalVariables(c("num_pat_num","num_Phenx", "start_date"))

#'@title Transforms an alphanumeric data frame to a numeric one
#'@description Transforms an alphanumeric data frame to a numeric one and returns it together with the look up tables for the patient id and phenx.
#'@param dbmart The alpha numeric dbmart that should be converted.
#'@param originalPhenxLookUP The dataframe containing a subset of the phenx and their numeric representations, which should be also assigned to the phenx in this dbmart, Optional
#'@returns A list of data frames containing the numeric datamart and the patient id and phenx look up tables.
#'
transformDbMartToNumeric<- function(dbmart, originalPhenxLookUP){
  `%>%` <- magrittr::`%>%`
  patient_num <- c(unique(dbmart$patient_num))
  patLookUp <- as.data.frame(patient_num)
  setDT(patLookUp)
  patLookUp[,num_pat_num := .I]
  patLookUp <- patLookUp %>% dplyr::mutate(num_pat_num = num_pat_num -1)
  dbmart_num <- dbmart %>% dplyr::left_join(patLookUp, by="patient_num")
  
    
  phenx <- c(unique(dbmart$phenx))
  phenxLookUp <- as.data.frame(phenx)
  setDT(phenxLookUp)
  
  if(!missing(originalPhenxLookUP)){
    print("Found original PhenxLookUp table, using the numeric representation for phenx stored in the original lookUp.")
    print("CAUTION: This might resolve in gaps between the assigned numeric values, if the original lookup contains phenx that are not in the current dbmart.")
    originalPhenx <- subset(originalPhenxLookUP, originalPhenxLookUP$phenx %in% phenxLookUp$phenx)
    largestPhenx <- max(originalPhenxLookUP$num_Phenx)
    phenxLookUp <- subset(phenxLookUp,!(phenxLookUp$phenx %in% originalPhenx$phenx))
    phenxLookUp[,num_Phenx := .I + largestPhenx]
    phenxLookUp <- rbind(originalPhenx,phenxLookUp)
  }else{
    phenxLookUp[,num_Phenx := .I]
  }

  
  dbmart_num <- dbmart_num %>% dplyr::left_join(phenxLookUp, by="phenx")
  dbmart_num <- dbmart_num %>% dplyr::select(num_pat_num, num_Phenx, start_date)
  dbmart_num <- dbmart_num[order(dbmart_num$num_pat_num,dbmart_num$start_date),]
  rownames(dbmart_num) = seq(length=nrow(dbmart_num))
  dbmart_num$start_date <- as.Date(dbmart_num$start_date)
  
  
  out <- list()
  out$patientLookUp <- patLookUp
  out$phenxLookUp <- phenxLookUp
  out$dbMart <- dbmart_num
  
  return (out)
}

#'@title Split a numeric dbMart in multiple chunks that can be sequenced
#'@description Split a numeric dbMart in multiple chunks that can be sequenced. This function should be used if the original numeric dbmart contains to many entries. 
#'  Either due to memory limitation or the that the number of sequences is larger than 2**31-1, which is the maximum of entries in an R vector.
#'@param dbmart_num The numeric dbmart that should be splitt in chunks.
#'@param includeCorSeq A boolean parameter, should be true if the corseq flag will be set to true during sequencing. Default is false!
#'@param buffer a additional buffer how many additional bytes should not be considered as available. Default is 10 MB. 
#'@returns A list of 2 list. The first list contains the dbmart chunks, the second one contains the look up tables to translate the chunk-patnum to one from the orignal dbmart.
#'
splitdbMartInChunks <-function(dbmart_num, includeCorSeq = FALSE, buffer = 10000000){
  `%>%` <- magrittr::`%>%`
  uniquePatients <-unique(dbmart_num$num_pat_num)
  entriesPerPatient <- count(dbmart_num$num_pat_num)
  dbmartEntrySizeInCPP <- 16
  sequenceSizeInCPP <- 16
  sequencesSizeInR <- 20
  if(includeCorSeq == TRUE){
    sequencesSizeInR <- sequencesSizeInR + 16
  }

  entriesPerPatient <- entriesPerPatient %>% dplyr::mutate(seqs = (freq*(freq-1))/2)
  entriesPerPatient$seqs <- as.numeric(entriesPerPatient$seqs)
  entriesPerPatient <- entriesPerPatient %>% dplyr::mutate(mem = freq*dbmartEntrySizeInCPP + seqs*sequenceSizeInCPP + seqs*sequencesSizeInR)
  
  availMem <- as.numeric(memuse::Sys.meminfo()$freeram)
  dbMartSize <- as.numeric(object.size(dbmart_num))
  availMem <- availMem - dbMartSize - buffer
  
  memInChunk <- 0
  firstEntryInChunk <- c(0)
  pat <-0
  seqCount <- 0
  maxSeqCount <- 2**31-1
  #calculate chunk sizes
  for(pat in uniquePatients){
    pat <- pat+1 #(pats are starting with 0 for c++ so add one to use it )
    memInChunk <- memInChunk + entriesPerPatient$mem[pat]
    seqCount <- seqCount + entriesPerPatient$seqs[pat]
    if(memInChunk >= availMem ||seqCount >= maxSeqCount){
      if(entriesPerPatient$mem[pat] > availMem || entriesPerPatient$seqs[pat]> maxSeqCount){
        stop(paste0("Cannont split the dbMart in adaptive chunks, the required memory for one Patient exceeds the available memory! To much memory required for patient: ", pat-1))
      }
      firstEntryInChunk <- c(firstEntryInChunk, pat-1)
      memInChunk <-entriesPerPatient$mem[pat]
      seqCount <- entriesPerPatient$mem[pat]
    }
  }
  
  # chop chop
  chunks<- list()
  chunks$chunks <-list()
  chunks$lookUps <- list()
  for(i in 1:length(firstEntryInChunk)){
    start <- firstEntryInChunk[i];
    if(i < length(firstEntryInChunk)){
      stop <-firstEntryInChunk[i+1]
    }else{
      stop <- nrow(dbmart_num)
    }

    chunk <- dbmart_num %>% dplyr::filter(dbmart_num$num_pat_num >= start & dbmart_num$num_pat_num < stop)

    patientsInChunk <- c(unique(chunk$num_pat_num))
    
    chunkLookUp <- as.data.frame(patientsInChunk)
    colnames(chunkLookUp) <- c("num_pat_num")
    setDT(chunkLookUp)
    
    chunkLookUp[,chunk_pat_num := .I]
    chunkLookUp <- chunkLookUp %>% dplyr::mutate(chunk_pat_num = chunk_pat_num -1)
    chunks$chunks <- append(chunks$chunks, list(chunk %>% dplyr::left_join(chunkLookUp, by="num_pat_num") 
                                                %>% dplyr::select(num_pat_num=chunk_pat_num, num_Phenx, start_date)))
    chunks$lookUps <- append(chunks$lookUps, list(chunkLookUp))
  }  
  return(chunks)
}
