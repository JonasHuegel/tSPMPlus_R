#' Transforms an alphanumeric data frame to a numeric one
#'
#'@description Transforms an alphanumeric data frame to a numeric one and returns it together with the look up tables for the patient id and phenx 
#'@param dbmart the alpha numeric dnmart that should be converted
#'@returns a list of data frames containing the numeric datamart and the patient id and phenx look up tables
#'
transformDbMartToNumeric<- function(dbmart){
  require(dplyr)
  require(tidyr)
  require(DT)
  require(data.table)
  patient_num <- c(unique(dbmart$patient_num))
  
  patLookUp <- as.data.frame(patient_num)
  
  setDT(patLookUp)
  
  patLookUp[,num_pat_num := .I]
  patLookUp <- patLookUp %>% mutate(num_pat_num = num_pat_num -1)
  
  dbmart_num <- dbmart %>% dplyr::left_join(patLookUp, by="patient_num")
  
  phenx <- c(unique(dbmart$phenx))
  
  phenxLookUp <- as.data.frame(phenx)
  
  setDT(phenxLookUp)
  
  phenxLookUp[,num_Phenx := .I]
  
  dbmart_num <- dbmart_num %>% dplyr::left_join(phenxLookUp, by="phenx")
  
  dbmart_num <- dbmart_num %>% select(num_pat_num, num_Phenx, start_date)
  
  dbmart_num <- dbmart_num[order(dbmart_num$num_pat_num,dbmart_num$start_date),]
  
  rownames(dbmart_num) = seq(length=nrow(dbmart_num))
  
  dbmart_num$start_date <- as.Date(dbmart_num$start_date)
  
  out <- list()
  out$patientLookUp <- patLookUp
  out$phenxLookUp <- phenxLookUp
  out$dbMart <- dbmart_num
  
  return (out)
}