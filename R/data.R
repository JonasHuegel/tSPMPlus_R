#' Data set of synthetic patient data 
#'
#' A dataset containing synthetic patient data and used in the vignets for the long covid identification. 
#' The data set is based on the synthetic data covid data from Synthea. 
#' (Walonoski J, Klaus S, Granger E, Hall D, Gregorowicz A, Neyarapally G, Watson A, Eastman J. Synthea™ Novel coronavirus (COVID-19) model and synthetic data set. Intelligence-Based Medicine. 2020 Nov;1:100007. https://doi.org/10.1016/j.ibmed.2020.100007)
#'
#' @format A data frame with 670941 rows and 4 variables:
#' \describe{
#'   \item{start_date}{the date when this entry was recorded}
#'   \item{patient_num}{the corresponding patient number}
#'   \item{phenx}{the identifier for the entry}
#'   \item{description}{the longer description for the entry}
#' }
"dbmart"