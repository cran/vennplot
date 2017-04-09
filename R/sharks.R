#' Data on human encounters with great white sharks.
#'
#' @format A dataset of logicals with 65 rows and 9 columns.  Each entry classifies whether or not the encounter (row) was of a given type (column):
#' \describe{
#'   \item{Male}{Male victim}
#'   \item{Female}{Female victim}
#'   \item{AM}{Morning encounter}
#'   \item{Australia}{Encounter in Australia}
#'   \item{USA}{Encounter in the United States}
#'   \item{Surfing}{Surfing incident}
#'   \item{Scuba}{Scuba-diving incident}
#'   \item{Fatality}{Whether or not there was a fatality}
#'   \item{Injury}{Whether or not there was an injury}
#' }
#' @source \url{http://sharkattackinfo.com/shark_attack_news_sas.html}.  Data collected by Professor Pierre-Jerome Bergeron, University of Ottawa.
#' @examples
#' combinations = sharks[,c(1,3:5,8)]
#' vennplot(combinations = combinations)
"sharks"
