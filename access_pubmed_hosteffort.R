

# ====================== Extract effort over time (annual PubMed publications) between 1950 and 2020 (70 year timespan) ===========================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2")

# full harmoniesd associations database (EID2, GMPD2, HP3)
dd = read.csv("./output/data_processed/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
dd = dd[ !is.na(dd$Database), ]
dd = dd[ !is.na(dd$ParasiteType), ]
dd = dd[ !is.na(dd$Host_Harmonised), ]



# ======================= function to extract for each species ====================

# dependencies
library(RISmed)

# function to scrape pubmed results
getPubMedYears = function(species_binomial){
  
  # print species name
  print(sprintf("Processing: %s", species_binomial))
  
  # create pubmed search term
  search_term = paste(strsplit(species_binomial," ")[[1]][1], "[TIAB]", "AND", strsplit(species_binomial," ")[[1]][2], "[TIAB]", sep=" ")
  
  # run searches: firstly to get number of publications, then set this as retmax to ensure all are returned 
  e = simpleError("test error")
  search1 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2020), error=function(e) e)
  if(class(search1)[1] == "simpleError"){ return(data.frame(Host = species_binomial, Year=NA, NumPubs=NA, Note="Lookup error")) }
  numpubs_total = RISmed::QueryCount(search1)
  if(numpubs_total == 0){ return(data.frame(Host = species_binomial, Year=2020, NumPubs=0, Note="No publications")) }
  
  # second search with numpubs_total set as max
  Sys.sleep(0.5)
  search2 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2020, retmax=numpubs_total), error=function(e) e)
  Sys.sleep(0.5)
  if(class(search2)[1] == "simpleError"){ return(data.frame(Host = species_binomial, Year=NA, NumPubs=NA, Note="Lookup error")) }
  
  # extract years of all pubs from pubmed
  years = YearPubmed(RISmed::EUtilsGet(search2))
  if(class(years)[1] == "simpleError"){ return(data.frame(Host = species_binomial, Year=NA, NumPubs=NA, Note="Lookup error")) }
  
  # otherwise summarise
  resx = as.data.frame(table(years)) %>%
    dplyr::rename("Year"=years, "NumPubs"=Freq) %>%
    dplyr::mutate(Host = species_binomial,
                  Year = as.numeric(as.vector(Year)))

  # sleep for 0.25 seconds to prevent over-requesting
  Sys.sleep(0.25)
  return(resx)
}

# ============== run scrape and append records to csv ============

# create filenames
output_loc = "./output/data_processed/host_effort/"
save_file = paste(output_loc, "PubMed_Hosts_PubsByYear_22102020.csv", sep="")

# for all host species get time series of publication effort
species_vec = unique(dd$Host_Harmonised)

# append each new query to csv
for(i in 90:length(species_vec)){

  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  resx = tryCatch(getPubMedYears(species_vec[i]))
  
  # initialise file on first iteration, and then append
  if(class(resx)[1] == "simpleError"){ next
  } else if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
}
