

# ====================== Extract effort over time (annual PubMed publications) between 1950 and 2020 (70 year timespan) ===========================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "RISmed")

# full harmoniesd associations database (EID2, GMPD2, HP3)
dd = read.csv("./output/data_processed/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
dd = dd[ !is.na(dd$Database), ]
dd = dd[ !is.na(dd$ParasiteType), ]
dd = dd[ !is.na(dd$Host_Harmonised), ]

# unique species separated by host synonyms
spp = dd[ !duplicated(dd$Host_Harmonised), ] %>%
  dplyr::arrange(desc(HostClass), HostOrder, HostFamily, Host_Harmonised) %>%
  tidyr::separate_rows(HostSynonyms, sep=",")



# ======================= function to extract for each species ====================

# dependencies
library(RISmed)

# function to scrape pubmed results
getPubMedYears = function(species_binomial){
  
  # print species name
  print(sprintf("Processing: %s", species_binomial))
  
  # create search term with all synonyms
  datx = spp[ spp$Host_Harmonised == species_binomial, ]
  sppx = str_trim(c(datx$Host_Harmonised[1], datx$HostSynonyms[ datx$HostSynonyms != ""]))
  sppx = strsplit(sppx, " ")
  sppx = sppx[ which(unlist(lapply(sppx, length)) < 3) ]
  createSearchTerm = function(x){ paste("(", x[1], " [TIAB] AND ", x[2], " [TIAB])", sep="") }
  search_term = paste(unlist(lapply(sppx, createSearchTerm)), collapse=" OR ")
  
  # create pubmed search term
  #search_term = paste(strsplit(species_binomial," ")[[1]][1], "[TIAB]", "AND", strsplit(species_binomial," ")[[1]][2], "[TIAB]", sep=" ")
  
  # run searches: firstly to get number of publications, then set this as retmax to ensure all are returned 
  e = simpleError("test error")
  search1 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2020), error=function(e) e)
  if(class(search1)[1] == "simpleError"){ return(data.frame(Year=NA, NumPubs=NA, Host = species_binomial, Note="Lookup error")) }
  numpubs_total = RISmed::QueryCount(search1)
  if(numpubs_total == 0){ return(data.frame(Year=2018, NumPubs=0, Host = species_binomial, Note="No publications")) }
  
  # second search with numpubs_total set as max
  Sys.sleep(0.5)
  search2 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2020, retmax=numpubs_total), error=function(e) e)
  Sys.sleep(0.5)
  if(class(search2)[1] == "simpleError"){ return(data.frame(Year=NA, NumPubs=NA, Host = species_binomial, Note="Lookup error")) }
  
  # extract years of all pubs from pubmed
  years = YearPubmed(RISmed::EUtilsGet(search2))
  if(class(years)[1] == "simpleError"){ return(data.frame(Year=NA, NumPubs=NA, Host = species_binomial, Note="Lookup error")) }
  
  # otherwise summarise
  resx = as.data.frame(table(years)) %>%
    dplyr::rename("Year"=years, "NumPubs"=Freq) %>%
    dplyr::mutate(Host = species_binomial,
                  Year = as.numeric(as.vector(Year)))
  resx$Note = ""

  # sleep for 0.25 seconds to prevent over-requesting
  Sys.sleep(0.25)
  return(resx)
}


# ============== run scrape and append records to csv ============

# create filenames
output_loc = "./output/data_processed/host_effort/"
save_file = paste(output_loc, "PubMed_Hosts_PubsByYear_22102020_syns_update.csv", sep="")

# for all host species get time series of publication effort
species_vec = unique(spp$Host_Harmonised)
species_vec = species_vec[ !species_vec %in% c("axis axis", "indicator indicator")]

# re-run for lookup error
rr = read.csv(paste(output_loc, "PubMed_Hosts_PubsByYear_22102020_syns.csv", sep=""), stringsAsFactors = FALSE)
species_vec = rr$Host[ rr$Note == "Lookup error" ]

# append each new query to csv
for(i in 1:length(species_vec)){

  # run query (5 attempts)
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  for(attempt in 1:5){
    
    resx = tryCatch(getPubMedYears(species_vec[i]), error=function(e) e)
    if(class(resx)[1] != "simpleError"){ 
      break
    } else{
      Sys.sleep(0.5)
      resx = tryCatch(getPubMedYears(species_vec[i]), error=function(e) e)
      }
    }
  
  # return lookup error if still throwing errors
  if(class(resx)[1] == "simpleError"){ 
    resx = data.frame(Year=NA, NumPubs=NA, Host = species_vec[i], Note="Lookup error")
  }
  
  # initialise file on first iteration, and then append
  if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
}

# combine
rr1 = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_22102020_syns.csv", stringsAsFactors = FALSE)
rr2 = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_22102020_syns_update.csv", stringsAsFactors = FALSE)
rr = rbind(rr1, rr2)
write.csv(rr, "./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_Final.csv", row.names = FALSE)
