


# ====================== Extract effort over time (annual PubMed publications) between 1950 and 2020 (70 year timespan) ===========================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "RISmed", "taxize", "Hmisc")

# unique host species separated by host synonyms
spp = read.csv("./data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Host, HostClass, HostOrder, HostFamily, HostSynonyms) %>%
  dplyr::filter(!duplicated(Host)) %>%
  dplyr::arrange(desc(HostClass), HostOrder, HostFamily, Host) %>%
  tidyr::separate_rows(HostSynonyms, sep=",")





# ======================= function to extract for each species ====================

# function to scrape pubmed results
getPubMedYears = function(species_binomial, taxid_search = FALSE, virus_related = FALSE){
  
  # print species name
  print(sprintf("Processing: %s", species_binomial))
  
  # create search term with all synonyms
  datx = spp[ spp$Host == species_binomial, ]
  sppx = str_trim(c(datx$Host[1], Hmisc::capitalize(datx$HostSynonyms[ datx$HostSynonyms != ""])))
  sppx = strsplit(sppx, " ")
  sppx = sppx[ which(unlist(lapply(sppx, length)) < 3) ]
  
  if (!taxid_search) {
    # Search by name(s)
    createSearchTerm = function(x){ paste("(", x[1], " [TIAB] AND ", x[2], " [TIAB])", sep="") }
    search_term = paste(unlist(lapply(sppx, createSearchTerm)), collapse=" OR ")
    search_db = "pubmed"
    if(virus_related){
      search_term = paste(search_term, "AND (virus [TIAB] OR viral [TIAB])")
    }
      
  } else {
    # Search by taxid
    taxonomy_id = taxize::get_uid(species_binomial, rank_query = "species", 
                                  ask = FALSE, messages = FALSE)
    
    if (attr(taxonomy_id, "multiple_matches"))
      warning("Multiple names matched for ", species_binomial, ", returning NA")
    
    if (attr(taxonomy_id, "pattern_match")) {
      # Don't want partial matches
      warning("Partial match for", species_binomial, ", returning NA")
      taxonomy_id = NA_character_
    }
    
    search_term = paste0("txid", taxonomy_id,"[Organism:exp]") # Papers mentioning this species or any subspecies
    search_db = "PMC"  # This query only works on Pubmed Central
  }
  
  
  # create pubmed search term
  #search_term = paste(strsplit(species_binomial," ")[[1]][1], "[TIAB]", "AND", strsplit(species_binomial," ")[[1]][2], "[TIAB]", sep=" ")
  
  # run searches: firstly to get number of publications, then set this as retmax to ensure all are returned 
  e = simpleError("test error")
  search1 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db=search_db, datetype='pdat', mindate=1930, maxdate=2020), error=function(e) e)
  if(class(search1)[1] == "simpleError"){ return(data.frame(Year=NA, NumPubs=NA, Host = species_binomial, Note="Lookup error")) }
  numpubs_total = RISmed::QueryCount(search1)
  if(numpubs_total == 0){ return(data.frame(Year=2018, NumPubs=0, Host = species_binomial, Note="No publications")) }
  
  # second search with numpubs_total set as max
  Sys.sleep(0.5)
  search2 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db=search_db, datetype='pdat', mindate=1930, maxdate=2020, retmax=numpubs_total), error=function(e) e)
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
output_loc = "./output/"
save_file = paste(output_loc, "PubMed_HostsEffort_PerYear_1930.csv", sep="")

# for all host species get time series of publication effort (axis axis == impossible to lookup)
species_vec = unique(spp$Host)
#species_vec = species_vec[ !species_vec %in% c("Axis axis", "Indicator indicator", "Canis lupus familiaris")]

# re-run for lookup error
# rr = read.csv(paste(output_loc, "PubMed_Hosts_PubsByYear_22102020_syns.csv", sep=""), stringsAsFactors = FALSE)
# species_vec = rr$Host[ rr$Note == "Lookup error" ]

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
# rr1 = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_22102020_syns.csv", stringsAsFactors = FALSE)
# rr2 = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_22102020_syns_update.csv", stringsAsFactors = FALSE)
# rr = rbind(rr1, rr2)
# write.csv(rr, "./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_Final.csv", row.names = FALSE)



# ============== run virus-related scrape and append records to csv ============

# create filenames
output_loc = "./output/"
save_file = paste(output_loc, "PubMed_HostsEffort_PerYear_VirusRelated.csv", sep="")

# for all host species get time series of publication effort
species_vec = unique(spp$Host)
#species_vec = species_vec[ !species_vec %in% c("Axis axis", "Indicator indicator", "Canis lupus familiaris")]

# re-run for lookup error
# rr = read.csv(paste(output_loc, "PubMed_Hosts_PubsByYear_22102020_syns.csv", sep=""), stringsAsFactors = FALSE)
# species_vec = rr$Host[ rr$Note == "Lookup error" ]

# append each new query to csv
for(i in 1:length(species_vec)){
  
  # run query (5 attempts)
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  for(attempt in 1:5){
    
    resx = tryCatch(getPubMedYears(species_vec[i], virus_related = TRUE), error=function(e) e)
    if(class(resx)[1] != "simpleError"){ 
      break
    } else{
      Sys.sleep(0.8)
      resx = tryCatch(getPubMedYears(species_vec[i], virus_related = TRUE), error=function(e) e)
    }
  }
  
  # return lookup error if still throwing errors
  if(class(resx)[1] == "simpleError"){ 
    resx = data.frame(Year=NA, NumPubs=NA, Host = species_vec[i], Note="Lookup error")
  }
  
  # add "virus related" field
  resx$VirusRelated = TRUE
  
  # initialise file on first iteration, and then append
  if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
}




# =========================== Create compiled data for other projects =========================

# hx = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_19302020.csv") %>%
#   dplyr::group_by(Host) %>%
#   dplyr::summarise(Pubs_All = sum(NumPubs))
# 
# hy = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_19302020.csv") %>%
#   dplyr::group_by(Host) %>%
#   dplyr::summarise(Pubs_VirusRelated = sum(NumPubs))
# 
# effort = left_join(hx, hy)
# effort = left_join(effort, spp[ !duplicated(spp$Host), c("Host", "HostOrder", "HostFamily")])
# write.csv(effort, "./output/host_effort/PubMed_HostCounts_Total_CLOVER.csv", row.names=FALSE)
# 
# ggplot(effort) + 
#   geom_boxplot(aes(HostOrder, Pubs_VirusRelated)) +
#   theme(axis.text.x = element_text(angle=90))










