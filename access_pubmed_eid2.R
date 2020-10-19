
# ============= Scrape NCBI databases for EID2 publication details inc. year ============

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")
pacman::p_load("RISmed", "dplyr", "magrittr")




# ================ load data and prepare for PubMed query =================

# eid2 for mammals
eid2 = read.csv("./data/host_pathogen_2020/data/hostparasitedb_raw/EID2/Wardehetal_2015_EID2/SpeciesInteractions_EID2.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Carrier.classification %in% c("Human", "Mammal", "Domestic", "Primate", "Rodent"))

# EID2 data with PMIDs (stored in 'Publications' column)
# 'database' field corresponds to exact name of db in NCBI databases (used to lookup relevant db)
pm = eid2[ eid2$Publications != "", -which(names(eid2) == "Sequences")] %>%
  tidyr::separate_rows(Publications, sep=";") %>%
  dplyr::rename("id" = Publications)
pm$Database = "pubmed"

# data with NCBI nucleotide IDs
# set database to nuccore
nuc = eid2[ eid2$Sequences != "", -which(names(eid2) == "Publications")] %>%
  tidyr::separate_rows(Sequences, sep=";") %>%
  dplyr::rename("id" = Sequences)
nuc$Database = "nuccore"

# combine into dataframe to query (147,000 unique ids)
query_df = rbind(pm, nuc)
eid_df = query_df # for later

# prep for scrape
query_df = query_df[ !duplicated(query_df$id), ]
query_df$record = 1:nrow(query_df)



# ================== extract PubMed publications for associations with PMIDs ==================

# function to scrape pubmed results
#' @param x 1:nrow(query_df), i.e. the nth row of query_df

searchNCBI = function(x){
  
  # print identifier
  print(sprintf("Processing: %s", query_df$id[x]))
  
  # run scrape (5 attempts and exit if successful)
  e = simpleError("test error")
  for(attempt in 1:5){
    search = tryCatch(RISmed::EUtilsGet(query_df$id[x], type="efetch", db=query_df$Database[x]), error=function(e) e)
    if(class(search)[1] != "simpleError"){ break }
    Sys.sleep(0.5)
  }
  
  # PubMed
  if(query_df$Database[x] == "pubmed"){
    
    # if not returned successfully
    if(class(search)[1] != "Medline"){
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = FALSE,
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI1 = NA,
        Date_NCBI2 = NA,
        PathName_NCBI = NA
      ) 
      Sys.sleep(0.5)
      return(res)
    }
    
    # if no records exist
    if(length(search@PMID) == 0){
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID = "no record",
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI1 = NA,
        Date_NCBI2 = NA,
        PathName_NCBI = NA
      ) 
    } else{
      # otherwise return records
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID = search@PMID,
        Year = search@YearPubmed,
        Journal = search@MedlineTA,
        Title = search@ArticleTitle,
        Date_NCBI1 = NA,
        Date_NCBI2 = NA,
        PathName_NCBI = NA
      )
    }
  }
  
  # Nucleotide
  if(query_df$Database[x] == "nuccore") {
    
    # error
    if(class(search)[1] == "simpleError"){  
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = FALSE,
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI1 = NA,
        Date_NCBI2 = NA,
        PathName_NCBI = NA
      ) 
      Sys.sleep(0.5)
      return(res)
    }
    
    # any vs no records
    if(length(as.vector(search[7])) == 0){
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI1= "no record",
        Date_NCBI2= "no record",
        PathName_NCBI = NA
      ) 
    } else{
      res = data.frame(
        record = query_df$record[x],
        id = query_df$id[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI1 = as.vector(search[7]),
        Date_NCBI2 = as.vector(search[8]),
        PathName_NCBI = as.vector(search[9]) 
        )
    }
  }
  
  # sleep for 0.25 seconds to prevent over-requesting on API (might be able to get away with shorter)
  Sys.sleep(0.5)
  return(res)
}


# ============== run scrape and append records to csv ============

# # create filenames
# output_loc = "./output/data_processed/pathogens/eid2_scrape/"
# save_file = paste(output_loc, "EID2_FullScrape_19102020.csv", sep="")
# 
# # append each new query to csv
# for(i in 1:nrow(query_df)){
# 
#   # run query
#   cat(paste(i, "...", sep=""))
#   e = simpleError("test error")
#   resx = tryCatch(searchNCBI(i))
# 
#   # initialise file on first iteration, and then append
#   if(class(resx)[1] == "simpleError"){ next
#   } else if(i == 1){
#     write.csv(resx, save_file, row.names=FALSE)
#   } else{
#     write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
#   }
# }


# ============== combine with EID2 records ==============

# combine and save
scrape = read.csv("./output/data_processed/pathogens/eid2_scrape/EID2_FullScrape_19102020.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(id = as.character(id))
eid_dfx = left_join(eid_df, scrape, by=c("id", "Database")) %>%
  filter(!is.na(Lookup_Successful))
write.csv(eid_dfx, "./output/data_processed/pathogens/EID2_Year_full.csv", row.names=FALSE)

