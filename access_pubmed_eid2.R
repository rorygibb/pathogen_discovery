
# ============= Scrape NCBI databases for EID2 publication details inc. year ============

# dependencies and basedir
pacman::p_load("RISmed", "dplyr", "magrittr")
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")




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

# combine into dataframe to query (176,000 associations)
query_df = rbind(pm, nuc)
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
  
  # return error result if still throwing error
  #if(class(search)[1] == "simpleError"){  
  if(class(search)[1] != "Medline"){
    res = data.frame(
      id = query_df$id[x],
      Database = query_df$Database[x],
      PMID = NA,
      Year = NA,
      Journal = NA,
      Title = NA,
      Date_NCBI = NA,
      PathName_NCBI = NA
    ) 
    Sys.sleep(0.5)
    return(res)
  }
  
  # PubMed
  if(query_df$Database[x] == "pubmed"){
    
    # if no records exist
    if(length(search@PMID) == 0){
      res = data.frame(
        id = query_df$id[x],
        Database = query_df$Database[x],
        PMID = "no record",
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI = NA,
        PathName_NCBI = NA
      ) 
    } else{
      # otherwise return records
      res = data.frame(
        id = query_df$id[x],
        Database = query_df$Database[x],
        PMID = search@PMID,
        Year = search@YearPubmed,
        Journal = search@MedlineTA,
        Title = search@ArticleTitle,
        Date_NCBI = NA,
        PathName_NCBI = NA
      )
    }
  }
  
  # Nucleotide
  if(query_df$Database[x] == "nuccore") {
    
    if(length(as.vector(search[7])) == 0){
      res = data.frame(
        id = query_df$id[x],
        Database = query_df$Database[x],
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI = "no record",
        PathName_NCBI = NA
      ) 
    } else{
      res = data.frame(
        id = query_df$id[x],
        Database = query_df$Database[x],
        PMID = NA,
        Year = NA,
        Journal = NA,
        Title = NA,
        Date_NCBI = as.vector(search[7]),
        PathName_NCBI = as.vector(search[9]) 
        )
    }
  }
  
  # sleep for 0.25 seconds to prevent over-requesting on API (might be able to get away with shorter)
  Sys.sleep(0.5)
  return(res)
}


# ============== run scrape and append records to csv ============

# create filenames
output_loc = "./output/data_processed/pathogens/eid2_scrape/"
save_file = paste(output_loc, "EID2_Scrape_16102020.csv", sep="")

# append each new query to csv
for(i in 1:nrow(query_df)){

  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  resx = tryCatch(searchNCBI(i))
  
  # create save file or append
  if(class(resx)[1] == "simpleError"){ next 
  } else if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE) # append
  }
}









# ================= OLD =======================

# # ============== run scrape in increments of 500, saving output after each ===============
# 
# # create filenames
# output_loc = "./output/data_processed/pathogens/"
# 
# # will take ~24 hours to do 176,000 (EID2 mammal records)
# # increments and output loc
# inc = seq(1, nrow(query_df), by=500)
# 
# # run scrape in increments of 500 (deal with errors, save progress)
# for(i in 1:length(inc)){
#   
#   # print
#   print(paste(inc[i], "...", sep=""))
#   
#   # specify records to search; if file already exists (e.g. running second time round) skip to next
#   if(i < length(inc)){ range = inc[i]:(inc[i+1]-1) }
#   if(i == length(inc)){ range = inc[i]:nrow(query_df) }
#   filename = paste(output_loc, "EID2Year_Mammals_", max(range), ".R", sep="")
#   if(file.exists(filename)){ next }
#   
#   # run scrape and print time after each 500 iterations
#   sx = Sys.time()
#   e = simpleError("test error")
#   resx = tryCatch(lapply(range, searchNCBI), error=function(e) e)
#   if(class(resx)[1] == "simpleError"){ next 
#   } else{
#     save(resx, file=filename)
#     }
#   ex = Sys.time()
#   print(ex-sx)
# }
# 
# 
# # ======================= second run: identify records where lookup was unsuccessful and try again ====================
# 
# # load and combine into df
# ff = list.files(output_loc, pattern=".R", full.names = TRUE)
# result = data.frame()
# for(i in 1:length(ff)){ 
#   load(ff[i]) 
#   resx = do.call(rbind.data.frame, resx)
#   result = rbind(result, resx)
# }
# 
# # remove duplicates and combine
# result = result[ !duplicated(result$id), ]
# eid = left_join(query_df, result, by=c("id"="id", "Database"="Database"))
# eid$record = 1:nrow(eid)
# #write.csv(eid, "./output/data_processed/EID2_scrape_v1_gaps.csv", row.names=FALSE)
# 
# 
# # ------------- extract gaps for mammals -------------
# 
# # NAs for mammals (not human)
# eid = eid[ eid$Carrier.classification != "Human", ]
# query_df = eid[ which(eid$Database == "pubmed" & is.na(eid$PMID) | eid$Database == "nuccore" & is.na(eid$Date_NCBI)), ]
# 
# # create filenames
# output_loc = "./output/data_processed/pathogens/update/"
# 
# # increments and output loc
# inc = seq(1, nrow(query_df), by=500)
# 
# # run scrape in increments of 500 (deal with errors, save progress)
# for(i in 1:length(inc)){
#   
#   # print
#   print(paste(inc[i], "...", sep=""))
#   
#   # specify records to search; if file already exists (e.g. running second time round) skip to next
#   if(i < length(inc)){ range = inc[i]:(inc[i+1]-1) }
#   if(i == length(inc)){ range = inc[i]:nrow(query_df) }
#   filename = paste(output_loc, "EID2Year_Mammals_update_", max(range), ".R", sep="")
#   if(file.exists(filename)){ next }
#   
#   # run scrape and print time after each 500 iterations
#   sx = Sys.time()
#   e = simpleError("test error")
#   resx = tryCatch(lapply(range, searchNCBI), error=function(e) e)
#   if(class(resx)[1] == "simpleError"){ next 
#   } else{
#     save(resx, file=filename)
#   }
#   ex = Sys.time()
#   print(ex-sx)
# }
# 
# # load and compile
# ff = list.files(output_loc, pattern=".R", full.names = TRUE)
# result = data.frame()
# for(i in 1:length(ff)){ 
#   load(ff[i]) 
#   resx = do.call(rbind.data.frame, resx)
#   result = rbind(result, resx)
# }
# 
# # subset to get records that still don't have lookup successfully
# eid1 = eid[! eid$id %in% query_df$id, ]
# eid2 = eid[ (eid$id %in% query_df$id) & (!eid$id %in% result$id), ] # still missing
# eid3 = left_join(eid[ (eid$id %in% query_df$id) & (eid$id %in% result$id), c(1:8, 15)], result[ !duplicated(result$id), ])
# 
# # fine except a few; lookup
# eid = rbind(eid1, eid3)
# 
# 
# 
# # -------------- run again to fill last gaps -----------------
# 
# # create filenames
# output_loc = "./output/data_processed/"
# query_df = eid2
# 
# # create save file 
# for(i in 1:nrow(query_df)){
# #for(i in 1:50){
#     
#   cat(paste(i, "...", sep=""))
#   e = simpleError("test error")
#   resx = tryCatch(searchNCBI(i))
#   if(class(resx)[1] == "simpleError"){ next 
#   } else if(i == 1){
#     write.csv(resx, paste(output_loc, "EID_scrape_16102020.csv", sep=""), row.names=FALSE)
#   } else{
#     write.table(resx, paste(output_loc, "EID_scrape_16102020.csv", sep=""), append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE) # append
#   }
# }
# 
# # load outcome
# rr = read.csv("./output/data_processed/EID_scrape_16102020.csv", stringsAsFactors = FALSE)
# 
# # combine with EID2
# eid2x = eid2[ , c(1:8, 15) ]
# eid2x = left_join(eid2x, rr[ !duplicated(rr$id), ])
# 
# # combine all
# eid_fin = rbind(eid, eid2x)
# write.csv(eid_fin, "./output/data_processed/pathogens/EID2_Year_PubMedscrape_16102020.csv", row.names=FALSE)
# 
# library(RISmed)
# library(dplyr)
# library(magrittr)
# 
# setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")


# 
# # increments and output loc
# query_df = eid2
# inc = seq(1, nrow(query_df), by=500)
# 
# # run scrape in increments of 500 (deal with errors, save progress)
# for(i in 1:length(inc)){
#   
#   # print
#   print(paste(inc[i], "...", sep=""))
#   
#   # specify records to search; if file already exists (e.g. running second time round) skip to next
#   if(i < length(inc)){ range = inc[i]:(inc[i+1]-1) }
#   if(i == length(inc)){ range = inc[i]:nrow(query_df) }
#   filename = paste(output_loc, "EID2Year_Mammals_update_", max(range), ".R", sep="")
#   if(file.exists(filename)){ next }
#   
#   # run scrape and print time after each 500 iterations
#   sx = Sys.time()
#   e = simpleError("test error")
#   resx = tryCatch(lapply(range, searchNCBI), error=function(e) e)
#   if(class(resx)[1] == "simpleError"){ next 
#   } else{
#     save(resx, file=filename)
#   }
#   ex = Sys.time()
#   print(ex-sx)
# }
# 
# # load and compile
# ff = list.files(output_loc, pattern=".R", full.names = TRUE)
# result = data.frame()
# for(i in 1:length(ff)){ 
#   load(ff[i]) 
#   resx = do.call(rbind.data.frame, resx)
#   result = rbind(result, resx)
# }
# 
