#
# Retreive dates for nucleotide records in EID2
#

library(readr)
library(dplyr)
library(stringr)
library(rentrez)
library(retry)

batch_size = 100
max_retries = 5  # Number of reconnection attempts

api_key <- readline("NCBI API key: ")
eid_data <- read_csv("data/eid2/EID2_NucleotideRecords.csv", col_types = cols(.default = "c"))

stopifnot(all(eid_data$Database == "nuccore")) # Code below assumes only this db needs to be checked
identifiers <- unique(eid_data$id)

pub_dates <- data.frame()

for (batch_start in seq(1, length(identifiers), by = batch_size)) {
  batch_stop = min(batch_start + batch_size, length(identifiers))
  cat(sprintf("%d-%d of %d\r", batch_start, batch_stop, length(identifiers)))
  
  retry({
    post_handle <- entrez_post(db = "nuccore", id = identifiers[batch_start:batch_stop], api_key = api_key)
    
    accessions <- entrez_fetch(db = "nuccore", web_history = post_handle, rettype = "acc", api_key = api_key)
    summary_data <- entrez_summary(db = "nuccore", web_history = post_handle, api_key = api_key)
  },
    when = ".",
    max_tries = max_retries,
    interval = 60
  )
  
  accessions <- str_remove(accessions, "\n$") %>%   # trailing \n
    str_split("\n", simplify = TRUE) %>% 
    as.vector()
  
  dates <- extract_from_esummary(summary_data, "createdate") %>% 
    unname()
  
  pub_dates <- rbind(pub_dates, 
                     data.frame(id = identifiers[batch_start:batch_stop],
                                accession = accessions,
                                creation_date = dates))
  
}


# Check for duplicates (all versions of a sequence return the same create_date)
pub_dates <- pub_dates %>% 
  rename(accession_version = .data$accession) %>% 
  mutate(accession = str_remove(.data$accession_version, "\\.[[:digit:]]$"))

if (n_distinct(pub_dates$accession) != nrow(pub_dates))
  warning("Some accessions are duplicated", call. = FALSE)

# Output
pub_dates %>% 
  select(.data$id, .data$accession_version, .data$accession, .data$creation_date) %>% 
  write_excel_csv("data/eid2/eid2_nucleotide_dates.csv")
