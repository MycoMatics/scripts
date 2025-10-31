# the input is a UTF8 cs file with 1 colomn named taxon and in that column taxon names that you want to fuzzy match to the GBIF bacbbone
# match taxa names to GBIF backbone 
library(rgbif)
library(dplyr)
library(tidyverse)
library(openxlsx)
taxon_list <- read_csv("D:/Givdieri/OneDrive - UGent/UGent-PC/Givdieri/Documents/Projects/doctoraat_Fagus_deadwood/data/data_identifications/ID_GBIF4.csv")
fetch_gbif_info <- function(taxon_name) {
  result <- tryCatch(
    {
      name_backbone(name = taxon_name)
    },
    error = function(e) {
      warning(paste("Error querying GBIF for taxon:", taxon_name, "-", e$message))
      return(NULL)
    }
  )
  
  # Handle unmatched taxa or errors
  if (is.null(result) || length(result) == 0 || is.null(result$scientificName)) {
    return(data.frame(
      taxon = taxon_name,
      scientificName = NA,
      latestName = NA,
      author_year = NA,
      matchType = "No match",
      status = NA,
      note = "Unmatched or error",
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract latest accepted name if the provided name is a synonym
  latest_name <- if (!is.null(result$acceptedUsageKey)) {
    accepted <- name_usage(key = result$acceptedUsageKey)
    if (!is.null(accepted) && !is.null(accepted$data)) {
      accepted$data$scientificName
    } else {
      NA
    }
  } else {
    result$scientificName
  }
  
  # Return the taxon information
  data.frame(
    taxon = taxon_name,
    scientificName = result$scientificName,
    latestName = latest_name,
    author_year = ifelse(!is.null(result$authorship), result$authorship, NA),
    matchType = ifelse(!is.null(result$matchType), result$matchType, NA),
    status = ifelse(!is.null(result$status), result$status, NA),
    note = if (!is.null(result$acceptedUsageKey)) "Synonym, resolved" else "Accepted name",
    stringsAsFactors = FALSE
  )
}

# Apply the function to the taxon list
taxon_info <- bind_rows(lapply(taxon_list$taxon, fetch_gbif_info))

# Merge the results with the original taxon list
taxon_complete <- taxon_list %>%
  left_join(taxon_info, by = "taxon")

write.xlsx(
  taxon_complete, 
  file = "taxon_with_authorship_synonymy4.xlsx"
)

