library(rvest)
library(stringr)
library(dplyr)
library(readr)
library(data.table)

# links and directories
base_url <- "https://data.ris.ripe.net/rrc00/"
OUTPUT_DIR <- "final_weighted_networks"
STATUS_LOG <- "network_status_log.csv"

start_month <- "2016.05"
end_month <- "2026.04"
MIN_NODES <- 10


# Scrape the website to find all expected months
main_page <- read_html(base_url)
links <- main_page %>% html_nodes("a") %>% html_attr("href")
monthly_folders <- links[str_detect(links, "^\\d{4}\\.\\d{2}/$")]
monthly_folders <- monthly_folders[str_remove(monthly_folders, "/") >= start_month & str_remove(monthly_folders, "/") <= end_month]

results <- data.frame(month = character(), status = character(), stringsAsFactors = FALSE)

#  Cross-reference website data with local saved files
for (folder in monthly_folders) {
  raw_yyyy_mm <- str_remove(folder, "/")
  num_yyyy_mm <- str_replace(raw_yyyy_mm, "\\.", "")
  yyyy_mm <- paste0(substr(num_yyyy_mm, 1, 4), "_", substr(num_yyyy_mm, 5, 6))
  
  file_path <- file.path(OUTPUT_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz"))
  
  if (!file.exists(file_path)) {
    # File is on the website but not saved locally
    status <- "Corrupted/Missing"
  } else {
    # File exists locally, fast read just the ASNs to check node count
    edges_df <- tryCatch({
      fread(file_path, select = c(1, 2), colClasses = c("numeric", "numeric"), showProgress = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(edges_df) || nrow(edges_df) == 0) {
      status <- "Corrupted/Missing"
    } else {
      num_nodes <- length(unique(c(edges_df[[1]], edges_df[[2]])))
      
      if (num_nodes < MIN_NODES) {
        status <- "MissingData"
      } else {
        status <- "AnalysisReady"
      }
    }
  }
  
  cat(sprintf("Checked %s: %s\n", yyyy_mm, status))
  results <- rbind(results, data.frame(month = yyyy_mm, status = status))
}

write_csv(results, STATUS_LOG)