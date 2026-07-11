library(rvest)
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(data.table)

base_url <- "https://data.ris.ripe.net/rrc00/"
dest_dir <- "bview_data"

# Reference directory (Read-only)
REFERENCE_DIR <- "processed_networks"

# Output directory (Write-only)
NEW_OUTPUT_DIR <- "processed_networks_fast_run"

LOG_FILE <- "processing_log_fast_run.txt"
METADATA_FILE <- "metadata_bview_files_fast_run.csv"

start_month <- "2016.05" 
end_month <- "2026.04"

if (!dir.exists(dest_dir)) dir.create(dest_dir)
if (!dir.exists(NEW_OUTPUT_DIR)) dir.create(NEW_OUTPUT_DIR)

writeLines("BGP Processing Log - Safe Mode Start", LOG_FILE)

log_msg <- function(msg) {
  cat(msg, "\n")
  cat(paste(Sys.time(), "-", msg, "\n"), file = LOG_FILE, append = TRUE)
}

existing_files <- list.files(dest_dir, full.names = TRUE)

if (length(existing_files) > 0) {
  target_limit <- str_replace(end_month, "\\.", "") 
  file_dates <- str_extract(basename(existing_files), "\\d{8}")
  file_yyyy_mm <- substr(file_dates, 1, 6)
  
  excess_indices <- !is.na(file_yyyy_mm) & file_yyyy_mm > target_limit
  if (any(excess_indices)) {
    file.remove(existing_files[excess_indices])
  }
}

metadata <- data.frame(month_directory=character(), downloaded_filename=character(), status=character())
main_page <- read_html(base_url)
links <- main_page %>% html_nodes("a") %>% html_attr("href")

monthly_folders <- links[str_detect(links, "^\\d{4}\\.\\d{2}/$")]
monthly_folders <- monthly_folders[str_remove(monthly_folders, "/") >= start_month & str_remove(monthly_folders, "/") <= end_month]

for (folder in monthly_folders) {
  raw_yyyy_mm <- str_remove(folder, "/")
  num_yyyy_mm <- str_replace(raw_yyyy_mm, "\\.", "") 
  yyyy_mm <- paste0(substr(num_yyyy_mm, 1, 4), "_", substr(num_yyyy_mm, 5, 6)) 
  
  file_in_reference <- file.path(REFERENCE_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz"))
  file_in_new <- file.path(NEW_OUTPUT_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz"))
  
  if (file.exists(file_in_reference) || file.exists(file_in_new)) {
    metadata <- bind_rows(metadata, data.frame(month_directory=raw_yyyy_mm, downloaded_filename="SKIPPED", status="Already processed"))
    next 
  }
  
  month_url <- paste0(base_url, folder)
  month_page <- tryCatch({ read_html(month_url) }, error = function(e) { NULL })
  
  if (is.null(month_page)) next
  
  file_links <- month_page %>% html_nodes("a") %>% html_attr("href")
  bview_files <- file_links[str_detect(file_links, "^bview\\.\\d{8}\\.\\d{4}\\.gz$")]
  
  if (length(bview_files) == 0) next
  
  selected_file <- ifelse(paste0("bview.", num_yyyy_mm, "01.0000.gz") %in% bview_files, 
                          paste0("bview.", num_yyyy_mm, "01.0000.gz"), 
                          bview_files[1])
  
  dest_path <- file.path(dest_dir, selected_file)
  
  if (!file.exists(dest_path)) {
    log_msg(sprintf("Downloading: %s...", selected_file))
    status_msg <- tryCatch({ download.file(paste0(month_url, selected_file), destfile = dest_path, mode = "wb", quiet = TRUE); "Downloaded" }, error=function(e) "Error")
  } else {
    status_msg <- "Already exists"
  }
  metadata <- bind_rows(metadata, data.frame(month_directory=raw_yyyy_mm, downloaded_filename=selected_file, status=status_msg))
}

write_csv(metadata, METADATA_FILE)

bview_files <- list.files(dest_dir, pattern = "^bview\\..*\\.gz$", full.names = TRUE)

if (length(bview_files) > 0) {
  
  awk_script_path <- tempfile(fileext = ".awk")
  writeLines('
  BEGIN { FS="|"; OFS="\t" }
  {
    n = split($7, asns, " ");
    prev = "";
    for (i=1; i<=n; i++) {
      curr = asns[i];
      if (curr ~ /^[0-9]+$/) {
        if (prev != "" && prev != curr) {
          if (prev + 0 < curr + 0) print prev, curr;
          else print curr, prev;
        }
        prev = curr;
      }
    }
  }', awk_script_path)
  
  for (filepath in bview_files) {
    filename <- basename(filepath)
    date_str <- str_split(filename, "\\.")[[1]][2]
    yyyy_mm <- paste0(substr(date_str, 1, 4), "_", substr(date_str, 5, 6))
    
    file_in_reference <- file.path(REFERENCE_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz"))
    file_in_new <- file.path(NEW_OUTPUT_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz"))
    
    if (file.exists(file_in_reference) || file.exists(file_in_new)) {
      next
    }
    
    log_msg(sprintf("Processing %s...", filename))
    
    tmp_txt <- tempfile(fileext = ".txt")
    cmd <- paste0("bgpdump -m ", shQuote(filepath), " 2>/dev/null | awk -f ", shQuote(awk_script_path), " > ", shQuote(tmp_txt))
    
    system(cmd, wait = TRUE)
    
    edges_dt <- tryCatch({
      fread(tmp_txt, col.names = c("asn1", "asn2"), colClasses = c("integer", "integer"), nThread = 4)
    }, error = function(e) data.table(asn1=integer(), asn2=integer()))
    
    file.remove(tmp_txt) 
    
    total_edges <- nrow(edges_dt)
    
    if (total_edges == 0) {
      log_msg(sprintf("Month: %s | ERROR: No paths extracted.", yyyy_mm))
      next
    }
    
    log_msg(sprintf("Month: %s | Raw edges extracted: %d", yyyy_mm, total_edges))
    
    weighted_edges <- edges_dt[, .(weight = .N), by = .(asn1, asn2)]
    unweighted_edges <- weighted_edges[, .(asn1, asn2)]
    
    fwrite(weighted_edges, file.path(NEW_OUTPUT_DIR, paste0("as_edges_", yyyy_mm, "_weighted.tsv.gz")), sep = "\t", compress = "gzip")
    fwrite(unweighted_edges, file.path(NEW_OUTPUT_DIR, paste0("as_edges_", yyyy_mm, "_unweighted.tsv.gz")), sep = "\t", compress = "gzip")
  }
  
  file.remove(awk_script_path)
  log_msg("Pipeline completed.")
}