library(readr)
library(dplyr)
library(stringr)


input_dir <- "final_weighted_networks"
output_detail_dir <- "edge_turnover_details"

# Create the output directory if it doesn't exist to store the monthly logs
if (!dir.exists(output_detail_dir)) {
  dir.create(output_detail_dir)
  cat(sprintf("Created directory: %s\n", output_detail_dir))
}

# Get files and sort them chronologically
network_files <- list.files(path = input_dir, pattern = "^as_edges_.*\\.tsv\\.gz$", full.names = TRUE)
network_files <- sort(network_files)

if (length(network_files) < 2) stop("Need at least 2 months of data to calculate turnover.")

cat(sprintf("Starting EDGE-LEVEL Churn Analysis for %d consecutive periods...\n", length(network_files) - 1))


file1 <- network_files[1]
month1_name <- str_extract(basename(file1), "\\d{4}_\\d{2}")

df1 <- read_tsv(file1, show_col_types = FALSE, progress = FALSE)

# Edges are already pre-sorted (min to max), so we use blazing-fast string concatenation
edges1 <- unique(paste0(df1$asn1, "_", df1$asn2))

# comparison loop
for (i in 2:length(network_files)) {
  
  file2 <- network_files[i]
  month2_name <- str_extract(basename(file2), "\\d{4}_\\d{2}")
  
  cat(sprintf("Tracking detailed edge shifts: %s -> %s... ", month1_name, month2_name))
  
  # Load Month 2
  df2 <- read_tsv(file2, show_col_types = FALSE, progress = FALSE)
  edges2 <- unique(paste0(df2$asn1, "_", df2$asn2))
  
  #  Isolate the exact edge identities using Set Theory
  persistent_edges   <- intersect(edges1, edges2)
  appearing_edges    <- setdiff(edges2, edges1)
  disappearing_edges <- setdiff(edges1, edges2)
  
  #  Compile them into a single, clean dataframe for this specific period
  period_df <- data.frame(
    edge_id = c(persistent_edges, appearing_edges, disappearing_edges),
    status = c(
      rep("persistent", length(persistent_edges)),
      rep("appearing", length(appearing_edges)),
      rep("disappearing", length(disappearing_edges))
    ),
    stringsAsFactors = FALSE
  )
  
  #  Save to the dedicated folder 
  period_name <- paste0("edge_tracking_", month1_name, "_to_", month2_name, ".csv")
  out_path <- file.path(output_detail_dir, period_name)
  
  write_csv(period_df, out_path)
  
  cat(sprintf("Saved %d edges.\n", nrow(period_df)))
  
  #  Carry over 2nd month as first one for the next iteration (Relay Handoff)
  month1_name <- month2_name
  edges1 <- edges2
}

