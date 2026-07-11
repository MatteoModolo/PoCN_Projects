library(readr)
library(dplyr)
library(stringr)
library(igraph)

# Direwctories
input_dir <- "final_weighted_networks"
output_dir <- "degree_distribution"

# Create the output folder if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)

network_files <- list.files(path = input_dir, 
                            pattern = "^as_edges_.*_weighted\\.tsv\\.gz$", 
                            full.names = TRUE)

if (length(network_files) == 0) {
  stop("No files found in the 'final_weighted_networks' directory ")
}

cat(sprintf("Extracting degree distributions for %d networks\n", length(network_files)))


for (file_path in network_files) {
  
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  # Define the exact output file name
  output_csv <- file.path(output_dir, paste0("degree_dist_", yyyy_mm, ".csv"))
  
  # If the file already exists skip to save time
  if (file.exists(output_csv)) {
    cat(sprintf(" %s already processed \n", yyyy_mm))
    next
  }
  
  # Read the network
  edges_df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
  
  # Node mapping
  all_asns <- unique(c(edges_df$asn1, edges_df$asn2))
  id1 <- match(edges_df$asn1, all_asns)
  id2 <- match(edges_df$asn2, all_asns)
  
  # Build gtraph
  edge_matrix <- cbind(id1, id2)
  g <- graph_from_edgelist(edge_matrix, directed = FALSE)
  
  # Calculate degree distro
  deg_probs <- degree_distribution(g, mode = "all")
  
  # create dataframe
  degree_df <- data.frame(
    degree = 0:(length(deg_probs) - 1),
    probability = deg_probs
  )
  
  # clean 0 prob degree and save
  degree_df_clean <- degree_df %>% filter(probability > 0)
  
  write_csv(degree_df_clean, output_csv)
  
  # Print progress with the maximum degree found to confirm it worked
  max_deg <- max(degree_df_clean$degree)
  cat(sprintf(" - Saved %s: Highest degree found = %d\n", yyyy_mm, max_deg))
}

