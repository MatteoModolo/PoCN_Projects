library(readr)
library(dplyr)
library(stringr)
library(igraph)
library(ggplot2)

# Directories
input_dir <- "final_weighted_networks"
output_csv <- "top_100_components_monthly.csv"
fragmented_csv <- "fragmented_networks_log.csv" 

network_files <- list.files(path = input_dir, 
                            pattern = "^as_edges_.*_weighted\\.tsv\\.gz$", 
                            full.names = TRUE)

if (length(network_files) == 0) {
  stop("No files found in the 'final_weighted_networks' directory")
}

results_list <- list()
fragmented_list <- list()
cat(sprintf("Scanning %d networks to extract top 100 components \n", length(network_files)))

# Compute the top components
for (file_path in network_files) {
  
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  edges_df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
  
  # FAST NODE MAPPING
  all_asns <- unique(c(edges_df$asn1, edges_df$asn2))
  id1 <- match(edges_df$asn1, all_asns)
  id2 <- match(edges_df$asn2, all_asns)
  
  edge_matrix <- cbind(id1, id2)
  g <- graph_from_edgelist(edge_matrix, directed = FALSE)
  
  # compute componentrs
  comp_sizes <- components(g)$csize
  
  # Compute total network size for comparison
  total_nodes <- vcount(g)
  
  # EXTRACT AND FORMAT TOP 100
  sorted_sizes <- sort(comp_sizes, decreasing = TRUE)
  
  # Change NA into 0 (networks should be fully connected, kept analysis to check for rare cases)
  top_100 <- numeric(100) # Creates exactly 100 zeros safely
  valid_comps <- min(length(sorted_sizes), 100)
  top_100[1:valid_comps] <- sorted_sizes[1:valid_comps]
  
  df_row <- as.data.frame(t(top_100))
  colnames(df_row) <- paste0("comp_", 1:100)
  
  results_list[[yyyy_mm]] <- bind_cols(month = yyyy_mm, total_nodes = total_nodes, df_row)
  
  
  cat(sprintf(" Processed %s: %d total nodes, LCC has %d nodes \n", yyyy_mm, total_nodes, top_100[1]))
  
  # FIX: Define lcc_size before using it in the check below
  lcc_size <- top_100[1] 
  
  if (lcc_size < total_nodes) {
    disconnected_nodes <- total_nodes - lcc_size
    total_components <- length(comp_sizes)
    
    fragmented_list[[yyyy_mm]] <- data.frame(
      month = yyyy_mm,
      total_nodes = total_nodes,
      lcc_size = lcc_size,
      disconnected_nodes = disconnected_nodes,
      total_components = total_components,
      second_largest_comp = ifelse(length(sorted_sizes) >= 2, sorted_sizes[2], 0)
    )
  }
}

#Save data
final_top_100_df <- bind_rows(results_list)
write_csv(final_top_100_df, output_csv)

# Check which networks are not fully connected
if (length(fragmented_list) > 0) {
  fragmented_df <- bind_rows(fragmented_list)
  fragmented_df <- fragmented_df %>% arrange(month)
  write_csv(fragmented_df, fragmented_csv)
  cat(sprintf(" Found %d networks that were not perfectly connected.\n", nrow(fragmented_df)))
  cat(sprintf("Diagnostic data for fragmented networks saved to: '%s'\n", fragmented_csv))
} else {
  cat("All networks perfectly connected (LCC == Total Nodes). No fragmented log generated.\n")
}

# Calculate the ratio of the Largest Connected Component to total nodes
final_top_100_df <- final_top_100_df %>%
  mutate(lcc_ratio = comp_1 / total_nodes)

# Extract first and last month for clean X-axis anchoring
first_month <- min(final_top_100_df$month)
last_month <- max(final_top_100_df$month)

# Generate the plot
lcc_plot <- ggplot(final_top_100_df, aes(x = month, y = lcc_ratio, group = 1)) +
  geom_line(color = "#984EA3", linewidth = 1) +
  geom_point(color = "#984EA3", size = 2, alpha = 0.8) +
  scale_x_discrete(breaks = c(first_month, last_month)) + 
  
  scale_y_continuous(limits = c(0.98, 1), labels = scales::percent_format(accuracy = 0.01)) + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Largest Connected Component (LCC) relative to Total Network Size",
    x = "Timeframe",
    y = "Percentage of Nodes in LCC"
  )

# Save the plot
plot_output <- "lcc_ratio_timeline.png"
ggsave(plot_output, plot = lcc_plot, width = 10, height = 5, dpi = 300, bg = "white")