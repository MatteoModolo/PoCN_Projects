library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

# dir and files
OUTPUT_DIR <- "final_weighted_networks"
METRICS_FILE <- "network_metrics_summary.csv"
PLOT_FILE <- "network_growth_timeline.png"

cat("Starting Metrics Extraction...\n")

# Find all network files in the directory
files <- list.files(OUTPUT_DIR, pattern = "^as_edges_.*_weighted\\.tsv\\.gz$", full.names = TRUE)

if (length(files) == 0) {
  stop("CRITICAL ERROR: No network files found in the target directory.")
}

metrics <- data.frame(month = character(), total_nodes = numeric(), total_edges = numeric(), stringsAsFactors = FALSE)

# Read each file and calculate size
for (file_path in files) {
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  edges_df <- tryCatch({
    read_tsv(file_path, col_select = c(1, 2), show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(edges_df) && nrow(edges_df) > 0) {
    num_edges <- nrow(edges_df)
    num_nodes <- length(unique(c(edges_df[[1]], edges_df[[2]])))
    
    cat(sprintf("Analyzed %s | Nodes: %d | Edges: %d\n", yyyy_mm, num_nodes, num_edges))
    metrics <- rbind(metrics, data.frame(month = yyyy_mm, total_nodes = num_nodes, total_edges = num_edges))
  } else {
    cat(sprintf("FAILED to read: %s (File might be corrupted)\n", yyyy_mm))
  }
}


# Mathematically sort chronologically
metrics <- metrics[order(metrics$month), ]

# Save the final statistical summary CSV
write_csv(metrics, METRICS_FILE)

# plot

first_month <- min(metrics$month)
last_month <- max(metrics$month)

# Plot nodes
p1 <- ggplot(metrics, aes(x = month, y = total_nodes, group = 1)) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_point(color = "#377EB8", size = 1.5, alpha = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Removes dates
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Macro-Structural Evolution of the Internet",
    subtitle = "Total Autonomous Systems (Nodes)",
    y = "Node Count"
  )

# Plot edges
p2 <- ggplot(metrics, aes(x = month, y = total_edges, group = 1)) +
  geom_line(color = "#E41A1C", linewidth = 1) +
  geom_point(color = "#E41A1C", size = 1.5, alpha = 0.8) +
  scale_x_discrete(breaks = c(first_month, last_month)) + # ONLY shows start and end dates
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
    subtitle = "Total BGP Routing Agreements (Edges)",
    x = "Timeframe",
    y = "Edge Count"
  )

# Combine the plots dynamically using patchwork
final_plot <- p1 / p2

# Save the plot
ggsave(PLOT_FILE, plot = final_plot, width = 10, height = 7, dpi = 300, bg = "white")