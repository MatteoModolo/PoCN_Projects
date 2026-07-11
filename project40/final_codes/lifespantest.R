library(readr)
library(dplyr)
library(stringr)
library(ggplot2)


input_dir <- "final_weighted_networks"
plot_output <- "top_100_lifespan_gantt_chart_clean.png"

# Get files and sort chronologically
network_files <- list.files(path = input_dir, pattern = "^as_edges_.*\\.tsv\\.gz$", full.names = TRUE)
network_files <- sort(network_files)

if (length(network_files) == 0) stop("No network files found in the target directory.")

results_list <- list()

for (file_path in network_files) {
  
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  # Load the raw edges
  edges_df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
  
  # Fast degree calculation
  all_nodes <- c(edges_df$asn1, edges_df$asn2)
  degree_table <- table(all_nodes)
  
  # Convert to dataframe, sort, and keep ONLY the Top 100 for this month
  month_top_100 <- data.frame(
    node = as.character(names(degree_table)),
    degree = as.integer(degree_table),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(degree)) %>%
    slice_head(n = 100) %>%
    mutate(month = yyyy_mm)
  
  results_list[[yyyy_mm]] <- month_top_100
  
  cat(sprintf(" Extracted Top 100 for %s\n", yyyy_mm))
}

# Combine all months into one master Top 100 occurrence ledger
top_100_occurrences <- bind_rows(results_list) %>% arrange(month)

unique_months <- sort(unique(top_100_occurrences$month))
first_month <- unique_months[1]
last_month <- unique_months[length(unique_months)]


# Calculate lifespans inside the Top 100 VIP Club
node_lifespans <- top_100_occurrences %>%
  group_by(node) %>%
  summarize(
    first_seen = min(month),
    last_seen = max(month),
    total_months = n()
  ) %>%
  mutate(
    status = case_when(
      first_seen == first_month & last_seen == last_month ~ "Founding Elite (Persistent)",
      first_seen > first_month & last_seen == last_month ~ "Rising Star (Emerged)",
      first_seen == first_month & last_seen < last_month ~ "Fallen Giant (Disappeared)",
      first_seen > first_month & last_seen < last_month ~ "Transient (Flashed & Faded)",
      TRUE ~ "Unknown"
    )
  )


# Convert the "YYYY_MM" strings to date objects
node_lifespans <- node_lifespans %>%
  mutate(
    start_date = as.Date(paste0(first_seen, "_01"), format="%Y_%m_%d"),
    end_date = as.Date(paste0(last_seen, "_01"), format="%Y_%m_%d")
  )

# Sort the nodes so the chart looks like a clean waterfall
node_lifespans <- node_lifespans %>%
  arrange(status, start_date, desc(total_months)) %>%
  mutate(node_fct = factor(as.character(node), levels = unique(as.character(node))))



# Define custom, distinct colors for the 4 states
status_colors <- c(
  "Founding Elite (Persistent)" = "#4DAF4A",  # Green
  "Rising Star (Emerged)"       = "#377EB8",  # Blue
  "Fallen Giant (Disappeared)"  = "#E41A1C",  # Red
  "Transient (Flashed & Faded)" = "#FF7F00"   # Orange
)

p <- ggplot(node_lifespans, aes(x = start_date, xend = end_date, y = node_fct, yend = node_fct, color = status)) +
  # Draw the lifespan blocks (slightly thinner than before to accommodate more nodes)
  geom_segment(linewidth = 2.0, alpha = 0.9) +
  
  scale_color_manual(values = status_colors) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  
  theme_minimal() +
  theme(
    # Completely remove the Y-axis text and ticks for a clean look
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Timeframe",
    y = "Unique Autonomous Systems (Sorted by Status)"
  )

# Save the plot
ggsave(plot_output, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
