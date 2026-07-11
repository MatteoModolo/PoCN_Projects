library(readr)
library(dplyr)
library(ggplot2)
library(scales) # Needed for clean percentage formatting on the plot


input_file <- "top_30_nodes_monthly.csv"
plot_output <- "rank_1_persistence_bar_chart.png"

# Load the data
tryCatch({
  top_nodes <- read_csv(input_file, show_col_types = FALSE)
}, error = function(e) {
  stop(sprintf("Error: Could not find '%s'. Make sure you are in the right directory.", input_file))
})

# Get the absolute total number of months in the dataset
total_months <- n_distinct(top_nodes$month)

cat(sprintf("Analyzing dataset containing %d total months...\n", total_months))


# Nodes that were rank #1 in at least one month
rank_1_nodes <- top_nodes %>%
  filter(rank == 1) %>%
  pull(node) %>%
  unique()

cat(sprintf("Found %d unique node(s) that achieved Rank 1 over the decade.\n", length(rank_1_nodes)))

# For those specific nodes, count how many total months they spent in the Top 10
persistence_df <- top_nodes %>%
  filter(node %in% rank_1_nodes & rank <= 10) %>%
  group_by(node) %>%
  summarise(months_in_top_10 = n(), .groups = 'drop') %>%

  mutate(
    fraction = months_in_top_10 / total_months,
    node_label = as.character(node)
  ) %>%
  # Sort from highest fraction to lowest
  arrange(desc(fraction))

# bar plot

p <- ggplot(persistence_df, aes(x = reorder(node_label, -fraction), y = fraction, fill = node_label)) +
  geom_col(alpha = 0.85, width = 0.6) +
  
  # Add the exact percentage number hovering just above each bar
  geom_text(aes(label = percent(fraction, accuracy = 0.1)), 
            vjust = -0.6, fontface = "bold", size = 4.5) +
  
  # Format Y axis as percentages, expand slightly so the text doesn't get cut off
  scale_y_continuous(labels = percent_format(), limits = c(0, max(persistence_df$fraction) * 1.1)) +
  
  # A pleasant, colorblind-friendly palette
  scale_fill_brewer(palette = "Set1") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "none", # Hide legend because the X-axis already has the node IDs
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  labs(
    title = "Top 10 Persistence of Rank #1 Nodes",
    subtitle = sprintf("Percentage of the decade the highest-degree nodes spent inside the Top 10 (Total months: %d)", total_months),
    x = "Autonomous System Number (Rank #1 Nodes)",
    y = "Fraction of Time Spent in Top 10"
  )

# Save the plot
ggsave(plot_output, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
cat(sprintf("Done! Bar plot saved to '%s'\n", plot_output))