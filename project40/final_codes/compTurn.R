library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# files and dir
input_dir <- "final_weighted_networks"
output_csv <- "network_turnover_analysis.csv"

# Get files and sort them chronologically
network_files <- list.files(path = input_dir, pattern = "^as_edges_.*\\.tsv\\.gz$", full.names = TRUE)
network_files <- sort(network_files)

if (length(network_files) < 2) stop("Need at least 2 months of data to calculate turnover ")

results_list <- list()
cat(sprintf("Starting Churn Analysis for %d consecutive periods \n", length(network_files) - 1))


# Read the very first month to start the loop
file1 <- network_files[1]
month1_name <- str_extract(basename(file1), "\\d{4}_\\d{2}")

df1 <- read_tsv(file1, show_col_types = FALSE, progress = FALSE)
nodes1 <- unique(c(df1$asn1, df1$asn2))
# Create a fast, unique string ID for every edge 
edges1 <- paste0(df1$asn1, "_", df1$asn2) 

# comparison loop, the second month is kept for the next itearation
for (i in 2:length(network_files)) {
  
  file2 <- network_files[i]
  month2_name <- str_extract(basename(file2), "\\d{4}_\\d{2}")
  
  cat(sprintf("Comparing %s to %s \n", month1_name, month2_name))
  
  # Load Month 2
  df2 <- read_tsv(file2, show_col_types = FALSE, progress = FALSE)
  nodes2 <- unique(c(df2$asn1, df2$asn2))
  edges2 <- paste0(df2$asn1, "_", df2$asn2)
  
  # checking persistent, appearing and disappearsing nodes
  persistent_nodes   <- length(intersect(nodes1, nodes2))
  appearing_nodes    <- length(setdiff(nodes2, nodes1))
  disappearing_nodes <- length(setdiff(nodes1, nodes2))
  union_nodes        <- length(union(nodes1, nodes2))
  
  # checking persistent, appearing and disappearsing edges
  persistent_edges   <- length(intersect(edges1, edges2))
  appearing_edges    <- length(setdiff(edges2, edges1))
  disappearing_edges <- length(setdiff(edges1, edges2))
  union_edges        <- length(union(edges1, edges2))
  
  # save results
  period_name <- paste0(month1_name, "_to_", month2_name)
  
  results_list[[period_name]] <- data.frame(
    period = period_name,
    month1 = month1_name,
    month2 = month2_name,
    
    # Node Stats
    nodes_persistent = persistent_nodes,
    nodes_appearing = appearing_nodes,
    nodes_disappearing = disappearing_nodes,
    nodes_jaccard = persistent_nodes / union_nodes,
    
    # Edge Stats
    edges_persistent = persistent_edges,
    edges_appearing = appearing_edges,
    edges_disappearing = disappearing_edges,
    edges_jaccard = persistent_edges / union_edges
  )
  
  # carry over 2nd month as first one for next iteration
  month1_name <- month2_name
  nodes1 <- nodes2
  edges1 <- edges2
}

# save datas
final_df <- bind_rows(results_list)
write_csv(final_df, output_csv)

# (reading data from file to not repeat the analysis)
churn_data <- read_csv("network_turnover_analysis.csv", show_col_types = FALSE)

# Reshape data for plotting edges
plot_data <- churn_data %>%
  select(month2, edges_persistent, edges_appearing, edges_disappearing) %>%
  # Make disappearing negative
  mutate(edges_disappearing = -edges_disappearing) %>%
  pivot_longer(cols = c(edges_persistent, edges_appearing, edges_disappearing), 
               names_to = "Status", values_to = "Count")

# Extract first and last month for the graph to be tidier
unique_months <- unique(plot_data$month2)
start_month <- unique_months[1]
end_month <- unique_months[length(unique_months)]

# Plot
c <- ggplot(plot_data, aes(x = month2, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("edges_appearing" = "blue", 
                               "edges_persistent" = "grey70", 
                               "edges_disappearing" = "red")) +
  scale_x_discrete(breaks = c(start_month, end_month)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11)) +
  labs(title = "Month-over-Month BGP Edge Turnover",
       x = "Timeframe",
       y = "Number of Edges",
       fill = "Edge Status")

print(c)