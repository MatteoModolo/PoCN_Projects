library(igraph)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mclust)    

louvain_dir <- "louvain_communities"
leiden_dir  <- "leiden_communities"

sbm_dir <- "sbm_mdl_blocks"
sbm_prefix <- "sbm_nodes_"
sbm_col <- "sbm_block_id"

# Find all months we have data for across Leiden
leiden_files <- list.files(leiden_dir, pattern = "^leiden_nodes_.*\\.csv$")
months <- stringr::str_extract(leiden_files, "\\d{4}_\\d{2}")

results <- list()

for (yyyy_mm in months) {
  
  # 1. File Paths
  f_louvain <- file.path(louvain_dir, paste0("louvain_nodes_", yyyy_mm, ".csv"))
  f_leiden  <- file.path(leiden_dir, paste0("leiden_nodes_", yyyy_mm, ".csv"))
  f_sbm     <- file.path(sbm_dir, paste0(sbm_prefix, yyyy_mm, ".csv"))
  
  # Skip if any model failed to generate a file for this month
  if (!file.exists(f_louvain) || !file.exists(f_leiden) || !file.exists(f_sbm)) next
  
  #  Load Data
  louvain <- read_csv(f_louvain, show_col_types = FALSE)
  leiden  <- read_csv(f_leiden, show_col_types = FALSE)
  sbm     <- read_csv(f_sbm, show_col_types = FALSE)
  
  # Merge data to ensure 1-to-1 node alignment
  merged_df <- leiden %>%
    rename(leiden_id = community_id) %>%
    inner_join(louvain %>% rename(louvain_id = community_id), by = "asn") %>%
    inner_join(sbm %>% rename(sbm_id = all_of(sbm_col)), by = "asn")
  
  # Calculate Community Counts
  n_leiden  <- length(unique(merged_df$leiden_id))
  n_louvain <- length(unique(merged_df$louvain_id))
  n_sbm     <- length(unique(merged_df$sbm_id))
  
  # Calculate divergence with built in mclust function
  ari_louvain <- adjustedRandIndex(merged_df$leiden_id, merged_df$louvain_id)
  ari_sbm     <- adjustedRandIndex(merged_df$leiden_id, merged_df$sbm_id)
  
  div_louvain <- 1 - ari_louvain
  div_sbm     <- 1 - ari_sbm
  
  # Store
  results[[yyyy_mm]] <- data.frame(
    month = yyyy_mm,
    Count_Leiden = n_leiden,
    Count_Louvain = n_louvain,
    Count_SBM = n_sbm,
    Divergence_Louvain = div_louvain,
    Divergence_SBM = div_sbm
  )
}

final_df <- bind_rows(results)

unique_months <- sort(unique(final_df$month))
start_month <- unique_months[1]
end_month <- unique_months[length(unique_months)]

# Pivot Community Counts
df_counts <- final_df %>%
  select(month, Leiden = Count_Leiden, Louvain = Count_Louvain, SBM = Count_SBM) %>%
  pivot_longer(cols = -month, names_to = "Model", values_to = "Community_Count")

# Pivot Divergence (Comparing against Leiden)
df_divergence <- final_df %>%
  select(month, Louvain = Divergence_Louvain, SBM = Divergence_SBM) %>%
  pivot_longer(cols = -month, names_to = "Model", values_to = "Divergence")

# Ensure consistent colors
model_colors <- c("Leiden" = "#E41A1C", "Louvain" = "#377EB8", "SBM" = "#4DAF4A")



# Plot of number of communities
p1 <- ggplot(df_counts, aes(x = month, y = Community_Count, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_x_discrete(breaks = c(start_month, end_month)) +
  scale_color_manual(values = model_colors) +
  
  coord_cartesian(ylim = c(0, 200)) + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Community Detection Models",
    subtitle = "Total Number of Communities / Blocks Detected ",
    y = "Number of Communities"
  )

# plot divergence from leiden
p2 <- ggplot(df_divergence, aes(x = month, y = Divergence, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_x_discrete(breaks = c(start_month, end_month)) +
  scale_color_manual(values = model_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    subtitle = "Structural Divergence from the Leiden Baseline",
    x = "Timeframe",
    y = "Divergence Score"
  )

# Combine the plots dynamically using patchwork
final_plot <- p1 / p2 + plot_layout(heights = c(1, 1))

print(final_plot)
ggsave("structural_divergence_analysis.png", plot = final_plot, width = 10, height = 7, dpi = 300, bg = "white")