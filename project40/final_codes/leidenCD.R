library(igraph)
library(readr)
library(dplyr)
library(stringr)

# Directories and files
input_dir <- "final_weighted_networks"
output_summary_csv <- "leiden_modularity_summary.csv"
comm_output_dir <- "leiden_communities"

# point to the Louvain directory strictly to have its layout coordinates to have same rapresentation
louvain_layout_dir <- "louvain_communities"

target_plot_month <- "2016_01" 

# Create a dedicated folder for the community assignment files and plots
if (!dir.exists(comm_output_dir)) dir.create(comm_output_dir)

network_files <- list.files(path = input_dir, pattern = "^as_edges_.*\\.tsv\\.gz$", full.names = TRUE)

if (length(network_files) == 0) stop("No network files found in the input directory.")

results_list <- list()
cat(sprintf("Starting Leiden community detection for %d months \n", length(network_files)))

# loop on networks
for (file_path in network_files) {
  
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  # Read the network edge list and build igraph object
  df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
  g <- graph_from_data_frame(df, directed = FALSE)
  
  # run leiden script using modularity so that is comparable with the louvain mnethod
  leiden_comm <- cluster_leiden(g, objective_function = "modularity")
  
  # extract global metrics
  mod_score <- modularity(g, membership(leiden_comm))
  num_communities <- length(leiden_comm)
  num_nodes <- vcount(g)
  num_edges <- ecount(g)
  
  cat(sprintf("[%s] Nodes: %5d | Communities: %4d | Modularity (Q): %.4f\n", 
              yyyy_mm, num_nodes, num_communities, mod_score))
  
  # Extract node assignments
  comm_df <- data.frame(
    asn = V(g)$name,
    community_id = as.integer(membership(leiden_comm))
  )
  
  # Save the same month exact community footprint
  comm_filename <- file.path(comm_output_dir, paste0("leiden_nodes_", yyyy_mm, ".csv"))
  write_csv(comm_df, comm_filename)
  
  # take saved layout and plot
  if (yyyy_mm == target_plot_month) {
    cat(sprintf(" Preparing plot for %s\n", yyyy_mm))
    
    # Check if the Louvain layout exists so we can reuse the exact coordinates
    layout_file <- file.path(louvain_layout_dir, paste0("layout_", yyyy_mm, ".rds"))
    
    if (file.exists(layout_file)) {
      graph_layout <- readRDS(layout_file)
    } else {
      set.seed(42)
      graph_layout <- layout_with_drl(g)
      saveRDS(graph_layout, file.path(comm_output_dir, paste0("layout_", yyyy_mm, ".rds")))
    }
    
    # Plot the graph using Leiden communities but  Louvain spatial coordinates
    plot_filename <- file.path(comm_output_dir, paste0("leiden_plot_", yyyy_mm, ".png"))
    png(plot_filename, width = 3000, height = 3000, res = 300)
    
    v_sizes <- log10(degree(g) + 1) * 1.2
    
    plot(leiden_comm, g, 
         layout = graph_layout,
         vertex.size = v_sizes,
         vertex.label = NA,
         edge.width = 0.05,
         edge.color = rgb(0.5, 0.5, 0.5, alpha = 0.2),
         mark.groups = NULL,
         main = paste("Leiden Community Detection - BGP Network", yyyy_mm)
    )
    
    dev.off()
    cat(sprintf("Community plot saved to %s \n", plot_filename))
  }
  
  # Store summary metrics 
  results_list[[yyyy_mm]] <- data.frame(
    month = yyyy_mm,
    total_nodes = num_nodes,
    total_edges = num_edges,
    leiden_communities = num_communities,
    leiden_modularity = mod_score
  )
}

# save data
final_summary <- bind_rows(results_list)
write_csv(final_summary, output_summary_csv)
