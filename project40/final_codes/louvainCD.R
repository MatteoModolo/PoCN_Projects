library(igraph)
library(readr)
library(dplyr)
library(stringr)

# directories and paths
input_dir <- "final_weighted_networks"
output_summary_csv <- "louvain_modularity_summary.csv"
comm_output_dir <- "louvain_communities"

# picked a random month to save structure and compare graphically the found communities
target_plot_month <- "2016_01" 

# Create a dedicated folder for the community assignment files and plots
if (!dir.exists(comm_output_dir)) dir.create(comm_output_dir)

network_files <- list.files(path = input_dir, pattern = "^as_edges_.*\\.tsv\\.gz$", full.names = TRUE)

if (length(network_files) == 0) stop("No network files found in the input directory ")

results_list <- list()
cat(sprintf("Starting Louvain community detection for %d months \n", length(network_files)))

# Community detection loop
for (file_path in network_files) {
  
  filename <- basename(file_path)
  yyyy_mm <- str_extract(filename, "\\d{4}_\\d{2}")
  
  # Read the network edge list
  df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
  
  # Build network from data
  g <- graph_from_data_frame(df, directed = FALSE)
  
  # Run built in louvain alg, the algorithm uses the weights in the data if present
  louvain_comm <- cluster_louvain(g)
  
  # Extract the global metrics
  mod_score <- modularity(louvain_comm)
  num_communities <- length(louvain_comm)
  num_nodes <- vcount(g)
  num_edges <- ecount(g)
  
  cat(sprintf("[%s] Nodes: %5d | Communities: %4d | Modularity (Q): %.4f\n", 
              yyyy_mm, num_nodes, num_communities, mod_score))
  
  # Extract assignment of nodes to the communities
  # V(g)$name gets the actual ASN. membership(louvain_comm) gets their community ID 
  comm_df <- data.frame(
    asn = V(g)$name,
    community_id = as.integer(membership(louvain_comm))
  )
  
  # Save this month's community footprint
  comm_filename <- file.path(comm_output_dir, paste0("louvain_nodes_", yyyy_mm, ".csv"))
  write_csv(comm_df, comm_filename)
  
  # saving layout and plot 
  
  if (yyyy_mm == target_plot_month) {
    cat(sprintf(" Generating layout for %s \n", yyyy_mm))
    
    # Use DrL layout  since its highly optimized for massive graphs and community structures
    set.seed(42) # Ensure the layout is reproducible
    graph_layout <- layout_with_drl(g)
    
    # Save the graph structure and layout to disk so Leiden/other methods can use the same coordinates
    saveRDS(g, file.path(comm_output_dir, paste0("graph_structure_", yyyy_mm, ".rds")))
    saveRDS(graph_layout, file.path(comm_output_dir, paste0("layout_", yyyy_mm, ".rds")))
    cat("Layout data saved \n")
    
    # Plot the graph
    plot_filename <- file.path(comm_output_dir, paste0("louvain_plot_", yyyy_mm, ".png"))
    png(plot_filename, width = 3000, height = 3000, res = 300)
    
    # Make hubs slightly larger, leaves smaller
    v_sizes <- log10(degree(g) + 1) * 1.2
    
    plot(louvain_comm, g, 
         layout = graph_layout,
         vertex.size = v_sizes,
         vertex.label = NA,           # Disable text labels 
         edge.width = 0.05,           # Ultra-thin edges to prevent mess
         edge.color = rgb(0.5, 0.5, 0.5, alpha = 0.2),
         mark.groups = NULL,          
         main = paste("Louvain Community Detection - BGP Network", yyyy_mm)
    )
    
    dev.off()
    cat(sprintf("Community plot saved to %s \n", plot_filename))
  }
  
  # Store summary metrics for the master file
  results_list[[yyyy_mm]] <- data.frame(
    month = yyyy_mm,
    total_nodes = num_nodes,
    total_edges = num_edges,
    louvain_communities = num_communities,
    louvain_modularity = mod_score
  )
}


final_summary <- bind_rows(results_list)
write_csv(final_summary, output_summary_csv)

