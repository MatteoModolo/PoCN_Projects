library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# union find
get_trajectory_uf <- function(n, rule, max_r = 1.5) {
  total_edges <- floor(n * max_r)
  
  parent <- 1:n
  comp_size <- rep(1, n)
  max_c <- 1
  
  # Vectorized edge pool
  u1 <- sample(1:n, total_edges, replace = TRUE)
  v1 <- sample(1:n, total_edges, replace = TRUE)
  
  if (rule %in% c("BF", "PR", "SR")) {
    u2 <- sample(1:n, total_edges, replace = TRUE)
    v2 <- sample(1:n, total_edges, replace = TRUE)
  }
  
  # Array to store the size of the giant component at every single step
  history <- numeric(total_edges)
  
  for (step in 1:total_edges) {
    
    # Path compression
    r_u1 <- u1[step]
    while(parent[r_u1] != r_u1) { parent[r_u1] <- parent[parent[r_u1]]; r_u1 <- parent[r_u1] }
    r_v1 <- v1[step]
    while(parent[r_v1] != r_v1) { parent[r_v1] <- parent[parent[r_v1]]; r_v1 <- parent[r_v1] }
    
    if (rule == "ER") {
      ru <- r_u1
      rv <- r_v1
    } else {
      r_u2 <- u2[step]
      while(parent[r_u2] != r_u2) { parent[r_u2] <- parent[parent[r_u2]]; r_u2 <- parent[r_u2] }
      r_v2 <- v2[step]
      while(parent[r_v2] != r_v2) { parent[r_v2] <- parent[parent[r_v2]]; r_v2 <- parent[r_v2] }
      
      s_u1 <- comp_size[r_u1]
      s_v1 <- comp_size[r_v1]
      s_u2 <- comp_size[r_u2]
      s_v2 <- comp_size[r_v2]
      
      # Model Rules
      choose_first <- TRUE
      if (rule == "BF") {
        choose_first <- (s_u1 == 1 && s_v1 == 1)
      } else if (rule == "PR") {
        choose_first <- ((s_u1 * s_v1) <= (s_u2 * s_v2))
      } else if (rule == "SR") {
        choose_first <- ((s_u1 + s_v1) <= (s_u2 + s_v2))
      }
      
      if (choose_first) {
        ru <- r_u1
        rv <- r_v1
      } else {
        ru <- r_u2
        rv <- r_v2
      }
    }
    
    # Union by size
    if (ru != rv) {
      if (comp_size[ru] < comp_size[rv]) {
        parent[ru] <- rv
        comp_size[rv] <- comp_size[rv] + comp_size[ru]
        if (comp_size[rv] > max_c) max_c <- comp_size[rv]
      } else {
        parent[rv] <- ru
        comp_size[ru] <- comp_size[ru] + comp_size[rv]
        if (comp_size[ru] > max_c) max_c <- comp_size[ru]
      }
    }
    
    # Record the maximum component size for this step
    history[step] <- max_c
  }
  
  return(history)
}


run_asymptotic_experiment <- function(n_values, rules = c("ER", "BF", "PR"), num_sims = 10) {
  
  results <- data.frame()
  
  cat("=== Extracting Phase Transition Boundaries ===\n")
  for (n in n_values) {
    cat(sprintf("Simulating for Network Size (n) = %d\n", n))
    
    threshold_0 <- sqrt(n)
    threshold_1 <- 0.5 * n
    
    for (rule in rules) {
      
      # Run simulation and save in matrix
      history_matrix <- do.call(cbind, lapply(1:num_sims, function(i) get_trajectory_uf(n, rule)))
      
      # Averaging values into Master Curve
      master_avg_history <- rowMeans(history_matrix)
      
      # ONE Interpolation per model for t0 and t1
      t0 <- which(master_avg_history >= threshold_0)[1] - 1
      t1 <- which(master_avg_history > threshold_1)[1]
      
      if (is.na(t0) || t0 < 1) t0 <- 1
      
      results <- rbind(results, data.frame(
        n = n,
        Rule = rule,
        t0_density = t0 / n,
        t1_density = t1 / n
      ))
    }
  }
  
  
  
  # Plot
  
  # Reshape data for plotting
  plot_data <- results %>%
    pivot_longer(cols = c(t0_density, t1_density), 
                 names_to = "Metric", 
                 values_to = "Density") %>%
    mutate(Metric = ifelse(Metric == "t0_density", "t0/n (Start)", "t1/n (End)"))
  
  # Base colors
  base_colors <- c("ER" = "black", "BF" = "blue", "PR" = "red", "SR" = "darkorange")
  
  p <- ggplot(plot_data, aes(x = n, y = Density, color = Rule, shape = Metric)) +
    
    # Plot the simulated points
    geom_point(size = 3.5, alpha = 0.9) +
    geom_line(aes(group = interaction(Rule, Metric)), alpha = 0.4, linewidth = 1) +
    
    
    scale_x_log10(breaks = n_values, labels = scales::comma) +
    scale_color_manual(values = base_colors) +
    scale_linetype_manual(values = c("t0/n (Start)" = "dotted", "t1/n (End)" = "dashed")) +
    
    labs(
      title = "Asymptotic Convergence of t0/n and t1/n",
      x = "Network Size (n) - Log Scale",
      y = "Edge Density (r = t/n)",
      color = "Model",
      shape = "Boundary",
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  print(p)
  ggsave("Convergence.png", plot = p, width = 10, height = 6, dpi = 300)
  
  return(list(raw_data = results))
}

# --- EXECUTION ---
set.seed(42)


n_test_values <- c(5000, 10000, 20000, 50000, 100000, 200000, 500000,1000000,2000000)

# Run multiple simulations per point
final_data <- run_asymptotic_experiment(n_test_values, rules = c("ER", "BF", "PR"), num_sims =50)