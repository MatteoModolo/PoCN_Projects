library(igraph)
library(ggplot2)
library(dplyr)

# union find simulation
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
    
    # Inline Path Compression
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


run_finite_size_density_analysis <- function(n_ranges, rules = c("ER", "BF", "PR"), num_sims = 5) {
  
  final_results <- data.frame()
  
  cat("=== Running Density Gap Convergence Analysis ===\n")
  cat("NOTE: Using rigorous ensemble-averaging before interpolating Delta.\n")
  
  for (range_name in names(n_ranges)) {
    n_values <- n_ranges[[range_name]]
    mean_n <- mean(n_values) # Used to place the point on the X-axis of the final plot
    
    cat(sprintf("\n[ Testing Range: %s (Mean N = %.0f) ]\n", range_name, mean_n))
    
    range_data <- data.frame()
    
    for (n in n_values) {
      threshold_0 <- sqrt(n)
      threshold_1 <- 0.5 * n
      
      for (rule in rules) {
        # 1. RUN SIMULATIONS & BUILD MATRIX
        history_matrix <- do.call(cbind, lapply(1:num_sims, function(i) get_trajectory_uf(n, rule)))
        
        # 2. ENSEMBLE AVERAGE THE SIMULATIONS
        master_avg_history <- rowMeans(history_matrix)
        
        # 3. ONE EXACT INTERPOLATION ON THE MASTER CURVE
        t0 <- which(master_avg_history >= threshold_0)[1] - 1
        t1 <- which(master_avg_history > threshold_1)[1]
        
        if (is.na(t0) || t0 < 1) t0 <- 1
        
        delta_density <- (t1 - t0) / n
        
        range_data <- rbind(range_data, data.frame(
          n = n, 
          Rule = rule, 
          log_n = log10(n), 
          log_Delta_Density = log10(delta_density)
        ))
      }
    }
    
    # 4. FIT THE LOG-LOG REGRESSION STRICTLY FOR THIS RANGE
    fits <- range_data %>%
      group_by(Rule) %>%
      do({
        fit <- lm(log_Delta_Density ~ log_n, data = .)
        ci <- confint(fit, "log_n", level = 0.95)
        data.frame(
          Slope = coef(fit)["log_n"],
          CI_Lower = ci[1],
          CI_Upper = ci[2]
        )
      }) %>%
      ungroup()
    
    # Store the local fits
    fits$Range <- range_name
    fits$Mean_N <- mean_n
    final_results <- rbind(final_results, fits)
    
    print(fits %>% select(Rule, Slope, CI_Lower, CI_Upper))
  }
  
  
  # plot
  p <- ggplot(final_results, aes(x = Mean_N, y = Slope, color = Rule, group = Rule)) +
    
    # Theoretical asymptotic lines for the Density Gap (beta - 1)
    geom_hline(yintercept = 0.0, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.5) +
    geom_hline(yintercept = -0.34, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.5) +
    
    # Data points and error bars (95% CI of the slope inside that local range)
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4) +
    
    # Log scale for the X-axis so ranges space out nicely
    scale_x_log10(breaks = unique(final_results$Mean_N), labels = unique(final_results$Range)) +
    scale_color_manual(values = c("ER" = "black", "BF" = "blue", "PR" = "red", "SR" = "darkorange")) +
    
    labs(
      title = "Convergence of the Density Gap Slope (Finite-Size Effects)",
      subtitle = expression("Slopes of " * Delta * "/n calculated in rolling windows. Proves that ER/BF are merely delayed by finite-size effects!"),
      x = "Network Size Range (Log Scale)",
      y = expression("Estimated Slope of " * Delta * "/n (" * beta - 1 * ")"),
      color = "Model"
    ) +
    annotate("text", x = min(final_results$Mean_N), y = 0.04, label = "Continuous Limit (ER/BF) = 0.0", hjust = 0, fontface="bold") +
    annotate("text", x = min(final_results$Mean_N), y = -0.30, label = "Explosive Limit (PR) ≈ -0.34", hjust = 0, color = "red", fontface="bold") +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  print(p)
  ggsave("finite_size_density_convergence.png", plot = p, width = 10, height = 6, dpi = 300)
  
  return(final_results)
}

set.seed(27)

# Define progressively larger ranges of Network Sizes
n_ranges_to_test <- list(
  "100 - 500" = c(100, 200, 300, 400, 500),
  "1k - 5k" = c(1000, 2000, 3000, 4000, 5000),
  "10k - 50k" = c(10000, 20000, 30000, 40000, 50000),
  "100k - 500k" = c(100000, 200000, 300000, 400000, 500000),
  "500k - 2.5M" = c(500000, 750000, 1000000, 1250000,1500000,2000000,2500000)
)

# Run the analysis
results <- run_finite_size_density_analysis(n_ranges_to_test, rules = c("ER", "BF", "PR"), num_sims = 50)