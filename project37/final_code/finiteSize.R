library(igraph)
library(ggplot2)
library(dplyr)

# Simulate networks
get_delta_uf <- function(n, rule, max_r = 1.5) {
  total_edges <- floor(n * max_r)
  
  parent <- 1:n
  comp_size <- rep(1, n)
  max_c <- 1
  
  t0 <- NA
  t1 <- NA
  threshold_0 <- sqrt(n)
  threshold_1 <- 0.5 * n
  
  u1 <- sample(1:n, total_edges, replace = TRUE)
  v1 <- sample(1:n, total_edges, replace = TRUE)
  
  if (rule %in% c("BF", "PR", "SR")) {
    u2 <- sample(1:n, total_edges, replace = TRUE)
    v2 <- sample(1:n, total_edges, replace = TRUE)
  }
  
  for (step in 1:total_edges) {
    
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
    
    if (is.na(t0) && max_c >= threshold_0) t0 <- step - 1
    if (is.na(t1) && max_c > threshold_1) {
      t1 <- step
      break # stop immediately once Delta is found
    }
  }
  
  return(t1 - t0)
}

# Function to do the analysis with different n ranges
run_finite_size_analysis <- function(n_ranges, rules = c("ER", "BF", "PR"), num_sims = 5) {
  
  final_results <- data.frame()
  
  for (range_name in names(n_ranges)) {
    n_values <- n_ranges[[range_name]]
    mean_n <- mean(n_values) # Used for X-axis placement
    
    cat(sprintf("\n=== Testing Range: %s (Mean N = %.0f) ===\n", range_name, mean_n))
    
    # Collect all deltas for this specific range
    range_data <- data.frame()
    for (n in n_values) {
      for (rule in rules) {
        deltas <- sapply(1:num_sims, function(i) get_delta_uf(n, rule))
        range_data <- rbind(range_data, data.frame(
          n = n, 
          Rule = rule, 
          log_n = log10(n), 
          log_Delta = log10(mean(deltas))
        ))
      }
    }
    
    # Fit the log-log regression independently for THIS RANGE
    fits <- range_data %>%
      group_by(Rule) %>%
      do({
        fit <- lm(log_Delta ~ log_n, data = .)
        ci <- confint(fit, "log_n", level = 0.95)
        data.frame(
          Slope = coef(fit)["log_n"],
          CI_Lower = ci[1],
          CI_Upper = ci[2]
        )
      }) %>%
      ungroup()
    
    # Store results
    fits$Range <- range_name
    fits$Mean_N <- mean_n
    final_results <- rbind(final_results, fits)
    
    print(fits %>% select(Rule, Slope, CI_Lower, CI_Upper))
  }

  
  # PLot convergence
  p <- ggplot(final_results, aes(x = Mean_N, y = Slope, color = Rule, group = Rule)) +
    # Theoretical asymptotic lines
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.5) +
    geom_hline(yintercept = 0.66, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.5) +
    
    # Data points and error bars representing the 95% CI of the slope in that range
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, linewidth = 0.8) +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    
    # Log scale for the X-axis 
    scale_x_log10(breaks = unique(final_results$Mean_N), labels = unique(final_results$Range)) +
    scale_color_manual(values = c("ER" = "black", "BF" = "blue", "PR" = "red", "SR" = "darkorange")) +
    
    labs(
      title = "Convergence of the Scaling Exponent (Finite-Size Effects)",
      x = "Network Size Range (Log Scale)",
      y = expression("Estimated Slope " * beta),
      color = "Model"
    ) +
    annotate("text", x = min(final_results$Mean_N), y = 1.02, label = "Theoretical Limit (ER/BF) = 1.0", hjust = 0) +
    annotate("text", x = min(final_results$Mean_N), y = 0.68, label = "Theoretical Limit (PR) ≈ 0.66", hjust = 0, color = "red") +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
    )
  
  print(p)
  ggsave("finite_size_effects_convergence.png", plot = p, width = 10, height = 6, dpi = 300)
  
  return(final_results)
}

set.seed(42)

# Define progressively larger ranges of Network Sizes
n_ranges_to_test <- list(
  "100 - 500" = c(100, 200, 300, 400, 500),
  "1k - 5k" = c(1000, 2000, 3000, 4000, 5000),
  "10k - 50k" = c(10000, 20000, 30000, 40000, 50000),
  "100k - 500k" = c(100000, 200000, 300000, 400000, 500000),
  "500k - 3M" = c(500000, 750000, 1000000, 1250000,1500000,2000000,2500000,3000000)
)

# Run the analysis 
results <- run_finite_size_analysis(n_ranges_to_test, rules = c("ER", "BF", "PR"), num_sims = 50)