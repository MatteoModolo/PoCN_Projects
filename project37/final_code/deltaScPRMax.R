library(igraph)
library(ggplot2)
library(gridExtra)
library(dplyr)

simulate_percolation_uf <- function(n, rule, max_r = 1.5, record_step = 50, save_trajectory = FALSE, file_name = NULL) {
  total_edges <- floor(n * max_r)
  
  # Initialize Union-Find
  parent <- 1:n
  comp_size <- rep(1, n)
  max_c <- 1
  
  t0 <- NA
  t1 <- NA
  threshold_0 <- sqrt(n)
  threshold_1 <- 0.5 * n
  
  # Vectorized random sampling for speed
  u1 <- sample(1:n, total_edges, replace = TRUE)
  v1 <- sample(1:n, total_edges, replace = TRUE)
  
  if (rule %in% c("BF", "PR", "PRmax", "SR")) {
    u2 <- sample(1:n, total_edges, replace = TRUE)
    v2 <- sample(1:n, total_edges, replace = TRUE)
  }
  
  # Pre-allocate tracking arrays
  num_records <- floor(total_edges / record_step)
  r_values <- numeric(num_records)
  c_n_values <- numeric(num_records)
  record_idx <- 1
  
  if (save_trajectory) history_max_c <- numeric(total_edges)
  
  # Simulation loop
  for (step in 1:total_edges) {
    
    # Find with path compression (Edge 1)
    r_u1 <- u1[step]
    while(parent[r_u1] != r_u1) { parent[r_u1] <- parent[parent[r_u1]]; r_u1 <- parent[r_u1] }
    r_v1 <- v1[step]
    while(parent[r_v1] != r_v1) { parent[r_v1] <- parent[parent[r_v1]]; r_v1 <- parent[r_v1] }
    
    if (rule == "ER") {
      ru <- r_u1
      rv <- r_v1
    } else {
      # Find with path compression (Edge 2)
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
      } else if (rule == "PRmax") {
        choose_first <- ((s_u1 * s_v1) >= (s_u2 * s_v2))
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
    
    # Update thresholds for the transition window
    if (is.na(t0) && max_c >= threshold_0) t0 <- step - 1
    if (is.na(t1) && max_c > threshold_1) t1 <- step
    
    # Recording
    if (step %% record_step == 0) {
      r_values[record_idx] <- step / n
      c_n_values[record_idx] <- max_c / n
      record_idx <- record_idx + 1
    }
    
    if (save_trajectory) history_max_c[step] <- max_c
  }
  
  if (save_trajectory && !is.null(file_name)) {
    write.csv(data.frame(step = 1:total_edges, max_c = history_max_c), file_name, row.names = FALSE)
  }
  
  return(list(
    curve_data = data.frame(r = r_values, C_n = c_n_values),
    delta = t1 - t0
  ))
}

# Experiment wrapper
run_scaling_experiment <- function(n_values, subset_plot_n, num_sims = 20, min_fit_n = 50000) {
  all_delta_results <- data.frame()
  all_curve_results <- data.frame()
  
  for (n in n_values) {
    cat(sprintf("Simulating n = %d\n", n))
    
    # ADDED PRmax to the rule iteration
    for (rule in c("ER", "BF", "PR", "PRmax", "SR")) {
      cat(sprintf("  Rule: %s, Runs: %d\n", rule, num_sims))
      
      results <- lapply(1:num_sims, function(i) {
        simulate_percolation_uf(n, rule, save_trajectory = FALSE)
      })
      
      deltas <- sapply(results, function(x) x$delta)
      mean_d <- mean(deltas, na.rm = TRUE)
      
      all_delta_results <- rbind(all_delta_results, data.frame(
        n = n, 
        Rule = rule, 
        Mean_Delta = mean_d, 
        log_n = log10(n),
        log_Delta = log10(mean_d),
        Delta_div_n = mean_d / n
      ))
      
      if (n %in% subset_plot_n) {
        curves <- do.call(rbind, lapply(results, function(x) x$curve_data))
        avg_curve <- aggregate(C_n ~ r, data = curves, FUN = mean)
        avg_curve$n <- n
        avg_curve$Rule <- rule
        all_curve_results <- rbind(all_curve_results, avg_curve)
      }
    }
  }
  
  # Extract slopes and 95% intervals STRICTLY for large N (N >= min_fit_n)
  fits <- all_delta_results %>%
    filter(n >= min_fit_n) %>%
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
  
  # Print only the fitting results to avoid confusion
  cat(sprintf("\n--- Delta Analysis Results (Fitting only for N >= %d) ---\n", min_fit_n))
  print(fits %>% select(Rule, Slope, CI_Lower, CI_Upper))
  cat("----------------------------------------------------------------\n")
  
  # Merge to create dynamic labels for the plots
  all_delta_results <- merge(all_delta_results, fits, by = "Rule")
  all_delta_results$Legend_Label <- sprintf("%s (Slope = %.2f, 95%% CI [%.2f, %.2f])", 
                                            all_delta_results$Rule, 
                                            all_delta_results$Slope,
                                            all_delta_results$CI_Lower,
                                            all_delta_results$CI_Upper)
  
  # Base colors mapping including PRmax
  rule_colors <- c("black", "blue", "red", "purple", "darkorange")
  
  # Log-Log Plot for scaling
  p_loglog <- ggplot(all_delta_results, aes(x = log_n, y = log_Delta, color = Legend_Label, fill = Legend_Label)) +
    geom_point(size = 4, alpha = 0.8) +
    # Draw regression bands ONLY for points >= min_fit_n
    geom_smooth(data = filter(all_delta_results, n >= min_fit_n), method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2, linetype = "dashed") +
    scale_color_manual(values = setNames(rule_colors, unique(all_delta_results$Legend_Label))) +
    scale_fill_manual(values = setNames(rule_colors, unique(all_delta_results$Legend_Label))) +
    labs(
      title = sprintf("Log-Log demonstration of Delta scaling (Fit N >= %d)", min_fit_n),
      subtitle = expression("ER/BF/PRmax (Extensive: "*beta%approx%1.0*") vs PR/SR (Sub-Linear: "*beta%approx%0.66*")"),
      x = expression(Log[10] * "(N) - Network Size"),
      y = expression(Log[10] * "(Delta) - Transition Window"),
      color = "Model and scaling exponent",
      fill = "Model and scaling exponent"
    ) +
    theme_bw() +
    theme(legend.position = "top", legend.direction = "vertical")
  
  # Classic Delta/n Plot
  p_delta_classic <- ggplot(all_delta_results, aes(x = n, y = Delta_div_n, color = Rule, group = Rule)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("ER" = "black", "BF" = "blue", "PR" = "red", "PRmax" = "purple", "SR" = "darkorange")) +
    labs(
      title = expression("Classic scaling of the transition window " * Delta / n),
      x = "Network Size (n)",
      y = expression(Delta / n)
    ) +
    theme_classic() +
    theme(legend.position = "top")
  
  # C/n Curves
  all_curve_results$n_factor <- as.factor(all_curve_results$n)
  
  p_curves <- ggplot(all_curve_results, aes(x = r, y = C_n, 
                                            color = Rule, 
                                            linetype = n_factor, 
                                            group = interaction(Rule, n_factor))) +
    geom_line(linewidth = 0.8) +
    coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1)) +
    scale_color_manual(values = c("ER" = "black", "BF" = "blue", "PR" = "red", "PRmax" = "purple", "SR" = "darkorange")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
    labs(
      title = "Percolation transition curves",
      x = "r (Edges / Nodes)",
      y = "C / n",
      color = "Model",
      linetype = "Size (n)"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.key.width = unit(1.5, "cm")
    )
  
  # Print plots to RStudio Viewer
  print(p_loglog)
  print(p_delta_classic)
  print(p_curves)
  
  # Save plots to disk
  cat("\nSaving plots as PNG images...\n")
  ggsave("plot_1_loglog_scaling.png", plot = p_loglog, width = 8, height = 6, dpi = 300)
  ggsave("plot_2_delta_classic.png", plot = p_delta_classic, width = 8, height = 6, dpi = 300)
  ggsave("plot_3_cn_curves.png", plot = p_curves, width = 9, height = 6, dpi = 300)
  cat("Plots saved successfully.\n")
  
  invisible(list(plot_loglog = p_loglog, plot_delta_classic = p_delta_classic, plot_curves = p_curves))
}

set.seed(42)

# Many more data points for the asymptotic regime
n_range <- c(500, 5000, 20000, 50000, 75000, 100000, 125000, 150000, 175000, 200000, 250000, 300000, 400000, 500000)
n_to_plot <- c(500, 50000, 100000, 500000)

# Run the experiment, regression is performed only for network sizes >= 100000
run_scaling_experiment(n_values = n_range, subset_plot_n = n_to_plot, num_sims = 5, min_fit_n = 100000)