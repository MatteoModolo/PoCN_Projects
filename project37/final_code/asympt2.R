library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# Union-Find simulator
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
  
  # Array for max component size
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
      
      # Model rules
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
    
    # Record max component size
    history[step] <- max_c
  }
  
  return(history)
}

run_asymptotic_experiment <- function(n_values, rules = c("ER", "BF", "PR"), num_sims = 10) {
  
  results <- data.frame()
  

  for (n in n_values) {
    cat(sprintf("Simulating n = %s\n", format(n, big.mark=",")))
    
    threshold_0 <- sqrt(n)
    threshold_1 <- 0.5 * n
    
    for (rule in rules) {
      
      # Run simulations
      history_matrix <- do.call(cbind, lapply(1:num_sims, function(i) get_trajectory_uf(n, rule)))
      
      # Average into Master Curve
      master_avg_history <- rowMeans(history_matrix)
      
      # Interpolate t0 and t1
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
  
  # Reshape data
  plot_data <- results %>%
    pivot_longer(cols = c(t0_density, t1_density), 
                 names_to = "Metric", 
                 values_to = "Density") %>%
    mutate(Metric = ifelse(Metric == "t0_density", "t0/n (Start)", "t1/n (End)"))
  
  
  
  fit_powerlaw_asymptote <- function(data_subset, metric_label) {
    # Initial estimates
    c_init <- tail(data_subset$Density, 1)
    beta_init <- -0.33
    a_init <- (data_subset$Density[1] - c_init) / (data_subset$n[1]^beta_init)
    
    # Try non-linear fit
    tryCatch({
      fit <- nls(Density ~ c + a * (n ^ beta), 
                 data = data_subset, 
                 start = list(c = c_init, a = a_init, beta = beta_init),
                 control = nls.control(maxiter = 1000, warnOnly = TRUE))
      
      # Get parameters
      s_fit <- summary(fit)
      c_est <- coef(fit)["c"]
      beta_est <- coef(fit)["beta"]
      se_c <- s_fit$parameters["c", "Std. Error"]
      
      data.frame(Metric = metric_label, 
                 Asymptote = c_est, 
                 Fitted_Beta = beta_est,
                 CI_Lower = c_est - 1.96 * se_c, 
                 CI_Upper = c_est + 1.96 * se_c)
      
    }, error = function(e) {
      # Fallback: Linear model (beta = -1/3)
      fit_lm <- lm(Density ~ I(n^(-1/3)), data = data_subset)
      ci_lm <- confint(fit_lm, "(Intercept)", level = 0.95)
      
      data.frame(Metric = metric_label, 
                 Asymptote = coef(fit_lm)["(Intercept)"], 
                 Fitted_Beta = -0.3333,
                 CI_Lower = ci_lm[1], 
                 CI_Upper = ci_lm[2])
    })
  }
  
  fits_t0 <- plot_data %>% filter(Metric == "t0/n (Start)") %>% group_by(Rule) %>% do(fit_powerlaw_asymptote(., "t0/n (Start)")) %>% ungroup()
  fits_t1 <- plot_data %>% filter(Metric == "t1/n (End)") %>% group_by(Rule) %>% do(fit_powerlaw_asymptote(., "t1/n (End)")) %>% ungroup()
  
  final_asymptotes <- rbind(fits_t0, fits_t1) %>% arrange(Rule, Metric)
  print(final_asymptotes %>% mutate(across(where(is.numeric), ~ round(., 4))))
  

  
  pr_data <- plot_data %>% filter(Rule == "PR")
  pr_data$is_t0 <- ifelse(pr_data$Metric == "t0/n (Start)", 1, 0)
  pr_data$is_t1 <- ifelse(pr_data$Metric == "t1/n (End)", 1, 0)
  
  # Initial estimates
  c_init <- mean(tail(pr_data$Density, 2))
  beta_init <- -0.33
  a0_init <- (pr_data$Density[pr_data$is_t0==1][1] - c_init) / (pr_data$n[pr_data$is_t0==1][1]^beta_init)
  a1_init <- (pr_data$Density[pr_data$is_t1==1][1] - c_init) / (pr_data$n[pr_data$is_t1==1][1]^beta_init)
  
  joint_pr_results <- tryCatch({
    fit <- nls(Density ~ c + is_t0 * a0 * (n ^ beta0) + is_t1 * a1 * (n ^ beta1),
               data = pr_data,
               start = list(c = c_init, a0 = a0_init, beta0 = beta_init, a1 = a1_init, beta1 = beta_init),
               control = nls.control(maxiter = 1000, warnOnly = TRUE))
    
    s_fit <- summary(fit)
    c_val <- coef(fit)["c"]
    se_val <- s_fit$parameters["c", "Std. Error"]
    
    list(
      asymptote = c_val,
      ci_low = c_val - 1.96 * se_val,
      ci_high = c_val + 1.96 * se_val,
      b0 = coef(fit)["beta0"],
      b1 = coef(fit)["beta1"],
      fallback = FALSE
    )
  }, error = function(e) {
    # Fallback: Linear model
    fit_lm <- lm(Density ~ 1 + Metric:I(n^(-1/3)), data = pr_data)
    ci_lm <- confint(fit_lm, "(Intercept)", level = 0.95)
    
    list(
      asymptote = coef(fit_lm)["(Intercept)"],
      ci_low = ci_lm[1],
      ci_high = ci_lm[2],
      b0 = -0.3333,
      b1 = -0.3333,
      fallback = TRUE
    )
  })
  
  pr_asymptote <- joint_pr_results$asymptote
  
  if(joint_pr_results$fallback) {
    cat("[!] NLS failed for PR. Using linear model.\n")
  }
  
  cat(sprintf("Joint Asymptote (r_c): %.5f\n", pr_asymptote))
  cat(sprintf("95%% CI: [%.5f, %.5f]\n", joint_pr_results$ci_low, joint_pr_results$ci_high))
  cat(sprintf("Fitted Beta (Start): %.4f, Fitted Beta (End): %.4f\n", joint_pr_results$b0, joint_pr_results$b1))
  
  
  # plots
  base_colors <- c("ER" = "black", "BF" = "blue", "PR" = "red", "SR" = "darkorange")
  
  p <- ggplot(plot_data, aes(x = n, y = Density, color = Rule, shape = Metric)) +
    
    # Plot data points
    geom_point(size = 3.5, alpha = 0.9) +
    geom_line(aes(group = interaction(Rule, Metric)), alpha = 0.4, linewidth = 1) +
    
    # Plot asymptotes
    geom_hline(data = final_asymptotes, 
               aes(yintercept = Asymptote, color = Rule, linetype = Metric), 
               linewidth = 1.2, alpha = 0.7) +
    
    # Add PR joint limit label
    annotate("text", x = max(n_values), y = pr_asymptote + 0.015, 
             label = sprintf("PR Joint r_c = %.4f", pr_asymptote), 
             color = "red", fontface = "bold", hjust = 1) +
    
    scale_x_log10(breaks = n_values, labels = scales::comma) +
    scale_color_manual(values = base_colors) +
    scale_linetype_manual(values = c("t0/n (Start)" = "dotted", "t1/n (End)" = "dashed")) +
    
    labs(
      title = "Asymptotic Convergence of t0/n and t1/n",
      x = "Network Size (n) - Log Scale",
      y = "Edge Density (r = t/n)",
      color = "Model",
      shape = "Boundary",
      linetype = "Asymptotic Limit"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  print(p)
  ggsave("asymptotic_convergence_extended.png", plot = p, width = 11, height = 7, dpi = 300)
  
  return(list(raw_data = plot_data, independent_asymptotes = final_asymptotes, joint_pr_asymptote = pr_asymptote))
}


set.seed(42)

# Extended ranges
n_test_values <- c(5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000)

# Run simulations
final_data <- run_asymptotic_experiment(n_test_values, rules = c("ER", "BF", "PR"), num_sims = 50)