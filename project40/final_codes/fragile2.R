library(readr)
library(dplyr)
library(igraph)
library(ggplot2)

# input
input_csv <- "ultimate_distribution_confrontation.csv"

df <- read_csv(input_csv, show_col_types = FALSE)

# Extract parameters for specific evolutionary milestones
target_years <- c("2002", "2016", "2024")
selected_months <- c()

for (y in target_years) {
  matched_month <- df %>% filter(grepl(paste0("^", y), month)) %>% pull(month) %>% first()
  if (!is.na(matched_month)) selected_months <- c(selected_months, matched_month)
}

cat("--------------------------------------------------\n")
cat("SELECTED NETWORK SNAPSHOTS FOR SIMULATION:\n")
print(selected_months)
cat("--------------------------------------------------\n\n")

# Configurations
N <- 100000                 # Number of nodes per simulation
n_realizations <- 50       # Number of realizations for smoothing
p_seq <- seq(1, 0, by = -0.05) 
network_types <- c("Empirical Power-Law", "Empirical Log-Normal", 
                   "Empirical Exponential", "Erdos-Renyi (Baseline)")

# Generator function now accepts 'er_avg_deg' to dynamically match the ER baseline
generate_empirical_graph <- function(type, N, p_xmin, p_alpha, p_ln_mu, p_ln_sig, p_exp_l, er_avg_deg = NULL) {
  
  if (type == "Empirical Power-Law") {
    degs <- round(p_xmin * (1 - runif(N))^(-1/(p_alpha - 1)))
    
  } else if (type == "Empirical Log-Normal") {
    F_xmin <- plnorm(p_xmin, meanlog = p_ln_mu, sdlog = p_ln_sig)
    u <- runif(N, min = F_xmin, max = 1)
    degs <- round(qlnorm(u, meanlog = p_ln_mu, sdlog = p_ln_sig))
    
  } else if (type == "Empirical Exponential") {
    degs <- round(rexp(N, rate = p_exp_l) + p_xmin)
    
  } else if (type == "Erdos-Renyi (Baseline)") {
    # MATHEMATICAL FIX: Match ER perfectly to the Power-Law's average degree
    er_p <- er_avg_deg / (N - 1)
    return(sample_gnp(n = N, p = er_p))
  }
  
  # Clean degree sequence (min degree 1, max N-1, even sum)
  degs[degs < 1] <- 1
  degs[degs >= N] <- N - 1 
  if (sum(degs) %% 2 != 0) degs[1] <- degs[1] + 1
  
  return(simplify(sample_degseq(degs, method = "simple")))
}

# Simulate percolation
simulate_percolation <- function(g, attack_type, p_seq) {
  N_orig <- vcount(g)
  lcc_sizes <- numeric(length(p_seq))
  
  # TARGETED ATTACK: Sort by degree descending (Remove Titans first)
  # RANDOM FAILURE: Random shuffle (Natural disasters)
  node_order <- if(attack_type == "targeted") order(degree(g), decreasing = TRUE) else sample(V(g))
  
  for (i in seq_along(p_seq)) {
    p <- p_seq[i]
    num_keep <- round(N_orig * p)
    
    if (num_keep == 0) lcc_sizes[i] <- 0
    else if (num_keep == N_orig) lcc_sizes[i] <- max(components(g)$csize) / N_orig
    else {
      # Keep nodes from the bottom of the removal list
      nodes_to_keep <- node_order[(N_orig - num_keep + 1):N_orig]
      lcc_sizes[i] <- max(components(induced_subgraph(g, nodes_to_keep))$csize) / N_orig
    }
  }
  return(lcc_sizes)
}


results <- data.frame()

for (target_month in selected_months) {
  cat(sprintf("Simulating topology parameters for: %s\n", target_month))
  
  row_data <- df %>% filter(month == target_month)
  p_xmin   <- row_data$x_min
  p_alpha  <- row_data$alpha
  p_ln_mu  <- row_data$LN_mu
  p_ln_sig <- row_data$LN_sigma
  p_exp_l  <- row_data$EXP_lambda
  
  for (r in 1:n_realizations) {
    # 1. Generate Power-Law FIRST so we can measure its actual average degree
    g_pl <- generate_empirical_graph("Empirical Power-Law", N, p_xmin, p_alpha, p_ln_mu, p_ln_sig, p_exp_l)
    current_avg_degree <- mean(degree(g_pl)) # Extract exact average degree
    
    # 2. Generate the rest, passing the exact average degree to the ER baseline
    graphs <- list(
      "Empirical Power-Law" = g_pl,
      "Empirical Log-Normal" = generate_empirical_graph("Empirical Log-Normal", N, p_xmin, p_alpha, p_ln_mu, p_ln_sig, p_exp_l),
      "Empirical Exponential" = generate_empirical_graph("Empirical Exponential", N, p_xmin, p_alpha, p_ln_mu, p_ln_sig, p_exp_l),
      "Erdos-Renyi (Baseline)" = generate_empirical_graph("Erdos-Renyi (Baseline)", N, p_xmin, p_alpha, p_ln_mu, p_ln_sig, p_exp_l, er_avg_deg = current_avg_degree)
    )
    
    # 3. Subject all graphs to both attack strategies
    for (type in network_types) {
      g <- graphs[[type]]
      
      lcc_random <- simulate_percolation(g, "random", p_seq)
      results <- rbind(results, data.frame(
        Snapshot = target_month, Type = type, Realization = r, 
        Attack = "Random", p = p_seq, LCC = lcc_random
      ))
      
      lcc_targeted <- simulate_percolation(g, "targeted", p_seq)
      results <- rbind(results, data.frame(
        Snapshot = target_month, Type = type, Realization = r, 
        Attack = "Targeted", p = p_seq, LCC = lcc_targeted
      ))
    }
  }
}

plot_data <- results %>%
  group_by(Snapshot, Type, Attack, p) %>%
  summarise(LCC_mean = mean(LCC), .groups = 'drop')


my_colors <- c("Erdos-Renyi (Baseline)" = "#E41A1C", 
               "Empirical Log-Normal" = "#4DAF4A", 
               "Empirical Power-Law" = "#377EB8", 
               "Empirical Exponential" = "#984EA3")

p_random <- ggplot(plot_data %>% filter(Attack == "Random"), 
                   aes(x = p, y = LCC_mean, color = Type, linetype = Type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = my_colors) +
  scale_x_reverse() + 
  facet_wrap(~ Snapshot, ncol = 3) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  labs(title = "Evolution of Random Failure Tolerance (2002 vs 2016 vs 2024)",
       
       x = "Fraction of Nodes Remaining (p)", y = "Normalized LCC Size")

ggsave("percolation_evolution_random.png", plot = p_random, width = 10, height = 5, bg = "white")

p_targeted <- ggplot(plot_data %>% filter(Attack == "Targeted"), 
                     aes(x = p, y = LCC_mean, color = Type, linetype = Type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = my_colors) +
  scale_x_reverse() +
  facet_wrap(~ Snapshot, ncol = 3) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  labs(title = "Evolution of Targeted Attack Vulnerability (2002 vs 2016 vs 2024)",
       
       x = "Fraction of Nodes Remaining (p)", y = "Normalized LCC Size")

ggsave("percolation_evolution_targeted.png", plot = p_targeted, width = 10, height = 5, bg = "white")
