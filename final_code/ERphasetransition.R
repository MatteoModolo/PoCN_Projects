library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(42)

# used parameters
n_values <- seq(200, 3000, by = 200) # Range of network sizes
r_values <- c(0.2, 0.4, 0.6, 0.8)    # Edge to node ratios
num_sims <- 300                      # Simulations per (n, r) pair

# generate union find networks
results_list <- list()
row_idx <- 1

for (r in r_values) {
  for (n in n_values) {
    m <- floor(r * n) # Number of edges
    lcc_sims <- numeric(num_sims)
    
    for (i in 1:num_sims) {
      # Create Erdős-Rényi graph G(n, m)
      g <- sample_gnm(n, m, directed = FALSE)
      
      # Find the size of the Largest Connected Component (LCC)
      lcc_sims[i] <- max(components(g)$csize)
    }
    
    # Store the averaged LCC
    results_list[[row_idx]] <- data.frame(
      n = n,
      r = as.factor(r),
      r_num = r,       
      avg_LCC = mean(lcc_sims)
    )
    row_idx <- row_idx + 1
  }
}

df <- bind_rows(results_list)

# fitting with 3 different hypthesis
df$pred_linear <- NA
df$pred_log <- NA

for (current_r in r_values) {
  sub_df <- df[df$r_num == current_r, ]
  
  # Fit both models globally for each r
  fit_linear <- lm(avg_LCC ~ n, data = sub_df)
  fit_log    <- lm(avg_LCC ~ log(n), data = sub_df)
  
  # Store predictions for plotting
  df$pred_linear[df$r_num == current_r] <- predict(fit_linear)
  df$pred_log[df$r_num == current_r]    <- predict(fit_log)
}

df_plot <- df %>%
  pivot_longer(
    cols = c("pred_linear", "pred_log"),
    names_to = "Model",
    values_to = "Prediction"
  ) %>%
  mutate(Model = recode(Model, 
                        "pred_linear" = "Linear Fit",
                        "pred_log" = "Logarithmic Fit"))

p <- ggplot(df_plot, aes(x = n)) +
  # Empirical points
  geom_point(aes(y = avg_LCC), color = "black", alpha = 0.5, size = 2.5) +
  # Model fit lines
  geom_line(aes(y = Prediction, color = Model, linetype = Model), linewidth = 1.2) +
  facet_wrap(~ r, scales = "free_y", labeller = label_bquote(r == .(as.character(r)))) +
  scale_color_manual(values = c("Linear Fit" = "#d95f02", 
                                "Logarithmic Fit" = "#1b9e77")) +
  scale_linetype_manual(values = c("Linear Fit" = "dashed", 
                                   "Logarithmic Fit" = "solid")) +
  labs(
    title = "The Phase Transition: Logarithmic vs Linear Scaling",
    subtitle = "Notice how the Logarithmic fit perfectly tracks r < 0.5, but completely fails for r > 0.5",
    x = "Network Size (n)",
    y = "Size of Largest Connected Component (LCC)",
    color = "Fitted Model",
    linetype = "Fitted Model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        strip.background = element_rect(fill = "#f0f0f0"),
        strip.text = element_text(face = "bold", size = 12))

print(p)