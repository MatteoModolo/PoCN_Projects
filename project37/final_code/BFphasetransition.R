library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)


# Parameters for the scan
scan_n_values <- seq(200, 10000, by = 200)
scan_r_values <- seq(0.1, 0.9, by = 0.05)
scan_sims <- 40

# Collect data for the scan
scan_results <- list()
idx <- 1
for (r in scan_r_values) {
  for (n in scan_n_values) {
    m <- floor(r * n)
    lccs <- numeric(scan_sims)
    for (i in 1:scan_sims) {
      g <- sample_gnm(n, m, directed = FALSE)
      lccs[i] <- max(components(g)$csize)
    }
    scan_results[[idx]] <- data.frame(r = r, n = n, avg_LCC = mean(lccs))
    idx <- idx + 1
  }
}
df_scan <- bind_rows(scan_results)

# Calculate Log-Likelihoods for Linear and Log models at EACH r
ll_table <- data.frame(r = scan_r_values, ll_log = NA, ll_lin = NA)

for (i in 1:nrow(ll_table)) {
  curr_r <- ll_table$r[i]
  sub_df <- df_scan[df_scan$r == curr_r, ]
  
  fit_lin <- lm(avg_LCC ~ n, data = sub_df)
  fit_log <- lm(avg_LCC ~ log(n), data = sub_df)
  
  ll_table$ll_lin[i] <- as.numeric(logLik(fit_lin))
  ll_table$ll_log[i] <- as.numeric(logLik(fit_log))
}

# Find the MLE for the Phase Transition point (r_c)
candidate_rc <- scan_r_values[scan_r_values > 0.2 & scan_r_values < 0.8]
mle_scores <- numeric(length(candidate_rc))

for (i in seq_along(candidate_rc)) {
  rc <- candidate_rc[i]
  # Assumption: Logarithmic before rc, Linear at and after rc
  score <- sum(ll_table$ll_log[ll_table$r < rc]) + sum(ll_table$ll_lin[ll_table$r >= rc])
  mle_scores[i] <- score
}

# Extract the best model statistics
best_rc <- candidate_rc[which.max(mle_scores)]
best_ll_pt <- max(mle_scores)

# Calculate Global Likelihoods for pure models
global_ll_lin <- sum(ll_table$ll_lin)
global_ll_log <- sum(ll_table$ll_log)

# Calculate Log-Likelihood Ratios (LLR) wrt Phase Transition Model
llr_lin <- global_ll_lin - best_ll_pt
llr_log <- global_ll_log - best_ll_pt

# Create a dataframe of the likelihood landscape for plotting
mle_results <- data.frame(
  Candidate_rc = candidate_rc, 
  Combined_Log_Likelihood = round(mle_scores, 2)
)

# --- CONSOLE OUTPUTS ---
cat("\n--- Maximum Likelihood Estimation Landscape ---\n")
print(mle_results)

cat(sprintf("\n>> MLE Estimated Critical Point (r_c): %.2f <<\n", best_rc))
cat(sprintf(">> Peak Log-Likelihood (Phase Transition Model): %.2f <<\n\n", best_ll_pt))

cat("--- Log-Likelihood Ratios (vs Phase Transition Hypothesis) ---\n")
cat("Note: A massive negative LLR mathematically rejects the pure models in favor of the PT model.\n")
cat(sprintf("Global LL (Strictly Linear):      %.2f\n", global_ll_lin))
cat(sprintf("Global LL (Strictly Logarithmic): %.2f\n", global_ll_log))
cat(sprintf("LLR (Linear vs PT): %.2f\n", llr_lin))
cat(sprintf("LLR (Log vs PT):    %.2f\n\n", llr_log))

# PLOT THE LIKELIHOOD LANDSCAPE
p_mle <- ggplot(mle_results, aes(x = Candidate_rc, y = Combined_Log_Likelihood)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(color = "purple", size = 3) +
  geom_vline(xintercept = best_rc, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Log-Likelihood Landscape for the ER Phase Transition",
    subtitle = sprintf("The peak indicates the most probable critical point (r_c = %.2f)", best_rc),
    x = "Candidate Critical Point (r_c)",
    y = "Combined Log-Likelihood Score"
  ) +
  theme_classic() +
  annotate("text", x = best_rc + 0.02, y = min(mle_scores), 
           label = "Maximum Likelihood", hjust = 0, angle = 90)

print(p_mle)



# Dynamically pick 2 points before and 2 points after the estimated r_c
plot_r_values <- c(best_rc - 0.15, best_rc - 0.05, best_rc + 0.05, best_rc + 0.15)
plot_r_values <- round(plot_r_values, 2) # Clean up floating points

plot_n_values <- seq(200, 3000, by = 200)
plot_sims <- 150

cat(sprintf("Selected r values for plotting: %s\n", paste(plot_r_values, collapse = ", ")))
cat("Running heavy simulations. Please wait...\n")

plot_results <- list()
idx <- 1
for (r in plot_r_values) {
  for (n in plot_n_values) {
    m <- floor(r * n)
    lccs <- numeric(plot_sims)
    for (i in 1:plot_sims) {
      g <- sample_gnm(n, m, directed = FALSE)
      lccs[i] <- max(components(g)$csize)
    }
    plot_results[[idx]] <- data.frame(n = n, r = as.factor(r), r_num = r, avg_LCC = mean(lccs))
    idx <- idx + 1
  }
}
df_plot_data <- bind_rows(plot_results)

# Fit models for plotting
df_plot_data$pred_linear <- NA
df_plot_data$pred_log <- NA

for (curr_r in plot_r_values) {
  sub_df <- df_plot_data[df_plot_data$r_num == curr_r, ]
  
  fit_linear <- lm(avg_LCC ~ n, data = sub_df)
  fit_log    <- lm(avg_LCC ~ log(n), data = sub_df)
  
  df_plot_data$pred_linear[df_plot_data$r_num == curr_r] <- predict(fit_linear)
  df_plot_data$pred_log[df_plot_data$r_num == curr_r]    <- predict(fit_log)
}

# Pivot data for ggplot
df_final <- df_plot_data %>%
  pivot_longer(
    cols = c("pred_linear", "pred_log"),
    names_to = "Model",
    values_to = "Prediction"
  ) %>%
  mutate(Model = recode(Model, 
                        "pred_linear" = "Linear Fit",
                        "pred_log" = "Logarithmic Fit"))

# Generate the faceted plot
p <- ggplot(df_final, aes(x = n)) +
  geom_point(aes(y = avg_LCC), color = "black", alpha = 0.5, size = 2.5) +
  geom_line(aes(y = Prediction, color = Model, linetype = Model), linewidth = 1.2) +
  facet_wrap(~ r, scales = "free_y", labeller = label_bquote(r == .(as.character(r)))) +
  scale_color_manual(values = c("Linear Fit" = "#d95f02", "Logarithmic Fit" = "#1b9e77")) +
  scale_linetype_manual(values = c("Linear Fit" = "dashed", "Logarithmic Fit" = "solid")) +
  labs(
    title = sprintf("Phase Transition Verified: Estimated r_c = %.2f", best_rc),
    subtitle = "Top row: Subcritical regime (r < r_c) | Bottom row: Supercritical regime (r > r_c)",
    x = "Network Size (n)",
    y = "Size of Largest Connected Component (LCC)",
    color = "Fitted Model",
    linetype = "Fitted Model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        strip.background = element_rect(fill = "#e5e5e5"),
        strip.text = element_text(face = "bold", size = 12))

print(p)