library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(42)

# set parameters
n_values <- seq(200, 30000, by = 200) # Range of network sizes
r_values <- c(0.2, 0.4, 0.6, 0.8)    # 2 before rc, 2 after rc
num_sims <- 50                     


# generate syntetic networks
results_list <- list()
row_idx <- 1

for (r in r_values) {
  for (n in n_values) {
    m <- floor(r * n) 
    lcc_sims <- numeric(num_sims)
    
    for (i in 1:num_sims) {
      g <- sample_gnm(n, m, directed = FALSE)
      lcc_sims[i] <- max(components(g)$csize)
    }
    
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

# Convert to log scale to find the polynomial scaling
df$log_n <- log(df$n)
df$log_LCC <- log(df$avg_LCC)

# fitting in log scale

# Global fitting
# Use as.factor(r) so each curve gets its own intercept but they share the same slope
fit_global <- lm(log_LCC ~ log_n + as.factor(r), data = df)
beta_global <- coef(fit_global)["log_n"]

# Two scaling model
df_sub <- df %>% filter(r_num < 0.5)
df_super <- df %>% filter(r_num >= 0.5)

fit_sub <- lm(log_LCC ~ log_n + as.factor(r), data = df_sub)
fit_super <- lm(log_LCC ~ log_n + as.factor(r), data = df_super)

beta_sub <- coef(fit_sub)["log_n"]
beta_super <- coef(fit_super)["log_n"]

# Extract Log-Likelihoods 
ll_global <- as.numeric(logLik(fit_global))
ll_pt <- as.numeric(logLik(fit_sub)) + as.numeric(logLik(fit_super))

# Generate predictions and transform back to normal scale using exp()
df$pred_global <- exp(predict(fit_global, newdata = df))
df$pred_pt <- NA
df$pred_pt[df$r_num < 0.5] <- exp(predict(fit_sub, newdata = df_sub))
df$pred_pt[df$r_num >= 0.5] <- exp(predict(fit_super, newdata = df_super))


# Print llr
cat("Note: LCC ~ n^beta\n")
cat(sprintf("Global Exponent:      beta = %.3f\n", beta_global))
cat(sprintf("Subcritical Exponent (r < 0.5):          beta = %.3f\n", beta_sub))
cat(sprintf("Supercritical Exponent (r > 0.5):        beta = %.3f\n", beta_super))


cat(sprintf("Global LL (One Power Fit):    %.2f\n", ll_global))
cat(sprintf("Global LL (Phase Transition): %.2f\n", ll_pt))

cat(sprintf("LLR (Global vs PT): %.2f\n\n", ll_global - ll_pt))


# Plot
df_plot <- df %>%
  pivot_longer(
    cols = c("pred_global", "pred_pt"),
    names_to = "Model",
    values_to = "Prediction"
  ) %>%
  mutate(Model = recode(Model, 
                        "pred_global" = sprintf("Global Power Fit (beta = %.2f)", beta_global),
                        "pred_pt" = sprintf("Phase Transition Fit (beta1 = %.2f, beta2 = %.2f)", beta_sub, beta_super)))

# Generate the plot
p <- ggplot(df_plot, aes(x = n)) +
  geom_point(aes(y = avg_LCC), color = "black", alpha = 0.5, size = 2.5) +
  geom_line(aes(y = Prediction, color = Model, linetype = Model), linewidth = 1.2) +
  facet_wrap(~ r, scales = "free_y", labeller = label_bquote(r == .(as.character(r)))) +
  scale_color_manual(values = setNames(c("#d95f02", "blue"), unique(df_plot$Model))) +
  scale_linetype_manual(values = setNames(c("dashed", "solid"), unique(df_plot$Model))) +
  labs(
    title = "Phase transition verification",
    x = "Network Size (n)",
    y = "Size of Largest Connected Component (LCC)",
    color = "Fitted Model",
    linetype = "Fitted Model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        strip.background = element_rect(fill = "#e5e5e5"),
        strip.text = element_text(face = "bold", size = 12))

print(p)