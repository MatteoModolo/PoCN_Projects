library(plotly)
library(dplyr)

# Fast Union-Find Simulator
get_delta_uf <- function(n, gamma, A) {
  max_edges <- floor(n * 1.5)
  parent <- 1:n
  comp_size <- rep(1, n)
  max_c <- 1
  
  t0 <- NA; t1 <- NA
  threshold_0 <- n^gamma
  threshold_1 <- A * n
  
  u1 <- sample(1:n, max_edges, replace = TRUE)
  v1 <- sample(1:n, max_edges, replace = TRUE)
  u2 <- sample(1:n, max_edges, replace = TRUE)
  v2 <- sample(1:n, max_edges, replace = TRUE)
  
  for (step in 1:max_edges) {
    r_u1 <- u1[step]; while(parent[r_u1] != r_u1) { parent[r_u1] <- parent[parent[r_u1]]; r_u1 <- parent[r_u1] }
    r_v1 <- v1[step]; while(parent[r_v1] != r_v1) { parent[r_v1] <- parent[parent[r_v1]]; r_v1 <- parent[r_v1] }
    r_u2 <- u2[step]; while(parent[r_u2] != r_u2) { parent[r_u2] <- parent[parent[r_u2]]; r_u2 <- parent[r_u2] }
    r_v2 <- v2[step]; while(parent[r_v2] != r_v2) { parent[r_v2] <- parent[parent[r_v2]]; r_v2 <- parent[r_v2] }
    
    # Product Rule
    if ((comp_size[r_u1] * comp_size[r_v1]) <= (comp_size[r_u2] * comp_size[r_v2])) {
      ru <- r_u1; rv <- r_v1
    } else {
      ru <- r_u2; rv <- r_v2
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
    
    if (is.na(t0) && max_c >= threshold_0) t0 <- step - 1
    if (max_c > threshold_1) {
      t1 <- step
      break 
    }
  }
  return(t1 - t0)
}

# 2. Negative Log-Likelihood Function
mle_nll <- function(params, gamma_data, beta_data) {
  mu <- params[1]
  lambda <- params[2]
  sigma <- params[3]
  
  if (sigma <= 0 || lambda <= 0) return(1e9)
  
  pred_beta <- (mu - gamma_data) / lambda
  -sum(dnorm(beta_data, mean = pred_beta, sd = sigma, log = TRUE))
}

# 3. 3D Experiment Wrapper (Generates Multiple Plots)
run_mle_beta_experiment <- function(num_sims = 3) {
  
  gamma_values <- seq(0, 0.6, by = 0.1)
  n_values <- c(2000, 5000, 10000, 20000) 
  A_values <- c(0.2, 0.4, 0.6, 0.8)
  
  grid <- expand.grid(gamma = gamma_values, n = n_values, A = A_values)
  grid$empirical_beta <- NA
  
  cat("=== Simulating and Extracting Empirical Beta (Local Fitting) ===\n")
  for (i in 1:nrow(grid)) {
    g_val <- grid$gamma[i]
    target_n <- grid$n[i]
    a_val <- grid$A[i]
    
    cat(sprintf("Target n = %5d, gamma = %.1f, A = %.1f ", target_n, g_val, a_val))
    
    # Define a local neighborhood around the target N (+/- 20%)
    # This guarantees the mean of the range is exactly target_n
    local_n_range <- round(target_n * c(0.8, 0.9, 1.0, 1.1, 1.2))
    
    mean_deltas <- numeric(length(local_n_range))
    
    for (k in seq_along(local_n_range)) {
      cat(".")
      current_n <- local_n_range[k]
      deltas <- numeric(num_sims)
      for (j in 1:num_sims) {
        deltas[j] <- get_delta_uf(current_n, g_val, a_val)
      }
      mean_deltas[k] <- mean(deltas)
    }
    cat("\n")
    
    # Calculate true local slope (beta) via log-log regression over the neighborhood
    fit <- lm(log(mean_deltas) ~ log(local_n_range))
    grid$empirical_beta[i] <- coef(fit)["log(local_n_range)"]
  }
  
  
  
  plots_list <- list()
  mle_results <- list()
  
  for (a_val in A_values) {
    # Subset data for the specific A value
    sub_grid <- grid %>% filter(A == a_val)
    
    # Fit MLE independently for this A
    mle_fit <- optim(
      par = c(mu = 1.0, lambda = 1.0, sigma = 0.1),
      fn = mle_nll,
      gamma_data = sub_grid$gamma,
      beta_data = sub_grid$empirical_beta,
      method = "L-BFGS-B",
      lower = c(0.1, 0.1, 1e-5),
      upper = c(5.0, 5.0, 1.0)
    )
    
    mu_est <- mle_fit$par["mu"]
    lambda_est <- mle_fit$par["lambda"]
    
    mle_results[[as.character(a_val)]] <- list(mu = mu_est, lambda = lambda_est)
    
    cat(sprintf("A = %.1f -> Extracted μ: %.4f, λ: %.4f\n", a_val, mu_est, lambda_est))
    
    # Generate the theoretical plane using extracted parameters
    teoria_matrix <- outer(n_values, gamma_values, FUN = function(n, g) {
      (mu_est - g) / lambda_est
    })
    
    # Create independent Plotly figure for this A
    fig <- plot_ly() %>%
      add_trace(
        data = sub_grid, 
        x = ~gamma, y = ~n, z = ~empirical_beta,
        type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = ~empirical_beta, colorscale = "Viridis", showscale = TRUE, line = list(color = "black", width = 1)),
        name = "Empirical Beta (Local Fit)"
      ) %>%
      add_surface(
        x = ~gamma_values, y = ~n_values, z = ~teoria_matrix,
        opacity = 0.4, colorscale = list(c(0, 1), c("lightgray", "gray")), showscale = FALSE,
        name = "Fitted Plane"
      ) %>%
      layout(
        title = sprintf("Threshold A = %.1f | Fitted: μ = %.3f, λ = %.3f", a_val, mu_est, lambda_est),
        scene = list(
          xaxis = list(title = "Gamma"),
          yaxis = list(title = "Target System Size (n)", type = "log"),
          zaxis = list(title = "Fitted Local Beta "),
          camera = list(eye = list(x = 1.6, y = -1.4, z = 1.1))
        )
      )
    
    # Print the plot to the viewer
    print(fig)
    plots_list[[as.character(a_val)]] <- fig
  }
  
  return(list(grid = grid, mle_results = mle_results, plots = plots_list))
}


set.seed(42)
mle_3d_data <- run_mle_beta_experiment(num_sims = 50)