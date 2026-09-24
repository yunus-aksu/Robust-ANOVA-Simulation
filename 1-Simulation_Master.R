# ==============================================================================
# Monte Carlo Simulation: One-Way ANOVA Alternatives Benchmark
# ------------------------------------------------------------------------------
# Author: Yunus Aksu
# Repository: https://github.com/yunus-aksu/Robust-ANOVA-Simulation
# Description: Benchmarking Classic ANOVA, Welch ANOVA, Brown-Forsythe, 
#              Kruskal-Wallis, 20% Trimmed Welch ANOVA, and Harrell-Davis 
#              Quantile (Median) with Percentile Bootstrap under 
#              Heteroscedasticity and Extreme Skewness.
# ==============================================================================

# --- 0. PACKAGE INSTALLATION AND SETUP ---
required_pkgs <- c("WRS2", "onewaytests", "doParallel", "foreach")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(WRS2)        # t1way (20% Trimmed Welch ANOVA)
library(onewaytests) # bf.test (Brown-Forsythe ANOVA)
library(doParallel)  
library(foreach)     

set.seed(2026)

# --- 1. HELPER FUNCTIONS ---

# Harrell-Davis Quantile Estimator (q = 0.50 -> Median)
hd_quantile <- function(x, q = 0.5) {
  x <- sort(x)
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- 1:n
  w <- pbeta(vec / n, m1, m2) - pbeta((vec - 1) / n, m1, m2)
  sum(w * x)
}

# Percentile Bootstrap Test for Harrell-Davis Median Differences
hd_boot_test <- function(data_df, nboot = 200) {
  groups <- split(data_df$y, data_df$group)
  K <- length(groups)
  obs_meds <- sapply(groups, hd_quantile, q = 0.5)
  obs_stat <- var(obs_meds)
  
  grand_hd <- hd_quantile(data_df$y, q = 0.5)
  centered_groups <- lapply(1:K, function(k) groups[[k]] - obs_meds[k] + grand_hd)
  
  boot_stats <- numeric(nboot)
  for (b in 1:nboot) {
    boot_g <- lapply(centered_groups, function(g) sample(g, length(g), replace = TRUE))
    boot_meds <- sapply(boot_g, hd_quantile, q = 0.5)
    boot_stats[b] <- var(boot_meds)
  }
  return(mean(boot_stats >= obs_stat))
}

# --- 2. DATA GENERATION FUNCTION (With Centering and Scaling) ---

generate_data <- function(n_vec, var_vec, dist_type = "normal", effect_size = c(0, 0, 0)) {
  K <- length(n_vec)
  y <- numeric(0)
  group <- factor(rep(1:K, times = n_vec))
  
  for (i in 1:K) {
    if (dist_type == "normal") {
      y_group <- rnorm(n_vec[i], mean = effect_size[i], sd = sqrt(var_vec[i]))
    } else if (dist_type == "lognormal") {
      raw_val <- rlnorm(n_vec[i], meanlog = 0, sdlog = 1)
      # Standardizing to mean = 0, sd = 1
      standardized_val <- (raw_val - exp(0.5)) / sqrt(exp(1) * (exp(1) - 1))
      y_group <- (standardized_val * sqrt(var_vec[i])) + effect_size[i]
    } else if (dist_type == "exponential") {
      raw_val <- rexp(n_vec[i], rate = 1)
      # Centering and scaling to mean = 0, sd = 1
      centered_val <- raw_val - 1
      y_group <- (centered_val * sqrt(var_vec[i])) + effect_size[i]
    }
    y <- c(y, y_group)
  }
  return(data.frame(y = y, group = group))
}

# --- 3. TEST EVALUATION FUNCTION ---

run_tests <- function(data) {
  p_vals <- numeric(6)
  
  # 1. Classic ANOVA
  p_vals[1] <- summary(aov(y ~ group, data = data))[[1]][1, "Pr(>F)"]
  
  # 2. Welch ANOVA
  p_vals[2] <- oneway.test(y ~ group, data = data, var.equal = FALSE)$p.value
  
  # 3. Brown-Forsythe ANOVA
  p_vals[3] <- onewaytests::bf.test(y ~ group, data = data, verbose = FALSE)$p.value
  
  # 4. Kruskal-Wallis Test
  p_vals[4] <- kruskal.test(y ~ group, data = data)$p.value
  
  # 5. 20% Trimmed Welch ANOVA
  p_vals[5] <- WRS2::t1way(y ~ group, data = data, tr = 0.2)$p.value
  
  # 6. Harrell-Davis Percentile Bootstrap Test
  p_vals[6] <- hd_boot_test(data, nboot = 200)
  
  names(p_vals) <- c("ANOVA", "Welch", "BF", "KW", "Trimmed", "HD_Boot")
  return(p_vals)
}

# --- 4. SIMULATION CONFIGURATION ---

n_scenarios   <- list("10-10-10" = c(10, 10, 10), "30-30-30" = c(30, 30, 30), "5-10-20" = c(5, 10, 20))
var_scenarios <- list("1-1-1" = c(1, 1, 1), "1-2-4" = c(1, 2, 4), "4-2-1" = c(4, 2, 1))
dist_scenarios <- c("normal", "lognormal", "exponential")

scenarios <- expand.grid(
  n_id = names(n_scenarios),
  var_id = names(var_scenarios),
  dist = dist_scenarios,
  stringsAsFactors = FALSE
)

iterations <- 10000
alpha <- 0.05

# Initialize Parallel Processing Backend
cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(cores)
registerDoParallel(cl)

cat(sprintf("Simulation initialized on %d parallel cores. Total scenarios: %d\n", cores, nrow(scenarios)))

# --- 5. PARALLEL MONTE CARLO EXECUTION LOOP ---

results_master <- foreach(
  i = 1:nrow(scenarios),
  .combine = rbind,
  .packages = c("WRS2", "onewaytests")
) %dopar% {
  
  current_n    <- n_scenarios[[scenarios$n_id[i]]]
  current_var  <- var_scenarios[[scenarios$var_id[i]]]
  current_dist <- scenarios$dist[i]
  
  # Storage Matrices: Iteration x 6 Methods
  p_h0 <- matrix(NA, nrow = iterations, ncol = 6)
  colnames(p_h0) <- c("ANOVA", "Welch", "BF", "KW", "Trimmed", "HD_Boot")
  
  # -----------------------------------------------------------------
  # STAGE 1: H0 Simulation (Empirical Type I Error & Size-Adjusted Cutoffs)
  # -----------------------------------------------------------------
  for (sim in 1:iterations) {
    sim_data_H0 <- generate_data(current_n, current_var, current_dist, effect_size = c(0, 0, 0))
    p_h0[sim, ] <- run_tests(sim_data_H0)
  }
  
  t1_rates <- colMeans(p_h0 < alpha, na.rm = TRUE)
  adj_cutoffs <- apply(p_h0, 2, quantile, probs = alpha, na.rm = TRUE)
  
  # -----------------------------------------------------------------
  # STAGE 2: H1 Simulation (Standardized Effect Sizes: Cohen's d)
  # -----------------------------------------------------------------
  effects <- c(0.2, 0.5, 0.8, 1.2)
  scenario_rows <- list()
  
  for (d in effects) {
    p_h1 <- matrix(NA, nrow = iterations, ncol = 6)
    colnames(p_h1) <- colnames(p_h0)
    
    # Applying location shift to Group 3 based on Cohen's d
    eff_vec <- c(0, 0, d)
    
    for (sim in 1:iterations) {
      sim_data_H1 <- generate_data(current_n, current_var, current_dist, effect_size = eff_vec)
      p_h1[sim, ] <- run_tests(sim_data_H1)
    }
    
    raw_pwr <- colMeans(p_h1 < alpha, na.rm = TRUE)
    adj_pwr <- colMeans(sweep(p_h1, 2, adj_cutoffs, "<"), na.rm = TRUE)
    
    for (m in colnames(p_h0)) {
      scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
        Distribution = current_dist,
        N_Design = scenarios$n_id[i],
        Var_Design = scenarios$var_id[i],
        Effect_Size_d = d,
        Method = m,
        Type1_Error = t1_rates[m],
        Raw_Power = raw_pwr[m],
        Size_Adjusted_Power = adj_pwr[m]
      )
    }
  }
  
  do.call(rbind, scenario_rows)
}

stopCluster(cl)

results_master$Type1_Error         <- sprintf("%.4f", results_master$Type1_Error)
results_master$Raw_Power           <- sprintf("%.4f", results_master$Raw_Power)
results_master$Size_Adjusted_Power <- sprintf("%.4f", results_master$Size_Adjusted_Power)

# --- 6. EXPORTING RESULTS ---

write.csv(results_master, "Master_Results.csv", row.names = FALSE)

# Summary Table for Type I Error
t1_summary <- unique(results_master[, c("Distribution", "N_Design", "Var_Design", "Method", "Type1_Error")])
write.csv(t1_summary, "Type1_Error_Summary.csv", row.names = FALSE)

cat("\nSimulation successfully completed! Results saved to 'Master_Results.csv' and 'Type1_Error_Summary.csv'.\n")
