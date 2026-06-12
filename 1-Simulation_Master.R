# ==============================================================================
# REQUIRED PACKAGES
# ==============================================================================
# Install if not available: install.packages(c("WRS2", "doParallel", "foreach", "onewaytests"))
library(WRS2)       # For Trimmed Welch ANOVA (t1way)
library(doParallel) # Parallel computing backend
library(foreach)    # Parallel loop
library(onewaytests)

Sys.setlocale("LC_ALL", "English")
# ==============================================================================
# 1. DATA GENERATION FUNCTION (Protected with Centering and Scaling)
# ==============================================================================
generate_data <- function(n_vec, var_vec, dist_type = "normal", effect_size = c(0,0,0)) {
  k <- length(n_vec)
  y <- numeric(0)
  group <- factor(rep(1:k, times = n_vec))
  
  for(i in 1:k) {
    if(dist_type == "normal") {
      y_group <- rnorm(n_vec[i], mean = effect_size[i], sd = sqrt(var_vec[i])) 
      
    } else if(dist_type == "lognormal") {
      raw_val <- rlnorm(n_vec[i], meanlog = 0, sdlog = 1)
      standardized_val <- (raw_val - exp(0.5)) / sqrt(exp(1)*(exp(1)-1))
      y_group <- (standardized_val * sqrt(var_vec[i])) + effect_size[i]
      
    } else if(dist_type == "exponential") {
      raw_val <- rexp(n_vec[i], rate = 1)
      centered_val <- raw_val - 1
      y_group <- (centered_val * sqrt(var_vec[i])) + effect_size[i]
    }
    y <- c(y, y_group)
  }
  return(data.frame(y = y, group = group))
}

# ==============================================================================
# 2. ACCELERATED TEST EXECUTION FUNCTION
# ==============================================================================
run_tests <- function(data) {
  p_vals <- numeric(5)
  
  # 1. Classical ANOVA
  fit_aov <- aov(y ~ group, data = data)
  p_vals[1] <- summary(fit_aov)[[1]][["Pr(>F)"]][1]
  
  # 2. Welch ANOVA
  p_vals[2] <- oneway.test(y ~ group, data = data, var.equal = FALSE)$p.value
  
  # 3. Brown-Forsythe ANOVA
  # Uses the bf.test function from the onewaytests package
  bf_out <- onewaytests::bf.test(y ~ group, data = data, verbose = FALSE)
  p_vals[3] <- bf_out$p.value
  
  # 4. Kruskal-Wallis
  p_vals[4] <- kruskal.test(y ~ group, data = data)$p.value
  
  # 5. Trimmed Welch ANOVA (20% trimmed)
  p_vals[5] <- t1way(y ~ group, data = data, tr = 0.2)$p.value
  
  names(p_vals) <- c("ANOVA", "Welch", "BF", "KW", "Trimmed")
  return(p_vals)
}

# ==============================================================================
# 3. SCENARIO MATRIX SETUP
# ==============================================================================
n_scenarios <- list(c(10, 10, 10), c(30, 30, 30), c(5, 10, 20))
var_scenarios <- list(c(1, 1, 1), c(1, 2, 4), c(4, 2, 1))
dist_scenarios <- c("normal", "lognormal", "exponential")

scenarios <- expand.grid(
  n_id = 1:length(n_scenarios),
  var_id = 1:length(var_scenarios),
  dist = dist_scenarios,
  stringsAsFactors = FALSE
)

# ==============================================================================
# 4. DUAL-ENGINE PARALLEL MONTE CARLO LOOP
# ==============================================================================
cores <- parallel::detectCores() - 1 
cl <- makeCluster(cores)
registerDoParallel(cl)

iterations <- 10000     # Ideal 10,000 iterations for the manuscript
alpha <- 0.05           # Nominal significance level
shift_value <- 1.2      # Effect size added to the 3rd group for Statistical Power

cat("Starting simulation. Total scenarios:", nrow(scenarios), "\n")

# foreach loop to combine dual results into a single data.frame
results_combined <- foreach(i = 1:nrow(scenarios), .combine = rbind, .packages = c("WRS2")) %dopar% {
  
  current_n <- n_scenarios[[scenarios$n_id[i]]]
  current_var <- var_scenarios[[scenarios$var_id[i]]]
  current_dist <- scenarios$dist[i]
  
  # Counters
  rej_H0 <- numeric(5) # Number of rejections for Type I Error
  rej_H1 <- numeric(5) # Number of rejections for Statistical Power
  
  for(sim in 1:iterations) {
    # -------------------------------------------------------------------
    # PHASE 1: TYPE I ERROR (No true difference)
    sim_data_H0 <- generate_data(current_n, current_var, current_dist, effect_size = c(0,0,0))
    p_vals_H0 <- run_tests(sim_data_H0)
    rej_H0 <- rej_H0 + (p_vals_H0 < alpha)
    
    # -------------------------------------------------------------------
    # PHASE 2: STATISTICAL POWER
    sim_data_H1 <- sim_data_H0
    # Shifting the mean of the 3rd group by 'shift_value' (True difference exists!)
    sim_data_H1$y[sim_data_H1$group == 3] <- sim_data_H1$y[sim_data_H1$group == 3] + shift_value
    p_vals_H1 <- run_tests(sim_data_H1)
    rej_H1 <- rej_H1 + (p_vals_H1 < alpha)
  }
  
  # Calculate empirical rates
  type1_rates <- rej_H0 / iterations
  power_rates <- rej_H1 / iterations
  
  # Return both T1 and Power results in a single row
  data.frame(
    N_Design = paste(current_n, collapse = "-"),
    Var_Design = paste(current_var, collapse = "-"),
    Distribution = current_dist,
    
    # Type I Error Results
    T1_ANOVA = type1_rates[1], T1_Welch = type1_rates[2], T1_BF = type1_rates[3], 
    T1_KW = type1_rates[4], T1_Trim = type1_rates[5],
    
    # Statistical Power Results
    PWR_ANOVA = power_rates[1], PWR_Welch = power_rates[2], PWR_BF = power_rates[3], 
    PWR_KW = power_rates[4], PWR_Trim = power_rates[5]
  )
}

stopCluster(cl)

# ==============================================================================
# 5. SEPARATING AND SAVING TABLES
# ==============================================================================
# Split the results into two separate tables in accordance with the manuscript format
results_Type1 <- results_combined[, c(1:3, 4:8)]
results_Power <- results_combined[, c(1:3, 9:13)]

# Print to console
cat("\n--- TYPE I ERROR TABLE ---\n")
print(head(results_Type1))

cat("\n--- STATISTICAL POWER TABLE ---\n")
print(head(results_Power))

# Save as CSV
write.csv(results_Type1, file = "Table_S1_Type1_Error.csv", row.names = FALSE)
write.csv(results_Power, file = "Table_S2_Statistical_Power.csv", row.names = FALSE)
write.csv(results_combined, file = "Simulation_Master_Results.csv", row.names = FALSE)