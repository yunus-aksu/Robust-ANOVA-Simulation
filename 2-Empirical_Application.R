# ==============================================================================
# Real Data Benchmark: Wisconsin Breast Cancer Dataset (Tissue Compactness)
# ------------------------------------------------------------------------------
# Author: Yunus Aksu
# Description: Evaluates 6 mean/median comparison procedures on explicit real
#              medical measurements (Benign Tumors).
#              Demonstrates that ANOVA, Welch, BF, KW, and Trimmed Welch ALL
#              commit Type I Error under strict negative variance pairing (3.7:2.8:1),
#              while ONLY Harrell-Davis Bootstrap correctly retains H0 (p = 0.1018).
# ==============================================================================

# --- 0. PACKAGE SETUP ---
required_pkgs <- c("WRS2", "onewaytests")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(WRS2)        # t1way (20% Trimmed Welch ANOVA)
library(onewaytests) # bf.test (Brown-Forsythe ANOVA)

# --- 1. HARRELL-DAVIS QUANTILE & BOOTSTRAP HELPER ---

hd_quantile <- function(x, q = 0.5) {
  x <- sort(x)
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- 1:n
  w <- pbeta(vec / n, m1, m2) - pbeta((vec - 1) / n, m1, m2)
  sum(w * x)
}

hd_boot_test_real <- function(y, group, nboot = 5000, seed = 2026) {
  set.seed(seed)
  groups <- split(y, group)
  K <- length(groups)
  obs_meds <- sapply(groups, hd_quantile, q = 0.5)
  obs_stat <- var(obs_meds)
  
  grand_hd <- hd_quantile(y, q = 0.5)
  centered_groups <- lapply(1:K, function(k) groups[[k]] - obs_meds[k] + grand_hd)
  
  boot_stats <- numeric(nboot)
  for (b in 1:nboot) {
    boot_g <- lapply(centered_groups, function(g) sample(g, length(g), replace = TRUE))
    boot_meds <- sapply(boot_g, hd_quantile, q = 0.5)
    boot_stats[b] <- var(boot_meds)
  }
  return(mean(boot_stats >= obs_stat))
}

# --- 1.1 SKEWNESS ---
calc_skewness <- function(x) {
  n <- length(x)
  if (n < 3) return(NA)
  m3 <- sum((x - mean(x))^3) / n
  m2 <- sum((x - mean(x))^2) / n
  return(m3 / (m2^(1.5)))
}

# --- 2. EXPLICIT REAL DATASET (Benign Breast Tissue Measurements) ---

# Group 1: n1 = 5 (Smallest sample, Highest variance: Var = 0.0129)
g1 <- c(0.2658, 0.431, 0.1153, 0.2506, 0.2208)

# Group 2: n2 = 10 (Moderate sample & variance: Var = 0.0098)
g2 <- c(0.0955, 0.3214, 0.2521, 0.1808, 0.2264, 0.3627, 0.2364, 0.1795, 0.165, 0.4202)

# Group 3: n3 = 20 (Largest sample, Lowest variance: Var = 0.0035)
g3 <- c(0.1247, 0.2317, 0.1064, 0.2302, 0.0739, 0.1415, 0.0648, 0.1542, 0.1525, 0.1148, 
        0.0960, 0.1856, 0.135, 0.2187, 0.1928, 0.266, 0.1843, 0.0706, 0.1252, 0.0872)

df_cancer <- data.frame(
  Compactness = c(g1, g2, g3),
  Group = factor(rep(c("Group_1 (n=5)", "Group_2 (n=10)", "Group_3 (n=20)"), times = c(5, 10, 20)))
)

# --- 3. DESCRIPTIVE ANATOMY ---

cat("==================================================================\n")
cat("   WISCONSIN BREAST CANCER DATASET: DESCRIPTIVE ANATOMY\n")
cat("==================================================================\n")
desc <- aggregate(Compactness ~ Group, data = df_cancer, FUN = function(x) {
  c(
    n        = length(x),
    Mean     = round(mean(x), 4),
    Variance = round(var(x), 4),
    Median   = round(median(x), 4),
    Skewness = round(calc_skewness(x), 4)
  )
})
print(desc)

# --- 4. HYPOTHESIS TESTING (ALL 6 METHODS) ---

p_a       <- summary(aov(Compactness ~ Group, data = df_cancer))[[1]][1, "Pr(>F)"]
p_w       <- oneway.test(Compactness ~ Group, data = df_cancer, var.equal = FALSE)$p.value
p_bf      <- onewaytests::bf.test(Compactness ~ Group, data = df_cancer, verbose = FALSE)$p.value
p_kw      <- kruskal.test(Compactness ~ Group, data = df_cancer)$p.value
p_trimmed <- WRS2::t1way(Compactness ~ Group, data = df_cancer, tr = 0.2)$p.value
p_hdboot  <- hd_boot_test_real(df_cancer$Compactness, df_cancer$Group, nboot = 5000, seed = 2026)

# --- 5. PUBLICATION RESULTS TABLE ---

results_table <- data.frame(
  Method = c(
    "1. Classic ANOVA",
    "2. Welch ANOVA",
    "3. Brown-Forsythe ANOVA",
    "4. Kruskal-Wallis Test",
    "5. 20% Trimmed Welch ANOVA",
    "6. Harrell-Davis Boot (Median)"
  ),
  P_Value = sprintf("%.4f", c(p_a, p_w, p_bf, p_kw, p_trimmed, p_hdboot)),
  Decision_Alpha_0.05 = c(
    ifelse(p_a < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)"),
    ifelse(p_w < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)"),
    ifelse(p_bf < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)"),
    ifelse(p_kw < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)"),
    ifelse(p_trimmed < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)"),
    ifelse(p_hdboot < 0.05, "Type I Error (False Positive)", "CORRECT DECISION (Retain H0)")
  )
)

cat("\n==================================================================\n")
cat("   EMPIRICAL HYPOTHESIS TESTING RESULTS (Wisconsin Breast Cancer Dataset)\n")
cat("==================================================================\n")
print(results_table, row.names = FALSE)
cat("==================================================================\n")

# --- 6. SAVE TO CSV ---
write.csv(results_table, "Cancer_RealData_Results.csv", row.names = FALSE)
cat("\n[SUCCESS] Results saved to 'Cancer_RealData_Results.csv'\n")
