# ==============================================================================
# Master R Script: Publication-Ready Figures for ANOVA Benchmark
# ------------------------------------------------------------------------------
# Author: Yunus Aksu
# Description: Generates ALL 4 publication figures (300 DPI) in a single run.
#
# Output Files:
#   1. Figure1_Type1_Error_Heatmap.png    (Bradley's Type I Error Heatmap)
#   2. Figure2_3x3_Size_Power_Matrix.png  (Size-Power Trade-off Matrix)
#   3. Figure3_Power_Comparison.png       (Spurious vs. Size-Adjusted Power)
#   4. Figure4_Decision_Matrix.png        (Practical Decision Grid)
# ==============================================================================

# --- 0. PACKAGE INSTALLATION & SETUP ---
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "ggrepel")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)

# Set Working Directory & Data Ingestion
master_results <- read.csv("Master_Results.csv", stringsAsFactors = FALSE)
t1_summary     <- read.csv("Type1_Error_Summary.csv", stringsAsFactors = FALSE)

# Method Name Mapping for Publication Display
method_labels <- c(
  "ANOVA"   = "Classic ANOVA",
  "Welch"   = "Welch ANOVA",
  "BF"      = "Brown-Forsythe",
  "KW"      = "Kruskal-Wallis",
  "Trimmed" = "20% Trimmed Welch",
  "HD_Boot" = "Harrell-Davis Boot"
)

master_results$Method_Name <- factor(method_labels[master_results$Method], levels = method_labels)
t1_summary$Method_Name     <- factor(method_labels[t1_summary$Method], levels = method_labels)

cat("======================================================================\n")
cat("[START] Generating All Publication Figures...\n")
cat("======================================================================\n")

# ==============================================================================
# FIGURE 1: TYPE I ERROR CONTROL HEATMAP - SMOOTH GRADIENT & CLEAN LEGEND
# ==============================================================================

library(ggplot2)
library(dplyr)
library(scales)

# 1. Data Formatting & Label Adjustments
t1_formatted <- t1_summary %>%
  mutate(
    Var_Label = gsub("-", ":", Var_Design),
    Scenario  = paste0("n: ", N_Design, "\nσ²: ", Var_Label),
    
    Distribution_Label = factor(
      case_when(
        Distribution == "normal"      ~ "Normal Distribution",
        Distribution == "exponential" ~ "Exponential Distribution",
        Distribution == "lognormal"   ~ "Lognormal Distribution"
      ),
      levels = c("Normal Distribution", "Exponential Distribution", "Lognormal Distribution")
    )
  )

fig1 <- ggplot(t1_formatted, aes(x = Method_Name, y = Scenario, fill = Type1_Error)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.4f", Type1_Error)), size = 3.0, color = "black", fontface = "bold") +
  
  scale_fill_gradientn(
    colors = c("#2c7bb6", "#2ca25f", "#fdae61", "#d7191c"), 
    values = rescale(c(0.015, 0.050, 0.085, 0.150)),      
    limits = c(0.015, 0.150),
    oob    = squish, 
    breaks = c(0.025, 0.050, 0.075, 0.100, 0.150),
    labels = c("< 0.025", "0.050", "0.075", "0.100", "≥ 0.150"), 
    name   = "Type I Error Rate\n(Nominal α = 0.05)"
  ) +
  facet_wrap(~ Distribution_Label, ncol = 3) +
  labs(
    title   = "Empirical Type I Error Rates Across Experimental Scenarios",
    x       = NULL,
    y       = "Design Scenario",
    caption = "Note: n denotes group sample sizes (n₁, n₂, n₃); σ² denotes group variance ratios (σ₁², σ₂², σ₃²).\nBradley's liberal robust interval is [0.025, 0.075]. Blue: Conservative (< 0.025); Green: Robust (0.025–0.075); Red: Severe Inflation (≥ 0.10)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, vjust = 1, face = "bold", color = "black"),
    axis.text.y      = element_text(size = 9, face = "bold", color = "black", lineheight = 0.9),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5), 
    plot.caption     = element_text(size = 9, color = "grey20", hjust = 0, margin = margin(t = 10)),
    legend.position  = "bottom",
    legend.key.width = unit(2.5, "cm")
  )

ggsave("Figure1_Type1_Error_Heatmap.png", plot = fig1, width = 13, height = 7.5, dpi = 300)

cat("[✔ 1/4] Figure 1 Exported: 'Figure1_Type1_Error_Heatmap.png'\n")
# ==============================================================================
# FIGURE 2: 3x3 MATRIX PLOT FOR SIZE-POWER TRADE-OFF
# ==============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. Data Preparation
df_fig2 <- master_results %>% filter(Effect_Size_d == 1.2)

scenarios_map <- c(
  "10-10-10_1-1-1" = "Homogeneous\n(n: 10-10-10, σ²: 1:1:1)",
  "5-10-20_1-2-4"  = "Positive Pairing\n(n: 5-10-20, σ²: 1:2:4)",
  "5-10-20_4-2-1"  = "Negative Pairing\n(n: 5-10-20, σ²: 4:2:1)"
)

df_fig2$Scen_Key <- paste0(df_fig2$N_Design, "_", df_fig2$Var_Design)
df_fig2 <- df_fig2 %>% filter(Scen_Key %in% names(scenarios_map))
df_fig2$Scenario <- factor(scenarios_map[df_fig2$Scen_Key], levels = scenarios_map)

dist_map_fig2 <- c(
  "normal"      = "Normal\nDistribution",
  "lognormal"   = "Lognormal\nDistribution",
  "exponential" = "Exponential\nDistribution"
)
df_fig2$Distribution_Label <- factor(dist_map_fig2[df_fig2$Distribution], levels = dist_map_fig2)

short_method_map <- c(
  "ANOVA"   = "ANOVA",
  "Welch"   = "Welch",
  "BF"      = "Brown-Forsythe",
  "KW"      = "Kruskal-Wallis",
  "Trimmed" = "Trimmed Welch",
  "HD_Boot" = "Harrell-Davis Boot"
)
df_fig2$Test <- factor(short_method_map[df_fig2$Method], levels = short_method_map)

fig2 <- ggplot(df_fig2, aes(x = Type1_Error, y = Size_Adjusted_Power, color = Test, shape = Test)) +
  # Bradley's Limit (0.025 - 0.075)
  annotate("rect", xmin = 0.025, xmax = 0.075, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "black") +
  geom_vline(xintercept = c(0.025, 0.075), linetype = "dotted", color = "gray30", linewidth = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_point(size = 3.5, stroke = 1.2, alpha = 0.9) +
  
  geom_text_repel(
    aes(label = Test), 
    size = 3, 
    box.padding = 0.5, 
    point.padding = 0.3, 
    max.overlaps = Inf, 
    show.legend = FALSE
  ) +
  
  facet_grid(Distribution_Label ~ Scenario) +
  
  scale_x_continuous(breaks = seq(0, 0.20, by = 0.05), limits = c(0, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.20), limits = c(0.1, 1.05)) +
  
  scale_color_manual(values = c(
    "ANOVA"              = "#B2182B",
    "Welch"              = "#2166AC",
    "Brown-Forsythe"     = "#1B7837",
    "Kruskal-Wallis"     = "#E69F00",
    "Trimmed Welch"      = "#762A83",
    "Harrell-Davis Boot" = "#008080"
  )) +
  scale_shape_manual(values = c(15, 17, 18, 16, 8, 3)) +
  
  labs(
    title   = "Type I Error and Size-Adjusted Power Performance Across Scenarios",
    x       = "Empirical Type I Error Rate (Nominal α = 0.05)",
    y       = "Size-Adjusted Statistical Power (1 - β)",
    caption = "Note: n denotes group sample sizes (n₁, n₂, n₃); σ² denotes group variance ratios (σ₁², σ₂², σ₃²).\nThe shaded region represents Bradley's liberal robust interval [0.025, 0.075] around nominal α = 0.05."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5), # Ortalanmış başlık
    plot.caption    = element_text(size = 9, color = "grey20", hjust = 0, margin = margin(t = 10)),
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

ggsave("Figure2_3x3_Size_Power_Matrix.png", plot = fig2, width = 13, height = 10, dpi = 300)

cat("[✔ 2/4] Figure 2 Exported: 'Figure2_3x3_Size_Power_Matrix.png'\n")


# ==============================================================================
# FIGURE 3: RAW VS. SIZE-ADJUSTED POWER GRID
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 1. Data Preparation
pwr_all <- master_results %>%
  group_by(Distribution, Method_Name, Effect_Size_d) %>%
  summarise(
    Mean_Raw_Power      = mean(Raw_Power),
    Mean_Size_Adj_Power = mean(Size_Adjusted_Power),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = c(Mean_Raw_Power, Mean_Size_Adj_Power),
    names_to  = "Power_Type",
    values_to = "Power_Value"
  ) %>%
  mutate(
    Power_Type_Label = ifelse(Power_Type == "Mean_Raw_Power", "Raw Power (Unadjusted)", "Size-Adjusted Power (Honest)"),
    Effect_Label     = factor(paste0("Cohen's d = ", Effect_Size_d), levels = c("Cohen's d = 0.2", "Cohen's d = 0.5", "Cohen's d = 0.8", "Cohen's d = 1.2")),
    Distribution_Label = factor(
      case_when(
        Distribution == "normal"      ~ "Normal Distribution",
        Distribution == "exponential" ~ "Exponential Distribution",
        Distribution == "lognormal"   ~ "Lognormal Distribution"
      ),
      levels = c("Normal Distribution", "Exponential Distribution", "Lognormal Distribution")
    )
  )

fig3 <- ggplot(pwr_all, aes(x = Method_Name, y = Power_Value, fill = Power_Type_Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_text(
    aes(label = sprintf("%.2f", Power_Value)),
    position = position_dodge(width = 0.8),
    vjust    = -0.4,
    size     = 3.2,
    fontface = "bold"
  ) +
  
  facet_grid(Effect_Label ~ Distribution_Label) +
  
  scale_fill_manual(values = c(
    "Raw Power (Unadjusted)"       = "#D95F02", 
    "Size-Adjusted Power (Honest)" = "#2B5C8F"  
  )) +
  
  scale_y_continuous(
    limits = c(0, 1.00), 
    breaks = seq(0, 1, by = 0.25), 
    labels = percent_format(),
    expand = c(0, 0) 
  ) +
  
  coord_cartesian(ylim = c(0, 1.05), clip = "off") +
  
  labs(
    title   = "Empirical Power Comparison: Raw vs. Size-Adjusted Power Across Distributions",
    x       = NULL,
    y       = "Empirical Statistical Power",
    fill    = "Power Metric:",
    caption = "Note: Size-Adjusted Power calibrates test statistics to empirical 95th percentile critical thresholds under H0.\nThis reveals the artificial/spurious power advantages of Kruskal-Wallis and Trimmed Welch tests under skewed distributions."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, face = "bold", color = "black", size = 8.5),
    strip.text      = element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 9.5),
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5), # Ortalanmış başlık
    plot.caption    = element_text(size = 9, color = "grey20", hjust = 0, margin = margin(t = 12)),
    panel.spacing.y = unit(0.8, "lines") # Satırlar arası dikey mesafeyi ferahlatma
  )

ggsave("Figure3_Power_Comparison.png", plot = fig3, width = 13.5, height = 12.5, dpi = 300)

cat("[✔ 3/4] Figure 3 Exported: 'Figure3_Power_Comparison.png'\n")

# ==============================================================================
# FIGURE 4: EVIDENCE-BASED DECISION MATRIX
# ==============================================================================

library(ggplot2)
library(dplyr)


decision_grid <- expand.grid(
  Sample_Design   = c("Equal\n(n: 10-10-10, 30-30-30)", "Unequal\n(n: 5-10-20)"),
  Variance_Design = c("Homogeneous\n(σ²: 1:1:1)", "Positive Pairing\n(σ²: 1:2:4)", "Negative Pairing\n(σ²: 4:2:1)"),
  Skewness        = c("Normal Distribution", "Exponential Distribution", "Lognormal Distribution"),
  stringsAsFactors = FALSE
)

decision_grid <- decision_grid %>%
  mutate(
    Recommended_Method = case_when(
      Skewness == "Normal Distribution" & Variance_Design == "Homogeneous\n(σ²: 1:1:1)" ~ "Classic ANOVA",
      Skewness == "Normal Distribution" & Variance_Design != "Homogeneous\n(σ²: 1:1:1)" ~ "Welch ANOVA",
      
      Skewness == "Exponential Distribution" & Variance_Design == "Homogeneous\n(σ²: 1:1:1)" ~ "Welch / Kruskal-Wallis",
      Skewness == "Exponential Distribution" & Variance_Design != "Homogeneous\n(σ²: 1:1:1)" ~ "Brown-Forsythe",
      
      Skewness == "Lognormal Distribution" & Sample_Design == "Equal\n(n: 10-10-10, 30-30-30)" ~ "Harrell-Davis Boot",
      Skewness == "Lognormal Distribution" & Sample_Design == "Unequal\n(n: 5-10-20)" & Variance_Design == "Negative Pairing\n(σ²: 4:2:1)" ~ "Harrell-Davis Boot",
      TRUE ~ "Brown-Forsythe /\nHarrell-Davis Boot"
    ),
    
    Skewness = factor(Skewness, levels = c("Normal Distribution", "Exponential Distribution", "Lognormal Distribution")),
    Variance_Design = factor(Variance_Design, levels = c("Homogeneous\n(σ²: 1:1:1)", "Positive Pairing\n(σ²: 1:2:4)", "Negative Pairing\n(σ²: 4:2:1)")),
    Sample_Design = factor(Sample_Design, levels = c("Equal\n(n: 10-10-10, 30-30-30)", "Unequal\n(n: 5-10-20)"))
  )

fig4 <- ggplot(decision_grid, aes(x = Sample_Design, y = Variance_Design, fill = Recommended_Method)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = Recommended_Method), fontface = "bold", size = 3.3, color = "black") +
  
  facet_wrap(~ Skewness, ncol = 3) +
  
  scale_fill_manual(
    values = c(
      "Classic ANOVA"            = "#A6CEE3", 
      "Welch ANOVA"              = "#1F78B4", 
      "Brown-Forsythe"           = "#B2DF8A", 
      "Welch / Kruskal-Wallis"   = "#FDBF6F", 
      "Harrell-Davis Boot"       = "#E31A1C", 
      "Brown-Forsythe /\nHarrell-Davis Boot" = "#1B9E77"  
    ),
    name = "Recommended Procedure:"
  ) +
  
  labs(
    title   = "Evidence-Based Decision Matrix for One-Way Location Procedures",
    x       = "Sample Size Allocation",
    y       = "Variance Structure",
    caption = "Note: n denotes group sample sizes (n₁, n₂, n₃); σ² denotes group variance ratios (σ₁², σ₂², σ₃²).\nDecision guidelines strictly prioritize Type I error preservation within Bradley's limits [0.025, 0.075], followed by Size-Adjusted Power optimization."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(face = "bold", size = 9.5, color = "black", lineheight = 0.9),
    axis.text.y     = element_text(face = "bold", size = 9.5, color = "black", lineheight = 0.9),
    axis.title      = element_text(face = "bold", size = 11),
    strip.text      = element_text(face = "bold", size = 11),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption    = element_text(size = 9, color = "grey20", hjust = 0, margin = margin(t = 12)),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 10)
  )

ggsave("Figure4_Decision_Matrix.png", plot = fig4, width = 13.5, height = 7.0, dpi = 300)

cat("[✔ 4/4] Figure 4 Exported: 'Figure4_Decision_Matrix.png'\n")

cat("======================================================================\n")
cat("[COMPLETE] ALL 4 FIGURES SUCCESSFULLY GENERATED AND SAVED AT 300 DPI!\n")
cat("======================================================================\n")
