# ==============================================================================
# FIGURE 1: 3x3 MATRIX PLOT FOR SIZE-POWER TRADE-OFF
# ==============================================================================

library(ggplot2)
library(ggrepel)

# Master Data Frame combining all data (with English Titles)
df_multi <- data.frame(
  Test = rep(rep(c("ANOVA", "Welch", "BF", "Kruskal-Wallis", "Trimmed"), 3), 3),
  
  # Distributions (Rows) - Long texts split into two lines with \n
  Distribution = factor(rep(c("1. Normal\nDistribution", 
                              "2. Lognormal\nDistribution", 
                              "3. Exponential\nDistribution"), each = 15)),
  
  # Scenarios (Columns)
  Scenario = factor(rep(rep(c("A) Homogeneous\n(n: 10-10-10, Var: 1-1-1)", 
                              "B) Positive Pairing\n(n: 5-10-20, Var: 1-2-4)", 
                              "C) Negative Pairing\n(n: 5-10-20, Var: 4-2-1)"), each = 5), 3)),
  
  Type1_Error = c(
    # --- NORMAL DISTRIBUTION ---
    0.0558, 0.0538, 0.0535, 0.0492, 0.0554,  # A (Homogeneous)
    0.0156, 0.0494, 0.0512, 0.0229, 0.0583,  # B (Positive Pairing)
    0.1501, 0.0590, 0.0589, 0.0890, 0.0879,  # C (Negative Pairing)
    
    # --- LOGNORMAL DISTRIBUTION ---
    0.0320, 0.0336, 0.0235, 0.0425, 0.0330,  # A (Homogeneous)
    0.0359, 0.0359, 0.0366, 0.0954, 0.0749,  # B (Positive Pairing)
    0.1241, 0.1523, 0.0879, 0.2150, 0.1712,  # C (Negative Pairing)
    
    # --- EXPONENTIAL DISTRIBUTION ---
    0.0456, 0.0472, 0.0384, 0.0485, 0.0476,  # A (Homogeneous)
    0.0271, 0.0416, 0.0445, 0.0523, 0.0565,  # B (Positive Pairing)
    0.1402, 0.1158, 0.0770, 0.1576, 0.1389   # C (Negative Pairing)
  ),
  
  Power = c(
    # --- NORMAL DISTRIBUTION ---
    0.7581, 0.7357, 0.7503, 0.7228, 0.6036,
    0.3342, 0.4046, 0.5141, 0.3476, 0.3081,
    0.6855, 0.4916, 0.3669, 0.6074, 0.3567,
    
    # --- LOGNORMAL DISTRIBUTION ---
    0.8232, 0.9077, 0.8019, 0.9803, 0.9743,
    0.5064, 0.7174, 0.6541, 0.8952, 0.7211,
    0.7834, 0.8112, 0.6010, 0.9509, 0.8775,
    
    # --- EXPONENTIAL DISTRIBUTION ---
    0.7647, 0.8282, 0.7528, 0.8909, 0.8013,
    0.3438, 0.5361, 0.5531, 0.5411, 0.3515,
    0.7201, 0.6607, 0.4762, 0.8275, 0.6385
  )
)

# Plotting the Graph (3x3 Matrix)
p = ggplot(df_multi, aes(x = Type1_Error, y = Power, color = Test, shape = Test)) +
  
  # Bradley's Robustness Interval (Gray Area) - xmax 0.22
  annotate("rect", xmin = 0.025, xmax = 0.075, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "black") +
  geom_vline(xintercept = c(0.025, 0.075), linetype = "dotted", color = "gray30", linewidth = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Points and Labels
  geom_point(size = 3.5, stroke = 1.2, alpha = 0.9) +
  geom_text_repel(aes(label = Test), size = 3, box.padding = 0.5, 
                  point.padding = 0.3, show.legend = FALSE, family = "serif") +
  
  # 3x3 MATRIX STRUCTURE
  facet_grid(Distribution ~ Scenario) +
  
  # Axis Settings
  scale_x_continuous(breaks = seq(0, 0.20, by = 0.05), limits = c(0, 0.22)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.20), limits = c(0.2, 1.0)) +
  
  # Colors and Shapes
  scale_color_manual(values = c("ANOVA" = "#B2182B", 
                                "Kruskal-Wallis" = "#E69F00",
                                "Welch" = "#2166AC", 
                                "BF" = "#1B7837", 
                                "Trimmed" = "#762A83")) +
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
  
  # English Axis Labels
  labs(
    x = "Empirical Type I Error Rate (\u03B1 = 0.05)",
    y = "Statistical Power (1 - \u03B2)"
  ) +
  
  # Academic Theme
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 1),
    # Font size reduced from 10 to 9 to fit perfectly into the boxes
    strip.text = element_text(face = "bold", size = 9, color = "black"), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggsave(filename = "Figure_1.tiff", 
       plot = p, 
       width = 10,          # Width (inches) - Ideal for 3 columns
       height = 8,          # Height (inches) - Ideal for 3 rows
       units = "in",        # Unit of measurement (inches)
       dpi = 600,           # High resolution requested by journals (300 or 600)
       compression = "lzw"  # Compression that reduces file size without degrading image quality
)
