R Code for: Type I Error and Size-Adjusted Power Behavior of One-Way ANOVA Alternatives Under Heteroscedasticity and Extreme Skewness: A Monte Carlo Study
Overview  
This repository contains the R scripts, datasets, and high-resolution visual outputs for the paper "Type I Error and Size-Adjusted Power Behavior of One-Way ANOVA Alternatives Under Heteroscedasticity and Extreme Skewness: A Monte Carlo Study".

The provided code allows researchers to fully replicate the 10,000-iteration Monte Carlo simulation (spanning 135 experimental conditions), reproduce the 4 publication-ready figures (300 DPI), and run the empirical application on the Wisconsin Breast Cancer (WDBC) dataset.

Evaluated Statistical Procedures  
The study benchmarks six location procedures across K = 3 independent groups under non-normality, heteroscedasticity, and unbalanced designs:

Classic One-Way ANOVA (F-test)  
Welch ANOVA  
Brown-Forsythe ANOVA (F^*)  
Kruskal-Wallis Test  
20% Trimmed Welch ANOVA  
Harrell-Davis Quantile Bootstrap Test (HD_Boot)  
Repository Structure and Files
--> 1-Simulation_Master.R          # Core parallel Monte Carlo simulation script (10,000 iterations)
--> 2-Empirical_Application.R         # Empirical validation on Wisconsin Breast Cancer dataset (WDBC)
--> 3-Data_Visualization.R             # Generates all 4 publication-ready figures (300 DPI)
--> S1_Type1_Error_Summary.csv   # Compiled Type I error rates across simulation scenarios
--> S2_Master_Results.csv        # Master dataset (Type I error, Raw power, Size-Adjusted power)
--> Cancer_RealData_Desc.csv     # Descriptive anatomy of clinical tissue compactness groups
--> Cancer_RealData_Results.csv   # Hypothesis decisions and p-values on clinical data
--> Figure1_Type1_Error_Heatmap.png   # Heatmap of Type I error control (Bradley's limits)
--> Figure2_3x3_Size_Power_Matrix.png # Size-power trade-off matrix
--> Figure3_Power_Comparison.png      # Spurious vs. Size-Adjusted power comparison
--> Figure4_Decision_Matrix.png       # Practical evidence-based decision matrix
Prerequisites
To run the replication scripts, ensure you have R (>= 4.0.0) installed along with the following packages:

# Required Packages
install.packages(c("ggplot2", "dplyr", "tidyr", "scales", "ggrepel",
                   "WRS2", "onewaytests", "doParallel", "foreach"))
How to Run the Replication
Clone the Repository:

git clone https://github.com/yunus-aksu/Robust-ANOVA-Simulation.git
cd Robust-ANOVA-Simulation
Run Monte Carlo Simulation: Open Simulation_Master.R in RStudio or run via terminal:

source("1-Simulation_Master.R")
Note: Automatically uses parallel processing cores to accelerate 10,000 iterations.

Generate Publication Figures: Run Master_Figures.R to reproduce Figures 1–4 at 300 DPI:

source("3-Data_Visualization.R")
Run Real Data Application: Run RealData_Benchmark.R to reproduce the clinical application on benign tissue compactness:

source("2-Empirical_Application.R")
Author & Contact
Yunus Aksu (Corresponding Author)
Email: yunusaksu011@gmail.com
ORCID: 0009-0003-0437-0121
