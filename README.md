# R Code for: On Robust One-Way Analysis of Variance: A Simulation Study

## Overview
This repository contains the R scripts and datasets used in the methodology article, *"On Robust One-Way Analysis of Variance: A Simulation Study"*. The provided codes allow researchers to fully replicate the 10,000-iteration Monte Carlo simulation, the 3x3 matrix visualizations, and the empirical application assessing the robustness and statistical power of five distinct mean-comparison tests (Classical ANOVA, Welch ANOVA, Brown-Forsythe, Kruskal-Wallis, and 20% Trimmed Welch ANOVA) under pathological data conditions (non-normality, unbalanced designs, and severe variance heterogeneity).

## Repository Structure and Files

### 1. R Scripts
* `Simulation_Master.R`: The core script that executes the dual-engine parallel Monte Carlo simulation. It calculates empirical Type I error rates and computes statistical power across 27 combined scenarios based on Bradley's liberal robustness criterion.
* `Empirical_Application.R`: An algorithmic script that searches and extracts a specific sub-sample from the historical Michelson-Morley speed of light dataset (`morley`). It empirically demonstrates the pathological breakdown of Classical ANOVA under negative variance pairing (4:2:1 ratio) in real-world data.
* `Figure_1_Plot.R`: The visualization script utilizing `ggplot2` and `ggrepel` to generate the publication-ready 3x3 matrix plot (Figure 1), illustrating the Size-Power trade-off across different variance structures.

### 2. Output Data (Supplementary Materials)
* `Table_S1_Type1_Error.csv`: The complete compiled results for empirical Type I error rates across all 27 simulated scenarios.
* `Table_S2_Statistical_Power.csv`: The complete compiled results for statistical power.
* `Figure_1.tiff`: The high-resolution (600 DPI) output of the trade-off plot.

## Prerequisites
To run the scripts seamlessly, ensure you have R installed along with the following packages:
- **Statistical Analysis:** `WRS2` (For Trimmed Welch ANOVA), `onewaytests` (For Brown-Forsythe test)
- **Parallel Computing:** `doParallel`, `foreach`
- **Visualization:** `ggplot2`, `ggrepel`

## How to Run the Replication
1. Clone this repository to your local machine or download the `.R` files directly.
2. Open the scripts in RStudio.
3. Install the required packages if you haven't already using `install.packages()`.
4. Source `Simulation_Master.R`. The parallel computing backend will automatically utilize your machine's available cores to accelerate the 10,000 iterations.
5. Run `Figure_1_Plot.R` to reproduce the high-resolution visualizations based on the generated data.
6. Run `Empirical_Application.R` to observe the real-world vulnerability of classical methods using the Michelson-Morley dataset.
