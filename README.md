# Partial Network Regression

This repository contains R and C++ code implementing the **Partial Network Regressive Model** for network reconstruction and community detection, based on the methodology described in:

> **Partial Network Regressive model: Individual-centered Network Reconstruction and Community Detection based on Signal Lasso**

## Overview

The method combines signal recovery with network structure detection, particularly useful for:
- Reconstructing unknown network connections from individual-centered Network
- Community detection in network structures

## Folder Structure

```
github/
├── SL1.R                 # Signal Lasso main function
├── SMC1.R                # Structure-based Momentum Clustering
├── RRCN.R               # Robust Regularized Composite Network
├── beta_OLS.R            # OLS estimation for regression coefficients
├── Y_theta.R            # Generate response variable
├── get_A_ER.R           # Generate Erdős-Rényi random graph
├── get_A_ER_dense.R      # Generate dense Erdős-Rényi random graph
├── get_y_lin.R          # Generate linear response variable
├── generate_Z_b_A.R      # Generate latent variables and coefficients
├── main.cpp            # C++ implementations for performance
├── main.R            # Main execution functions
├── simulation/        # Simulation studies
│   ├── simulation_ER.R
│   ├── simulation_ER_dense.R
│   ├── simulation_ER_sensitivity.R
│   ├── simulation_group.R
│   └── simulation_SBM.R
└── real data/        # Real-world data analysis
    ├── SSE180/           # Shanghai Stock Exchange 180
    └── urban crime/        # Urban crime data
```

## Simulation

The `simulation/` folder contains various simulation studies to validate the method:

- **simulation_ER.R**: Sparse Erdős-Rényi random graph model simulation
- **simulation_ER_dense.R**: Dense Erdős-Rényi random graph model simulation
- **simulation_ER_sensitivity.R**: Sensitivity analysis
- **simulation_group.R**: Group-based simulation
- **simulation_SBM.R**: Stochastic block model simulation

## Real Data Analysis

The `real data/` folder contains two real-world applications:

### 1. Shanghai Stock Exchange 180 (SSE180)
Analysis of stock market data for 180 listed companies, reconstructing the stock network based on shareholder relationships.

### 2. Urban Crime
Analysis of urban crime data across 138 census tracts, detecting communities in crime patterns.

## Reference

If you use this code in your research, please cite:

> Partial Network Regressive model: Individual-centered Network Reconstruction and Community Detection based on Signal Lasso
