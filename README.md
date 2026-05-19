# Extreme Risk Modeling with Copulas

Quantitative risk modeling project using Extreme Value Theory (EVT), copulas and Monte Carlo simulation in R.

## Project Overview

This project focuses on modeling extreme operational losses and dependency structures between two portfolios using advanced statistical and actuarial techniques.

The study includes:
- Extreme Value Theory (EVT)
- GEV and GPD distributions
- Threshold selection diagnostics
- Copula modeling
- Tail dependence analysis
- Monte Carlo simulation

---

## Main Objectives

- Model extreme losses for two risk series
- Estimate 99.5% quantiles
- Compute 200-year return losses
- Capture tail dependence between portfolios
- Simulate joint extreme scenarios

---

## Methodologies

### Marginal Modeling
- Normal distribution benchmark
- GEV (Generalized Extreme Value)
- GPD (Generalized Pareto Distribution)
- Peaks Over Threshold (POT)

### Dependence Modeling
- Gaussian Copula
- Clayton Copula
- Frank Copula
- Gumbel Copula

### Model Selection
- AIC/BIC comparison
- Mean Excess Plot
- Stability diagnostics
- Goodness-of-fit testing

---

## Final Model

### Selected Dependence Structure
- Gumbel Copula

### Selected Tail Model
- GPD (Generalized Pareto Distribution)

The final model captures:
- heavy-tail behavior,
- upper tail dependence,
- joint extreme loss scenarios.

---

## Technologies Used

- R
- ismev
- evd
- copula
- MASS
- Monte Carlo Simulation

---

## Key Features

- EVT calibration
- Tail risk estimation
- Threshold optimization
- Copula calibration
- Monte Carlo simulation
- Extreme dependence modeling
- Diagnostic visualizations

---

## Results

- 99.5% extreme quantile estimation
- 200-year return level estimation
- Tail dependence analysis
- Joint portfolio loss simulation

---

## Project Structure

```text
.
├── projet_grands_risques.R
├── Rapport_GRVE.pdf
└── README.md
```

---

## Authors

CY Tech – Applied Mathematics & Finance
