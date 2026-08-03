# Avellaneda & Lee (2008) Replication

## Overview
Two linked notebooks on the EURO STOXX 50 universe. The first studies covariance estimation and statistical factor models in their own right; the second builds on that PCA machinery to implement and stress-test the statistical arbitrage strategy from Avellaneda & Lee's "Statistical Arbitrage in the U.S. Equities Market" (2008).

## Methods
- Ledoit-Wolf shrinkage estimators (constant-correlation and single-factor targets)
- PCA-based factor models and eigenvalue-spectrum denoising
- Ornstein-Uhlenbeck calibration of idiosyncratic residuals
- Modified s-score trading signal with an explicit signal-to-execution lag
- Rolling, out-of-sample backtests with transaction costs

## Notebooks
1. **`01_covariance_estimation_and_factor_models.ipynb`** — sample covariance vs. shrinkage estimators, evaluated through out-of-sample minimum-variance and mean-variance portfolios (turnover, exposure, cost impact); PCA factor stability and GICS sector interpretation; shrinkage vs. PCA denoising on the same eigenvalue spectrum.
2. **`02_statistical_arbitrage.ipynb`** — the full statistical arbitrage pipeline: PCA factor model, O-U calibration, s-score construction, rolling backtest with transaction costs, sensitivity to the entry threshold and to the number of factors, a pre/post-2020 regime split, and a volume-adjusted "trading time" variant.

## Results
- Out-of-sample comparison of sample, constant-correlation, and single-factor covariance estimators across turnover, gross exposure, condition number, and net-of-cost Sharpe ratio.
- Stability of PCA factor loadings over time (cosine similarity) and their interpretation against GICS sectors.
- Statistical arbitrage strategy performance, gross and net of transaction costs, benchmarked against a vol-equalized equal-weight portfolio.
- Sensitivity of the strategy to the s-score entry threshold and to the number of PCA factors, including a pre/post-2020 regime comparison.

## How to run
Open either notebook and run all cells; both expect the repository root (this folder) as the working directory, since data is read from `data/` via relative paths. Requires `numpy`, `pandas`, and `matplotlib`.

## Structure
- `data/`: EURO STOXX 50 prices, volumes, and ticker/sector reference data
- `utilities/`: shared covariance, PCA, shrinkage, portfolio optimization, statistical arbitrage, and backtest modules
- `01_covariance_estimation_and_factor_models.ipynb`, `02_statistical_arbitrage.ipynb`: the two notebooks described above
