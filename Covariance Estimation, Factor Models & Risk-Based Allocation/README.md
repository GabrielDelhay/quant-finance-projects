# Risk-Based Allocation Strategies

## Overview
Mean-variance optimization is highly sensitive to estimation error in expected returns. This project studies two allocation approaches that sidestep the problem by only using the covariance matrix — Equal Risk Contribution (ERC) and Hierarchical Risk Parity (HRP) — on the EURO STOXX 50 universe, with monthly rebalancing and a 2-year rolling estimation window.

## Methods
- Equal Risk Contribution via constrained numerical optimization
- Hierarchical Risk Parity: correlation-distance clustering, dendrogram-based and recursive-bisection allocation traversals, market-mode detoning
- Ledoit-Wolf covariance shrinkage (constant-correlation and single-factor targets) as alternative inputs to both allocators
- Cluster stability measured via the Rand index across rebalances
- Drift-aware turnover accounting and transaction-cost-adjusted backtests

## Notebook
**`risk_based_allocation.ipynb`**:
- **Part I — ERC**: closed-form sanity checks (ERC reduces to inverse-volatility under a diagonal covariance), then an out-of-sample backtest across the three covariance estimators.
- **Part II — HRP**: detoning intuition (dendrogram + quasi-diagonalised heatmap), a generic parametric implementation, cluster stability over time, and sensitivity of the resulting weights to the linkage method and allocation traversal.
- **Part III — Connecting the dots**: ERC vs. HRP head-to-head — robustness to the covariance estimator, performance/vol/drawdown/turnover, and the impact of transaction costs on the turnover gap between the two.

## Results
- HRP shows materially higher turnover than ERC — the gap survives across all three covariance estimators and shows up as a real net-of-cost performance drag once transaction costs are applied.
- Detoning improves cluster stability (higher Rand index) at every cluster count tested.
- The choice of linkage method and allocation traversal has a larger effect on HRP portfolio concentration than the choice of covariance estimator.

## How to run
Open `risk_based_allocation.ipynb` and run all cells; it expects this folder as the working directory, since data is read from `data/` via relative paths. Requires `numpy`, `pandas`, `matplotlib`, `scipy`, `scikit-learn`, and `joblib`.

## Structure
- `data/`: EURO STOXX 50 daily prices
- `utilities/`: covariance, shrinkage, PCA/detoning, hierarchical clustering, HRP, portfolio optimization, and backtest modules
- `risk_based_allocation.ipynb`: the notebook described above
