# Implied Volatility Surface — Variance Risk Premium

## Overview
Does the S&P 500 implied volatility surface predict future realized volatility beyond what the return history alone provides? This project builds the daily IV surface from raw SPX options chains and tests it against a battery of historical-return and hybrid econometric models, evaluated strictly out-of-sample.

## Methods
- Black-Scholes IV inversion (Brent's method) with a per-(date, expiry) implied discount rate and forward price extracted via put-call parity regression (Brenner & Galai, 1986), rather than a fixed risk-free rate
- Constant-maturity (30-day) ATM IV, skew, and convexity via variance-space interpolation between bracketing expiries — the CBOE/VIX methodology
- Historical benchmarks: historical mean, HAR (Corsi, 2009), EGARCH(1,1)
- IV-based models: nested OLS (ATM IV → skew → convexity → VIX–ATM spread), Mincer-Zarnowitz unbiasedness test
- VAR/VECM robustness checks: stationarity, Johansen cointegration, Granger causality
- Panel fixed/random effects across expiries (Hausman test) to motivate a full-surface factor model
- PCA on the IV surface (level/slope/curvature factors)
- Hybrid HAR-IV and HAR-PCA models testing whether IV/surface factors add information beyond the return history
- Expanding-window out-of-sample evaluation (RMSE, MAE, QLIKE, Diebold-Mariano vs. HAR), with an explicit robustness check excluding the February 2018 "Volmageddon" spike
- A simple variance-risk-premium trading strategy built from the best OOS model

## Notebooks
1. **`Features_Engineering.ipynb`**: builds the cleaned options dataset, the IV surface, and the daily ATM IV / skew / convexity / realized-volatility time series, saved to `time_series.pkl`.
2. **`Econometrics_Models.ipynb`**: loads the pickled features and runs every model above, from the historical benchmarks through the OOS horse race and the trading application.

## Results
- IV carries real predictive content over a pure return-history baseline, but the OOS ranking is sensitive to the estimation window and the February 2018 volatility spike — see the Volmageddon robustness section for the honest read.
- The variance risk premium (IV² − realized variance) is positive on average, consistent with a standard insurance-premium interpretation of index options.
- Full numbers depend on re-running the notebooks; see `Econometrics_Models.ipynb` section 7 for the current OOS comparison table.

## How to run
Run `Features_Engineering.ipynb` first (from this folder; it reads `data/*.csv` and writes `time_series.pkl` — the IV computation over the full dataset takes several minutes), then `Econometrics_Models.ipynb`. Requires `numpy`, `pandas`, `scipy`, `statsmodels`, `linearmodels`, `arch`, `scikit-learn`, `yfinance`, and `matplotlib`.

## Structure
- `data/`: daily SPX options chains (one CSV per trading date)
- `utilities/`: IV computation, data cleaning, time-series/panel construction, OLS, and fixed/random-effects helpers
- `Features_Engineering.ipynb`, `Econometrics_Models.ipynb`: the two notebooks described above
