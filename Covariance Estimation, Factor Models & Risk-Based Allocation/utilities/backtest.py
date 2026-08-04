import pandas as pd


def _align_portfolios_and_returns(
    portfolios: pd.DataFrame,
    returns: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align a portfolio schedule with the return matrix on the common investable universe.

    The portfolio decided at a rebalance date's close is traded the next business day
    (T+1) and only starts earning P&L from T+2: positions are lagged by two trading
    days relative to the raw rebalance schedule.

    Parameters:
        portfolios (pd.DataFrame): Portfolio weights indexed by rebalance date.
        returns (pd.DataFrame): Asset returns indexed by trading date.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: The aligned portfolio and return matrices.
    """
    if portfolios.empty or portfolios.shape[1] == 0:
        raise ValueError("portfolios must contain at least one asset")

    overlapping_assets = portfolios.columns.intersection(returns.columns)
    if overlapping_assets.empty:
        raise ValueError("portfolios and returns must share at least one asset")

    portfolios = portfolios.set_axis(pd.to_datetime(portfolios.index), axis=0)
    aligned_portfolios = (
        portfolios.loc[:, overlapping_assets]
        .reindex(returns.index)
        .ffill()
        .shift(2)
        .dropna(axis=0, how="all")
    )
    if aligned_portfolios.empty:
        raise ValueError(
            "portfolios and returns must overlap on at least one investable date"
        )

    aligned_returns = returns.loc[aligned_portfolios.index, overlapping_assets]

    return aligned_portfolios.fillna(0.0), aligned_returns


def _drift_aware_turnover(
    portfolios: pd.DataFrame,
    returns: pd.DataFrame,
    overlapping_assets: pd.Index,
) -> pd.Series:
    """
    One-way turnover at each rebalance date, comparing the new target weights against
    the *actual* pre-rebalance weights -- the previous target drifted forward by the
    realized returns of the holding period, not frozen at their old values.

    On the first rebalance the portfolio opens from cash, so the actual weights are
    zero and the full target book counts as turnover.

    Parameters:
        portfolios (pd.DataFrame): Target portfolio weights indexed by rebalance date.
        returns (pd.DataFrame): Asset returns indexed by trading date.
        overlapping_assets (pd.Index): Assets common to both portfolios and returns.

    Returns:
        pd.Series: One-way turnover at each rebalance date.
    """
    portfolios = portfolios.set_axis(pd.to_datetime(portfolios.index), axis=0)
    reb_dates = portfolios.index.intersection(returns.index)
    turnover = pd.Series(index=reb_dates, dtype=float)

    prev_weights = None
    prev_date = None
    for date in reb_dates:
        target = portfolios.loc[date, overlapping_assets].fillna(0.0)

        if prev_weights is None:
            turnover[date] = target.abs().sum()
        else:
            window_returns = (
                returns.loc[prev_date:date, overlapping_assets].iloc[1:].fillna(0.0)
            )
            cum_return = (1.0 + window_returns).prod() - 1.0
            drifted = prev_weights * (1.0 + cum_return)
            denom = drifted.sum()
            actual = drifted / denom if denom != 0 else drifted
            turnover[date] = (target - actual).abs().sum()

        prev_weights = target
        prev_date = date

    return turnover


def portfolio_returns(
    portfolios: pd.DataFrame,
    returns: pd.DataFrame,
    transaction_costs: float = 0.0,
) -> pd.Series:
    """
    Compute the time series of realized portfolio returns from a rebalance schedule.

    Parameters:
        portfolios (pd.DataFrame): DataFrame where each column represents the weights of a
            portfolio over time (dates as index).
        returns (pd.DataFrame): DataFrame of asset returns (dates as index, assets as columns).
        transaction_costs (float): All-in transaction cost rate applied to one-way turnover
            (e.g. 0.001 = 10 bps). Defaults to 0.0.

    Returns:
        pd.Series: Series of realized portfolio returns over time.
    """
    aligned_portfolios, aligned_returns = _align_portfolios_and_returns(
        portfolios=portfolios,
        returns=returns,
    )

    gross_returns = aligned_portfolios.multiply(aligned_returns).sum(axis=1)

    if transaction_costs != 0.0:
        overlapping_assets = portfolios.columns.intersection(returns.columns)
        turnover_at_rebalance = _drift_aware_turnover(
            portfolios=portfolios,
            returns=returns,
            overlapping_assets=overlapping_assets,
        )
        # Turnover is realized on the same day the new weights start earning P&L.
        aligned_turnover = (
            turnover_at_rebalance.reindex(returns.index)
            .shift(2)
            .reindex(aligned_portfolios.index)
            .fillna(0.0)
        )
        gross_returns -= aligned_turnover * transaction_costs

    return gross_returns


def backtest(
    portfolios: pd.DataFrame,
    returns: pd.DataFrame,
    transaction_costs: float = 0.0,
) -> pd.Series:
    """
    Compute cumulative portfolio performance from a rebalance schedule.

    Parameters:
        portfolios (pd.DataFrame): DataFrame where each column represents the weights of a
            portfolio over time (dates as index).
        returns (pd.DataFrame): DataFrame of asset returns (dates as index, assets as columns).
        transaction_costs (float): All-in transaction cost rate applied to one-way turnover
            (e.g. 0.001 = 10 bps). Defaults to 0.0.

    Returns:
        pd.Series: Series representing the cumulative returns of the portfolios over time.
    """
    return (
        portfolio_returns(
            portfolios=portfolios,
            returns=returns,
            transaction_costs=transaction_costs,
        )
        .add(1)
        .cumprod()
    )
