function price = EuropeanKOCall_ClosedFormula(F0, K, KO, B, T, sigma)
    % EuropeanKOCall_ClosedFormula computes the price of an U&O european
    % barrier call option with closed formula

    % INPUT:
    % F0:    forward price
    % K:     strike
    % KO:    Up & Out Barrier level (checked only at maturity)
    % B:     discount factor
    % T:     time-to-maturity
    % sigma: volatility

    % If the strike is higher than the barrier, option value is zero
    if K >= KO
        price = 0;
        return;
    end

    % Call Vanilla with strike K 
    d1_K = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
    d2_K = d1_K - sigma * sqrt(T);
    Call_K = B * (F0 * normcdf(d1_K) - K * normcdf(d2_K));

    % Call Vanilla with strike KO 
    d1_KO = (log(F0 / KO) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
    d2_KO = d1_KO - sigma * sqrt(T);
    Call_KO = B * (F0 * normcdf(d1_KO) - KO * normcdf(d2_KO));

    % Digital option (Cash-or-Nothing) in KO
    Digital_KO = B * normcdf(d2_KO);

    % Final price
    price = Call_K - Call_KO - (KO - K) * Digital_KO;
    
end