function delta = DeltaKO(F0, K, KO, B, T, sigma, d)
% Let's program the function in the same logic of the Vega one
% Add a bump in S0 since the derivation is in S0
S0 = F0*B*exp(d*T);
h = S0/100;
F0_up = (S0+h)*exp(-d*T)/B;
F0_down = (S0-h)*exp(-d*T)/B;

% Just use the closed formula since we can

up_value = EuropeanKOCall_ClosedFormula(F0_up, K, KO, B, T, sigma);
down_value = EuropeanKOCall_ClosedFormula(F0_down, K, KO, B, T, sigma);
% Using centered differences we get a much better convergence
delta = (up_value - down_value)/(2*h);
return