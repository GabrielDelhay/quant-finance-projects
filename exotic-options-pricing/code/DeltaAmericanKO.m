function delta = DeltaAmericanKO(F0, K, KO, B, T, sigma, q)
% Let's program the functino in the same style as we did for Vega

S0 = F0*B*exp(q*T);
h = S0/1000;
%let's bump S0 to get the actual derivative in S0
F0_up = (S0+h)*exp(-q*T)/B;
F0_down = (S0-h)*exp(-q*T)/B;

% Just use the closed formula since we can
up_value = EuropeanOptionAmericanBarrier(F0_up, K, KO, B, T, sigma, q);
down_value = EuropeanOptionAmericanBarrier(F0_down, K, KO, B, T, sigma, q);
% Using centered differences we get a much better convergence
delta = (up_value - down_value)/(2*h);
return