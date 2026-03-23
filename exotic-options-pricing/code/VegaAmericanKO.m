function Vega = VegaAmericanKO(F0, K, KO, B, T, sigma, q)
% Let's program the functino in the same style as we did for Vega

h = sigma/100;
% Just use the closed formula since we can
% using centered differences we get a much better convergence
up_value = EuropeanOptionAmericanBarrier(F0, K, KO, B, T, sigma+h, q);
down_value = EuropeanOptionAmericanBarrier(F0, K, KO, B, T, sigma-h, q);

Vega = (up_value - down_value)/(2*h*100);% /100 to have the result in euro
return