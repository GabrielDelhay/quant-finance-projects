function optionPrice=EuropeanOptionKOMC(F0,K,KO,B,T,sigma,N)
% EuropeanOptionKOCRRcomputes the price of an U&O european barrier call
% option with MC simulation
%
% INPUT:
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     either number of time steps (knots for CRR tree)
%        or number of simulations in MC   
% flag:  1 call, -1 put
rng(1); %let's fix the seed for consistency reason
Z = randn(N, 1); %generation of std normal random variables 
ST = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z); 
for i=1:N
    if ST(i)<KO
        payoff(i)= max(ST(i) - K, 0); % Call payoff at expiry
    else
        payoff(i)= 0;
    end
end

optionPrice = B * mean(payoff);  %actualizing and taking the mean

end % function EuropeanOptionMC

