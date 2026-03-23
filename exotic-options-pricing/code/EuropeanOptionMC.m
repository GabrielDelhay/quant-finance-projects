function optionPrice=EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
% Option Price with different pricing methods
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
rng(1);
Z = randn(N, 1); %generation of std normal random variables 

ST = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z); 

if flag == 1
    payoff = max(ST - K, 0); % Call payoff at expiry
else
    payoff = max(K - ST, 0); % Put payoff at expiry
end

optionPrice = B * mean(payoff);  %actualizing and taking the mean

end % function EuropeanOptionMC

