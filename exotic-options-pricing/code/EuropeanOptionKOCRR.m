function optionPrice = EuropeanOptionKOCRR(F0, K, KO, B, T, sigma, N)
% EuropeanOptionKOCRRcomputes the price of an U&O european barrier call
% option with CRR approach

% INPUT:
% F0:    forward price
% K:     strike
% KO:    Up & Out Barrier level (checked only at maturity)
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of time steps
    
dt = T/N;

u = exp(sigma * sqrt(dt)); % Up factor
d = exp(-sigma * sqrt(dt)); % Down factor
p = (1 - d) / (u - d); % Risk-neutral probability

% Initialize forward prices at maturity
FT = zeros(N+1, 1);
for i = 0:N
    FT(i+1) = F0 * d^i * u^(N-i);
end      

% Calculate option values at maturity (FT = ST)
V = max(0, FT - K);
% Check knock-out barrier
V(FT>KO)=0;

l = N;
% Backward induction
for j = N:-1:1
    V = p * V(1:l) + (1 - p) * V(2:l+1); % Update option values
    l = l-1;
end 

optionPrice= B * V;

end % function EuropeanOptionKOCRR

