function Bermudan_price = BermudanOptionPrice(F0, K, T, B, sigma, q, M)
% this function computes the price of a bermudan call using a CRR approach    
% INPUT:
    % F0:    forward price
    % B:     discount factor
    % K:     strike
    % T:     time-to-maturity
    % sigma: volatility
    % M:     either number of time steps (knots for CRR tree)
  
% risk free rate
r = -log(B) / T;

% exercices times
dt = T / M;
Exer_times = [ceil(M/4), ceil(M/2), ceil(3*M/4)];

% parameters (forward measure)
u = exp(sigma * sqrt(dt));
d = 1 / u;
p = (1 - d) / (u - d); 

% Final nodes
j = 0:M; % up-moves
F_T = F0 * (u.^j) .* (d.^(M - j));

% Payoff 
V = max(F_T - K, 0);

% Backward Induction 
for i = (M - 1):-1:0
    
    % Continuation Value 
    V = p * V(2:end) + (1 - p) * V(1:end-1);
    
    % Exerices times right
    if ismember(i, Exer_times)
        t_i = i * dt; 
        
        % Forward price
        j_current = 0:i;
        F_t = F0 * (u.^j_current) .* (d.^(i - j_current));
        
        % Spot Price
        S_t = F_t * exp(-(r - q) * (T - t_i));
        
        % Payoff 
        P_t_T = exp(-r * (T - t_i));
        Payoff_esercizio = max(S_t - K, 0) / P_t_T;
        
        % Bermudan condition
        V = max(V, Payoff_esercizio);
    end
end

% Discount
Bermudan_price = B * V(1);

end