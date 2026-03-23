function optionPrice = EuropeanOptionCRR(F0, K, B, T, sigma, N, flag)
    % INPUT:
    % F0:    forward price
    % B:     discount factor
    % K:     strike
    % T:     time-to-maturity
    % sigma: volatility
    % N:     either number of time steps (knots for CRR tree)
    %        or number of simulations in MC   
    % flag:  1 call, -1 put

    % Initial Parameters
    dt = T / N;                  % size of the single time step
    u = exp(sigma * sqrt(dt));   % up factor
    d = 1 / u;                   % down factor
    p = (1 - d) / (u - d);       % risk-neutral probability
    
    % Initialize the vector that will store the option values
    V = zeros(N + 1, 1);
    
    % Payoff calculation at maturity (Time T, Step N)
    % i represents the number of "Up" movements made
    for i = 0:N
        % The final price depends on how many "Up" (i) and "Down" (N-i) movements occurred
        ST = F0 * (u^(N - i)) * (d^i);
        
        % Payoff calculation
        if flag == 1
            V(i + 1) = max(ST - K, 0); % Call Payoff
        else
            V(i + 1) = max(K - ST, 0); % Put Payoff
        end
    end
    
    % Backward Induction (Moving backwards in time)
    % step starts from N-1 and goes down to 0 (the initial node)
    for step = N-1 : -1 : 0
        
        % At each time step, we have 'step + 1' nodes
        % i always represents the number of "Up" movements accumulated so far
        for i = 0:step
            
            % The value in case of an upward move (Up) is at the next node (i+2 in the array)
            V_up = V(i + 1);
            
            % The value in case of a downward move (Down) is at the current node (i+1 in the array)
            V_down = V(i + 2);
            
            % The option value at the current node is the expected value
            V(i + 1) = p * V_up + (1 - p) * V_down;
            
        end
    end
    
    % Discounting to time zero
    optionPrice = B * V(1);
end