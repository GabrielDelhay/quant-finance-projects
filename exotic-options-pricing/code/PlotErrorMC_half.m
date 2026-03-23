function [M_half, stdEstim] = PlotErrorMC_half(F0, K, B, T, sigma)
    % PLOTERRORMC Computes and plots the Monte Carlo standard error decay.
    
    % INPUT:
    % F0:    forward price
    % K:     strike
    % B:     discount factor
    % T:     time-to-maturity
    % sigma: volatility
    
    % Parameters
    flag = 1;             % Call option
    m = 1:20;             % Let's simulate the half
    M_half = 2.^m;         % Simulations
    stdEstim = zeros(1, length(M_half)); 
    tol = 0.0001 * F0;    % 1 basis point tolerance
    
    % Monte Carlo Loop
    for i = 1:length(M_half)
        %rng(1); %let's fix the seed for explainability advantages
        % Generate random shocks and terminal prices
        Z_up = randn(M_half(i), 1);
        Z_down = -Z_up;
        ST_up = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z_up);
        ST_down = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z_down);

        % Discounted payoffs
        
        payoffs_up= B * max(flag * (ST_up - K), 0);
        payoffs_down= B * max(flag * (ST_down - K), 0);
        payoffs = (payoffs_up+payoffs_down)./2;
        % Standard Error: std / sqrt(M)
        %let's multiply the M_half to get to the real number of symmetric observation
        stdEstim(i) = std(payoffs) / sqrt(M_half(i)); 
    end
    
    % Visualization
    % Plot the calculated Standard Error
    loglog(M_half, stdEstim, '-s', 'LineWidth', 1.2);
    hold on;
    
    % Add the 1/sqrt(M) reference line
    % We anchor it to the first point: Error(M) = Error(1) * sqrt(M1) / sqrt(M)
    refLine = stdEstim(10) * sqrt(M_half(10)) ./ sqrt(M_half);
    loglog(M_half, refLine, '--k', 'LineWidth', 1.5);

    % Formatting
    title('ErrorMC VS ErrorMC Antithetic');
    xlabel('Number of Simulations (M)');
    ylabel('Standard Error');
    legend('MC Standard Error', '1/sqrt(M) Reference', '1/2 bp Tolerance', 'MC Antithetic error', 'Location', 'best');
    grid on;
end
