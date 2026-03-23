function [M_vec, stdEstim] = PlotErrorMC(F0, K, B, T, sigma)
    % PLOTERRORMC Computes and plots the Monte Carlo standard error decay.
    
    % INPUT:
    % F0:    forward price
    % B:     discount factor
    % K:     strike
    % T:     time-to-maturity
    % sigma: volatility
    
    % OUTPUT:
    %   M_vec    : Row vector containing the number of steps used (2^1, ..., 2^10).
    %   errorCRR : Row vector containing the absolute error for each M.

    % Parameters
    flag = 1;             % Call option
    m = 1:20;             % Exponents from 1 to 20
    M_vec = 2.^m;         % Simulations: 2^1, 2^2, ..., 2^20
    stdEstim = zeros(1, length(M_vec)); 
    bp = 1/2 * 0.0001;    % 1 basis point tolerance
    
    
    % Monte Carlo Loop
    for i = 1:length(M_vec)
        M = M_vec(i);
        %rng(3); %let's fix the seed for explainability advantages
        % Generate random shocks and terminal prices
        Z = randn(M, 1);
        ST = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z);
        
        % Discounted payoffs
        payoffs = B * max(flag * (ST - K), 0);
        
        % Standard Error: std / sqrt(M)
        stdEstim(i) = std(payoffs) / sqrt(M); 
    end
    
    % Visualization
    figure;
    % Plot the calculated Standard Error
    loglog(M_vec, stdEstim, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    hold on;
    
    % Add the 1/sqrt(M) reference line
    % We anchor it to the first point: Error(M) = Error(1) * sqrt(M1) / sqrt(M)
    refLine = stdEstim(20) * sqrt(M_vec(20)) ./ sqrt(M_vec);
    loglog(M_vec, refLine, '--k', 'LineWidth', 1.5);
    
    % Add the tolerance line
    yline(bp, '--r', '1/2 bp Tolerance');
    
    % Formatting
    title('Monte Carlo Convergence');
    xlabel('Number of Simulations (M)');
    ylabel('Standard Error');
    legend('MC Standard Error', '1/sqrt(M) Reference');
    grid on;
    
    % M Selection
    % Find first M below tolerance (ignoring zero results at start)
    idx = find(stdEstim <= bp & stdEstim > 0, 1, 'first');
    
    fprintf('\n--- MC ANALYSIS ---\n');
    if ~isempty(idx)
        fprintf('Tolerance reached at M = 2^%d (%d)\n', m(idx), M_vec(idx));
    else
        fprintf('Tolerance not reached within M = 2^20\n');
    end
end