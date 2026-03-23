function [M, errorCRR] = PlotErrorCRR(F0, K, B, T, sigma)
    % PLOTERRORCRR Computes and plots the CRR method error as the number of steps M varies.
    
    % INPUT:
    %   F0, K, B, T, sigma : Standard option pricing parameters.
    
    % OUTPUT:
    %   M        : Row vector containing the number of steps used (2^1, ..., 2^10).
    %   errorCRR : Row vector containing the absolute error for each M.

    % Define the exponents from 1 to 10
    m = 1:10;
    
    % Generate the row vector M as powers of 2 (2^1, 2^2, ..., 2^10)
    M = 2.^m;
    
    % Pre-allocate the row vector for error results
    errorCRR = zeros(1, length(M));
    
    % Calculate the "exact value" using the closed formula (Call option, flag=1)
    exact = EuropeanOptionClosed(F0, K, B, T, sigma, 1);
    
    % Loop to calculate the error for each value of M
    for i = 1:length(M)
        % Price the option using the CRR tree with M(i) steps
        approx = EuropeanOptionCRR(F0, K, B, T, sigma, M(i), 1);
        
        % Absolute error calculation
        errorCRR(i) = abs(approx - exact);
    end
    
    % Visualization
    figure('Name', 'Convergence of CRR Model');
    
    % Create log-log plot to visualize convergence rate
    loglog(M, errorCRR, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    hold on;

    % Add reference line for O(1/M) convergence
    loglog(M, (errorCRR(10)*M(10))./M, '--k', 'LineWidth', 1.5);
    
    % Tolerance line (1 bp of the Forward as per Hint 2)
    bp =1/2 * 0.0001;
    yline(bp, '--r', '1/2 bp Tolerance', 'LabelHorizontalAlignment', 'left');
    
    % Formatting the plot
    title('CRR Tree Error Convergence');
    xlabel('Number of Steps (M)');
    ylabel('Absolute Error');
    legend('CRR Error', 'O(1/M) Reference', '1/2 bp Tolerance');
    grid on;
    
    % Console Output for M selection
    valid_idx = find(errorCRR <= bp, 1);
    fprintf('\n--- CRR ERROR ANALYSIS ---\n');
    if ~isempty(valid_idx)
        fprintf('Criteria satisfied! Selected M: %d (2^%d)\n', M(valid_idx), m(valid_idx));
        fprintf('Error reached: %.6f\n', errorCRR(valid_idx));
    else
        fprintf('Warning: 1 bp tolerance (%.4f) NOT reached within M = 1024.\n', bp);
    end
end