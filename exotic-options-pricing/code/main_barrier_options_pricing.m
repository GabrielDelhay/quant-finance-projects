% =========================================================================
%  RUN_BARRIER_OPTIONS_PRICING.m
%
%  Pricing and Greeks analysis of European and American Barrier Options
%  using Closed-Form, CRR Binomial Tree, and Monte Carlo methods.
%
%  Topics covered:
%    - European option pricing (Closed-Form / CRR / Monte Carlo)
%    - Convergence analysis: CRR steps and MC simulations
%    - European Up-and-Out (KO) Call option pricing and Vega
%    - American barrier vs. European barrier comparison (Delta & Vega)
%    - Variance reduction via Antithetic Variables
%    - Bermudan option pricing vs. European as a function of dividend yield
%
%  Author : GabrielDelhay
% =========================================================================

clearvars; close all; clc;

%% -------------------------------------------------------------------------
%  1. MODEL PARAMETERS
% --------------------------------------------------------------------------

S0    = 1;       % Initial underlying price
K     = 1.1;     % Strike price
r     = 0.025;   % Risk-free rate (continuous compounding)
TTM   = 1/3;     % Time to maturity (in years)
sigma = 0.212;   % Implied volatility
d     = 0.02;    % Continuous dividend yield
flag  = 1;       % Option type: +1 = Call, -1 = Put

% Derived quantities
B  = exp(-r * TTM);               % Discount factor
F0 = S0 * exp(-d * TTM) / B;     % Forward price

%% -------------------------------------------------------------------------
%  2. EUROPEAN OPTION PRICING
% --------------------------------------------------------------------------

fprintf('\n--- European Option Pricing ---\n');

% Pricing mode: 1 = Closed-Form | 2 = CRR Binomial Tree | 3 = Monte Carlo
pricingMode = 3;
M = 10000;   % Number of MC simulations (or CRR steps)

OptionPrice = EuropeanOptionPrice(F0, K, B, TTM, sigma, pricingMode, M, flag);
fprintf('Option Price (mode=%d): %.6f\n', pricingMode, OptionPrice);

%% -------------------------------------------------------------------------
%  3. CONVERGENCE ANALYSIS
% --------------------------------------------------------------------------

fprintf('\n--- Convergence Analysis ---\n');

% CRR: pricing error as a function of number of steps
[nCRR, errCRR] = PlotErrorCRR(F0, K, B, TTM, sigma);

% Monte Carlo: std of estimator as a function of number of simulations
[nMC, stdEstim] = PlotErrorMC(F0, K, B, TTM, sigma);

%% -------------------------------------------------------------------------
%  4. UP-AND-OUT (KO) EUROPEAN CALL OPTION
% --------------------------------------------------------------------------

fprintf('\n--- KO Option Pricing ---\n');

KO = 1.4;   % Barrier level

% Closed-form reference price
Call_KO_True = EuropeanKOCall_ClosedFormula(F0, K, KO, B, TTM, sigma);
fprintf('KO Call (Closed-Form) : %.6f\n', Call_KO_True);

% CRR approximation
M_CRR = 2^7;
Call_KO_CRR = EuropeanOptionKOCRR(F0, K, KO, B, TTM, sigma, M_CRR);
fprintf('KO Call (CRR, M=%d)   : %.6f\n', M_CRR, Call_KO_CRR);

% Monte Carlo approximation
M_MC = 2^20;
Call_KO_MC = EuropeanOptionKOMC(F0, K, KO, B, TTM, sigma, M_MC);
fprintf('KO Call (MC, M=%d)    : %.6f\n', M_MC, Call_KO_MC);

%% -------------------------------------------------------------------------
%  5. VEGA OF THE KO OPTION — THREE METHODS
% --------------------------------------------------------------------------

fprintf('\n--- KO Option Vega ---\n');

S0_vector  = 0.65:0.01:1.45;
F0_vector  = S0_vector .* exp(-d * TTM) ./ B;

M_CRR_vega = 2^10;   % Higher resolution for CRR Vega

vega_exact   = zeros(size(S0_vector));
vega_num_crr = zeros(size(S0_vector));
vega_num_mc  = zeros(size(S0_vector));

for i = 1:length(S0_vector)
    vega_exact(i)   = VegaKO(F0_vector(i), K, KO, B, TTM, sigma, M,         3);
    vega_num_crr(i) = VegaKO(F0_vector(i), K, KO, B, TTM, sigma, M_CRR_vega, 1);
    vega_num_mc(i)  = VegaKO(F0_vector(i), K, KO, B, TTM, sigma, M_MC,      2);
end

figure;
plot(S0_vector, vega_exact,   'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
plot(S0_vector, vega_num_crr, 'r-',  'LineWidth', 1.5);
plot(S0_vector, vega_num_mc,  'g-',  'LineWidth', 1.5);
xline(KO, '--k', sprintf('Barrier (%.1f)', KO), ...
    'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
xlabel('Underlying Price S_0 (€)');
ylabel('Vega (€)');
title('Vega of Up-and-Out European Call');
legend('Exact (Analytical)', 'Numerical (CRR)', 'Numerical (MC)', 'Location', 'best');
grid on;

%% -------------------------------------------------------------------------
%  6. AMERICAN BARRIER — COMPARISON WITH EUROPEAN BARRIER
% --------------------------------------------------------------------------

fprintf('\n--- American vs. European Barrier ---\n');

S0_vector = 0.65:0.001:1.45;
F0_vector = S0_vector .* exp(-d * TTM) ./ B;

% Price at reference point
Call_American_KO = EuropeanOptionAmericanBarrier(F0, K, KO, B, TTM, sigma, d);
fprintf('American KO Call Price: %.6f\n', Call_American_KO);

% --- Delta comparison ---
Eur_delta      = zeros(length(F0_vector), 1);
American_delta = zeros(length(F0_vector), 1);

for i = 1:length(F0_vector)
    Eur_delta(i)      = DeltaKO(F0_vector(i), K, KO, B, TTM, sigma, d);
    American_delta(i) = DeltaAmericanKO(F0_vector(i), K, KO, B, TTM, sigma, d);
end

figure;
plot(S0_vector, Eur_delta,      '-', 'LineWidth', 1.5);
hold on;
plot(S0_vector, American_delta, '-', 'LineWidth', 1.5);
xline(KO, '--k', sprintf('Barrier (%.1f)', KO), 'LineWidth', 1.2);
xlabel('Underlying Price S_0 (€)');
ylabel('\Delta');
title('Delta: European vs. American Barrier');
legend('European Barrier', 'American Barrier', 'Location', 'best');
grid on;

% --- Vega comparison ---
Eur_vega      = zeros(length(F0_vector), 1);
American_vega = zeros(length(F0_vector), 1);

for i = 1:length(F0_vector)
    Eur_vega(i)      = VegaKO(F0_vector(i), K, KO, B, TTM, sigma, M, 3);
    American_vega(i) = VegaAmericanKO(F0_vector(i), K, KO, B, TTM, sigma, d);
end

figure;
plot(S0_vector, Eur_vega,      '-', 'LineWidth', 1.5);
hold on;
plot(S0_vector, American_vega, '-', 'LineWidth', 1.5);
xline(KO, '--k', sprintf('Barrier (%.1f)', KO), 'LineWidth', 1.2);
xlabel('Underlying Price S_0 (€)');
ylabel('Vega (€)');
title('Vega: European vs. American Barrier');
legend('European Barrier', 'American Barrier', 'Location', 'best');
grid on;

%% -------------------------------------------------------------------------
%  7. VARIANCE REDUCTION — ANTITHETIC VARIABLES
% --------------------------------------------------------------------------

fprintf('\n--- Antithetic Variables Variance Reduction ---\n');

figure;
[M_vec,      stdEstim]      = PlotErrorMC(F0, K, B, TTM, sigma);
hold on;
[M_vec_half, stdEstim_half] = PlotErrorMC_half(F0, K, B, TTM, sigma);
legend('Standard MC', 'Antithetic Variables', 'Location', 'best');
title('MC Std. Error: Standard vs. Antithetic Variables');
xlabel('Number of Simulations M');
ylabel('Std. Error of Estimator');
grid on;

%% -------------------------------------------------------------------------
%  8. BERMUDAN OPTION PRICING
% --------------------------------------------------------------------------

fprintf('\n--- Bermudan Option ---\n');

Bermudan = BermudanOptionPrice(F0, K, TTM, B, sigma, d, M);
fprintf('Bermudan Option Price : %.6f\n', Bermudan);

%% -------------------------------------------------------------------------
%  9. BERMUDAN VS. EUROPEAN — SENSITIVITY TO DIVIDEND YIELD
% --------------------------------------------------------------------------

fprintf('\n--- Bermudan vs. European as a function of dividend yield ---\n');

q_vector         = linspace(0, 0.05, 100);
Bermudan_vector  = zeros(1, length(q_vector));
European_vector  = zeros(1, length(q_vector));

for i = 1:length(q_vector)
    F0_q = S0 * exp(-q_vector(i) * TTM) / B;
    Bermudan_vector(i) = BermudanOptionPrice(F0_q, K, TTM, B, sigma, q_vector(i), M);
    European_vector(i) = EuropeanOptionClosed(F0_q, K, B, TTM, sigma, 1);
end

figure;
plot(q_vector, Bermudan_vector, '-',  'LineWidth', 1.5);
hold on;
plot(q_vector, European_vector, '--', 'LineWidth', 1.5);
xlabel('Dividend Yield q');
ylabel('Option Price (€)');
title('Bermudan vs. European Call Price as a Function of Dividend Yield');
legend('Bermudan', 'European', 'Location', 'best');
grid on;

fprintf('\nDone.\n');




