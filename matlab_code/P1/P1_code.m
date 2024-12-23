clc; clear; close all;

% SETAR(2,1,1) Parameter Estimation with Threshold using Prediction Error Method

% Number of observations
n = [500, 2e3, 5e3, 10e3, 50e3, 100e3];

% Generate time series data for a SETAR(2,1,1) model
% Define true parameters for simulation
true_threshold = 0.5; % True threshold for regime change
a1 = 4; b1 = 0.7;     % Parameters for regime 1
a2 = -4; b2 = 0.7;    % Parameters for regime 2
sigma = 1;            % Noise standard deviation

theta0 = [a1, b1, a2, b2, true_threshold];

for i = 1:length(n)

    % Initialize time series
    X = zeros(n(i), 1);
    rng(1);  % For reproducibility
    
    % Simulate the SETAR(2,1,1) process
    for t = 2:n(i)
        if X(t-1) <= true_threshold
            X(t) = a1 + b1 * X(t-1) + sigma * randn;
        else
            X(t) = a2 + b2 * X(t-1) + sigma * randn;
        end
    end
    
    % Objective function for parameter estimation
    % theta = [a1, b1, a2, b2, threshold]
    objective_function = @(theta) setar_rss_with_threshold(theta, X);
    
    % Set up a grid of initial guesses
    initial_guesses = [
        3.9, 0.6, -3.9, 0.6, 0.4;
        4.5, 0.7, -4.5, 0.7, 0.6;
        4, 0.7, -4, 0.7, 0.5;
    ];
    
    best_params = [];
    best_rss = Inf;
    
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'TolX', 1e-8, 'TolFun', 1e-8);
    
    for ig = 1:size(initial_guesses, 1)
        initial_guess = initial_guesses(ig, :);
        
        % Run fminunc from each initial guess
        problem = createOptimProblem('fminunc', 'objective', objective_function, 'x0', initial_guess, 'options', options);
        ms = MultiStart('UseParallel', true);
        [estimated_params, current_rss] = run(ms, problem, 50); 
        
        if current_rss < best_rss
            best_rss = current_rss;
            best_params = estimated_params;
        end
    end
    
    % Output the best result
    fprintf('Best estimated parameters: a1 = %.4f, b1 = %.4f, a2 = %.4f, b2 = %.4f, threshold = %.4f\n', ...
        best_params(1), best_params(2), best_params(3), best_params(4), best_params(5));

    
    % Collect data
    param_matrix{i} = best_params;
    RSS_matrix(i) = sum(abs(theta0 - param_matrix{i}));

end

%%

fontSize = 13;
fontWeight = 'bold';
customBlue = [24/255, 54/255, 104/255];

figure;
plot(n, RSS_matrix, 'LineWidth', 2, 'Color', customBlue);

xticks([0, 2e4, 4e4, 6e4, 8e4, 10e4]);
xlabel('Number of Observations (n)', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('Residual Sum of Squares (RSS)', 'FontWeight', fontWeight, 'FontSize', fontSize);
title('RSS vs Sample Size');
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize);
box on;

%% Function to compute residual sum of squares (RSS) for SETAR(2,1,1) with threshold
function rss = setar_rss_with_threshold(theta, X)
    a1 = theta(1);
    b1 = theta(2);
    a2 = theta(3);
    b2 = theta(4);
    threshold = theta(5);
    
    n = length(X);
    residuals = zeros(n, 1);
    
    % Calculate residuals based on the SETAR(2,1,1) model with estimated threshold
    for t = 2:n
        if X(t-1) <= threshold
            % Regime 1
            X_pred = a1 + b1 * X(t-1);
        else
            % Regime 2
            X_pred = a2 + b2 * X(t-1);
        end
        residuals(t) = X(t) - X_pred;
    end
    
    % Compute the sum of squared residuals
    rss = sum(residuals.^2);
end

%%

current_file = 'P1_code.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX2/matlab_code/P1';
copyfile(current_file, new_directory);

