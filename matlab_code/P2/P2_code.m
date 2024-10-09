clc; clear; close all;

% SETAR(2,1,1) Parameter Estimation for Contour Plots

% Define the grid for the two parameters
p1_values = linspace(3.5, 4.5, 100);  % Example grid for p1 (a1 -> theta_1)
p2_values = linspace(0.6, 0.8, 100);  % Example grid for p2 (b1 -> theta_2)

% Generate the grid
[p1_grid, p2_grid] = meshgrid(p1_values, p2_values);

% Initialize the loss function grid (Q_N)
Q_N = zeros(size(p1_grid));

% Define the true parameters
true_threshold = 0.5;  % Threshold for regime change
a1 = 4; b1 = 0.7;     % True parameters for regime 1
a2 = -4; b2 = 0.7;    % True parameters for regime 2
sigma = 1;            % Noise standard deviation

% Generate synthetic time series data
n = 10e3;  % Total number of observations
X = zeros(n, 1);
rng(1);  % For reproducibility

% Simulate the SETAR(2,1,1) process
for t = 2:n
    if X(t-1) <= true_threshold
        X(t) = a1 + b1 * X(t-1) + sigma * randn;
    else
        X(t) = a2 + b2 * X(t-1) + sigma * randn;
    end
end

% Define the loss function calculation (Q_N)
% This function computes the sum of squared residuals for a given p1, p2
for i = 1:length(p1_values)
    for j = 1:length(p2_values)
        p1 = p1_grid(i, j);
        p2 = p2_grid(i, j);
        theta = [p1, p2, a2, b2, true_threshold];  % Keeping a2, b2, threshold fixed
        
        % Calculate residual sum of squares (RSS) for this combination of p1, p2
        Q_N(i, j) = setar_rss_with_threshold(theta, X);
    end
end

% Subsets of the data to be used for contour plots
subsets = {1:10e3, 1:5e3, 1:1e3, 1:1e2, 5001:5500, 5001:5050};

% Create a large figure window
figure('Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);  % This makes the figure fullscreen

% Create a 3x2 tiled layout for more compact subplots
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
t.TileSpacing = 'compact';  % Remove extra space between plots
t.Padding = 'compact';      % Remove padding around the layout

for s = 1:length(subsets)
    subset = subsets{s};
    
    %  Q_N subsets
    Q_N_subset = zeros(size(p1_grid));
    for i = 1:length(p1_values)
        for j = 1:length(p2_values)
            p1 = p1_grid(i, j);
            p2 = p2_grid(i, j);
            theta = [p1, p2, a2, b2, true_threshold];  % Keeping a2, b2, threshold fixed
            Q_N_subset(i, j) = setar_rss_with_threshold(theta, X(subset));
        end
    end
    
    % Create subplot in 3x2 layout
    nexttile;
    contourf(p1_grid, p2_grid, Q_N_subset, 20, 'LineColor', 'none');
    colorbar;
    
    % Label the axes with theta_1 and theta_2
    xlabel('\theta_1', 'FontSize', 12, 'FontWeight', 'bold');  % For parameter p1 (a1)
    ylabel('\theta_2', 'FontSize', 12, 'FontWeight', 'bold');  % For parameter p2 (b1)
    
    title(sprintf('Subset %d:%d', subset(1), subset(end)), 'FontSize', 12);
    
    % Plot the specific point (0.7, 4) in each subplot
    hold on;
    plot(4, 0.7, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % Red circle for the point
    hold off;
end

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

current_file = 'P2_code.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX2/matlab_code/P2';
copyfile(current_file, new_directory);
