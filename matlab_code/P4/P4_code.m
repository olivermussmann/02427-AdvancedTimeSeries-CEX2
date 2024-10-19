clc; clear; close all;

% Parameters
a_true          = 0.4;          % True value of 'a'
sigma_e         = 1;            % Standard deviation of measurement noise (e_t)
num_simulations = 20;           % Number of simulations
num_steps       = 1000;         % Number of time steps per simulation
t_vec           = 1:num_steps;

% State space model matrices (nonlinear in 'a')
f = @(Z) [Z(2) * Z(1); Z(2)];   % State transition function
h = @(Z) Z(1);                  % Observation function

% Jacobians for EKF
F_jacobian = @(Z) [Z(2), Z(1); 0, 1];   % Derivative of f(Z)
H_jacobian = [1, 0];                    % Derivative of h(Z)

% Noise covariance matrices
R = sigma_e^2;                          % Measurement noise covariance

% Define scenarios
scenarios = [0.5, 1, 10; -0.5, 1, 10; 0.5, 1, 1; -0.5, 1, 1; ...
             0.5, 10, 10; -0.5, 10, 10; 0.5, 10, 1; -0.5, 10, 1];

% Loop through each figure to group scenarios in pairs
for fig_idx = 1:4
    % Create new figure
    figure('Position', [700, 300, 900, 400]);
    
    for subplot_idx = 1:2
        % Determine the scenario index
        scenario_idx = (fig_idx - 1) * 2 + subplot_idx;
        
        % Extract parameters for the current scenario
        a0 = scenarios(scenario_idx, 1);
        Ra = scenarios(scenario_idx, 2);
        sigma_v = scenarios(scenario_idx, 3);

        % Store estimates and variance of a
        a_est = zeros(num_simulations, num_steps);
        var_a_est = zeros(num_simulations, num_steps);

        % Perform 20 simulations
        for sim = 1:num_simulations
            % Initialize variables
            Z = [0; a0];    % Initial state and parameter guess
            P = eye(2);     % Initial covariance estimate
            P(2, 2) = Ra;   % Set Ra as initial variance for parameter a
            x_true = 0;     % True state initialization
            
            % Storage for results
            x_estimates = zeros(num_steps, 1);
            a_estimates = zeros(num_steps, 1);
            var_estimates = zeros(num_steps, 1);

            for t = 1:num_steps
                % True system evolution
                v_t = sigma_v * randn;
                e_t = sigma_e * randn;
                x_true = a_true * x_true + v_t;
                y_t = x_true + e_t;

                % Prediction step
                Z_pred = f(Z);              % Predicted state
                F = F_jacobian(Z);          % Jacobian of f(Z)
                Q = [sigma_v^2, 0; 0, 0];   % Process noise covariance
                P_pred = F * P * F' + Q;    % Predicted covariance

                % Update step
                H = H_jacobian;  % Jacobian of h(Z)
                K = P_pred * H' / (H * P_pred * H' + R);    % Kalman gain
                Z = Z_pred + K * (y_t - h(Z_pred));         % Updated state estimate
                P = (eye(2) - K * H) * P_pred;              % Updated covariance estimate

                % Store results
                x_estimates(t) = Z(1);          % Estimated state
                a_estimates(t) = Z(2);          % Estimated parameter a
                var_estimates(t) = P(2, 2);     % Estimated variance of a
            end

            a_est(sim, :) = a_estimates;
            var_a_est(sim, :) = var_estimates;
        end

        % Plot for the current scenario
        subplot(1, 2, subplot_idx);
        hold on;

        % Plot the first 'a_est' with legend name, and the rest without
        plot(t_vec, a_est(1, :), 'Color', [1, 0, 0, 0.5], 'LineWidth', 1, 'DisplayName', '$\hat{\theta}$');
        for i = 2:num_simulations
            plot(t_vec, a_est(i, :), 'Color', [1, 0, 0, 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        % Plot the first 'var_a_est' with legend name, and the rest without
        plot(t_vec, var_a_est(1, :), 'Color', [0, 0, 1, 0.5], 'LineWidth', 1, 'DisplayName', 'Var($\hat{\theta}$)');
        for i = 2:num_simulations
            plot(t_vec, var_a_est(i, :), 'Color', [0, 0, 1, 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        % Plot true value of 'a'
        plot(t_vec, a_true * ones(1, length(t_vec)), 'Color', 'k', 'LineWidth', 3, 'DisplayName', 'True a');

        % Labels, grid, and legend
        ylabel('\theta', 'FontSize', 10, 'FontWeight', 'bold');
        xlabel('Time step', 'FontSize', 10, 'FontWeight', 'bold');
        title(['Scenario ', num2str(scenario_idx)], 'FontSize', 12, 'FontWeight', 'bold');
        ylim([-1.5, 1.5])
        grid on;
        box on;
        legend('show', 'FontSize', 10, 'FontWeight', 'bold', 'Location', 'best', 'interpreter', 'latex');
        hold off;
    end
end

current_file = 'P4_code.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX2/matlab_code/P4';
copyfile(current_file, new_directory);