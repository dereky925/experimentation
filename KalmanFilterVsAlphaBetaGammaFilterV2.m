% Maneuvering Target Tracking Simulation
% Comparing Extended Kalman Filter vs Alpha-Beta-Gamma Filter
% with noisy measurements and changing bias

clear;
clc;
close all;

% Simulation parameters
dt = 0.1;                    % Time step (s)
total_time = 60;             % Total simulation time (s)
time = 0:dt:total_time;      % Time vector
N = length(time);            % Number of time steps

% Target initial state [x, y, vx, vy, ax, ay]
X0 = [0; 0; 10; 5; 0; 0];

% Process noise parameters (acceleration variance)
sigma_a = 2;                 % Standard deviation of acceleration noise (m/s²)
Q_continuous = [0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 sigma_a^2 0;
                0 0 0 0 0 sigma_a^2];

% Measurement noise parameters
sigma_pos = 15;              % Position measurement noise (m)
bias_change_rate = 0.05;     % How quickly the bias can change
max_bias = sigma_pos;        % Maximum bias magnitude

% Generate true target trajectory
X_true = zeros(6, N);
X_true(:, 1) = X0;

% Initialize changing bias (starts at 0)
bias_x = zeros(1, N);
bias_y = zeros(1, N);

% Maneuver parameters
maneuver_times = [20, 40];   % When maneuvers occur
maneuver_acc = [5, -3; 3, -4]; % [ax1, ax2; ay1, ay2]

for k = 2:N
    % State transition matrix
    F = [1 0 dt 0 0.5*dt^2 0;
         0 1 0 dt 0 0.5*dt^2;
         0 0 1 0 dt 0;
         0 0 0 1 0 dt;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
    
    % Process noise matrix (discrete time)
    % Simplified calculation for computational efficiency
    G = [dt^2/2 0; 
         0 dt^2/2; 
         dt 0; 
         0 dt; 
         1 0; 
         0 1];
    Q = G * diag([sigma_a^2, sigma_a^2]) * G';
    
    % Add maneuvers
    if time(k-1) >= maneuver_times(1) && time(k-1) < maneuver_times(2)
        X_true(5:6, k-1) = maneuver_acc(:, 1);
    elseif time(k-1) >= maneuver_times(2)
        X_true(5:6, k-1) = maneuver_acc(:, 2);
    end
    
    % Generate process noise
    w = sqrt(Q) * randn(6, 1);
    
    % Update true state
    X_true(:, k) = F * X_true(:, k-1) + w;
    
    % Update bias with random walk
    bias_x(k) = bias_x(k-1) + bias_change_rate * randn * dt;
    bias_y(k) = bias_y(k-1) + bias_change_rate * randn * dt;
    
    % Limit bias magnitude
    bias_x(k) = min(max(bias_x(k), -max_bias), max_bias);
    bias_y(k) = min(max(bias_y(k), -max_bias), max_bias);
end

% Generate noisy measurements (only position)
Z = zeros(2, N);
for k = 1:N
    % Add white noise and bias to position measurements
    Z(:, k) = X_true(1:2, k) + sigma_pos * randn(2, 1) + [bias_x(k); bias_y(k)];
end

% Extended Kalman Filter Implementation
X_ekf = zeros(6, N);
P_ekf = zeros(6, 6, N);

% Initialize EKF state and covariance
X_ekf(:, 1) = [Z(:, 1); 0; 0; 0; 0]; % Initialize with first measurement
P_ekf(:, :, 1) = diag([sigma_pos^2, sigma_pos^2, 100, 100, 10, 10]); % Initial uncertainty

% EKF parameters
% Process noise covariance accounting for unknown maneuvers
Q_ekf = Q * 2; % Increase process noise to account for maneuvers

% Measurement model matrices
H = [1 0 0 0 0 0;
     0 1 0 0 0 0]; % Measurement matrix (only position)
R_ekf = sigma_pos^2 * eye(2) * 2; % Measurement noise (increased to account for bias)

for k = 2:N
    % State transition matrix
    F = [1 0 dt 0 0.5*dt^2 0;
         0 1 0 dt 0 0.5*dt^2;
         0 0 1 0 dt 0;
         0 0 0 1 0 dt;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
    
    % Prediction step
    X_ekf_pred = F * X_ekf(:, k-1);
    P_ekf_pred = F * P_ekf(:, :, k-1) * F' + Q_ekf;
    
    % Update step
    y = Z(:, k) - H * X_ekf_pred; % Innovation
    S = H * P_ekf_pred * H' + R_ekf; % Innovation covariance
    K = P_ekf_pred * H' / S; % Kalman gain
    
    X_ekf(:, k) = X_ekf_pred + K * y;
    P_ekf(:, :, k) = (eye(6) - K * H) * P_ekf_pred;
end

% Alpha-Beta-Gamma Filter Implementation
X_abg = zeros(6, N);
X_abg(:, 1) = [Z(:, 1); 0; 0; 0; 0]; % Initialize with first measurement

% Alpha-Beta-Gamma parameters
% These parameters are tuned for this specific scenario
alpha = 0.7;   % Position gain
beta = 0.4;    % Velocity gain
gamma = 0.1;   % Acceleration gain

for k = 2:N
    % Prediction step
    x_pred = X_abg(1, k-1) + X_abg(3, k-1) * dt + 0.5 * X_abg(5, k-1) * dt^2;
    y_pred = X_abg(2, k-1) + X_abg(4, k-1) * dt + 0.5 * X_abg(6, k-1) * dt^2;
    vx_pred = X_abg(3, k-1) + X_abg(5, k-1) * dt;
    vy_pred = X_abg(4, k-1) + X_abg(6, k-1) * dt;
    ax_pred = X_abg(5, k-1);
    ay_pred = X_abg(6, k-1);
    
    % Calculate residuals
    r_x = Z(1, k) - x_pred;
    r_y = Z(2, k) - y_pred;
    
    % Update step
    X_abg(1, k) = x_pred + alpha * r_x;
    X_abg(2, k) = y_pred + alpha * r_y;
    X_abg(3, k) = vx_pred + (beta * r_x) / dt;
    X_abg(4, k) = vy_pred + (beta * r_y) / dt;
    X_abg(5, k) = ax_pred + (gamma * r_x) / (0.5 * dt^2);
    X_abg(6, k) = ay_pred + (gamma * r_y) / (0.5 * dt^2);
end

% Calculate Errors
% Position errors
pos_error_ekf = sqrt((X_true(1,:) - X_ekf(1,:)).^2 + (X_true(2,:) - X_ekf(2,:)).^2);
pos_error_abg = sqrt((X_true(1,:) - X_abg(1,:)).^2 + (X_true(2,:) - X_abg(2,:)).^2);

% Velocity errors
vel_error_ekf = sqrt((X_true(3,:) - X_ekf(3,:)).^2 + (X_true(4,:) - X_ekf(4,:)).^2);
vel_error_abg = sqrt((X_true(3,:) - X_abg(3,:)).^2 + (X_true(4,:) - X_abg(4,:)).^2);

% Acceleration errors
acc_error_ekf = sqrt((X_true(5,:) - X_ekf(5,:)).^2 + (X_true(6,:) - X_ekf(6,:)).^2);
acc_error_abg = sqrt((X_true(5,:) - X_abg(5,:)).^2 + (X_true(6,:) - X_abg(6,:)).^2);

% Root Mean Square Errors (RMSE)
rmse_pos_ekf = sqrt(mean(pos_error_ekf.^2));
rmse_pos_abg = sqrt(mean(pos_error_abg.^2));
rmse_vel_ekf = sqrt(mean(vel_error_ekf.^2));
rmse_vel_abg = sqrt(mean(vel_error_abg.^2));
rmse_acc_ekf = sqrt(mean(acc_error_ekf.^2));
rmse_acc_abg = sqrt(mean(acc_error_abg.^2));

% Plot Results
% Trajectory plot
figure(1);
plot(X_true(1,:), X_true(2,:), 'k-', 'LineWidth', 2);
hold on;
plot(Z(1,:), Z(2,:), 'r.', 'MarkerSize', 3);
plot(X_ekf(1,:), X_ekf(2,:), 'b-', 'LineWidth', 1.5);
plot(X_abg(1,:), X_abg(2,:), 'g-', 'LineWidth', 1.5);
grid on;
legend('True Trajectory', 'Noisy Measurements', 'EKF Estimate', 'ABG Estimate');
title('Target Trajectory and Estimates');
xlabel('X Position (m)');
ylabel('Y Position (m)');

% Position error plot
figure(2);
subplot(3,1,1);
plot(time, pos_error_ekf, 'b-', time, pos_error_abg, 'g-');
grid on;
title('Position Error');
legend('EKF', 'ABG');
xlabel('Time (s)');
ylabel('Error (m)');

% Velocity error plot
subplot(3,1,2);
plot(time, vel_error_ekf, 'b-', time, vel_error_abg, 'g-');
grid on;
title('Velocity Error');
legend('EKF', 'ABG');
xlabel('Time (s)');
ylabel('Error (m/s)');

% Acceleration error plot
subplot(3,1,3);
plot(time, acc_error_ekf, 'b-', time, acc_error_abg, 'g-');
grid on;
title('Acceleration Error');
legend('EKF', 'ABG');
xlabel('Time (s)');
ylabel('Error (m/s²)');

% Plot measurements with bias
figure(3);
subplot(2,1,1);
plot(time, X_true(1,:), 'k-', time, Z(1,:), 'r.', time, bias_x, 'g-');
grid on;
title('X Position Measurements with Bias');
legend('True Position', 'Noisy Measurements', 'Bias');
xlabel('Time (s)');
ylabel('Position (m)');

subplot(2,1,2);
plot(time, X_true(2,:), 'k-', time, Z(2,:), 'r.', time, bias_y, 'g-');
grid on;
title('Y Position Measurements with Bias');
legend('True Position', 'Noisy Measurements', 'Bias');
xlabel('Time (s)');
ylabel('Position (m)');

% Summary of Results
fprintf('===== Performance Comparison =====\n');
fprintf('Position RMSE (m):\n  EKF: %.2f\n  ABG: %.2f\n', rmse_pos_ekf, rmse_pos_abg);
fprintf('Velocity RMSE (m/s):\n  EKF: %.2f\n  ABG: %.2f\n', rmse_vel_ekf, rmse_vel_abg);
fprintf('Acceleration RMSE (m/s²):\n  EKF: %.2f\n  ABG: %.2f\n', rmse_acc_ekf, rmse_acc_abg);

if rmse_pos_ekf < rmse_pos_abg
    fprintf('\nEKF performs better for position estimation by %.1f%\n', ...
        100*(rmse_pos_abg-rmse_pos_ekf)/rmse_pos_abg);
else
    fprintf('\nABG performs better for position estimation by %.1f%\n', ...
        100*(rmse_pos_ekf-rmse_pos_abg)/rmse_pos_ekf);
end

if rmse_vel_ekf < rmse_vel_abg
    fprintf('EKF performs better for velocity estimation by %.1f%\n', ...
        100*(rmse_vel_abg-rmse_vel_ekf)/rmse_vel_abg);
else
    fprintf('ABG performs better for velocity estimation by %.1f%\n', ...
        100*(rmse_vel_ekf-rmse_vel_abg)/rmse_vel_ekf);
end

if rmse_acc_ekf < rmse_acc_abg
    fprintf('EKF performs better for acceleration estimation by %.1f%\n', ...
        100*(rmse_acc_abg-rmse_acc_ekf)/rmse_acc_abg);
else
    fprintf('ABG performs better for acceleration estimation by %.1f%\n', ...
        100*(rmse_acc_ekf-rmse_acc_abg)/rmse_acc_ekf);
end