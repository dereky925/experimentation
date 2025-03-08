% Simulation Parameters
clear; clc;close all
dt = 0.005;              % Time step (s) ~200 Hz update rate
T = 10;                  % Total simulation time (s)
N = floor(T/dt);         % Number of simulation steps
time = (0:N-1)' * dt;    % Time vector

% True Target Trajectory (1D motion)
% Define a maneuvering target with piecewise acceleration:
% - 0 <= t < 2 s: steady motion (acc = 0)
% - 2 <= t < 3 s: rapid acceleration (acc = 20 m/s^2)
% - 3 <= t < 5 s: constant velocity (acc = 0)
% - 5 <= t < 6 s: deceleration (acc = -15 m/s^2)
% - 6 <= t <= 10 s: oscillatory maneuver
acc_true = zeros(N,1);
for k = 1:N
    t = time(k);
    if t < 2
        acc_true(k) = 0;
    elseif t < 3
        acc_true(k) = 20;
    elseif t < 5
        acc_true(k) = 0;
    elseif t < 6
        acc_true(k) = -15;
    else
        acc_true(k) = 5 * sin(2*pi*0.5*t);
    end
end

% Integrate acceleration to get velocity and position
vel_true = zeros(N,1);
pos_true = zeros(N,1);
for k = 2:N
    vel_true(k) = vel_true(k-1) + acc_true(k-1)*dt;
    pos_true(k) = pos_true(k-1) + vel_true(k-1)*dt + 0.5*acc_true(k-1)*dt^2;
end

% Simulated Measurements
% Measurement model: z = true position + bias + noise
% Bias is modeled as a linearly drifting offset.
b0 = 0;                      % Initial bias
bias_slope = 0.05;           % Bias drift (m/s)
bias_true = b0 + bias_slope * time;  % True bias over time

% White noise: Assume sigma = 0.1 m (adjust based on sensor specs)
sigma_noise = 0.1;
noise = sigma_noise * randn(N,1);

% Position measurements
z = pos_true + bias_true + noise;

% Extended Kalman Filter (EKF) Implementation
% Augmented state: x = [position; velocity; acceleration; bias]
% State transition (constant acceleration, bias assumed nearly constant):
%   x(k+1) = F*x(k) + w
F = [1 dt 0.5*dt^2 0;
     0 1  dt       0;
     0 0  1        0;
     0 0  0        1];
 
% Process noise covariance (allow changes in acceleration and bias)
q_acc = 10;        % Process noise variance for acceleration
q_bias = 0.001;    % Process noise variance for bias
Q = diag([0, 0, q_acc, q_bias]);

% Measurement model: z = H*x + v, where H = [1 0 0 1]
H = [1 0 0 1];
R = sigma_noise^2; % Measurement noise variance

% Initialization
ekf_est = zeros(4, N);
ekf_est(:,1) = [z(1); 0; 0; 0];  % Initialize using first measurement; others 0
P_k = eye(4);                   % Initial state covariance

% Run the EKF
for k = 2:N
    % Prediction step
    x_pred = F * ekf_est(:,k-1);
    P_pred = F * P_k * F' + Q;
    
    % Measurement update
    z_pred = H * x_pred;
    y_k = z(k) - z_pred;  % Innovation (residual)
    S = H * P_pred * H' + R;  % Innovation covariance
    K = P_pred * H' / S;      % Kalman gain
    
    % Update state and covariance
    x_upd = x_pred + K * y_k;
    P_k = (eye(4) - K * H) * P_pred;
    
    ekf_est(:,k) = x_upd;
end

% Alpha-Beta-Gamma (α–β–γ) Filter Implementation
% State: [position; velocity; acceleration]
% This filter does NOT account for bias.
abg_est = zeros(3, N);
abg_est(:,1) = [z(1); 0; 0];

% Tuning parameters (adjust based on performance needs)
alpha = 0.85;
beta = 0.005;
gamma = 0.0001;

% Run the α–β–γ filter
for k = 2:N
    % Prediction step (assume constant acceleration)
    x_pred_abg = abg_est(1,k-1) + abg_est(2,k-1)*dt + 0.5 * abg_est(3,k-1)*dt^2;
    v_pred_abg = abg_est(2,k-1) + abg_est(3,k-1)*dt;
    a_pred_abg = abg_est(3,k-1);
    
    % Compute measurement residual
    r = z(k) - x_pred_abg;
    
    % Update step
    abg_est(1,k) = x_pred_abg + alpha * r;
    abg_est(2,k) = v_pred_abg + (beta/dt) * r;
    abg_est(3,k) = a_pred_abg + (2*gamma/dt^2) * r;
end

% Performance Evaluation: Compute RMSE for each state
rmse_ekf = zeros(3,1);
rmse_abg = zeros(3,1);
% EKF (ignore the bias state for RMSE of physical states)
rmse_ekf(1) = sqrt(mean((ekf_est(1,:) - pos_true').^2));
rmse_ekf(2) = sqrt(mean((ekf_est(2,:) - vel_true').^2));
rmse_ekf(3) = sqrt(mean((ekf_est(3,:) - acc_true').^2));
% α–β–γ filter
rmse_abg(1) = sqrt(mean((abg_est(1,:) - pos_true').^2));
rmse_abg(2) = sqrt(mean((abg_est(2,:) - vel_true').^2));
rmse_abg(3) = sqrt(mean((abg_est(3,:) - acc_true').^2));

fprintf('EKF RMSE:\n  Position = %f\n  Velocity = %f\n  Acceleration = %f\n', rmse_ekf(1), rmse_ekf(2), rmse_ekf(3));
fprintf('Alpha-Beta-Gamma RMSE:\n  Position = %f\n  Velocity = %f\n  Acceleration = %f\n', rmse_abg(1), rmse_abg(2), rmse_abg(3));

% Plot Results

% Plot position
figure('Color','white','Position',[0 0 1000 1000]);
subplot(3,1,1);
hold on; grid minor;
plot(time, pos_true, 'k', 'LineWidth',1.5); hold on;
plot(time, ekf_est(1,:), 'b--', 'LineWidth',1);
plot(time, abg_est(1,:), 'r:', 'LineWidth',1);
title('Position');
legend('True', 'EKF', 'α–β–γ');
xlabel('Time (s)'); ylabel('Position (m)');
ax = gca;
ax.FontSize = 20;

% Plot velocity
subplot(3,1,2);
hold on; grid minor;
plot(time, vel_true, 'k', 'LineWidth',1.5); hold on;
plot(time, ekf_est(2,:), 'b--', 'LineWidth',1);
plot(time, abg_est(2,:), 'r:', 'LineWidth',1);
title('Velocity');
legend('True', 'EKF', 'α–β–γ');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
ax = gca;
ax.FontSize = 20;

% Plot acceleration
subplot(3,1,3);
hold on; grid minor;
plot(time, acc_true, 'k', 'LineWidth',1.5); hold on;
plot(time, ekf_est(3,:), 'b--', 'LineWidth',1);
plot(time, abg_est(3,:), 'r:', 'LineWidth',1);
title('Acceleration');
legend('True', 'EKF', 'α–β–γ');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
ax = gca;
ax.FontSize = 20;

% Plot bias (EKF bias estimate vs. true bias)
figure('Color','white','Position',[2000 0 1000 1000]);
hold on; grid minor;
plot(time, bias_true, 'k', 'LineWidth',1.5); hold on;
plot(time, ekf_est(4,:), 'b--', 'LineWidth',1);
title('Measurement Bias');
legend('True Bias', 'EKF Bias Estimate');
xlabel('Time (s)'); ylabel('Bias (m)');