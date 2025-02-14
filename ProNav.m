% PN Guidance Simulation for Missile Interception with Configurable Target Maneuver
% This script simulates a missile intercepting a target using Proportional
% Navigation (PN) guidance. The missile is assumed to have constant speed,
% and its lateral acceleration is computed based on the LOS rate. The PN
% law commands an acceleration proportional to the rate of change of the
% line-of-sight (LOS) angle between the missile and target.
%
% The simulation works as follows:
% 1. The missile and target positions are initialized. The missile starts
%    at the origin with an initial heading, while the target is set to move
%    with a constant velocity.
% 2. At every time step, the LOS angle (the angle between the missile-to-target
%    line and the horizontal) is calculated. The LOS rate is estimated using
%    finite differences.
% 3. Using the PN law (a_c = N*V_m*los_rate), the missile computes the required
%    lateral acceleration command.
% 4. The missile’s heading is updated based on this acceleration (with its speed
%    kept constant), and its new position is calculated.
% 5. Meanwhile, the target moves along a predetermined path and performs a
%    configurable maneuver (e.g., a constant turn) during a specified time window.
% 6. The simulation logs key parameters and plots the missile and target trajectories,
%    LOS angle, LOS rate, and commanded acceleration.
%
% Additionally, the start positions of the interceptor and target are marked 
% (blue star for interceptor start, red star for target start), and the computed 
% intercept point (when miss distance < 0.5 m) is marked with a yellow star having 
% black edges.
%
% This script provides hands-on insight into the workings of PN guidance in missile guidance,
% and allows you to study the effect of a maneuvering target.

clear; close all; clc;

% Simulation Parameters
dt = 0.001;           % Time step [s]
T_end = 120;         % Simulation end time [s]
time = 0:dt:T_end;   % Time vector

% Missile Parameters
V_m = 400;           % Missile speed [m/s]
N = 3;               % Navigation constant

% Initial missile state: position and heading
missile.pos = [0; 0];         % [x; y] initial position [m]
missile.theta = deg2rad(45);  % Initial heading angle [rad]

% Target Parameters
V_t = 250;                   % Target speed [m/s]
target.heading = deg2rad(0); % Target initial heading [rad] (moving horizontally)
target.pos = [5000; 5000];   % Initial target position [m]

% Configurable Target Maneuver Parameters
target.maneuver = true;         % Set to true to enable maneuver
target.maneuver_start = 30;     % Time (s) when maneuver begins
target.maneuver_end = 60;       % Time (s) when maneuver ends
target.turn_rate = deg2rad(-5);  % Target turn rate during maneuver [rad/s]

% Preallocate arrays for logging data
numSteps = length(time);
missile_pos_log = zeros(2, numSteps);
target_pos_log = zeros(2, numSteps);
los_log         = zeros(1, numSteps);
los_rate_log    = zeros(1, numSteps);
acc_cmd_log     = zeros(1, numSteps);
missile_theta_log = zeros(1, numSteps);

% Main Simulation Loop
prev_los = 0; % Initialize previous LOS angle

for k = 1:numSteps
    % Log current states
    missile_pos_log(:, k) = missile.pos;
    target_pos_log(:, k)  = target.pos;
    missile_theta_log(k)  = missile.theta;
    
    % Compute the line-of-sight (LOS) angle from missile to target
    relative_pos = target.pos - missile.pos;
    los = atan2(relative_pos(2), relative_pos(1));
    los_log(k) = los;
    
    % Estimate LOS rate using finite differences (with angle wrapping)
    if k == 1
        los_rate = 0;  % No previous value at t=0
    else
        dlos = wrapToPi(los - prev_los);  % Ensure continuity in angle
        los_rate = dlos / dt;
    end
    los_rate_log(k) = los_rate;
    
    % PN Guidance Law: Compute commanded lateral acceleration
    a_cmd = N * V_m * los_rate;
    acc_cmd_log(k) = a_cmd;
    
    % Update missile state:
    % The missile's heading rate is given by theta_dot = a_cmd/V_m since speed is constant.
    missile.theta = missile.theta + (a_cmd / V_m) * dt;
    
    % Update missile position based on new heading
    missile.pos = missile.pos + V_m * [cos(missile.theta); sin(missile.theta)] * dt;
    
    % Update target state:
    % If target maneuver is enabled and current time is within the maneuver window,
    % update the target heading based on the turn rate.
    if target.maneuver && (time(k) >= target.maneuver_start) && (time(k) <= target.maneuver_end)
        target.heading = target.heading + target.turn_rate * dt;
    end
    % Update target position with its (possibly updated) heading
    target.pos = target.pos + V_t * [cos(target.heading); sin(target.heading)] * dt;
    
    % Update previous LOS angle for next iteration
    prev_los = los;
    
    % Check intercept condition: if miss distance < 0.5 m, exit the loop
    miss_distance = norm(missile.pos - target.pos);
    if miss_distance < 0.5
        intercept_index = k;
        break;
    end
end

% If simulation ends without intercept, set intercept index to last index
if ~exist('intercept_index','var')
    intercept_index = numSteps;
end
intercept_point = missile_pos_log(:, intercept_index);

% Plot missile and target trajectories with markers for start positions and intercept point
figure('Position',[0 0 1000 1000],'Color','white');
plot(missile_pos_log(1, 1:intercept_index), missile_pos_log(2, 1:intercept_index), 'b-', 'LineWidth', 2);
hold on;
plot(target_pos_log(1, 1:intercept_index), target_pos_log(2, 1:intercept_index), 'r--', 'LineWidth', 2);
ax = gca;
ax.FontSize = 20;

% Mark starting positions
plot(missile_pos_log(1,1), missile_pos_log(2,1), 'b*', 'MarkerSize', 10, 'LineWidth', 2); % Interceptor start
plot(target_pos_log(1,1), target_pos_log(2,1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);   % Target start

% Mark interception point with a yellow star with black edges
plot(intercept_point(1), intercept_point(2), 'p', 'MarkerSize', 20, ...
    'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

xlabel('X Position [m]');
ylabel('Y Position [m]');
title('Missile and Target Trajectories using PN Guidance');
legend('Missile Path', 'Target Path', 'Interceptor Start', 'Target Start', 'Intercept Point');
grid on;
axis equal;

% Plot LOS Angle and LOS Rate over time
figure('Position',[1000 0 1000 1000],'Color','white');

subplot(2,1,1);
plot(time(1:intercept_index), rad2deg(los_log(1:intercept_index)), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('LOS Angle [deg]');
title('Line-of-Sight Angle vs. Time');
grid on;
ax = gca;
ax.FontSize = 20;

subplot(2,1,2);
plot(time(1:intercept_index), rad2deg(los_rate_log(1:intercept_index)), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('LOS Rate [deg/s]');
title('Line-of-Sight Rate vs. Time');
grid on;
ax = gca;
ax.FontSize = 20;

% Plot commanded acceleration over time
figure('Position',[2000 0 1000 1000],'Color','white');
plot(time(1:intercept_index), acc_cmd_log(1:intercept_index), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Commanded Acceleration [m/s^2]');
title('Commanded Lateral Acceleration vs. Time');
grid on;
ax = gca;
ax.FontSize = 20;



