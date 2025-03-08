clear; close all; clc;
% Define simulation parameters
fs = 300;             % Sampling rate in Hz
T = 1/fs;             % Sampling period (seconds)
t = 0:T:2;            % Time vector from 0 to 2 seconds
f = 1;                % Frequency of the sine wave in Hz

% Generate the true position signal as a sine wave and its derivative (true velocity)
pos_true = sin(2*pi*f*t);               % True position signal (sine wave)
vel_true = 2*pi*f*cos(2*pi*f*t);          % True velocity (analytical derivative of sine)

% Add white Gaussian noise to simulate measurement noise
noise_amplitude = 0.2;                   
pos_noisy = pos_true + noise_amplitude*randn(size(t));  % Noisy position signal

% -------------------------------------------------------------------------
% Manual Savitzky-Golay Derivative Estimation (using local polynomial fitting)
% This algorithm performs a local polynomial fit over a sliding window
% and uses the first-order coefficient as an estimate of the derivative.
order = 3;             % Order of the polynomial (3rd order)
framelen = 21;         % Frame length (number of points in the sliding window; must be odd)
half_window = (framelen-1)/2;  % Number of points on each side of the center point
vel_sg = zeros(size(pos_noisy));  % Preallocate array for estimated velocity

for n = half_window+1 : length(pos_noisy)-half_window
    % Create a window of time values centered at zero for polynomial fitting.
    % Multiply by T to convert index differences to time differences.
    x_window = ((-half_window:half_window)*T)'; 
    
    % Extract the segment of the noisy signal corresponding to the current window.
    y_window = pos_noisy(n-half_window:n+half_window)';
    
    % Construct the design (Vandermonde) matrix for a 3rd order polynomial:
    % The columns represent 1, x, x^2, and x^3.
    A = [ones(framelen,1), x_window, x_window.^2, x_window.^3];
    
    % Solve the least squares problem A*a = y_window to obtain polynomial coefficients.
    % The solution vector "a" contains the coefficients [a0; a1; a2; a3].
    a = A\y_window;  
    
    % The derivative at the center of the window (x=0) is the first derivative,
    % which is simply the coefficient a(2).
    vel_sg(n) = a(2);  
end

% -------------------------------------------------------------------------
% Manual implementation of a 2nd order Butterworth Low-Pass Filter using the bilinear transform
% This section converts an analog Butterworth filter to a digital filter.
% The filter is applied in a causal (forward-only) manner, which is suitable
% for real-time processing where new data arrives sample-by-sample.

% Set the desired cutoff frequency (in Hz) for the Butterworth filter.
cutoff = 1.5;                     

% Pre-warp the cutoff frequency using the tangent function to compensate
% for the non-linear frequency mapping introduced by the bilinear transform.
K = tan(pi*cutoff/fs);

% Compute the normalization factor. This factor normalizes the filter coefficients
% so that the digital filter's gain at DC is unity.
norm_factor = 1 + sqrt(2)*K + K^2;

% Calculate the numerator coefficients (b coefficients) for the filter.
% These coefficients originate from the analog Butterworth transfer function.
b0 = K^2 / norm_factor;  % Coefficient for the current input sample.
b1 = 2*b0;               % Coefficient for the previous input sample (2*b0).
b2 = b0;                 % Coefficient for the input sample before the previous one (same as b0).

% Calculate the denominator coefficients (a coefficients) for the filter.
% In the difference equation, a0 is normalized to 1.
a1 = 2*(K^2 - 1) / norm_factor;    % Coefficient for the first previous output sample.
a2 = (1 - sqrt(2)*K + K^2) / norm_factor;  % Coefficient for the second previous output sample.

% Pre-allocate the output vector for the filtered signal.
pos_butter = zeros(size(pos_noisy));

% Apply the digital Butterworth filter in a causal manner.
% For each sample, use the current and past input values and past output values.
for n = 1:length(pos_noisy)
    if n == 1
        % For the very first sample, there are no previous inputs or outputs.
        pos_butter(n) = b0*pos_noisy(n);
    elseif n == 2
        % For the second sample, only one previous input and one previous output are available.
        pos_butter(n) = b0*pos_noisy(n) + b1*pos_noisy(n-1) - a1*pos_butter(n-1);
    else
        % For subsequent samples, use the current sample, the two previous inputs,
        % and the two previous outputs.
        pos_butter(n) = b0*pos_noisy(n) + b1*pos_noisy(n-1) + b2*pos_noisy(n-2) ...
                        - a1*pos_butter(n-1) - a2*pos_butter(n-2);
    end
end

% Compute the derivative (velocity) from the filtered position signal using finite difference.
vel_butter = diff(pos_butter) / T;
% Append the last value to keep the velocity vector the same length as the input.
vel_butter = [vel_butter, vel_butter(end)];

% -------------------------------------------------------------------------
% Manual First Order Filter Implementation (RC Low-Pass Filter)
% This filter mimics a simple RC circuit where the output is computed as a weighted
% combination of the current input and the previous output.
cutoff = 0.5;           % Set a new cutoff frequency for the first order filter (in Hz)
RC = 1/(2*pi*cutoff);   % Calculate the RC time constant from the cutoff frequency
alpha_RC = T/(RC+T);    % Compute the filter coefficient (alpha), which determines the
                        % weight of the current input relative to the previous output.
pos_first = zeros(size(pos_noisy));  % Preallocate the filtered output vector.

for n = 1:length(pos_noisy)
    if n == 1
        % For the first sample, simply initialize with the first noisy measurement.
        pos_first(n) = pos_noisy(n);
    else
        % For subsequent samples, compute the filtered output as:
        % pos_first(n) = (alpha_RC * current input) + ((1 - alpha_RC) * previous output)
        pos_first(n) = alpha_RC*pos_noisy(n) + (1-alpha_RC)*pos_first(n-1);
    end
end

% Compute the velocity from the first order filtered signal using finite difference.
vel_first = diff(pos_first)/T;
vel_first = [vel_first, vel_first(end)];

% -------------------------------------------------------------------------
% Manual Alpha-Beta Filter Implementation
% The alpha-beta filter is a simple predictive filter that estimates both position and velocity.
% It uses a prediction step (based on the previous state) and a correction step (based on the new measurement).
alpha_ab = 0.85;  % Gain factor for the position correction
beta_ab = 0.01;  % Gain factor for the velocity correction

xhat = zeros(size(pos_noisy));  % Estimated (filtered) position
vhat = zeros(size(pos_noisy));  % Estimated velocity

% Initialize the filter with the first measurement and assume an initial velocity of zero.
xhat(1) = pos_noisy(1);
vhat(1) = 0;

for k = 2:length(pos_noisy)
    % Prediction step:
    % Predict the new position using the previous position and velocity.
    xhat_pred = xhat(k-1) + T*vhat(k-1);
    
    % Compute the residual (difference) between the actual measurement and the predicted position.
    residual = pos_noisy(k) - xhat_pred;
    
    % Correction step:
    % Update the position estimate by adding a fraction (alpha_ab) of the residual.
    xhat(k) = xhat_pred + alpha_ab*residual;
    
    % Update the velocity estimate by adding a term proportional to the residual.
    % The factor beta_ab/T adjusts the correction based on the time step.
    vhat(k) = vhat(k-1) + (beta_ab/T)*residual;
end

% -------------------------------------------------------------------------
% Plotting the results in a 1600x1600 pixels figure with 3 rows and 2 columns of subplots.
% The first subplot (spanning the entire first row) shows the position signals,
% and the remaining four subplots show the velocity estimates from different filtering algorithms.

figure('Position',[100,100,1600,1600],'Color','w');

% Plot the position signal (spanning the first row across two subplot spaces).
subplot(3,2,[1,2])
plot(t, pos_noisy, 'b', 'LineWidth', 1.5); hold on;
plot(t, pos_true, 'r', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize',20);
ylabel('Position', 'FontSize',20);
legend('Noisy','True','FontSize',20);
set(gca,'FontSize',20);
grid on;

% Plot the velocity estimated by the Savitzky-Golay filter.
subplot(3,2,3)
plot(t, vel_true, 'r', 'LineWidth', 2); hold on;
plot(t, vel_sg, 'k', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Velocity', 'FontSize',20);
legend('True','SG','FontSize',20);
set(gca,'FontSize',20);
grid on;

% Plot the velocity estimated by the Butterworth filter.
subplot(3,2,4)
plot(t, vel_true, 'r', 'LineWidth', 2); hold on;
plot(t, vel_butter, 'g', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Velocity', 'FontSize',20);
legend('True','Butterworth','FontSize',20);
set(gca,'FontSize',20);
grid on;

% Plot the velocity estimated by the First Order filter.
subplot(3,2,5)
plot(t, vel_true, 'r', 'LineWidth', 2); hold on;
plot(t, vel_first, 'm', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Velocity', 'FontSize',20);
legend('True','First Order','FontSize',20);
set(gca,'FontSize',20);
grid on;

% Plot the velocity estimated by the Alpha-Beta filter.
subplot(3,2,6)
plot(t, vel_true, 'r', 'LineWidth', 2); hold on;
plot(t, vhat, 'c', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Velocity', 'FontSize',20);
legend('True','Alpha-Beta','FontSize',20);
set(gca,'FontSize',20);
grid on;