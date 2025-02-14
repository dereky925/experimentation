%% P Controller
clear;clc;close all

% Define time and reference signal
t = 0:0.01:10;        % Time vector
ref = ones(size(t));  % Reference signal (step input)

% Define the plant (simple first-order system)
% Plant: G(s) = 1 / (s + 1)
num = 1;
den = [1 1];
plant = tf(num, den);

% Define Proportional controller
Kp = 2;  % Proportional gain
controller = Kp;  % Proportional control only

% Closed-loop system
CL = feedback(controller*plant, 1);

% Simulate
[y, t, x] = lsim(CL, ref, t);

% Plot response
figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
plot(t, ref, 'r--', 'LineWidth', 2); 
plot(t, y, 'b-', 'LineWidth', 2);
title('Proportional Control (P)');
xlabel('Time (s)');
ylabel('Output');
legend('Reference', 'Output');
ax=gca;
ax.FontSize = 20;

%% PI Controller
clear;clc; close all

% Define time and reference signal
t = 0:0.01:10;        % Time vector
ref = ones(size(t));  % Reference signal (step input)

% Define the plant (simple first-order system)
% Plant: G(s) = 1 / (s + 1)
num = 1;
den = [1 1];
plant = tf(num, den);

% Define PI controller: C(s) = Kp + Ki/s
Kp = 2;
Ki = 1;  % Integral gain
controller = Kp + tf([Ki], [1 0]);  % Kp + Ki/s

% Closed-loop system
CL = feedback(controller*plant, 1);

% Simulate
[y, t, x] = lsim(CL, ref, t);

% Plot response
figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
plot(t, ref, 'r--', 'LineWidth',2);
plot(t, y, 'b-', 'LineWidth', 2);
title('Proportional-Integral Control (PI)');
xlabel('Time (s)');
ylabel('Output');
legend('Reference', 'Output');
grid on;
ax=gca;
ax.FontSize = 20;

%% PID Controller

% Define the plant (simple first-order system)
% Plant: G(s) = 1 / (s + 1)
num = 1;
den = [1 1];
plant = tf(num, den);

% Define PID controller: C(s) = Kp + Ki/s + Kd*s
Kp = 2;  % Proportional gain
Ki = 1;  % Integral gain
Kd = 0.5;  % Derivative gain
controller = Kp + tf([Ki], [1 0]) + tf([Kd 0], [1]);  % PID: Kp + Ki/s + Kd*s

% Closed-loop system
CL = feedback(controller*plant, 1);

% Simulate
[y, t, x] = lsim(CL, ref, t);

% Plot response
figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
plot(t, ref, 'r--', 'LineWidth', 2); hold on;
plot(t, y, 'b-', 'LineWidth', 2);
title('Proportional-Integral-Derivative Control (PID)');
xlabel('Time (s)');
ylabel('Output');
legend('Reference', 'Output');
ax=gca;
ax.FontSize = 20;

%% P Control with steady state error
clear;clc;close all

% Define a new plant: G(s) = 1 / (2s + 1)
num = 1;
den = [2 1];
plant = tf(num, den);

% Proportional controller only
Kp = 1;
controller = Kp;

% Closed-loop system
CL = feedback(controller*plant, 1);

% Simulate
[y, t, x] = lsim(CL, ref, t);

% Plot response
figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
plot(t, ref, 'r--', 'LineWidth', 2);
plot(t, y, 'b-', 'LineWidth', 2);
title('Proportional Control: Steady-State Error Example');
xlabel('Time (s)');
ylabel('Output');
legend('Reference', 'Output');
grid on;
ax=gca;
ax.FontSize = 20;

%% Controller comparison
clear;clc;close all

% Define the original plant
num = 1;
den = [1 1];
plant = tf(num, den);

% Define controllers
Kp = 2;
Ki = 1;
Kd = 0.5;
P_controller = Kp;
PI_controller = Kp + tf([Ki], [1 0]);
PID_controller = Kp + tf([Ki], [1 0]) + tf([Kd 0], [1]);

% Closed-loop systems
P_CL = feedback(P_controller*plant, 1);
PI_CL = feedback(PI_controller*plant, 1);
PID_CL = feedback(PID_controller*plant, 1);

% Simulate
[yP, ~] = lsim(P_CL, ref, t);
[yPI, ~] = lsim(PI_CL, ref, t);
[yPID, ~] = lsim(PID_CL, ref, t);

% Plot all responses
figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
plot(t, ref, 'r--', 'LineWidth', 2); 
plot(t, yP, 'b-', 'LineWidth', 2, 'DisplayName', 'P Controller');
plot(t, yPI, 'g-', 'LineWidth', 2, 'DisplayName', 'PI Controller');
plot(t, yPID, 'm-', 'LineWidth', 2, 'DisplayName', 'PID Controller');
title('Comparison of P, PI, and PID Controllers');
xlabel('Time (s)');
ylabel('Output');
legend('Reference', 'P Controller', 'PI Controller', 'PID Controller');
grid on;
ax=gca;
ax.FontSize = 20;

