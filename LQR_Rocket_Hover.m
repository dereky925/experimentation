%{
**LQR tries to drive states (x) to 0 while minimizing control input (u)**

This self-contained MATLAB script demonstrates how to use LQR (Linear Quadratic 
Regulator) control for a simplified 1D rocket near a hover condition, but now 
we will test multiple Q and R values to see the impact on performance.

-------------------------------------------------------------------------------
SYSTEM DESCRIPTION (Same as before):
-------------------------------------------------------------------------------
1) 1D vertical rocket, hovering around altitude h_des.
2) State: x = [ e_h; e_v ], 
   where e_h = altitude error from h_des, e_v = velocity error from 0.
3) Control input: u = (Thrust - mg).
4) Linearized dynamics:
     dot(e_h) = e_v
     dot(e_v) = (1/m)*u
   => A = [0  1; 0  0], B = [0; 1/m].

-------------------------------------------------------------------------------
MODIFICATION: WE WILL VARY Q & R, THEN PLOT RESULTS ON THE SAME FIGURE.
%}

clear; clc; close all;

% 1) Rocket / Simulation Parameters
m = 1000;         % kg, rocket mass
g = 9.81;         % m/s^2, gravity
h_des = 100;      % desired hover altitude (m)
tspan = [0 60];   % simulate 60 seconds
x0 = [10; 1];     % initial error: e_h=10m above h_des, e_v=0

% System matrices
A = [0  1;
     0  0];
B = [0;
     1/m];

% 2) Define multiple (Q,R) sets to try
% We'll store them in a struct array or cell array
% label -> used in legend
% Q Penalizes deviations from desired state
% R penalizes excessive control effort

combos = {
    struct('Q', diag([10,1]),   'R', 0.1,    'label','Q=[10,1], R=0.1'), 
    struct('Q', diag([100,1]),  'R', 0.1,    'label','Q=[100,1], R=0.1'),
    struct('Q', diag([1000,1]),  'R', 0.1,      'label','Q=[1000,1], R=0.1')
};

% We'll store solutions for each (Q,R) set
results = cell(numel(combos),1);

% 3) Loop over each (Q,R) set, solve LQR, simulate, store
for k = 1:numel(combos)
    thisQ = combos{k}.Q;
    thisR = combos{k}.R;

    % Solve for LQR gain
    % Solving the Riccati equation by hand is very hard for larger than a 2x2
    % Recommend just using lqr()
    K = lqr(A,B,thisQ,thisR);

    % Simulate closed-loop
    [tSol, xSol] = ode45(@(t,x) rocket_closedloop(t, x, A, B, K), tspan, x0);

    % Extract e_h, e_v
    e_h_sol = xSol(:,1);
    e_v_sol = xSol(:,2);

    % Convert to actual altitude & velocity
    hSol = h_des + e_h_sol;
    vSol = e_v_sol;

    % Store in our cell array for plotting
    results{k} = struct( ...
        't', tSol, ...
        'hSol', hSol, ...
        'vSol', vSol, ...
        'e_h', e_h_sol, ...
        'e_v', e_v_sol, ...
        'label', combos{k}.label ...
    );
end

% 4) Plot: Actual Altitude & Velocity in one figure
figure('Name','Altitude & Velocity','Position',[0 50 1200 800]);

% Altitude
subplot(2,1,1); hold on; grid minor;
for k=1:numel(results)
    plot(results{k}.t, results{k}.hSol, 'LineWidth',2, 'DisplayName', results{k}.label);
end
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Rocket Altitude vs Time (Multiple Q,R)');
legend('Location','best');

% Velocity
subplot(2,1,2); hold on; grid minor;
for k=1:numel(results)
    plot(results{k}.t, results{k}.vSol, 'LineWidth',2, 'DisplayName', results{k}.label);
end
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('Rocket Velocity vs Time (Multiple Q,R)');
legend('Location','best');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

% 5) Plot: Altitude Error & Velocity Error in another figure
figure('Name','Error States','Position',[1500 50 1200 800]);

% Altitude Error
subplot(2,1,1); hold on; grid minor;
for k=1:numel(results)
    plot(results{k}.t, results{k}.e_h, 'LineWidth',2, 'DisplayName', results{k}.label);
end
xlabel('Time (s)'); ylabel('Altitude Error (m)');
title('Altitude Error vs Time (Multiple Q,R)');
legend('Location','best');

% Velocity Error
subplot(2,1,2); hold on; grid minor;
for k=1:numel(results)
    plot(results{k}.t, results{k}.e_v, 'LineWidth',2, 'DisplayName', results{k}.label);
end
xlabel('Time (s)'); ylabel('Velocity Error (m/s)');
title('Velocity Error vs Time (Multiple Q,R)');
legend('Location','best');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);


%{
 rocket_closedloop:
 Defines the closed-loop dynamics: dot(x) = (A - B*K)*x
 x=[e_h; e_v], control u=-K*x 
 => no need to explicitly compute u
%}
function dxdt = rocket_closedloop(~, x, A, B, K)
    dxdt = (A - B*K)*x;
end