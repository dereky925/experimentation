clc; clear; close all;

% Non-Linear Time Invariant systems can be represented with state space
% dot = Ax + Bu

% Non-LTI just means system is non-linear (ex. drag = F_drag = 0.5 * C_d * A * rho * v*abs(v);)
% Or the A or B matrix is not constant (ex. mass changes throughout flight)

% Methods like transfer functions and bode plots do not work

% Define constants
g       = 9.81;            % m/s^2, gravity
rho     = 1.225;           % kg/m^3, sea-level air density
C_d     = 0.3;             % Drag coefficient
A       = pi*(5^2);        % Reference area (say, 10 m diameter rocket)
T       = 25e6;             % Thrust (N), constant until burnout
mdot    = 500;             % Fuel burn rate (kg/s)
t_burn  = 120;             % Burn time (s)

% Initial conditions
h0 = 0;        % altitude (m)
v0 = 0;        % velocity (m/s)
m0 = 2.0e5;    % total mass (kg) at liftoff (including fuel)

% State vector: x = [h; v; m]
x0 = [h0; v0; m0];

% Simulation time
tspan = [0 300];  % 5 minutes total

% ODE solver
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t, X] = ode45(@(t,x) rocketDynamics(t, x, g, rho, C_d, A, ...
                                     T, mdot, t_burn), ...
               tspan, x0, options);

% Extract solutions
h = X(:,1);
v = X(:,2);
m = X(:,3);

% Plot altitude, velocity, mass vs. time
figure('Position',[50,50,1200,600],'Color','white');
subplot(3,1,1)
plot(t, h/1000, 'LineWidth',2);
xlabel('Time (s)'); ylabel('Altitude (km)'); grid on;
title('Rocket Altitude vs. Time');

subplot(3,1,2)
plot(t, v, 'LineWidth',2);
xlabel('Time (s)'); ylabel('Velocity (m/s)'); grid on;
title('Rocket Velocity vs. Time');

subplot(3,1,3)
plot(t, m/1000, 'LineWidth',2);
xlabel('Time (s)'); ylabel('Mass (tonnes)'); grid on;
title('Rocket Mass vs. Time');


function dxdt = rocketDynamics(t, x, g, rho, C_d, A, Thrust, mdot, t_burn)
    % x = [h; v; m]
    h = x(1);
    v = x(2);
    m = x(3);
    
    % Thrust schedule
    if t < t_burn
        T = Thrust;
        fuelFlow = mdot;
    else
        T = 0;
        fuelFlow = 0;
    end
    
    % Drag force (nonlinear in velocity)
    F_drag = 0.5 * C_d * A * rho * v*abs(v);
    
    % Equations of motion
    dhdt = v;
    dvdt = (T - F_drag) / m - g;   % Net acceleration
    dmdt = -fuelFlow;              % Mass is decreasing
    
    dxdt = [dhdt; dvdt; dmdt];
end