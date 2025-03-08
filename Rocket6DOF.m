% Clear workspace and close figures
% https://en.wikipedia.org/wiki/Starship_flight_test_6
clear; clc; close all;

% --------------------- Time & Integration Setup ---------------------
dt = 0.01;            % time step [s]
tf = 15*60;           % final time [s] (10 minutes)
time = 0:dt:tf;       % time vector
N = length(time);

% --------------------- Stage-Specific Parameters ---------------------
% --- Booster (Stage 1) ---
stage1.m_dry   = 2.0e5;        % 200,000 kg dry mass
stage1.m_prop  = 2.3e6;        % 2,300,000 kg propellant
stage1.m0      = stage1.m_dry + stage1.m_prop;  % total 2.5e6
stage1.T_mag   = 87.5e6;       % [N] engine thrust (booster)
stage1.Isp     = 330;          % [s] specific impulse (sea-level)
stage1.t_burn  = 150;          % [s] burn time
stage1.g0      = 9.81;         % [m/s^2]
stage1.A_ref   = pi*(9/2)^2;   % [m^2], diameter = 9 m
stage1.Cd      = 0.3;          % drag coefficient

% --- Ship (Stage 2) ---
stage2.m0     = 2.5e6;         % [kg] stage-2 mass at booster separation
stage2.T_mag  = 12.25e6;       % [N] engine thrust (ship)
stage2.Isp    = 350;           % [s] vacuum Isp
stage2.t_burn = 8*60;           % [s] burn time (ship)
stage2.g0     = 9.81;          % [m/s^2]
stage2.A_ref  = stage1.A_ref;  % reuse same reference area
stage2.Cd     = stage1.Cd;     

% --------------------- Stage Separation Time ---------------------
stageSepTime = 150;  % [s]

% --------------------- Combined (Pre-Separation) Parameters ---------------------
paramsCombined.g        = stage1.g0;
paramsCombined.burnTime = stage1.t_burn;  
paramsCombined.T_engine = stage1.T_mag;
paramsCombined.mdot     = stage1.m_prop / stage1.t_burn; 
paramsCombined.Cd       = stage1.Cd;
paramsCombined.A        = stage1.A_ref;
paramsCombined.C_damp   = 1000;

I_booster = diag([1e6, 2e6, 3e6]);
I_ship    = diag([1e6, 2e6, 3e6]);
paramsCombined.I        = I_booster + I_ship;
paramsCombined.rho0     = 1.225;  
paramsCombined.H        = 8000;   

% -- Multiple Maneuvers for the Combined Stage --
%   We define arrays of pitch intervals and thruster intervals.
%   You can add as many as you like.

paramsCombined.maneuvers.pitch = [
    struct('start', 30,  'end', 40,  'M',  5000),   % big pitch moment
    % struct('start', 100, 'end', 110, 'M',  -3000),  % negative pitch moment
];

% Example thruster intervals for the combined rocket:
paramsCombined.maneuvers.thruster = [
    % struct('start', 60, 'end', 90, 'T', 5e7)  % override thrust for 30s
    % You can add more intervals if desired
];

% --------------------- Ship (Post-Separation) Parameters ---------------------
paramsShip.g        = stage2.g0;
paramsShip.burnTime = stage2.t_burn;
paramsShip.T_engine = stage2.T_mag;
paramsShip.mdot     = stage2.m0 / stage2.t_burn;
paramsShip.Cd       = stage2.Cd;
paramsShip.A        = stage2.A_ref;
paramsShip.C_damp   = 1000;
paramsShip.I        = I_ship;  
paramsShip.rho0     = 1.225;
paramsShip.H        = 8000;

% Multiple maneuvers for the Ship
% (Currently empty, but you can add them like so:)
paramsShip.maneuvers.pitch = [
    struct('start', 150, 'end', 165, 'M', -3000), 
    struct('start', 300, 'end', 315, 'M', 3000),
    struct('start', 350, 'end', 370, 'M', -3000)
];
paramsShip.maneuvers.thruster = [
    struct('start', 150, 'end', 8*60, 'T', 9e6)
];

% --------------------- Booster (Post-Separation) Parameters ---------------------
paramsBoosterPost.g        = stage1.g0;
paramsBoosterPost.burnTime = 0;     % no post-sep burn
paramsBoosterPost.T_engine = 0;     % no thrust
paramsBoosterPost.mdot     = 0;     % no mass flow
paramsBoosterPost.C_damp   = 1000;
paramsBoosterPost.I        = I_ship;  % same as ship
paramsBoosterPost.Cd       = stage2.Cd;
paramsBoosterPost.A        = stage2.A_ref;
paramsBoosterPost.rho0     = 1.225;
paramsBoosterPost.H        = 8000;

% Multiple maneuvers for the Booster after separation
paramsBoosterPost.maneuvers.pitch = [
    % struct('start', 200, 'end', 210, 'M', 1000),
    % struct('start', 300, 'end', 305, 'M', -500)
];
paramsBoosterPost.maneuvers.thruster = [
    % struct('start', 180, 'end', 190, 'T', 3e6)
];

% --------------------- Initial Conditions ---------------------
x0 = 0; y0 = 0; z0 = 0;
vx0 = 0; vy0 = 0; vz0 = 0;
phi0 = 0; theta0 = 0; psi0 = 0;
p0 = 0; q0 = 0; r0 = 0;
massCombined = stage1.m0 + stage2.m0;
state0_combined = [x0; y0; z0; vx0; vy0; vz0; phi0; theta0; psi0; p0; q0; r0; massCombined];

% --------------------- Simulation: Combined Phase ---------------------
indexSep = find(time >= stageSepTime, 1);
state_history_combined = zeros(13, indexSep);
state_history_combined(:,1) = state0_combined;
currentState = state0_combined;

for k = 1:indexSep-1
    t = time(k);
    currentState = rk4_step(@(t,s) dynamics(t, s, paramsCombined), t, currentState, dt);
    state_history_combined(:, k+1) = currentState;
end

% At separation, set ship mass
shipState0 = state_history_combined(:, end);
shipState0(13) = stage2.m0;

% Booster mass = same as ship's (just for demonstration)
boosterState0 = state_history_combined(:, end);
boosterState0(13) = stage2.m0;  

% --------------------- Simulation: Post-Separation ---------------------
% Ship
numShipSteps = N - indexSep + 1;
state_history_ship = zeros(13, numShipSteps);
state_history_ship(:,1) = shipState0;
currentShipState = shipState0;
for k = 1:(numShipSteps-1)
    t = time(indexSep + k - 1);
    currentShipState = rk4_step(@(t,s) dynamics(t, s, paramsShip), t, currentShipState, dt);
    state_history_ship(:, k+1) = currentShipState;
end

% Booster
numBoosterSteps = N - indexSep + 1;
state_history_booster = zeros(13, numBoosterSteps);
state_history_booster(:,1) = boosterState0;
currentBoosterState = boosterState0;
for k = 1:(numBoosterSteps-1)
    t = time(indexSep + k - 1);
    currentBoosterState = rk4_step(@(t,s) dynamics(t, s, paramsBoosterPost), t, currentBoosterState, dt);
    state_history_booster(:, k+1) = currentBoosterState;
end

% --------------------- Merge for Plotting ---------------------
full_state_ship    = [state_history_combined, state_history_ship(:,2:end)];
full_state_booster = [state_history_combined, state_history_booster(:,2:end)];
full_time = time;

% --------------------- Plotting ---------------------
figure('Color','w','Position',[0 0 1100 1100]);
plot3(full_state_booster(1,:)./1000, ...
      full_state_booster(2,:)./1000, ...
      full_state_booster(3,:)./1000, 'b','LineWidth',2);
hold on;
plot3(full_state_ship(1,:)./1000, ...
      full_state_ship(2,:)./1000, ...
      full_state_ship(3,:)./1000, 'r--','LineWidth',2);
xlabel('x [km]'); ylabel('y [km]'); zlabel('Altitude [km]');
title('Rocket Trajectories: Booster (blue) and Ship (red)');
grid minor; view(3);
ax = gca; ax.FontSize = 20;

figure('Color','w','Position',[800 0 1100 1100]);
speed_booster = sqrt(sum(full_state_booster(4:6,:).^2,1));
speed_ship    = sqrt(sum(full_state_ship(4:6,:).^2,1));
plot(full_time, speed_booster./1000, 'b','LineWidth',2);
hold on;
plot(full_time, speed_ship./1000, 'r--','LineWidth',2);
xlabel('Time [s]'); ylabel('Speed [km/s]');
title('Speed vs. Time');
grid minor;
ax = gca; ax.FontSize = 20;

figure('Color','w','Position',[1500 0 1100 1100]);
subplot(3,1,1);
plot(full_time, rad2deg(full_state_booster(7,:)), 'b','LineWidth',2);
hold on;
plot(full_time, rad2deg(full_state_ship(7,:)), 'r','LineWidth',2);
ylabel('Roll [deg]');
grid minor; ax = gca; ax.FontSize = 20;

subplot(3,1,2);
plot(full_time, rad2deg(full_state_booster(8,:)), 'b','LineWidth',2);
hold on;
plot(full_time, rad2deg(full_state_ship(8,:)), 'r','LineWidth',2);
ylabel('Pitch [deg]');
grid minor; ax = gca; ax.FontSize = 20;

subplot(3,1,3);
plot(full_time, rad2deg(full_state_booster(9,:)), 'b','LineWidth',2);
hold on;
plot(full_time, rad2deg(full_state_ship(9,:)), 'r','LineWidth',2);
xlabel('Time [s]'); ylabel('Yaw [deg]');
title('Euler Angles vs. Time');
grid minor;
sgtitle('Angles');
ax = gca; ax.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions

function s_next = rk4_step(dyn, t, s, dt)
    k1 = dyn(t, s);
    k2 = dyn(t + dt/2, s + dt/2 * k1);
    k3 = dyn(t + dt/2, s + dt/2 * k2);
    k4 = dyn(t + dt, s + dt * k3);
    s_next = s + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function ds = dynamics(t, s, params)
    % s = [x; y; z; vx; vy; vz; phi; theta; psi; p; q; r; m]
    pos      = s(1:3);
    vel      = s(4:6);
    phi      = s(7);
    theta    = s(8);
    psi      = s(9);
    angRates = s(10:12);
    m        = s(13);
    
    % Gravity
    F_grav = [0; 0; -m * params.g];
    
    % 1) Base Thrust from engine if within burnTime
    T_current = 0;
    if t <= params.burnTime
        T_current = params.T_engine;
    end
    
    % 2) Check if there are any thruster intervals that override
    if isfield(params, 'maneuvers') && isfield(params.maneuvers, 'thruster')
        thrusterList = params.maneuvers.thruster;
        for i = 1:length(thrusterList)
            if (t >= thrusterList(i).start) && (t <= thrusterList(i).end)
                % Override the thrust with the interval's value
                T_current = thrusterList(i).T;
            end
        end
    end
    
    % Thrust in body z-axis
    F_thrust_body = [0; 0; T_current];
    
    % Rotation matrix (Z-Y-X)
    R = euler_to_rotation(phi, theta, psi);
    F_thrust = R * F_thrust_body;
    
    % Aerodynamic drag
    V   = norm(vel);
    alt = pos(3);
    rho = params.rho0 * exp(-max(alt,0)/params.H);
    if V > 0
        F_drag = -0.5 * rho * V^2 * params.Cd * params.A * (vel / V);
    else
        F_drag = [0; 0; 0];
    end
    
    % Net force & acceleration
    F_total = F_grav + F_thrust + F_drag;
    accel   = F_total / m;
    
    % Moments: damping + pitch intervals
    M_aero    = -params.C_damp * angRates;
    M_pitch   = 0;  % accumulate from all pitch intervals
    if isfield(params, 'maneuvers') && isfield(params.maneuvers, 'pitch')
        pitchList = params.maneuvers.pitch;
        for i = 1:length(pitchList)
            if (t >= pitchList(i).start) && (t <= pitchList(i).end)
                M_pitch = M_pitch + pitchList(i).M;
            end
        end
    end
    
    M_control = [0; M_pitch; 0];
    M_total   = M_aero + M_control;
    
    % Angular acceleration
    omega     = angRates;
    omega_dot = params.I \ (M_total - cross(omega, params.I * omega));
    
    % Euler angle rates (Z-Y-X)
    p = omega(1); q = omega(2); r = omega(3);
    phi_dot   = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r;
    theta_dot = cos(phi)*q - sin(phi)*r;
    psi_dot   = (sin(phi)/cos(theta))*q + (cos(phi)/cos(theta))*r;
    
    % Mass flow
    m_dot = 0;
    if (t <= params.burnTime) && (m > 0)
        m_dot = -params.mdot;
    end
    
    ds = zeros(13,1);
    ds(1:3)   = vel;
    ds(4:6)   = accel;
    ds(7)     = phi_dot;
    ds(8)     = theta_dot;
    ds(9)     = psi_dot;
    ds(10:12) = omega_dot;
    ds(13)    = m_dot;
end

function R = euler_to_rotation(phi, theta, psi)
    cphi   = cos(phi);   sphi   = sin(phi);
    ctheta = cos(theta); stheta = sin(theta);
    cpsi   = cos(psi);   spsi   = sin(psi);
    
    R = [ cpsi*ctheta,  cpsi*stheta*sphi - spsi*cphi,   cpsi*stheta*cphi + spsi*sphi;
          spsi*ctheta,  spsi*stheta*sphi + cpsi*cphi,   spsi*stheta*cphi - cpsi*sphi;
         -stheta,       ctheta*sphi,                    ctheta*cphi ];
end