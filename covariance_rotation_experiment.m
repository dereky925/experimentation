clear; clc; close all;

% Define the mean vector for 3 dimensions
mu = [1 1 1];

% Define a 3x3 covariance matrix
Sigma = [1 1.5 0.8;
         1.5 3 1.2;
         0.8 1.2 2];

Sigma = [10 0 0;
         0 1 0;
         0 0 1];

% Rotate covariance by X degrees with pre and post multiplication
rotx = 45;
rotz = 90;

% Rotate about y axis
Sigma2 = R2(deg2rad(rotx)) * Sigma * R2(deg2rad(rotx))';

% Rotate about z axis
Sigma2 = R3(deg2rad(rotz)) * Sigma2 * R3(deg2rad(rotz))';

% Set random seed for reproducibility
rng('default')

% Generate 1000 samples from a 3D multivariate normal distribution
R = mvnrnd(mu, Sigma, 10000);
R2 = mvnrnd(mu, Sigma2, 10000);

% Plot in 3D
figA = figure(1);
hold on
grid minor
scatter3(R(:,1), R(:,2), R(:,3), '+')
scatter3(R2(:,1), R2(:,2), R2(:,3), '+')
figA.Position = [10 10 1000 1000];

% Set axis labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot of Multivariate Normal Samples');
view(45, 30);  % Adjusts to a 3D viewing angle with azimuth 45° and elevation 30°
ax = gca;   
ax.FontSize = 16;
axis equal

% Display mean of each dimension to verify
disp('Mean of each dimension:')
mean(R)