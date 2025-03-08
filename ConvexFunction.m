% In a convex function, for any two points, the function lies below the line (or chord) connecting them
% Clear the workspace and close any open figures
clear; close all; clc;

% Define the function f(x) = x^2
f = @(x) x.^2;

% Define the range for x
x_min = -2.5;
x_max = 1.5;
x = linspace(x_min, x_max, 400);

% Compute the function values
y = f(x);

% Choose two points on the function (for example, x1 = -2 and x2 = 1)
x1 = -2;
x2 = 1;
y1 = f(x1);
y2 = f(x2);

% Create a vector for the chord between (x1,y1) and (x2,y2)
% For any x between x1 and x2, the chord value is given by linear interpolation
x_chord = linspace(x1, x2, 200);
y_chord = ((x_chord - x1)/(x2 - x1)) * y2 + ((x2 - x_chord)/(x2 - x1)) * y1;

% Plot the function
figure('Color','w','Position',[100 100 800 600]);
plot(x, y, 'b-', 'LineWidth',2); hold on;
grid on; 


% Plot the chord between the two chosen points
plot(x_chord, y_chord, 'r--', 'LineWidth',2);

% Highlight the chosen points on the function
plot(x1, y1, 'ko', 'MarkerSize',8, 'MarkerFaceColor','k');
plot(x2, y2, 'ko', 'MarkerSize',8, 'MarkerFaceColor','k');

% Add labels and title
xlabel('x', 'FontSize',14);
ylabel('f(x) = x^2', 'FontSize',14);
title('Demonstration of Convexity: f(x) lies below the chord connecting any two points', 'FontSize',16);

% Add legend to the plot
legend('f(x)=x^2','Chord between (x1,y1) and (x2,y2)', 'Selected Points', 'Location','Best');

% Annotate the selected points
text(x1, y1, sprintf('  (%.1f, %.1f)', x1, y1), 'FontSize',12);
text(x2, y2, sprintf('  (%.1f, %.1f)', x2, y2), 'FontSize',12);

hold off;
ax = gca;
ax.FontSize = 20;