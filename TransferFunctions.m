%% Transfer functions

clear;clc;close all

% Scaling the numerator directly scales the output
num1 = 1;
num2 = 2;
den = [1 1];

G1 = tf(num1, den);
G2 = tf(num2, den);

figure('Color','white','Position',[0 0 1200 1000]);
hold on; grid minor;
step(G1, 'b', G2, 'r--'); 
legend('G(s) = 1/(s+1)', 'G(s) = 2/(s+1)');
title('Effect of Scaling Numerator');
grid minor;

%%

clear;clc;close all

num1 = 1;
num2 = [1 2];  % Introduces a zero at s = -2
den = [1 1];

G1 = tf(num1, den);
G2 = tf(num2, den);

figure('Color','white','Position',[0 0 1200 1000]);
step(G1, 'b', G2, 'r--'); 
legend('G(s) = 1/(s+1)', 'G(s) = (s+2)/(s+1)');
title('Effect of Adding a Zero in the Numerator');
grid minor;

%%
clear;clc;close all

num1 = 1;
num2 = [1 3 2];  % Quadratic numerator
den = [1 1 1];

G1 = tf(num1, den);
G2 = tf(num2, den);


figure('Color','white','Position',[0 0 1200 1000]);
step(G1, 'b', G2, 'r--'); 
legend('G(s) = 1/(s+1)', 'G(s) = (s^2 + 3s + 2)/(s+1)');
title('Effect of Quadratic Numerator');
grid minor;

