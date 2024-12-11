%% Exercise 2
clearvars;
close all;
clc;

%% 1
% Coefficients of the Z-transform
num = [1 0 0];                     % Numerator: z^2
den = [2 -1 3];                    % Denominator: 2z^2 - z + 3

% Poles and zeros
[z, p, k] = tf2zp(num, den);

% Plot the pole-zero diagram
figure;
zplane(z, p);
title('Pole-Zero Diagram of X_1(z)');
xlabel('Real Part');
ylabel('Imaginary Part');
grid on;

% ROC for anti-causal
disp('Region of Convergence (ROC):');
disp('ROC: |z| > r_max = 1.225, ');

% Residues
[res, pk, ~] = residuez(num, den);

% Inverse Z-Transform (using residues and poles)
syms z n;
X1_z = sum(res ./ (1 - pk * z^(-1)));

x_n = iztrans(X1_z, z, n);

n_values = -20:20;
x_n_numeric = double(subs(x_n, n, n_values));

% Plot the inverse Z-transform
figure;
stem(x_n_numeric, x_n_numeric, 'filled');
title('Inverse Z-Transform x_1[n]');
xlabel('n');
ylabel('x[n]');
grid on;

%% 2

num = [0 0 1];
den = conv([1 -0.9], conv([1 0.65], [1 0.7 0.7]));  % Denominator: expanded form

% Poles and zeros
[z, p, k] = tf2zp(num, den);

% Plot the pole-zero diagram
figure;
zplane(z, p);
title('Pole-Zero Diagram of X_2(z)');
xlabel('Real Part');
ylabel('Imaginary Part');
grid on;

% ROC for non-causal signals (|z| < r_min)
disp('Region of Convergence (ROC):');
disp('ROC: |z| < r_min = 0.65');

% Residues
[res, pk, ~] = residuez(num, den);

% Inverse Z-Transform (using residues and poles)
syms z n;
X2_z = sum(res ./ (1 - pk * z^(-1)));

x_n = iztrans(X2_z, z, n);

n_values = -20:20;
x_n_numeric = double(subs(x_n, n, n_values));

% Plot the inverse Z-Transform
figure;
stem(n_values, x_n_numeric, 'filled');
title('Inverse Z-Transform x_2[n]');
xlabel('n');
ylabel('x[n]');
grid on;

%% 3

num = [1 2];                      % Numerator: 1 + 2z^(-1)
den = [1 -0.8 0.3];               % Denominator: 1 - 0.8z^(-1) + 0.3z^(-2)

% Poles and zeros
[z, p, k] = tf2zp(num, den);

% Plot the pole-zero diagram
figure;
zplane(z, p);
title('Pole-Zero Diagram of X_3(z)');
xlabel('Real Part');
ylabel('Imaginary Part');
grid on;

% ROC for causal signals (|z| > r_max)
disp('Region of Convergence (ROC):');
disp('ROC: |z| > r_max = 0.5477');

% Residues
[res, pk, ~] = residuez(num, den);

% Inverse Z-Transform (using residues and poles)
syms z n;
X3_z = sum(res ./ (1 - pk * z^(-1)));

x_n = iztrans(X3_z, z, n);

n_values = -20:0;
x_n_numeric = double(subs(x_n, n, n_values));

% Plot the inverse Z-transform
figure;
stem(n_values, x_n_numeric, 'filled');
title('Inverse Z-Transform x_3[n]');
xlabel('n');
ylabel('x[n]');
grid on;