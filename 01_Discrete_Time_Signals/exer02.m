%% Exercise 02

% Anonymous function for unit step
u = @(n) (n >= 0)*1.0;  

% Define range of n
n = -10:20;  

% Define the sampled signal x[n]
x_n = @(n) (0.75.^n) .* sin(2*pi*n/5) .* (u(n - 7) - u(n + 7));

% Plot original signal x[n]
figure;
subplot(2,1,1);
stem(n, x_n(n), 'filled');
title('Original Signal x[n]');
xlabel('n');
ylabel('x[n]');
grid on;

% Define the transformation xb[n] = 2x[n + 1] + x[n - 3] - x[2n] + x[5 - n]

x_b = @(n) 2*x_n(n + 1) + x_n(n - 3) - x_n(2*n) + x_n(5 - n); 

% Plot transformed signal xb[n]
subplot(2,1,2);
stem(n, x_b(n), 'filled');
title('Transformed Signal xb[n]');
xlabel('n');
ylabel('xb[n]');
grid on;

%% Energy and Power Calculation for Transformed Signal xb[n]
energy_xb = sum(abs(x_b(n)).^2);  % Energy of the transformed signal xb[n]
power_xb = 0;                     % Energy signals have null power

% Display the energy and power of the transformed signal
disp(['Energy of xb[n]: ', num2str(energy_xb)]);
disp(['Power of xb[n]: ', num2str(power_xb)]);