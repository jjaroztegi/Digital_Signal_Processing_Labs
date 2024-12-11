%% Exercise 1

% Define the parameters
fs = 2150;
T_s = 1/fs;

% Original frequencies
f1 = 1250;
f2 = 1775;
f3 = 775;

% 1. Calculation of the Period
% F = f1/fs = 1250/2150 = 25/43 --> N=43;
% F = f2/fs = 1775/2150 = 71/86 --> N=86;
% F = f3/fs = 775/2150 = 31/86 --> N=86;
% mcm (43, 86, 86) = 86
N = 86;

% 2. Determine discrete frecuencies
F1 = f1/fs;
F2 = f2/fs;
F3 = f3/fs;

% 3. Aliasing
% Aliasing occurs when a frequency in the analog signal exceeds the Nyquist
% frequency, which in this case is fs/2 = 1075Hz
% So, the frecuency components 1250 and 1775 will be aliased.

% 4. Alias frecuencies
F1_alias = F1 - round(F1);
F2_alias = F2 - round(F2);
F3_alias = F3 - round(F3);
% Display the results
disp(['F1 original: ', num2str(F1), ', F1 alias: ', num2str(F1_alias)]);
disp(['F2 original: ', num2str(F2), ', F2 alias: ', num2str(F2_alias)]);
disp(['F3 original: ', num2str(F3), ', F3 alias: ', num2str(F3_alias)]);

% Time vector for plotting one full period
t = 0:0.00001:(N-1)*T_s;

% Define the analog signal x(t)
x_t = @(t) sin(2*pi*f1*t) + cos(2*pi*f2*t + pi/2) - sin(2*pi*f3*t);

% Sampled signal x[n]
n = 0:N-1;
x_n = @(n) sin(2*pi*F1*n) + cos(2*pi*F2*n + pi/2) - sin(2*pi*F3*n);
x_n_alias = @(n) sin(2*pi*F1_alias*n) + cos(2*pi*F2_alias*n + pi/2) ...
                 - sin(2*pi*F3_alias*n);

% Plotting the continuous and sampled signals
figure;
plot(t, x_t(t), 'b', 'DisplayName', 'x(t)');
hold on;
stem(n*T_s, x_n(n), 'r', 'DisplayName', 'x[n]', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
stem(n*T_s, x_n_alias(n), 'g--', 'DisplayName', 'x[n] alias freq', 'MarkerFaceColor', 'r');
xlabel('Time (seconds)');
ylabel('Amplitude');
legend;
title('Continuous and Sampled Signals');
