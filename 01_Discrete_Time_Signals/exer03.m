%% Exercise 03

% Define the continuous-time signal parameters
f0 = 50;
t_max = 0.12;
         
% Define symbolic variable for time
syms t;
syms x(t);

x(t) = sin(2*pi*f0*t);


% Define sampling parameters
fs = 550;               % Discrete sampling frequency (550 Hz)
F = f0/fs;
n_vals = 0:t_max*fs;    % Discrete time samples from 0 to 120ms

% Sample the continuous-time signal at discrete time points
x_n = sin(2*pi*F*n_vals);

% Plot both signals on the same graph
figure;
subplot(2,1,1);
hold on;
fplot(x(t), [0 t_max]);
stem(n_vals/fs, x_n, 'r', 'filled'); 
title('x(t) and Sampled Signal x[n] at 550Hz');
xlabel('t (s)');
ylabel('Amplitude');
legend('Continuous-Time Signal x(t)', 'Sampled Signal x[n]');
grid on;

% Period
[num, den] = rat(F);
disp(['Period of x[n]: ', num2str(den)]);

%%
fs = 80;
F = f0/fs;
n_vals = 0:t_max*fs;

% Sample the continuous-time signal
x_n = sin(2*pi*F*n_vals);
% Plot both signals
subplot(2,1,2);
hold on;
fplot(x(t), [0 t_max]);
stem(n_vals/fs, x_n, 'r', 'filled'); 
title('x(t) and Sampled Signal x[n] at 80Hz');
xlabel('t (s)');
ylabel('Amplitude');
legend('Continuous-Time Signal x(t)', 'Sampled Signal x[n]');
grid on;

% Period
[num, den] = rat(F);
disp(['Period of x[n]: ', num2str(den)]);