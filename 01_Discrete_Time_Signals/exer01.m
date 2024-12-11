%% Exercise 01

% Define constants
a = 0.8162 + 0.4288i;

% Anonymous function for unit step
u = @(n) (n >= 0)*1.0;


% Define the signal x[n] = a^n * u[n+3]
x_n = @(n) (a.^n) .* u(n-3);

% Range of n
n = -3:50;

% Plot Real, Imaginary, Magnitude, and Phase

figure;

% Real part
subplot(2,2,1);
stem(n, real(x_n(n)), 'filled');
title('Real Part of x[n]');
xlabel('n');
ylabel('Re\{x[n]\}');
grid on;

% Imaginary part
subplot(2,2,2);
stem(n, imag(x_n(n)), 'filled');
title('Imaginary Part of x[n]');
xlabel('n');
ylabel('Im\{x[n]\}');
grid on;

% Magnitude
subplot(2,2,3);
stem(n, abs(x_n(n)), 'filled');
title('Magnitude of x[n]');
xlabel('n');
ylabel('|x[n]|');
grid on;

% Phase
subplot(2,2,4);
stem(n, angle(x_n(n)), 'filled');
title('Phase of x[n]');
xlabel('n');
ylabel('Phase(x[n])');
grid on;

%% Determine Discrete Frequency and Exponential Decay Rate

alpha = abs(a);               % Exponential decay rate
freq_d = angle(a)/(2*pi);     % Discrete frequency

% Sampling frequency and analog frequency calculation
Fs = 125000;                    % Sampling frequency
f_analog = freq_d * Fs ;        % Analog frequency

% Display results
disp(['Exponential decay rate (alpha): ', num2str(alpha)]);
disp(['Discrete frequency (freq_d): ', num2str(freq_d)]);
disp(['Analog frequency in Hz: ', num2str(f_analog)]);

