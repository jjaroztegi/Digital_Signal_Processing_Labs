%% Exercise 02

% Frequencies of the two signals
fc1 = 126;
fc2 = 722;

% Periods of the two signals
T1 = 1 / fc1;
T2 = 1 / fc2;
T = 3 * max(T1, T2);

% We need fc1/fm1 = fc2/fm2
% Set sampling rate
fm1 = 20 * fc1;  % 2520 Hz
fm2 = (fc2 / fc1) * fm1;  % 14440 Hz

% Continuous signals
syms t;
x1(t) = 5 * sin(2 * pi * fc1 * t);
x2(t) = 5 * sin(2 * pi * fc2 * t);

% Range of n
n1 = 0:T*fm1;
n2 = 0:T*fm2;

% Sampled signals
x1_n = 5 * sin(2 * pi * fc1 / fm1 * n1);
x2_n = 5 * sin(2 * pi * fc2 / fm2 * n2);

% Plot the sampled signals
figure;
hold on;

% Plot continuous-time signals
fplot(x1, [0 T], 'LineWidth', 1.5);
fplot(x2, [0 T], 'LineWidth', 1.5);

% Plot sampled signals
stem(n1/fm1, x1_n, 'r', 'filled');
stem(n2/fm2, x2_n, 'b', 'filled');

% Labels, title, and legend
xlabel('Time (s)');
ylabel('Amplitude');
title('Continuous and Sampled Signals');
legend('x_1(t)', 'x_2(t)', 'x_1[n]', 'x_2[n]');
grid on;

%% Plot vectors x1[n] and x2[n] along a sample number axis
figure;
subplot(2, 1, 1);
stem(n1, x1_n, 'r', 'filled');
title('Sampled Signal x_1[n]');
xlabel('Sample number');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
stem(n2, x2_n, 'k', 'filled');
title('Sampled Signal x_2[n]');
xlabel('Sample number');
ylabel('Amplitude');
grid on;

%{
Continuous and Sampled Signals: The first plot shows the continuous signals
alongside their sampled versions. You can see how well the sampling
captures the signal depending on the rate.

As both signals have been sampled at their respective sampling rates,
aliasing should not be an issue if the sampling frequency is well above the
Nyquist rate. Even though the two signals have different frequencies, they
can still have similar sample points depending on how their sampling rates
are related.
%}