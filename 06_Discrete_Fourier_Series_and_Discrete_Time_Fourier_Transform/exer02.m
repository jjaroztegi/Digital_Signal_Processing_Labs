clearvars;
close all;

%% Exercise 2 a)
load din.dat;

N = 2048;
fs = 4096;
x = transpose(din);
n = 0:length(din)-1;

% Frequency vector for DTFT (centered around zero)
nn = (-N/2:N/2-1) * fs / N;  % Frequency vector in Hz
X = zeros(1, N);  % Preallocate X for DTFT coefficients

% DTFT
for i = 1:N
    F = (i - 1) / N;  % Normalized frequency
    X(i) = sum(x .* exp(-1i * 2 * pi * F * n));
end

% Normalize and convert to dB
Xs = 20 * log10(abs(X) / N);  % Magnitude in dB, normalized by N

% Plot DTFT results
figure;
stem(nn, [Xs(N/2:N) Xs(1:N/2-1)], 'filled', 'b');  % Plot DTFT (blue)
hold on;
title('20·log|X(F)|');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;

%% Exercise 2 b)

% Initialize the reconstructed signal
reconstructed_signal = zeros(size(x));  % Same size as the original signal

% Compute the time-domain signal using the IDTFT
for t = 1:length(reconstructed_signal)
    for k = 1:N
        F = (k - 1) / N;  % Normalized frequency
        reconstructed_signal(t) = reconstructed_signal(t) + (X(k) / N) * exp(1i * 2 * pi * F * (t - 1));
    end
end

% Convert to real
reconstructed_signal = real(reconstructed_signal);

% Plot original and reconstructed signals
t = (0:length(din)-1) / fs;  % Time vector for original signal

figure;
plot(t, din, 'b', 'LineWidth', 1.5);  % Plot original signal in blue
hold on;
plot(t, reconstructed_signal, 'r--', 'LineWidth', 1.5);  % Plot reconstructed signal in red dashed line
hold off;

% Add labels and title
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal vs Reconstructed Signal');
legend({'Original Signal', 'Reconstructed Signal'});
grid on;

%% Exercise 2 c)

% The ability of the IDTFT to accurately reconstruct the original signal
% depends on various factors, including the sampling frequency, the choice
% of Fourier coefficients, and the inherent characteristics of the signal.
% If the reconstructed signal closely matches the original, it demonstrates
% the effectiveness of the Fourier transform in capturing the essential
% features of the signal. On the other hand, noticeable discrepancies may
% indicate the need for adjustments in sampling rates, windowing
% techniques, or the number of points used in the transformation to enhance
% reconstruction fidelity.