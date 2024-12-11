clearvars;
close all;
clc;

% Radar parameters
fs = 1.023e6;               % Sampling frequency (Hz)
pulse_length = 1023;        % Pulse length (samples)
pulse_duration = 1e-3;      % Pulse duration (seconds)
N = 20;                     % Number of chirp pulses
duty_cycle = 0.02;          % Duty cycle
alpha = 0.1;                % Channel loss factor
k = 20;                     % Delay in samples

% Generate noise
s = cacode(6, 1);           % CA code for PRN 6
s = 2 * s - 1;              % Map '1' to -1 and '0' to 1

% Generate signal
pulsed_signal = zeros(1, N * (pulse_length + round((1 - duty_cycle) * pulse_length))); 
for i = 0:N-1
    start_index = i * (pulse_length + round((1 - duty_cycle) * pulse_length)) + 1;
    pulsed_signal(start_index:start_index + pulse_length - 1) = s;
end

% Simulate received signal y[n] (without noise)
y = zeros(1, length(pulsed_signal));
y(k + 1:k + pulse_length) = alpha * s;  % Pulse delayed by k

% Step 1: Cross-correlation of y[n] and s[n]
[r, lags] = xcorr(y, s);

% Estimate the delay from the peak value
[~, max_idx] = max(abs(r));         % Find index of peak
estimated_delay = lags(max_idx);

fprintf('Estimated delay without noise: %d samples\n', estimated_delay);

% Step 2: Add Noise
SNR_dB = -25;

% Signal power: ignore zero values
signal_power = var(alpha * s);  % Power of the clean signal (scaled by alpha)
noise_power = signal_power / (10^(SNR_dB / 10));
noise = sqrt(noise_power) * randn(1, pulse_length);  % Generate AWGN

y_noisy = y;
y_noisy(k + 1:k + pulse_length) = y(k + 1:k + pulse_length) + noise;

% Cross-correlation but with noise signal
[r_noisy, lags_noisy] = xcorr(y_noisy, s);  
[~, max_idx_noisy] = max(abs(r_noisy));

% Mean and standard deviation
num_trials = 100;
delays = zeros(1, num_trials);
max_lag = pulse_length;

for i = 1:num_trials
    % Generate new noise
    noise = sqrt(noise_power) * randn(1, pulse_length);  % Generate AWGN
    
    % Add noise
    y_noisy = y;
    y_noisy(k + 1:k + pulse_length) = y(k + 1:k + pulse_length) + noise;
    
    [r_noisy, lags_noisy] = xcorr(y_noisy, s, max_lag);
    
    [~, max_idx_noisy] = max(r_noisy);
    delays(i) = lags_noisy(max_idx_noisy);
end

mean_delay = mean(delays);
std_delay = std(delays);

fprintf('Mean estimated delay with noise: %.2f samples\n', mean_delay);
fprintf('Standard deviation of estimated delay: %.2f samples\n', std_delay);

% Step 3: Plot SNR - standard deviation
SNRs = -43:1:-10;
std_devs = zeros(size(SNRs));

for idx = 1:length(SNRs)
    SNR_dB = SNRs(idx);
    noise_power = signal_power / (10^(SNR_dB / 10));
    
    delays = zeros(1, num_trials);
    
    for i = 1:num_trials
        % Add noise
        noise = sqrt(noise_power) * randn(size(y));
        y_noisy = y + noise;
        % Cross-correlation
        [r_noisy, lags_noisy] = xcorr(y_noisy, s);
        [~, max_idx_noisy] = max(abs(r_noisy));
        delays(i) = lags_noisy(max_idx_noisy);
    end
    
    % Standard deviation of the delays
    std_devs(idx) = std(delays);
end

% Plot standard deviation vs SNR
figure;
plot(SNRs, std_devs, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Standard Deviation of Delay Estimates');
title('Standard Deviation of Delay Estimates vs. SNR');
grid on;