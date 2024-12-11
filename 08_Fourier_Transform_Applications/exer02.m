clearvars;
close all;
clc;

% Load the audio file
[y, fm] = audioread('audio_file.wav');
[num_samples, num_channels] = size(y);


% Parameters
window_size = 8192;
num_windows = floor(length(y)/window_size);
signal_percent = 0.999;

% Initialize arrays
compressed_signal = zeros(size(y));
num_coeffs_kept = zeros(num_windows, 1);
compression_ratios = zeros(num_windows, 1);

% Process each channel
for ch = 1:num_channels
    % Process each window in current channel
    for i = 1:num_windows
        % Extract window
        start_idx = (i-1)*window_size + 1;
        end_idx = i*window_size;
        window_data = y(start_idx:end_idx, ch);
        
        % Calculate FFT
        window_fft = fft(window_data);
        
        % Calculate power spectrum
        power_spectrum = abs(window_fft).^2;
        total_power = sum(power_spectrum);
        
        % Sort coefficients by magnitude
        [sorted_power, sort_idx] = sort(power_spectrum, 'descend');
        cumulative_power = cumsum(sorted_power);
        normalized_power = cumulative_power / total_power;
        
        % Find number of coefficients needed for 99.9% energy
        num_coeffs = find(normalized_power >= signal_percent, 1, 'first');
        num_coeffs_kept(i, ch) = num_coeffs;
        
        % Create new FFT array with only significant coefficients
        window_fft_compressed = zeros(size(window_fft));
        significant_indices = sort_idx(1:num_coeffs);
        window_fft_compressed(significant_indices) = window_fft(significant_indices);
        
        % Reconstruct window
        window_reconstructed = real(ifft(window_fft_compressed));
        compressed_signal(start_idx:end_idx, ch) = window_reconstructed;
        
        % Calculate compression ratio for this window
        compression_ratios(i, ch) = window_size / (num_coeffs * 2);
    end
end


% Plot results
figure(1);

% Plot a segment of the signal (first 3 windows)
plot_length = min(3*window_size, num_samples);
t = (1:plot_length)/fm;

% Left channel
subplot(2,2,1);
plot(t, y(1:plot_length, 1));
title('Original Signal - Left Channel');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2,2,2);
plot(t, y(1:plot_length, 2));
title('Original Signal - Right Channel');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2,2,3);
plot(t, compressed_signal(1:plot_length, 1));
title('Compressed Signal - Left Channel');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2,2,4);
plot(t, compressed_signal(1:plot_length, 2));
title('Compressed Signal - Right Channel');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Display compression statistics for both channels
fprintf('Left Channel:\n');
fprintf('Average coefficients kept: %.2f out of %d\n', ...
    mean(num_coeffs_kept(:,1)), window_size);
fprintf('Average compression ratio: %.2f\n\n', mean(compression_ratios(:,1)));

fprintf('Right Channel:\n');
fprintf('Average coefficients kept: %.2f out of %d\n', ...
    mean(num_coeffs_kept(:,2)), window_size);
fprintf('Average compression ratio: %.2f\n', mean(compression_ratios(:,2)));

% Plot spectrum of a single window for both channels
window_idx = 1;
start_idx = (window_idx-1)*window_size + 1;
end_idx = window_idx*window_size;
f = (0:window_size-1)*(fm/window_size);

figure(3);
% Left channel spectrum
subplot(2,2,1);
plot(f, abs(fft(y(start_idx:end_idx, 1))));
title('Original Spectrum - Left Channel');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,2);
plot(f, abs(fft(y(start_idx:end_idx, 2))));
title('Original Spectrum - Right Channel');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,3);
plot(f, abs(fft(compressed_signal(start_idx:end_idx, 1))));
title('Compressed Spectrum - Left Channel');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,4);
plot(f, abs(fft(compressed_signal(start_idx:end_idx, 2))));
title('Compressed Spectrum - Right Channel');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;